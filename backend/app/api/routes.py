from fastapi import APIRouter, WebSocket, WebSocketDisconnect, HTTPException
from fastapi.responses import Response, StreamingResponse
from typing import Optional
from app.models.molecular import (
    StructureRequest, StructureResponse, ChatMessage, ChatResponse
)
from app.services.molecular_service import MolecularService
from app.services.ai_service import AIService
from app.services.structure_generator import StructureGenerator
from app.services.reaction_service import ReactionService
from app.services.pubchem_service import PubChemService
from app.tasks.molecular_tasks import (
    task_smiles_to_structure,
    task_process_chat,
    task_generate_structure,
    task_optimize_structure,
    task_calculate_properties,
    task_validate_structure,
)
from app.celery_app import celery_app
from sse_starlette.sse import EventSourceResponse
import json
import redis
import asyncio
import os
from dotenv import load_dotenv
from rdkit import Chem
from io import BytesIO

load_dotenv()

router = APIRouter()
molecular_service = MolecularService()
ai_service = AIService()
structure_generator = StructureGenerator()
reaction_service = ReactionService()
pubchem_service = PubChemService()

@router.post("/structure/from-smiles", response_model=StructureResponse)
async def structure_from_smiles(request: StructureRequest, async_task: bool = False):
    """Convert SMILES to 3D structure. Set async_task=true for long-running operations."""
    if not request.smiles:
        raise HTTPException(status_code=400, detail="SMILES string is required")
    
    if async_task:
        # Use Celery for async processing
        task = task_smiles_to_structure.delay(request.smiles)
        return {"task_id": task.id, "status": "processing", "message": "Task submitted"}
    
    # Synchronous processing for quick operations
    structure = molecular_service.smiles_to_structure(request.smiles)
    if not structure:
        raise HTTPException(status_code=400, detail="Invalid SMILES string")
    
    properties = molecular_service.calculate_properties(request.smiles)
    
    return StructureResponse(
        structure=structure,
        molecular_weight=properties.get("molecular_weight"),
        formula=properties.get("formula"),
        properties=properties
    )

@router.post("/structure/to-smiles")
async def structure_to_smiles(structure: StructureRequest):
    """Convert structure to SMILES."""
    if not structure.structure:
        raise ValueError("Structure is required")
    
    smiles = molecular_service.structure_to_smiles(structure.structure)
    if not smiles:
        raise ValueError("Could not convert structure to SMILES")
    
    return {"smiles": smiles}

@router.post("/structure/properties")
async def get_properties(request: StructureRequest, async_task: bool = False):
    """Get properties of a structure. Set async_task=true for large molecules."""
    if async_task:
        # Use Celery for async property calculation
        structure_dict = request.structure.dict() if request.structure else None
        task = task_calculate_properties.delay(request.smiles, structure_dict)
        return {"task_id": task.id, "status": "processing", "message": "Property calculation task submitted"}
    
    # Synchronous processing
    if request.smiles:
        properties = molecular_service.calculate_properties(request.smiles)
    elif request.structure:
        smiles = molecular_service.structure_to_smiles(request.structure)
        if smiles:
            properties = molecular_service.calculate_properties(smiles)
        else:
            raise HTTPException(status_code=400, detail="Could not convert structure to SMILES")
    else:
        raise HTTPException(status_code=400, detail="SMILES or structure is required")
    
    return properties

@router.post("/structure/validate")
async def validate_structure(structure: StructureRequest, async_task: bool = False):
    """Validate a molecular structure. Set async_task=true for complex structures."""
    if not structure.structure:
        raise HTTPException(status_code=400, detail="Structure is required")
    
    if async_task:
        # Use Celery for async validation
        task = task_validate_structure.delay(structure.structure.dict())
        return {"task_id": task.id, "status": "processing", "message": "Validation task submitted"}
    
    # Synchronous processing
    is_valid, error = molecular_service.validate_structure(structure.structure)
    return {"is_valid": is_valid, "error": error}

@router.post("/chat")
async def chat(message: ChatMessage, async_task: bool = True):
    """Process chat message with AI. Uses async task only for AI API calls."""
    # Check if this is a simple request that can be handled instantly
    message_lower = message.message.lower().strip()
    common_molecules = ['water', 'h2o', 'methane', 'ch4', 'ethane', 'ethanol', 'benzene', 'aspirin', 'caffeine']
    
    # For common molecule requests, handle synchronously (instant)
    if any(mol in message_lower for mol in common_molecules) and any(word in message_lower for word in ['create', 'generate', 'make', 'show', 'build']):
        response = ai_service.process_chat_message(message)
        # If it's a simple molecule creation, return immediately
        if response.structure and not response.response.startswith("I can help"):
            return response.dict()
    
    # For AI-powered responses, use async task
    if async_task:
        structure_dict = message.structure.dict() if message.structure else None
        task = task_process_chat.delay(message.message, structure_dict)
        return {"task_id": task.id, "status": "processing", "message": "Chat task submitted"}
    
    # Synchronous processing (fallback)
    response = ai_service.process_chat_message(message)
    return response.dict()

@router.post("/structure/generate")
async def generate_structure(description: dict, async_task: bool = True):
    """
    Generate molecular structure from natural language description.
    Priority: 1) PubChem (FIRST), 2) Common molecules, 3) AI (last).
    """
    description_text = description.get('description', '')
    if not description_text:
        raise HTTPException(status_code=400, detail="Description is required")
    
    # PRIORITY 1: Try PubChem FIRST (before anything else)
    import re
    molecule_name = re.sub(r'\b(create|generate|make|show|build|draw|display)\b', '', description_text, flags=re.IGNORECASE).strip()
    molecule_name = molecule_name.strip('.,!?').strip()
    
    if molecule_name and len(molecule_name) > 2:
        try:
            pubchem_result = pubchem_service.search_by_name(molecule_name)
            if pubchem_result and pubchem_result.get('smiles'):
                structure = molecular_service.smiles_to_structure(pubchem_result['smiles'])
                if structure:
                    return {
                        "structure": structure.dict(),
                        "source": "pubchem",
                        "cid": pubchem_result.get('cid'),
                        "pubchem_name": pubchem_result.get('name')
                    }
        except Exception as e:
            print(f"PubChem search failed: {e}")
            # Continue to next priority
    
    # PRIORITY 2: Try common molecules and PubChem via ai_service (which also checks PubChem)
    structure = ai_service.generate_structure_from_text(description_text)
    if structure:
        # If we got here and PubChem didn't work above, it's from common molecules
        return {"structure": structure.dict(), "source": "fallback"}
    
    # PRIORITY 3: Use async task for AI-powered generation (only if PubChem and common molecules failed)
    if async_task:
        task = task_generate_structure.delay(description_text)
        return {"task_id": task.id, "status": "processing", "message": "Searching PubChem and AI..."}
    
    # Synchronous AI generation (fallback)
    structure = structure_generator.generate_from_text(description_text)
    if not structure:
        raise HTTPException(status_code=400, detail="Could not generate structure from description")
    
    return {"structure": structure.dict(), "source": "ai"}

@router.post("/structure/optimize")
async def optimize_structure(structure: StructureRequest, async_task: bool = True):
    """Optimize molecular structure. Uses async task by default for CPU-intensive operations."""
    if not structure.structure:
        raise HTTPException(status_code=400, detail="Structure is required")
    
    if async_task:
        # Use Celery for async optimization
        task = task_optimize_structure.delay(structure.structure.dict())
        return {"task_id": task.id, "status": "processing", "message": "Optimization task submitted"}
    
    # Synchronous processing (fallback)
    optimized = structure_generator.optimize_structure(structure.structure)
    return {"structure": optimized}

@router.get("/tasks/{task_id}")
async def get_task_status(task_id: str):
    """Get the status of a Celery task."""
    task = celery_app.AsyncResult(task_id)
    
    if task.state == "PENDING":
        response = {
            "task_id": task_id,
            "status": "pending",
            "state": task.state,
        }
    elif task.state == "PROGRESS":
        response = {
            "task_id": task_id,
            "status": "processing",
            "state": task.state,
            "current": task.info.get("current", 0),
            "total": task.info.get("total", 1),
            "message": task.info.get("message", ""),
        }
    elif task.state == "SUCCESS":
        response = {
            "task_id": task_id,
            "status": "completed",
            "state": task.state,
            "result": task.result,
        }
    else:  # FAILURE or other states
        response = {
            "task_id": task_id,
            "status": "failed",
            "state": task.state,
            "error": str(task.info) if task.info else "Unknown error",
        }
    
    return response


@router.get("/tasks/{task_id}/stream")
async def stream_task_progress(task_id: str):
    """Stream task progress updates via Server-Sent Events."""
    redis_url = os.getenv("CELERY_BROKER_URL", "redis://localhost:6379/0")
    redis_client = redis.from_url(redis_url, decode_responses=True)
    pubsub = redis_client.pubsub()
    
    # Subscribe to task progress channel
    channel = f"task_progress:{task_id}"
    pubsub.subscribe(channel)
    
    async def event_generator():
        try:
            # Send initial connection message
            yield {
                "event": "connected",
                "data": json.dumps({"task_id": task_id, "status": "connected"})
            }
            
            # Also check current task status
            task = celery_app.AsyncResult(task_id)
            if task.state == "SUCCESS":
                yield {
                    "event": "progress",
                    "data": json.dumps({
                        "task_id": task_id,
                        "status": "completed",
                        "state": task.state,
                        "result": task.result,
                    })
                }
                return
            elif task.state == "FAILURE":
                yield {
                    "event": "progress",
                    "data": json.dumps({
                        "task_id": task_id,
                        "status": "failed",
                        "state": task.state,
                        "error": str(task.info) if task.info else "Unknown error",
                    })
                }
                return
            
            # Listen for progress updates
            while True:
                message = pubsub.get_message(timeout=1.0)
                if message and message["type"] == "message":
                    try:
                        data = json.loads(message["data"])
                        yield {
                            "event": "progress",
                            "data": json.dumps(data)
                        }
                        # If task is complete or failed, close the stream
                        if data.get("status") in ["completed", "failed"]:
                            break
                    except json.JSONDecodeError:
                        continue
                
                # Also check task status periodically as fallback
                task = celery_app.AsyncResult(task_id)
                if task.state == "SUCCESS":
                    yield {
                        "event": "progress",
                        "data": json.dumps({
                            "task_id": task_id,
                            "status": "completed",
                            "state": task.state,
                            "result": task.result,
                        })
                    }
                    break
                elif task.state == "FAILURE":
                    yield {
                        "event": "progress",
                        "data": json.dumps({
                            "task_id": task_id,
                            "status": "failed",
                            "state": task.state,
                            "error": str(task.info) if task.info else "Unknown error",
                        })
                    }
                    break
                
                await asyncio.sleep(0.1)
        except asyncio.CancelledError:
            pass
        finally:
            pubsub.unsubscribe(channel)
            pubsub.close()
            redis_client.close()
    
    return EventSourceResponse(event_generator())


@router.websocket("/ws/chat")
async def websocket_chat(websocket: WebSocket):
    """WebSocket endpoint for real-time chat."""
    await websocket.accept()
    try:
        while True:
            data = await websocket.receive_text()
            message_data = json.loads(data)
            message = ChatMessage(**message_data)
            
            # For WebSocket, use synchronous processing for immediate response
            response = ai_service.process_chat_message(message)
            await websocket.send_json(response.dict())
    except WebSocketDisconnect:
        pass


@router.post("/structure/export/mol")
async def export_mol(structure: StructureRequest):
    """Export structure as MOL file."""
    if not structure.structure:
        raise HTTPException(status_code=400, detail="Structure is required")
    
    try:
        smiles = molecular_service.structure_to_smiles(structure.structure)
        if not smiles:
            raise HTTPException(status_code=400, detail="Could not convert structure to SMILES")
        
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            raise HTTPException(status_code=400, detail="Invalid SMILES")
        
        # Add hydrogens and generate 3D coordinates
        mol = Chem.AddHs(mol)
        from rdkit.Chem import AllChem
        AllChem.EmbedMolecule(mol, randomSeed=42)
        AllChem.MMFFOptimizeMolecule(mol)
        
        mol_block = Chem.MolToMolBlock(mol)
        
        return Response(
            content=mol_block,
            media_type="chemical/x-mdl-molfile",
            headers={
                "Content-Disposition": f'attachment; filename="molecule.mol"'
            }
        )
    except Exception as e:
        raise HTTPException(status_code=500, detail=f"Export error: {str(e)}")


@router.post("/structure/export/sdf")
async def export_sdf(structure: StructureRequest):
    """Export structure as SDF file."""
    if not structure.structure:
        raise HTTPException(status_code=400, detail="Structure is required")
    
    try:
        smiles = molecular_service.structure_to_smiles(structure.structure)
        if not smiles:
            raise HTTPException(status_code=400, detail="Could not convert structure to SMILES")
        
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            raise HTTPException(status_code=400, detail="Invalid SMILES")
        
        # Add hydrogens and generate 3D coordinates
        mol = Chem.AddHs(mol)
        from rdkit.Chem import AllChem
        AllChem.EmbedMolecule(mol, randomSeed=42)
        AllChem.MMFFOptimizeMolecule(mol)
        
        sdf_block = Chem.MolToMolBlock(mol)
        # SDF format is similar to MOL but with $$$$ separator
        sdf_content = sdf_block + "\n$$$$\n"
        
        return Response(
            content=sdf_content,
            media_type="chemical/x-mdl-sdfile",
            headers={
                "Content-Disposition": f'attachment; filename="molecule.sdf"'
            }
        )
    except Exception as e:
        raise HTTPException(status_code=500, detail=f"Export error: {str(e)}")


@router.post("/structure/export/json")
async def export_json(structure: StructureRequest):
    """Export structure as JSON file."""
    if not structure.structure:
        raise HTTPException(status_code=400, detail="Structure is required")
    
    try:
        structure_dict = structure.structure.dict()
        json_content = json.dumps(structure_dict, indent=2)
        
        return Response(
            content=json_content,
            media_type="application/json",
            headers={
                "Content-Disposition": f'attachment; filename="molecule.json"'
            }
        )
    except Exception as e:
        raise HTTPException(status_code=500, detail=f"Export error: {str(e)}")



@router.post("/reaction/combine")
async def combine_molecules(request: dict):
    """
    Combine two molecules and predict reaction products.
    Expects: { "reactant1": { "smiles": "...", "structure": {...} }, "reactant2": { "smiles": "...", "structure": {...} } }
    """
    try:
        reactant1 = request.get('reactant1', {})
        reactant2 = request.get('reactant2', {})
        
        # Get SMILES - try from smiles field first, then convert from structure
        smiles1 = reactant1.get('smiles')
        smiles2 = reactant2.get('smiles')
        
        # If SMILES not provided, try to convert from structure
        if not smiles1 and reactant1.get('structure'):
            from app.models.molecular import MolecularStructure
            struct1 = MolecularStructure(**reactant1['structure'])
            smiles1 = molecular_service.structure_to_smiles(struct1)
        
        if not smiles2 and reactant2.get('structure'):
            from app.models.molecular import MolecularStructure
            struct2 = MolecularStructure(**reactant2['structure'])
            smiles2 = molecular_service.structure_to_smiles(struct2)
        
        if not smiles1 or not smiles2:
            raise HTTPException(
                status_code=400, 
                detail=f"Both reactants must have SMILES strings. Reactant1: {'missing' if not smiles1 else 'ok'}, Reactant2: {'missing' if not smiles2 else 'ok'}"
            )
        
        # Validate and sanitize SMILES before calling reaction service
        from rdkit import Chem
        
        # Function to sanitize SMILES (remove explicit hydrogens, canonicalize)
        def sanitize_smiles(smiles: str) -> Optional[str]:
            """Sanitize SMILES by removing explicit hydrogens and canonicalizing."""
            try:
                # Try to parse the SMILES as-is first
                mol = Chem.MolFromSmiles(smiles)
                if mol:
                    # Remove explicit hydrogens and get canonical SMILES
                    mol = Chem.RemoveHs(mol)
                    canonical = Chem.MolToSmiles(mol, canonical=True)
                    return canonical
                
                # If parsing fails, try aggressive cleaning
                # This handles malformed SMILES with explicit hydrogens
                cleaned = smiles
                # Remove all [H] patterns
                import re
                cleaned = re.sub(r'\[H\]', '', cleaned)
                cleaned = re.sub(r'\(\[H\]\)', '', cleaned)
                # Fix aromatic ring notation
                cleaned = cleaned.replace('C1=', 'c1=').replace('=C1', '=c1')
                cleaned = cleaned.replace('C1', 'c1').replace('=1', '1')
                # Remove empty parentheses
                cleaned = re.sub(r'\(\)', '', cleaned)
                
                # Try parsing cleaned version
                if cleaned and cleaned != smiles:
                    mol = Chem.MolFromSmiles(cleaned)
                    if mol:
                        mol = Chem.RemoveHs(mol)
                        return Chem.MolToSmiles(mol, canonical=True)
                
                # Last resort: if structure is available, try converting from structure
                return None
            except Exception as e:
                print(f"Error sanitizing SMILES {smiles}: {e}")
                return None
        
        # Sanitize both SMILES
        sanitized_smiles1 = sanitize_smiles(smiles1)
        sanitized_smiles2 = sanitize_smiles(smiles2)
        
        if not sanitized_smiles1:
            raise HTTPException(status_code=400, detail=f"Invalid SMILES for reactant 1: {smiles1}")
        if not sanitized_smiles2:
            raise HTTPException(status_code=400, detail=f"Invalid SMILES for reactant 2: {smiles2}")
        
        # Use sanitized SMILES for reaction
        smiles1 = sanitized_smiles1
        smiles2 = sanitized_smiles2
        
        # Predict reaction
        reaction_result = reaction_service.predict_reaction(smiles1, smiles2)
        
        if not reaction_result.get('success'):
            raise HTTPException(status_code=400, detail=reaction_result.get('error', 'Reaction prediction failed'))
        
        # Simulate collision
        reactant1_struct = None
        reactant2_struct = None
        
        # Try to get structures for collision animation
        try:
            if reactant1.get('structure'):
                from app.models.molecular import MolecularStructure
                try:
                    reactant1_struct = MolecularStructure(**reactant1['structure'])
                except:
                    reactant1_struct = molecular_service.smiles_to_structure(smiles1)
            else:
                reactant1_struct = molecular_service.smiles_to_structure(smiles1)
            
            if reactant2.get('structure'):
                from app.models.molecular import MolecularStructure
                try:
                    reactant2_struct = MolecularStructure(**reactant2['structure'])
                except:
                    reactant2_struct = molecular_service.smiles_to_structure(smiles2)
            else:
                reactant2_struct = molecular_service.smiles_to_structure(smiles2)
        except Exception as e:
            print(f"Warning: Could not create structures for collision animation: {e}")
            # Continue without collision data - reaction still works
        
        
        if reactant1_struct and reactant2_struct:
            collision_data = reaction_service.simulate_collision(reactant1_struct, reactant2_struct)
            reaction_result['collision'] = collision_data
        
        return reaction_result
    except Exception as e:
        raise HTTPException(status_code=500, detail=f"Reaction error: {str(e)}")


@router.get("/structure/search/pubchem")
async def search_pubchem(name: str = None, cid: int = None):
    """
    Search PubChem database by compound name or CID.
    Returns compound data including SMILES.
    """
    try:
        if cid:
            result = pubchem_service.get_compound_by_cid(cid)
        elif name:
            result = pubchem_service.search_by_name(name)
        else:
            raise HTTPException(status_code=400, detail="Either 'name' or 'cid' parameter is required")
        
        if not result:
            raise HTTPException(status_code=404, detail="Compound not found in PubChem")
        
        # Convert SMILES to structure if available
        structure = None
        if result.get('smiles'):
            structure = molecular_service.smiles_to_structure(result['smiles'])
        
        return {
            "source": "pubchem",
            "cid": result.get('cid'),
            "name": result.get('name'),
            "smiles": result.get('smiles'),
            "formula": result.get('formula'),
            "molecular_weight": result.get('molecular_weight'),
            "iupac_name": result.get('iupac_name'),
            "structure": structure.dict() if structure else None
        }
    except HTTPException:
        raise
    except Exception as e:
        raise HTTPException(status_code=500, detail=f"PubChem search error: {str(e)}")
