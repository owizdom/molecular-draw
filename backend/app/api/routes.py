from fastapi import APIRouter, WebSocket, WebSocketDisconnect, HTTPException
from fastapi.responses import Response, StreamingResponse
from app.models.molecular import (
    StructureRequest, StructureResponse, ChatMessage, ChatResponse
)
from app.services.molecular_service import MolecularService
from app.services.ai_service import AIService
from app.services.structure_generator import StructureGenerator
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
    """Generate structure from natural language description. Uses async task only for AI calls."""
    desc = description.get("description", "")
    if not desc:
        raise HTTPException(status_code=400, detail="Description is required")
    
    # Try fast fallback first (common molecules)
    desc_lower = desc.lower().strip()
    common_molecules = ['water', 'h2o', 'methane', 'ch4', 'ethane', 'ethanol', 'benzene', 'aspirin', 'caffeine']
    
    if any(mol in desc_lower for mol in common_molecules):
        # Try instant generation
        structure = structure_generator.generate_from_text(desc)
        if structure:
            return {"structure": structure}
    
    # For complex requests, use async task
    if async_task:
        task = task_generate_structure.delay(desc)
        return {"task_id": task.id, "status": "processing", "message": "Structure generation task submitted"}
    
    # Synchronous processing (fallback)
    structure = structure_generator.generate_from_text(desc)
    if not structure:
        raise HTTPException(status_code=500, detail="Could not generate structure from description")
    
    return {"structure": structure}

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

