from app.celery_app import celery_app
from app.services.molecular_service import MolecularService
from app.services.ai_service import AIService
from app.services.structure_generator import StructureGenerator
from app.models.molecular import MolecularStructure, ChatMessage, ChatResponse
from typing import Optional, Dict, Any
import redis
import json
import os
from dotenv import load_dotenv

load_dotenv()

molecular_service = MolecularService()
ai_service = AIService()
structure_generator = StructureGenerator()

def _publish_progress(task_id: str, status: str, message: str = "", current: int = 0, total: int = 1, data: Optional[Dict] = None):
    """Publish progress update to Redis pub/sub."""
    try:
        redis_url = os.getenv("CELERY_BROKER_URL", "redis://localhost:6379/0")
        redis_client = redis.from_url(redis_url, decode_responses=True)
        channel = f"task_progress:{task_id}"
        progress_data = {
            "task_id": task_id,
            "status": status,
            "message": message,
            "current": current,
            "total": total,
        }
        if data:
            progress_data.update(data)
        redis_client.publish(channel, json.dumps(progress_data))
        redis_client.close()
    except Exception as e:
        print(f"Error publishing progress: {e}")


@celery_app.task(name="tasks.smiles_to_structure", bind=True)
def task_smiles_to_structure(self, smiles: str) -> Dict[str, Any]:
    """Convert SMILES to 3D structure asynchronously."""
    try:
        self.update_state(state="PROGRESS", meta={"current": 0, "total": 3, "message": "Parsing SMILES string..."})
        _publish_progress(self.request.id, "processing", "Parsing SMILES string...", 0, 3)
        
        structure = molecular_service.smiles_to_structure(smiles)
        if not structure:
            _publish_progress(self.request.id, "failed", "Invalid SMILES string", 0, 3)
            return {"error": "Invalid SMILES string", "structure": None}
        
        self.update_state(state="PROGRESS", meta={"current": 1, "total": 3, "message": "Generating 3D coordinates..."})
        _publish_progress(self.request.id, "processing", "Generating 3D coordinates...", 1, 3)
        
        self.update_state(state="PROGRESS", meta={"current": 2, "total": 3, "message": "Calculating properties..."})
        _publish_progress(self.request.id, "processing", "Calculating properties...", 2, 3)
        
        properties = molecular_service.calculate_properties(smiles)
        
        result = {
            "structure": structure.dict() if structure else None,
            "molecular_weight": properties.get("molecular_weight"),
            "formula": properties.get("formula"),
            "properties": properties,
        }
        
        _publish_progress(self.request.id, "completed", "Structure generated successfully", 3, 3, result)
        return result
    except Exception as e:
        _publish_progress(self.request.id, "failed", f"Error: {str(e)}", 0, 3)
        return {"error": str(e), "structure": None}


@celery_app.task(name="tasks.process_chat", bind=True)
def task_process_chat(self, message: str, structure: Optional[Dict[str, Any]] = None) -> Dict[str, Any]:
    """Process chat message with AI asynchronously."""
    try:
        self.update_state(state="PROGRESS", meta={"current": 0, "total": 2, "message": "Processing your message..."})
        _publish_progress(self.request.id, "processing", "Processing your message...", 0, 2)
        
        mol_structure = None
        if structure:
            mol_structure = MolecularStructure(**structure)
        
        chat_message = ChatMessage(message=message, structure=mol_structure)
        
        self.update_state(state="PROGRESS", meta={"current": 1, "total": 2, "message": "Generating AI response..."})
        _publish_progress(self.request.id, "processing", "Generating AI response...", 1, 2)
        
        response = ai_service.process_chat_message(chat_message)
        
        result = {
            "response": response.response,
            "structure": response.structure.dict() if response.structure else None,
            "suggestions": response.suggestions,
        }
        
        _publish_progress(self.request.id, "completed", "Response generated", 2, 2, result)
        return result
    except Exception as e:
        _publish_progress(self.request.id, "failed", f"Error: {str(e)}", 0, 2)
        return {"error": str(e), "response": None, "structure": None}


@celery_app.task(name="tasks.generate_structure", bind=True)
def task_generate_structure(self, description: str) -> Dict[str, Any]:
    """Generate structure from natural language description asynchronously."""
    try:
        self.update_state(state="PROGRESS", meta={"current": 0, "total": 4, "message": "Analyzing description..."})
        _publish_progress(self.request.id, "processing", "Analyzing description...", 0, 4)
        
        self.update_state(state="PROGRESS", meta={"current": 1, "total": 4, "message": "Searching web for molecule..."})
        _publish_progress(self.request.id, "processing", "Searching web for molecule...", 1, 4)
        
        self.update_state(state="PROGRESS", meta={"current": 2, "total": 4, "message": "Generating SMILES notation..."})
        _publish_progress(self.request.id, "processing", "Generating SMILES notation...", 2, 4)
        
        structure = structure_generator.generate_from_text(description)
        if not structure:
            _publish_progress(self.request.id, "failed", "Could not generate structure from description", 0, 4)
            return {"error": "Could not generate structure from description", "structure": None}
        
        self.update_state(state="PROGRESS", meta={"current": 3, "total": 4, "message": "Creating 3D structure..."})
        _publish_progress(self.request.id, "processing", "Creating 3D structure...", 3, 4)
        
        result = {"structure": structure.dict()}
        _publish_progress(self.request.id, "completed", "Structure generated successfully", 4, 4, result)
        return result
    except Exception as e:
        _publish_progress(self.request.id, "failed", f"Error: {str(e)}", 0, 3)
        return {"error": str(e), "structure": None}


@celery_app.task(name="tasks.optimize_structure", bind=True)
def task_optimize_structure(self, structure: Dict[str, Any]) -> Dict[str, Any]:
    """Optimize molecular structure asynchronously."""
    try:
        self.update_state(state="PROGRESS", meta={"current": 0, "total": 2, "message": "Preparing structure..."})
        _publish_progress(self.request.id, "processing", "Preparing structure...", 0, 2)
        
        mol_structure = MolecularStructure(**structure)
        
        self.update_state(state="PROGRESS", meta={"current": 1, "total": 2, "message": "Optimizing geometry..."})
        _publish_progress(self.request.id, "processing", "Optimizing geometry...", 1, 2)
        
        optimized = structure_generator.optimize_structure(mol_structure)
        
        if not optimized:
            _publish_progress(self.request.id, "failed", "Optimization failed", 0, 2)
            return {"error": "Optimization failed", "structure": None}
        
        result = {"structure": optimized.dict()}
        _publish_progress(self.request.id, "completed", "Structure optimized successfully", 2, 2, result)
        return result
    except Exception as e:
        _publish_progress(self.request.id, "failed", f"Error: {str(e)}", 0, 2)
        return {"error": str(e), "structure": None}


@celery_app.task(name="tasks.calculate_properties", bind=True)
def task_calculate_properties(self, smiles: Optional[str] = None, structure: Optional[Dict[str, Any]] = None) -> Dict[str, Any]:
    """Calculate molecular properties asynchronously."""
    try:
        self.update_state(state="PROGRESS", meta={"current": 0, "total": 2, "message": "Preparing calculation..."})
        _publish_progress(self.request.id, "processing", "Preparing calculation...", 0, 2)
        
        if smiles:
            self.update_state(state="PROGRESS", meta={"current": 1, "total": 2, "message": "Calculating properties..."})
            _publish_progress(self.request.id, "processing", "Calculating properties...", 1, 2)
            properties = molecular_service.calculate_properties(smiles)
        elif structure:
            mol_structure = MolecularStructure(**structure)
            smiles = molecular_service.structure_to_smiles(mol_structure)
            if smiles:
                self.update_state(state="PROGRESS", meta={"current": 1, "total": 2, "message": "Calculating properties..."})
                _publish_progress(self.request.id, "processing", "Calculating properties...", 1, 2)
                properties = molecular_service.calculate_properties(smiles)
            else:
                _publish_progress(self.request.id, "failed", "Could not convert structure to SMILES", 0, 2)
                return {"error": "Could not convert structure to SMILES"}
        else:
            _publish_progress(self.request.id, "failed", "SMILES or structure is required", 0, 2)
            return {"error": "SMILES or structure is required"}
        
        _publish_progress(self.request.id, "completed", "Properties calculated", 2, 2, properties)
        return properties
    except Exception as e:
        _publish_progress(self.request.id, "failed", f"Error: {str(e)}", 0, 2)
        return {"error": str(e)}


@celery_app.task(name="tasks.validate_structure", bind=True)
def task_validate_structure(self, structure: Dict[str, Any]) -> Dict[str, Any]:
    """Validate molecular structure asynchronously."""
    try:
        self.update_state(state="PROGRESS", meta={"current": 0, "total": 2, "message": "Validating structure..."})
        _publish_progress(self.request.id, "processing", "Validating structure...", 0, 2)
        
        mol_structure = MolecularStructure(**structure)
        is_valid, error = molecular_service.validate_structure(mol_structure)
        
        result = {"is_valid": is_valid, "error": error}
        status = "completed" if is_valid else "failed"
        message = "Structure is valid" if is_valid else f"Validation failed: {error}"
        _publish_progress(self.request.id, status, message, 2, 2, result)
        return result
    except Exception as e:
        _publish_progress(self.request.id, "failed", f"Error: {str(e)}", 0, 2)
        return {"is_valid": False, "error": str(e)}

