from typing import Optional
from app.models.molecular import MolecularStructure
from app.services.molecular_service import MolecularService
from app.services.ai_service import AIService

class StructureGenerator:
    def __init__(self):
        self.ai_service = AIService()
        self.molecular_service = MolecularService()
    
    def generate_from_text(self, description: str) -> Optional[MolecularStructure]:
        """Generate structure from natural language using AI."""
        return self.ai_service.generate_structure_from_text(description)
    
    def optimize_structure(self, structure: MolecularStructure) -> Optional[MolecularStructure]:
        """Optimize molecular geometry."""
        smiles = self.molecular_service.structure_to_smiles(structure)
        if smiles:
            # Re-generate with optimization
            return self.molecular_service.smiles_to_structure(smiles)
        return structure
    
    def suggest_modifications(self, structure: MolecularStructure, suggestion: str) -> Optional[MolecularStructure]:
        """Apply AI-suggested modifications to structure."""
        smiles = self.molecular_service.structure_to_smiles(structure)
        if not smiles:
            return None
        
        prompt = f"Modify this molecule (SMILES: {smiles}) according to: {suggestion}. Provide the new SMILES."
        new_structure = self.ai_service.generate_structure_from_text(prompt)
        return new_structure

