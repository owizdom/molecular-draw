from pydantic import BaseModel
from typing import List, Optional

class Atom(BaseModel):
    element: str
    x: float
    y: float
    z: float
    id: Optional[int] = None

class Bond(BaseModel):
    atom1_id: int
    atom2_id: int
    bond_type: int  # 1=single, 2=double, 3=triple

class MolecularStructure(BaseModel):
    atoms: List[Atom]
    bonds: List[Bond]
    smiles: Optional[str] = None

class StructureRequest(BaseModel):
    smiles: Optional[str] = None
    structure: Optional[MolecularStructure] = None

class StructureResponse(BaseModel):
    structure: MolecularStructure
    molecular_weight: Optional[float] = None
    formula: Optional[str] = None
    properties: Optional[dict] = None

class ChatMessage(BaseModel):
    message: str
    structure: Optional[MolecularStructure] = None

class ChatResponse(BaseModel):
    response: str
    structure: Optional[MolecularStructure] = None
    suggestions: Optional[List[str]] = None

