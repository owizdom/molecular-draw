from rdkit import Chem
from rdkit.Chem import AllChem, Descriptors
from typing import Optional, List, Tuple
from app.models.molecular import MolecularStructure, Atom, Bond

class MolecularService:
    @staticmethod
    def smiles_to_structure(smiles: str) -> Optional[MolecularStructure]:
        """Convert SMILES string to molecular structure with 3D coordinates."""
        try:
            mol = Chem.MolFromSmiles(smiles)
            if mol is None:
                return None
            
            # Generate 3D coordinates
            mol = Chem.AddHs(mol)
            AllChem.EmbedMolecule(mol, randomSeed=42)
            AllChem.MMFFOptimizeMolecule(mol)
            
            # Extract atoms
            atoms = []
            conf = mol.GetConformer()
            for i, atom in enumerate(mol.GetAtoms()):
                pos = conf.GetAtomPosition(i)
                atoms.append(Atom(
                    element=atom.GetSymbol(),
                    x=float(pos.x),
                    y=float(pos.y),
                    z=float(pos.z),
                    id=i
                ))
            
            # Extract bonds
            bonds = []
            for bond in mol.GetBonds():
                bond_type_float = bond.GetBondTypeAsDouble()
                # Convert to int (1=single, 2=double, 3=triple)
                # Aromatic bonds (1.5) are treated as single bonds
                bond_type_int = int(round(bond_type_float))
                bonds.append(Bond(
                    atom1_id=bond.GetBeginAtomIdx(),
                    atom2_id=bond.GetEndAtomIdx(),
                    bond_type=bond_type_int
                ))
            
            return MolecularStructure(
                atoms=atoms,
                bonds=bonds,
                smiles=smiles
            )
        except Exception as e:
            print(f"Error converting SMILES: {e}")
            return None
    
    @staticmethod
    def structure_to_smiles(structure: MolecularStructure) -> Optional[str]:
        """Convert molecular structure to SMILES string."""
        try:
            # If structure has SMILES already stored, use it (most reliable)
            if hasattr(structure, 'smiles') and structure.smiles:
                # Validate and return canonical version
                mol = Chem.MolFromSmiles(structure.smiles)
                if mol:
                    mol = Chem.RemoveHs(mol)
                    return Chem.MolToSmiles(mol, canonical=True)
            
            # Otherwise, reconstruct from atoms and bonds
            mol = Chem.RWMol()
            
            # Add atoms - skip hydrogen atoms (they'll be added implicitly)
            atom_map = {}
            for idx, atom in enumerate(structure.atoms):
                # Skip explicit hydrogens - RDKit will add them automatically
                if atom.element.upper() == 'H':
                    continue
                rdkit_atom = Chem.Atom(atom.element)
                rdkit_idx = mol.AddAtom(rdkit_atom)
                # Map by ID if available, otherwise by index
                atom_key = atom.id if atom.id is not None else idx
                atom_map[atom_key] = rdkit_idx
            
            # Add bonds (only between non-hydrogen atoms)
            for bond in structure.bonds:
                atom1_idx = atom_map.get(bond.atom1_id)
                atom2_idx = atom_map.get(bond.atom2_id)
                # Only add bond if both atoms are in the map (i.e., not hydrogens)
                if atom1_idx is not None and atom2_idx is not None:
                    try:
                        mol.AddBond(
                            atom1_idx,
                            atom2_idx,
                            Chem.BondType(int(bond.bond_type))
                        )
                    except:
                        # If bond type is invalid, default to single bond
                        mol.AddBond(atom1_idx, atom2_idx, Chem.BondType.SINGLE)
            
            # Sanitize the molecule
            try:
                Chem.SanitizeMol(mol)
            except Exception as e:
                print(f"Sanitization warning: {e}")
                # Try to fix common issues
                try:
                    # Set aromaticity
                    Chem.SetAromaticity(mol)
                except:
                    pass
            
            # Get canonical SMILES (without explicit hydrogens)
            mol = Chem.RemoveHs(mol)
            smiles = Chem.MolToSmiles(mol, canonical=True)
            return smiles
        except Exception as e:
            print(f"Error converting structure to SMILES: {e}")
            return None
    
    @staticmethod
    def calculate_properties(smiles: str) -> dict:
        """Calculate molecular properties from SMILES."""
        try:
            mol = Chem.MolFromSmiles(smiles)
            if mol is None:
                return {}
            
            return {
                "molecular_weight": Descriptors.MolWt(mol),
                "formula": Chem.rdMolDescriptors.CalcMolFormula(mol),
                "num_atoms": mol.GetNumAtoms(),
                "num_bonds": mol.GetNumBonds(),
                "num_rings": Chem.rdMolDescriptors.CalcNumRings(mol),
            }
        except Exception as e:
            print(f"Error calculating properties: {e}")
            return {}
    
    @staticmethod
    def validate_structure(structure: MolecularStructure) -> Tuple[bool, Optional[str]]:
        """Validate molecular structure."""
        try:
            smiles = MolecularService.structure_to_smiles(structure)
            if smiles is None:
                return False, "Invalid structure: Could not convert to SMILES"
            
            mol = Chem.MolFromSmiles(smiles)
            if mol is None:
                return False, "Invalid structure: SMILES is not valid"
            
            return True, None
        except Exception as e:
            return False, f"Validation error: {str(e)}"

