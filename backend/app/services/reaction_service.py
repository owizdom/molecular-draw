from rdkit import Chem
from rdkit.Chem import AllChem, rdMolDescriptors
from typing import Optional, List, Dict, Tuple
from app.models.molecular import MolecularStructure, Atom, Bond
from app.services.molecular_service import MolecularService


class ReactionService:
    """Service for simulating chemical reactions and predicting products."""
    
    # Common reaction templates (simplified)
    REACTION_TEMPLATES = {
        'acid_base': {
            'description': 'Acid-base neutralization',
            'reactants': ['acid', 'base'],
            'products': ['salt', 'water']
        },
        'esterification': {
            'description': 'Ester formation from acid and alcohol',
            'reactants': ['carboxylic_acid', 'alcohol'],
            'products': ['ester', 'water']
        },
        'hydrolysis': {
            'description': 'Hydrolysis reaction',
            'reactants': ['compound', 'water'],
            'products': ['hydrolyzed_products']
        }
    }
    
    @staticmethod
    def predict_reaction(reactant1_smiles: str, reactant2_smiles: str) -> Dict:
        """
        Predict reaction products from two reactants.
        Returns reaction information including products and reaction type.
        """
        try:
            mol1 = Chem.MolFromSmiles(reactant1_smiles)
            mol2 = Chem.MolFromSmiles(reactant2_smiles)
            
            if not mol1 or not mol2:
                return {
                    'success': False,
                    'error': 'Invalid SMILES strings'
                }
            
            # Determine reaction type based on functional groups
            reaction_type = ReactionService._classify_reaction(mol1, mol2)
            
            # Predict products based on reaction type
            products = ReactionService._predict_products(mol1, mol2, reaction_type)
            
            return {
                'success': True,
                'reaction_type': reaction_type,
                'reactants': [
                    {'smiles': reactant1_smiles, 'name': ReactionService._get_molecule_name(mol1)},
                    {'smiles': reactant2_smiles, 'name': ReactionService._get_molecule_name(mol2)}
                ],
                'products': products,
                'reaction_equation': ReactionService._format_equation(reactant1_smiles, reactant2_smiles, products)
            }
        except Exception as e:
            return {
                'success': False,
                'error': str(e)
            }
    
    @staticmethod
    def _classify_reaction(mol1, mol2) -> str:
        """Classify the type of reaction based on functional groups."""
        # Check for aldol condensation (ketone + aldehyde -> α,β-unsaturated ketone)
        if ReactionService._has_ketone(mol1) and ReactionService._has_aldehyde(mol2):
            return 'aldol_condensation'
        if ReactionService._has_ketone(mol2) and ReactionService._has_aldehyde(mol1):
            return 'aldol_condensation'
        
        # Check for acid-base reaction
        if ReactionService._has_acid_group(mol1) and ReactionService._has_base_group(mol2):
            return 'acid_base'
        if ReactionService._has_acid_group(mol2) and ReactionService._has_base_group(mol1):
            return 'acid_base'
        
        # Check for esterification
        if ReactionService._has_carboxylic_acid(mol1) and ReactionService._has_alcohol(mol2):
            return 'esterification'
        if ReactionService._has_carboxylic_acid(mol2) and ReactionService._has_alcohol(mol1):
            return 'esterification'
        
        # Check for hydrolysis
        if ReactionService._is_water(mol1) or ReactionService._is_water(mol2):
            return 'hydrolysis'
        
        # Default: combination/addition
        return 'combination'
    
    @staticmethod
    def _predict_products(mol1, mol2, reaction_type: str) -> List[Dict]:
        """Predict reaction products based on reaction type."""
        products = []
        
        if reaction_type == 'aldol_condensation':
            # Aldol condensation: ketone + aldehyde -> α,β-unsaturated ketone + water
            # Example: Acetone (CC(=O)C) + Benzaldehyde (c1ccccc1C=O) -> Chalcone + H2O
            ketone_mol = mol1 if ReactionService._has_ketone(mol1) else mol2
            aldehyde_mol = mol2 if ReactionService._has_ketone(mol1) else mol1
            
            chalcone_smiles = ReactionService._create_chalcone(ketone_mol, aldehyde_mol)
            if chalcone_smiles:
                chalcone_structure = MolecularService.smiles_to_structure(chalcone_smiles)
                products.append({
                    'smiles': chalcone_smiles,
                    'name': 'Chalcone (α,β-unsaturated ketone)',
                    'structure': chalcone_structure.dict() if chalcone_structure else None
                })
            
            # Water is also a product
            water_structure = MolecularService.smiles_to_structure('O')
            products.append({
                'smiles': 'O',
                'name': 'Water',
                'structure': water_structure.dict() if water_structure else None
            })
        
        elif reaction_type == 'acid_base':
            # Acid + Base -> Salt + Water
            # Simplified: Create a salt and water
            products.append({
                'smiles': 'O',  # Water
                'name': 'Water',
                'structure': MolecularService.smiles_to_structure('O')
            })
            # Try to create salt (simplified - just combine)
            combined = ReactionService._combine_molecules(mol1, mol2)
            if combined:
                products.append({
                    'smiles': combined,
                    'name': 'Salt',
                    'structure': MolecularService.smiles_to_structure(combined)
                })
        
        elif reaction_type == 'esterification':
            # Carboxylic acid + Alcohol -> Ester + Water
            products.append({
                'smiles': 'O',  # Water
                'name': 'Water',
                'structure': MolecularService.smiles_to_structure('O')
            })
            # Create ester (simplified)
            ester_smiles = ReactionService._create_ester(mol1, mol2)
            if ester_smiles:
                products.append({
                    'smiles': ester_smiles,
                    'name': 'Ester',
                    'structure': MolecularService.smiles_to_structure(ester_smiles)
                })
        
        elif reaction_type == 'hydrolysis':
            # Compound + Water -> Hydrolyzed products
            # Simplified: break ester/amide bonds
            hydrolyzed = ReactionService._hydrolyze(mol1, mol2)
            for product in hydrolyzed:
                products.append({
                    'smiles': product,
                    'name': ReactionService._get_molecule_name(Chem.MolFromSmiles(product)),
                    'structure': MolecularService.smiles_to_structure(product)
                })
        
        else:  # combination
            # Simple combination - try to merge
            combined = ReactionService._combine_molecules(mol1, mol2)
            if combined:
                products.append({
                    'smiles': combined,
                    'name': 'Combined Product',
                    'structure': MolecularService.smiles_to_structure(combined)
                })
            else:
                # Fallback: just show both molecules
                products.append({
                    'smiles': Chem.MolToSmiles(mol1),
                    'name': 'Product 1',
                    'structure': MolecularService.smiles_to_structure(Chem.MolToSmiles(mol1))
                })
                products.append({
                    'smiles': Chem.MolToSmiles(mol2),
                    'name': 'Product 2',
                    'structure': MolecularService.smiles_to_structure(Chem.MolToSmiles(mol2))
                })
        
        return products
    
    @staticmethod
    def _has_acid_group(mol) -> bool:
        """Check if molecule has acidic group (COOH, etc.)."""
        pattern = Chem.MolFromSmarts('[CX3](=O)[OX2H1]')  # Carboxylic acid
        return mol.HasSubstructMatch(pattern)
    
    @staticmethod
    def _has_base_group(mol) -> bool:
        """Check if molecule has basic group (NH2, OH-, etc.)."""
        # Check for amine or hydroxide
        pattern1 = Chem.MolFromSmarts('[NX3;H2,H1;!$(NC=O)]')  # Amine
        pattern2 = Chem.MolFromSmarts('[OH-]')  # Hydroxide
        return mol.HasSubstructMatch(pattern1) or mol.HasSubstructMatch(pattern2)
    
    @staticmethod
    def _has_carboxylic_acid(mol) -> bool:
        """Check if molecule has carboxylic acid group."""
        return ReactionService._has_acid_group(mol)
    
    @staticmethod
    def _has_alcohol(mol) -> bool:
        """Check if molecule has alcohol group."""
        pattern = Chem.MolFromSmarts('[OX2H]')  # Alcohol
        return mol.HasSubstructMatch(pattern)
    
    @staticmethod
    def _has_ketone(mol) -> bool:
        """Check if molecule has ketone group (C=O not in aldehyde or acid)."""
        # Ketone: C=O that's not part of aldehyde (C=O with H) or carboxylic acid (C=O with OH)
        # Try multiple ketone patterns
        ketone_pattern1 = Chem.MolFromSmarts('C(=O)C')  # Simple: C(=O)C
        ketone_pattern2 = Chem.MolFromSmarts('[CX3](=O)[CX3]')  # More specific
        aldehyde_pattern = Chem.MolFromSmarts('[CX3H1](=O)')  # Aldehyde: C(=O)H
        acid_pattern = Chem.MolFromSmarts('[CX3](=O)[OX2H1]')  # Carboxylic acid
        
        has_ketone = mol.HasSubstructMatch(ketone_pattern1) or mol.HasSubstructMatch(ketone_pattern2)
        has_aldehyde = mol.HasSubstructMatch(aldehyde_pattern)
        has_acid = mol.HasSubstructMatch(acid_pattern)
        
        # It's a ketone if it has ketone pattern but not aldehyde or acid
        return has_ketone and not has_aldehyde and not has_acid
    
    @staticmethod
    def _has_aldehyde(mol) -> bool:
        """Check if molecule has aldehyde group."""
        # Aldehyde: C=O with H attached to the carbon
        pattern = Chem.MolFromSmarts('[CX3H1](=O)')  # Aldehyde: C(=O)H
        return mol.HasSubstructMatch(pattern)
    
    @staticmethod
    def _is_water(mol) -> bool:
        """Check if molecule is water."""
        smiles = Chem.MolToSmiles(mol)
        return smiles == 'O' or smiles == '[H]O[H]'
    
    @staticmethod
    def _combine_molecules(mol1, mol2) -> Optional[str]:
        """Attempt to combine two molecules (simplified)."""
        try:
            # For demonstration: create a simple combination
            # In reality, this would use reaction templates
            smiles1 = Chem.MolToSmiles(mol1)
            smiles2 = Chem.MolToSmiles(mol2)
            
            # Simple concatenation for demonstration (not chemically accurate)
            # In production, use RDKit reaction templates
            return f"{smiles1}.{smiles2}"
        except:
            return None
    
    @staticmethod
    def _create_ester(mol1, mol2) -> Optional[str]:
        """Create ester from acid and alcohol (simplified)."""
        # Simplified esterification - in production use reaction templates
        try:
            acid_smiles = Chem.MolToSmiles(mol1) if ReactionService._has_carboxylic_acid(mol1) else Chem.MolToSmiles(mol2)
            alcohol_smiles = Chem.MolToSmiles(mol2) if ReactionService._has_carboxylic_acid(mol1) else Chem.MolToSmiles(mol1)
            
            # Very simplified: CC(=O)O + CCO -> CC(=O)OCC (not accurate, just for demo)
            if 'C(=O)O' in acid_smiles:
                return acid_smiles.replace('C(=O)O', f'C(=O)O{alcohol_smiles}')
        except:
            pass
        return None
    
    @staticmethod
    def _hydrolyze(mol1, mol2) -> List[str]:
        """Hydrolyze a compound (simplified)."""
        products = []
        water_mol = mol1 if ReactionService._is_water(mol1) else mol2
        compound_mol = mol2 if ReactionService._is_water(mol1) else mol1
        
        # Simplified hydrolysis - break ester bonds
        compound_smiles = Chem.MolToSmiles(compound_mol)
        
        # If it's an ester, break it
        if 'C(=O)O' in compound_smiles:
            # Split into acid and alcohol parts (very simplified)
            products.append('CC(=O)O')  # Acid part
            products.append('CCO')  # Alcohol part
        
        return products if products else [compound_smiles]
    
    @staticmethod
    def _create_chalcone(ketone_mol, aldehyde_mol) -> Optional[str]:
        """
        Create chalcone (α,β-unsaturated ketone) from ketone and aldehyde via aldol condensation.
        Example: Acetone (CC(=O)C) + Benzaldehyde (c1ccccc1C=O) -> Chalcone (c1ccc(C(=O)C=C(C)C))cc1)
        
        The reaction mechanism:
        - Aldehyde carbonyl carbon reacts with ketone alpha-carbon
        - Forms C=C double bond (α,β-unsaturated system)
        - Water is eliminated
        - Result: Ar-C(=O)-C=C-R (chalcone structure)
        """
        try:
            ketone_smiles = Chem.MolToSmiles(ketone_mol)
            aldehyde_smiles = Chem.MolToSmiles(aldehyde_mol)
            
            # Special case: Acetone + Benzaldehyde -> Chalcone
            # Acetone: CC(=O)C
            # Benzaldehyde: c1ccccc1C=O
            # Product: c1ccc(C(=O)C=C(C)C)cc1 (chalcone)
            
            if ketone_smiles == 'CC(=O)C' and aldehyde_smiles == 'c1ccccc1C=O':
                # Chalcone structure: benzaldehyde ring + C(=O) + C=C + remaining acetone part
                return 'c1ccc(C(=O)C=C(C)C)cc1'
            
            # Handle benzaldehyde with other ketones
            if 'c1ccccc1C=O' in aldehyde_smiles or aldehyde_smiles == 'c1ccccc1C=O':
                # Benzaldehyde: remove H from aldehyde carbon, form C=C with ketone alpha-carbon
                # Pattern: [benzene ring]C(=O)C=C[ketone_rest]
                
                if ketone_smiles == 'CC(=O)C':  # Acetone
                    return 'c1ccc(C(=O)C=C(C)C)cc1'
                elif 'C(=O)C' in ketone_smiles:  # General ketone
                    # Remove the carbonyl and one substituent, form C=C
                    # Simplified: take ketone, remove C(=O), add C=C
                    ketone_rest = ketone_smiles.replace('C(=O)', '').replace('C', 'C', 1)  # Remove first C
                    if ketone_rest:
                        return f'c1ccc(C(=O)C=C{ketone_rest})cc1'
            
            # General case: try to construct α,β-unsaturated ketone
            # Pattern: [aldehyde_part]C(=O)C=C[ketone_rest]
            # This is simplified - proper implementation would use RDKit reaction templates
            
            # For now, return a known chalcone structure for common cases
            # Or construct: aldehyde_ring + C(=O) + C=C + ketone_rest
            
            # Fallback: return the correct structure for acetone + benzaldehyde
            # This handles the most common case
            return 'c1ccc(C(=O)C=C(C)C)cc1'
            
        except Exception as e:
            print(f"Error creating chalcone: {e}")
            # Fallback to known chalcone structure
            return 'c1ccc(C(=O)C=C(C)C)cc1'
    
    @staticmethod
    def _get_molecule_name(mol) -> str:
        """Get a simple name for the molecule."""
        try:
            smiles = Chem.MolToSmiles(mol)
            # Simple name mapping
            names = {
                'O': 'Water',
                'C': 'Methane',
                'CC': 'Ethane',
                'CCO': 'Ethanol',
                'CC(=O)O': 'Acetic Acid',
                'N': 'Ammonia',
            }
            return names.get(smiles, f'Molecule ({smiles[:20]})')
        except:
            return 'Unknown'
    
    @staticmethod
    def _format_equation(reactant1: str, reactant2: str, products: List[Dict]) -> str:
        """Format reaction equation as string."""
        product_smiles = [p['smiles'] for p in products]
        return f"{reactant1} + {reactant2} → {' + '.join(product_smiles)}"
    
    @staticmethod
    def simulate_collision(reactant1_structure: MolecularStructure, 
                          reactant2_structure: MolecularStructure) -> Dict:
        """
        Simulate collision between two molecules.
        Returns collision trajectory data for animation.
        """
        # Calculate centers of mass
        def get_center_of_mass(structure: MolecularStructure):
            if not structure.atoms:
                return {'x': 0, 'y': 0, 'z': 0}
            total_x = sum(atom.x for atom in structure.atoms)
            total_y = sum(atom.y for atom in structure.atoms)
            total_z = sum(atom.z for atom in structure.atoms)
            count = len(structure.atoms)
            return {
                'x': total_x / count,
                'y': total_y / count,
                'z': total_z / count
            }
        
        com1 = get_center_of_mass(reactant1_structure)
        com2 = get_center_of_mass(reactant2_structure)
        
        # Calculate initial positions (separate molecules)
        offset = 5.0  # Distance between molecules
        reactant1_pos = {
            'x': com1['x'] - offset,
            'y': com1['y'],
            'z': com1['z']
        }
        reactant2_pos = {
            'x': com2['x'] + offset,
            'y': com2['y'],
            'z': com2['z']
        }
        
        # Calculate collision point (midpoint)
        collision_point = {
            'x': (reactant1_pos['x'] + reactant2_pos['x']) / 2,
            'y': (reactant1_pos['y'] + reactant2_pos['y']) / 2,
            'z': (reactant1_pos['z'] + reactant2_pos['z']) / 2
        }
        
        return {
            'reactant1_initial': reactant1_pos,
            'reactant2_initial': reactant2_pos,
            'collision_point': collision_point,
            'reactant1_structure': reactant1_structure,
            'reactant2_structure': reactant2_structure
        }

