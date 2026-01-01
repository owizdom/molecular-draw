import requests
from typing import Optional, Dict, List
import time

class PubChemService:
    """Service for searching PubChem database."""
    
    BASE_URL = "https://pubchem.ncbi.nlm.nih.gov/rest/pug"
    REQUEST_DELAY = 0.1  # Delay between requests to respect rate limits
    
    def __init__(self):
        self.last_request_time = 0
    
    def _rate_limit(self):
        """Respect PubChem rate limits."""
        current_time = time.time()
        time_since_last = current_time - self.last_request_time
        if time_since_last < self.REQUEST_DELAY:
            time.sleep(self.REQUEST_DELAY - time_since_last)
        self.last_request_time = time.time()
    
    def search_by_name(self, name: str) -> Optional[Dict]:
        """
        Search PubChem by compound name.
        Returns compound data including CID and SMILES.
        """
        try:
            self._rate_limit()
            
            # Clean the name
            clean_name = name.strip()
            
            # Try to get CID first
            cid_url = f"{self.BASE_URL}/compound/name/{clean_name}/cids/JSON"
            response = requests.get(cid_url, timeout=5)
            
            if response.status_code != 200:
                return None
            
            cid_data = response.json()
            if 'IdentifierList' not in cid_data or not cid_data['IdentifierList']['CID']:
                return None
            
            cids = cid_data['IdentifierList']['CID']
            if not cids:
                return None
            
            # Use the first CID
            cid = cids[0]
            
            # Get compound properties - try CanonicalSMILES first, fallback to ConnectivitySMILES
            self._rate_limit()
            props_url = f"{self.BASE_URL}/compound/cid/{cid}/property/CanonicalSMILES,ConnectivitySMILES,MolecularFormula,MolecularWeight,IUPACName/JSON"
            props_response = requests.get(props_url, timeout=5)
            
            if props_response.status_code != 200:
                return None
            
            props_data = props_response.json()
            if 'PropertyTable' not in props_data or not props_data['PropertyTable']['Properties']:
                return None
            
            properties = props_data['PropertyTable']['Properties'][0]
            
            # Use CanonicalSMILES if available, otherwise fallback to ConnectivitySMILES
            smiles = properties.get('CanonicalSMILES') or properties.get('ConnectivitySMILES')
            
            return {
                'cid': cid,
                'smiles': smiles,
                'formula': properties.get('MolecularFormula'),
                'molecular_weight': properties.get('MolecularWeight'),
                'iupac_name': properties.get('IUPACName'),
                'name': clean_name
            }
        except requests.exceptions.Timeout:
            print(f"PubChem timeout for: {name}")
            return None
        except requests.exceptions.RequestException as e:
            print(f"PubChem request error for {name}: {e}")
            return None
        except Exception as e:
            print(f"PubChem error for {name}: {e}")
            return None
    
    def get_compound_by_cid(self, cid: int) -> Optional[Dict]:
        """Get compound details by PubChem CID."""
        try:
            self._rate_limit()
            
            props_url = f"{self.BASE_URL}/compound/cid/{cid}/property/CanonicalSMILES,ConnectivitySMILES,MolecularFormula,MolecularWeight,IUPACName/JSON"
            response = requests.get(props_url, timeout=5)
            
            if response.status_code != 200:
                return None
            
            props_data = response.json()
            if 'PropertyTable' not in props_data or not props_data['PropertyTable']['Properties']:
                return None
            
            properties = props_data['PropertyTable']['Properties'][0]
            
            # Use CanonicalSMILES if available, otherwise fallback to ConnectivitySMILES
            smiles = properties.get('CanonicalSMILES') or properties.get('ConnectivitySMILES')
            
            return {
                'cid': cid,
                'smiles': smiles,
                'formula': properties.get('MolecularFormula'),
                'molecular_weight': properties.get('MolecularWeight'),
                'iupac_name': properties.get('IUPACName')
            }
        except Exception as e:
            print(f"PubChem CID lookup error for {cid}: {e}")
            return None
    
    def get_smiles_from_pubchem(self, name: str) -> Optional[str]:
        """
        Extract SMILES notation from PubChem by compound name.
        Returns SMILES string or None if not found.
        """
        result = self.search_by_name(name)
        if result and result.get('smiles'):
            return result['smiles']
        return None
    
    def search_multiple(self, names: List[str]) -> Dict[str, Optional[str]]:
        """
        Search multiple compound names in parallel (with rate limiting).
        Returns dict mapping names to SMILES.
        """
        results = {}
        for name in names:
            smiles = self.get_smiles_from_pubchem(name)
            results[name] = smiles
        return results

