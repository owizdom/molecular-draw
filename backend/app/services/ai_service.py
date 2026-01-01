import os
from typing import Optional
from app.models.molecular import MolecularStructure, ChatMessage, ChatResponse
from app.services.molecular_service import MolecularService
from app.services.pubchem_service import PubChemService
from dotenv import load_dotenv
import re
import requests
import time

load_dotenv()

# Try to import OpenAI (optional)
try:
    from openai import OpenAI
    OPENAI_AVAILABLE = True
except ImportError:
    OPENAI_AVAILABLE = False

# Try to import Google Gemini (optional)
try:
    import google.genai as genai
    GEMINI_AVAILABLE = True
except ImportError:
    try:
        # Fallback to old package
        import google.generativeai as genai
        GEMINI_AVAILABLE = True
    except ImportError:
        GEMINI_AVAILABLE = False

# Try to import Hugging Face (optional)
try:
    from huggingface_hub import InferenceClient
    HF_AVAILABLE = True
except ImportError:
    HF_AVAILABLE = False

class AIService:
    def __init__(self):
        # Initialize PubChem service
        self.pubchem_service = PubChemService()
        
        # Initialize OpenAI if available
        self.openai_client = None
        if OPENAI_AVAILABLE:
            api_key = os.getenv("OPENAI_API_KEY")
            if api_key:
                self.openai_client = OpenAI(api_key=api_key)
        
        # Initialize Google Gemini if available (FREE TIER)
        self.gemini_client = None
        self.gemini_api_key = None
        if GEMINI_AVAILABLE:
            api_key = os.getenv("GEMINI_API_KEY")
            if api_key:
                self.gemini_api_key = api_key
                try:
                    # Use old google.generativeai package (most reliable)
                    import google.generativeai as genai
                    genai.configure(api_key=api_key)
                    # Use gemini-pro (most compatible)
                    self.gemini_client = genai.GenerativeModel('gemini-pro')
                    self.gemini_use_old = True
                except Exception as e:
                    print(f"Gemini configuration error: {e}")
                    # Try new package as fallback
                    try:
                        import google.genai as genai_new
                        client = genai_new.Client(api_key=api_key)
                        model = client.models.get('gemini-pro')
                        self.gemini_client = model
                        self.gemini_use_new = True
                    except:
                        pass
        
        # Initialize Hugging Face if available (FREE TIER)
        self.hf_client = None
        if HF_AVAILABLE:
            api_key = os.getenv("HF_API_KEY")  # Optional, works without key for some models
            self.hf_client = InferenceClient(api_key=api_key) if api_key else InferenceClient()
        
        # Determine which provider to use
        self.provider = self._determine_provider()
        print(f"AI Service initialized with provider: {self.provider}")
    
    def _determine_provider(self) -> str:
        """Determine which AI provider to use based on available keys."""
        if self.openai_client:
            return "openai"
        elif self.gemini_client:
            return "gemini"
        elif self.hf_client:
            return "huggingface"
        else:
            return "fallback"
    
    def _get_chemistry_system_prompt(self) -> str:
        """Get the chemistry-focused system prompt."""
        return """You are an expert chemistry AI assistant specialized in molecular modeling and drug design. 
You help users create, modify, and understand molecular structures. You can:
- Generate molecular structures from natural language descriptions
- Suggest modifications to existing structures
- Explain chemical properties and interactions
- Validate chemical structures
- Provide synthesis pathway suggestions

Always respond in a helpful, educational manner. When suggesting structures, provide SMILES notation when possible."""
    
    def _extract_smiles_from_text(self, text: str) -> Optional[str]:
        """Extract SMILES string from text."""
        # Clean text first
        text = text.strip()
        
        # Common words to exclude
        exclude_words = {'notation', 'string', 'code', 'format', 'is', 'the', 'a', 'an', 'for', 'of', 'create', 'generate'}
        
        # Try various patterns (most specific first)
        patterns = [
            r'SMILES[:\s=]+([A-Za-z0-9@\[\]()=+\-\\/]{2,})',
            r'Canonical\s+SMILES[:\s=]+([A-Za-z0-9@\[\]()=+\-\\/]{2,})',
            r'([A-Za-z0-9@\[\]()=+\-\\/]{3,})',  # General SMILES pattern - must be at least 3 chars
        ]
        
        for pattern in patterns:
            matches = re.findall(pattern, text, re.IGNORECASE)
            for match in matches:
                smiles = match.strip()
                # Validate it looks like SMILES
                # Must have at least 2 chars, contain letters, and not be common words
                if (len(smiles) >= 2 and 
                    any(c.isalpha() for c in smiles) and 
                    smiles.lower() not in exclude_words and
                    not smiles.lower().startswith(('the ', 'is ', 'for ', 'of '))):
                    # Additional validation: should start with a letter, number, or bracket
                    if smiles[0].isalnum() or smiles[0] in '[]()':
                        # Check if it contains valid SMILES characters
                        valid_chars = sum(1 for c in smiles if c.isalnum() or c in '@[]()=+-\\/')
                        if valid_chars >= len(smiles) * 0.7:  # 70% valid characters
                            return smiles
        
        # If no pattern match, try to extract from the whole text if it looks like SMILES
        # Remove common words and see if what's left is SMILES
        cleaned = re.sub(r'\b(SMILES|notation|string|code|format|is|the|a|an|for|of|create|generate)\b', '', text, flags=re.IGNORECASE)
        cleaned = cleaned.strip().strip('.,:;!?')
        
        # Check if cleaned text is mostly SMILES characters
        if len(cleaned) >= 2:
            valid_chars = sum(1 for c in cleaned if c.isalnum() or c in '@[]()=+-\\/')
            if valid_chars >= len(cleaned) * 0.8 and any(c.isalpha() for c in cleaned):
                # Make sure it doesn't start with excluded words
                first_word = cleaned.split()[0].lower() if cleaned.split() else ''
                if first_word not in exclude_words:
                    return cleaned
        
        return None
    
    def _get_common_molecules(self) -> dict:
        """Get a dictionary of common molecules for fallback."""
        return {
            # Basic molecules
            'water': 'O',
            'h2o': 'O',
            'ammonia': 'N',
            'nh3': 'N',
            'carbon dioxide': 'O=C=O',
            'co2': 'O=C=O',
            'hydrogen peroxide': 'OO',
            'h2o2': 'OO',
            
            # Alkanes
            'methane': 'C',
            'ch4': 'C',
            'ethane': 'CC',
            'c2h6': 'CC',
            'propane': 'CCC',
            'c3h8': 'CCC',
            'butane': 'CCCC',
            'c4h10': 'CCCC',
            'pentane': 'CCCCC',
            'c5h12': 'CCCCC',
            'hexane': 'CCCCCC',
            'c6h14': 'CCCCCC',
            'heptane': 'CCCCCCC',
            'c7h16': 'CCCCCCC',
            'octane': 'CCCCCCCC',
            'c8h18': 'CCCCCCCC',
            'nonane': 'CCCCCCCCC',
            'c9h20': 'CCCCCCCCC',
            'decane': 'CCCCCCCCCC',
            'c10h22': 'CCCCCCCCCC',
            'neopentane': 'CC(C)(C)C',
            'isobutanol': 'CC(C)CO',
            
            # Alkenes
            'ethene': 'C=C',
            'ethylene': 'C=C',
            'c2h4': 'C=C',
            'propene': 'CC=C',
            'propylene': 'CC=C',
            'c3h6': 'CC=C',
            '1-butene': 'C=CCC',
            'c4h8': 'C=CCC',  # Default to 1-butene
            '2-butene': 'CC=CC',
            'styrene': 'c1ccc(cc1)C=C',
            'c8h8': 'c1ccc(cc1)C=C',
            
            # Alkynes
            'ethyne': 'C#C',
            'acetylene': 'C#C',
            'c2h2': 'C#C',
            'propyne': 'CC#C',
            'c3h4': 'CC#C',
            
            # Cycloalkanes
            'cyclopropane': 'C1CC1',
            'c3h6': 'C1CC1',  # Note: conflicts with propene, but cyclopropane is less common
            'cyclobutane': 'C1CCC1',
            'c4h8': 'C1CCC1',  # Note: conflicts with butene
            'cyclopentane': 'C1CCCC1',
            'c5h10': 'C1CCCC1',
            'cyclohexane': 'C1CCCCC1',
            'c6h12': 'C1CCCCC1',
            
            # Alcohols
            'methanol': 'CO',
            'ch4o': 'CO',
            'ethanol': 'CCO',
            'c2h6o': 'CCO',
            'c2h5oh': 'CCO',
            'isopropanol': 'CC(C)O',
            'isopropyl alcohol': 'CC(C)O',
            'c3h8o': 'CC(C)O',
            'tert-butanol': 'CC(C)(C)O',
            'tert butanol': 'CC(C)(C)O',
            't-butanol': 'CC(C)(C)O',
            '2-methylpropan-2-ol': 'CC(C)(C)O',
            'c4h10o': 'CC(C)(C)O',
            'propylene glycol': 'OCC(O)CO',
            'c3h8o2': 'OCC(O)CO',
            'glycerol': 'OCC(O)CO',
            'glycerin': 'OCC(O)CO',
            'c3h8o3': 'OCC(O)CO',
            'ethylene glycol': 'OCCO',
            'c2h6o2': 'OCCO',
            
            # Aldehydes and Ketones
            'formaldehyde': 'C=O',
            'ch2o': 'C=O',
            'acetaldehyde': 'CC=O',
            'c2h4o': 'CC=O',
            'acetone': 'CC(=O)C',
            'c3h6o': 'CC(=O)C',
            'butanone': 'CCC(=O)C',
            'c4h8o': 'CCC(=O)C',
            'benzaldehyde': 'c1ccccc1C=O',
            'c7h6o': 'c1ccccc1C=O',
            
            # Carboxylic Acids
            'formic acid': 'C(=O)O',
            'ch2o2': 'C(=O)O',
            'acetic acid': 'CC(=O)O',
            'c2h4o2': 'CC(=O)O',
            'propionic acid': 'CCC(=O)O',
            'c3h6o2': 'CCC(=O)O',
            'butyric acid': 'CCCC(=O)O',
            'c4h8o2': 'CCCC(=O)O',
            'valeric acid': 'CCCCC(=O)O',
            'c5h10o2': 'CCCCC(=O)O',
            'caproic acid': 'CCCCCC(=O)O',
            'c6h12o2': 'CCCCCC(=O)O',
            'lactic acid': 'CC(O)C(=O)O',
            'c3h6o3': 'CC(O)C(=O)O',
            'citric acid': 'C(C(=O)O)C(CC(=O)O)(C(=O)O)O',
            'c6h8o7': 'C(C(=O)O)C(CC(=O)O)(C(=O)O)O',
            'succinic acid': 'O=C(O)CC(=O)O',
            'c4h6o4': 'O=C(O)CC(=O)O',
            'malonic acid': 'C(C(=O)O)C(=O)O',
            'c3h4o4': 'C(C(=O)O)C(=O)O',
            'benzoic acid': 'c1ccccc1C(=O)O',
            'c7h6o2': 'c1ccccc1C(=O)O',
            'terephthalic acid': 'OC(=O)c1ccc(cc1)C(=O)O',
            'c8h6o4': 'OC(=O)c1ccc(cc1)C(=O)O',
            
            # Ethers
            'dimethyl ether': 'COC',
            'c2h6o': 'COC',  # Note: conflicts with ethanol
            'diethyl ether': 'CCOCC',
            'c4h10o': 'CCOCC',  # Note: conflicts with tert-butanol
            'tetrahydrofuran': 'C1CCOC1',
            'thf': 'C1CCOC1',
            'c4h8o': 'C1CCOC1',  # Note: conflicts with butanone
            'tetrahydropyran': 'C1CCOCCC1',
            'c5h10o': 'C1CCOCCC1',
            'tert-butyl methyl ether': 'COC(C)(C)C',
            'tbme': 'COC(C)(C)C',
            'c5h12o': 'COC(C)(C)C',
            
            # Aromatic Heterocycles
            'furan': 'c1ccco1',
            'c4h4o': 'c1ccco1',
            'pyridine': 'c1ccncc1',
            'c5h5n': 'c1ccncc1',
            'pyrrole': 'c1ccc[nH]1',
            'c4h5n': 'c1ccc[nH]1',
            'piperidine': 'C1CCNCC1',
            'c5h11n': 'C1CCNCC1',
            'indole': 'c1ccc2c(c1)[nH]cc2',
            'c8h7n': 'c1ccc2c(c1)[nH]cc2',
            'quinoline': 'c1ccc2cccnc2c1',
            'c9h7n': 'c1ccc2cccnc2c1',
            
            # Aromatic Compounds
            'benzene': 'c1ccccc1',
            'c6h6': 'c1ccccc1',
            'toluene': 'Cc1ccccc1',
            'c7h8': 'Cc1ccccc1',
            'xylene': 'Cc1cccc(C)c1',
            'o-xylene': 'Cc1cccc(C)c1',
            'c8h10': 'Cc1cccc(C)c1',
            'naphthalene': 'c1cccc2ccccc12',
            'c10h8': 'c1cccc2ccccc12',
            'pyrene': 'c1cc2cccc3c2c4c1cccc4cc3',
            'c16h10': 'c1cc2cccc3c2c4c1cccc4cc3',
            'phenol': 'c1ccccc1O',
            'c6h6o': 'c1ccccc1O',
            'catechol': 'c1ccc(c(c1)O)O',
            'c6h6o2': 'c1ccc(c(c1)O)O',
            'anisole': 'COc1ccccc1',
            'c7h8o': 'COc1ccccc1',
            
            # Amines
            'methylamine': 'CN',
            'ch5n': 'CN',
            'ethylamine': 'CCN',
            'c2h7n': 'CCN',
            'propylamine': 'CCCN',
            'c3h9n': 'CCCN',
            'butylamine': 'CCCCN',
            'c4h11n': 'CCCCN',
            'aniline': 'c1ccccc1N',
            'c6h7n': 'c1ccccc1N',
            'phenethylamine': 'c1ccc(cc1)CCN',
            'c8h11n': 'c1ccc(cc1)CCN',
            'nitrobenzene': 'c1ccc(cc1)[N+](=O)[O-]',
            'c6h5no2': 'c1ccc(cc1)[N+](=O)[O-]',
            
            # Nitriles
            'hydrogen cyanide': 'C#N',
            'hcn': 'C#N',
            'chn': 'C#N',
            'acetonitrile': 'CC#N',
            'c2h3n': 'CC#N',
            
            # Pharmaceuticals
            'aspirin': 'CC(=O)OC1=CC=CC=C1C(=O)O',
            'c9h8o4': 'CC(=O)OC1=CC=CC=C1C(=O)O',
            'caffeine': 'CN1C=NC2=C1C(=O)N(C(=O)N2C)C',
            'c8h10n4o2': 'CN1C=NC2=C1C(=O)N(C(=O)N2C)C',
            'melatonin': 'CC(=O)NCCC1=CNc2c1cc(OC)cc2',
            'c13h16n2o2': 'CC(=O)NCCC1=CNc2c1cc(OC)cc2',
            'benzocaine': 'COc1ccc(cc1)CC(=O)N',
            'c9h11no2': 'COc1ccc(cc1)CC(=O)N',
            'ibuprofen': 'CC(C)CC1=CC=C(C=C1)C(C)C(=O)O',
            'c13h18o2': 'CC(C)CC1=CC=C(C=C1)C(C)C(=O)O',
            'acetaminophen': 'CC(=O)NC1=CC=C(O)C=C1',
            'paracetamol': 'CC(=O)NC1=CC=C(O)C=C1',
            'c8h9no2': 'CC(=O)NC1=CC=C(O)C=C1',
            'diazepam': 'CN(C)C(=O)CN=C(c1ccccc1)c2ccccc2',
            'c16h13cln2o': 'CN(C)C(=O)CN=C(c1ccccc1)c2ccccc2',
            'nicotine': 'CN1CCC[C@H]1c2cccnc2',
            'c10h14n2': 'CN1CCC[C@H]1c2cccnc2',
            
            # Sugars
            'glucose': 'C([C@@H]1[C@H]([C@@H]([C@H](C(O1)O)O)O)O)O',
            'c6h12o6': 'C([C@@H]1[C@H]([C@@H]([C@H](C(O1)O)O)O)O)O',
            'glyceraldehyde': 'OC[C@H](O)C=O',
            'c3h6o3': 'OC[C@H](O)C=O',
            
            # Halogenated Compounds
            'methyl chloride': 'CCl',
            'ch3cl': 'CCl',
            'chloroethene': 'C=CCl',
            'vinyl chloride': 'C=CCl',
            'c2h3cl': 'C=CCl',
            '1,2-dichloroethane': 'ClCCCl',
            'c2h4cl2': 'ClCCCl',
            'carbon tetrachloride': 'C(Cl)(Cl)(Cl)Cl',
            'ccl4': 'C(Cl)(Cl)(Cl)Cl',
            'methyl bromide': 'CBr',
            'ch3br': 'CBr',
            'ethyl fluoride': 'CCF',
            'c2h5f': 'CCF',
            'benzoyl chloride': 'O=C(c1ccccc1)Cl',
            'c7h5clo': 'O=C(c1ccccc1)Cl',
            
            # Other
            'sulfur hexafluoride': 'FS(F)(F)(F)(F)F',
            'sf6': 'FS(F)(F)(F)(F)F',
            'phthalide': 'O=C1OCc2ccccc12',
            'c8h6o2': 'O=C1OCc2ccccc12',
            'heptazine': 'C1=NC2=NC=NC3=NC=NC(=N1)N23',
            'c6h3n7': 'C1=NC2=NC=NC3=NC=NC(=N1)N23',
            'ethane-1,1-dithiol': 'SC(S)C',
            'c2h6s2': 'SC(S)C',
            'pentaerythritol': 'C(C(CO)(CO)CO)(CO)O',
        }
    
    def _fallback_generate_smiles(self, description: str) -> Optional[str]:
        """Fallback method using common molecules dictionary."""
        description_lower = description.lower().strip()
        
        # Check common molecules
        common = self._get_common_molecules()
        if description_lower in common:
            return common[description_lower]
        
        # Try to find partial matches
        for key, smiles in common.items():
            if key in description_lower or description_lower in key:
                return smiles
        
        # Try to extract if already a SMILES
        smiles_match = re.search(r'([A-Za-z0-9@\[\]()=+\-\\/]{2,})', description)
        if smiles_match:
            potential_smiles = smiles_match.group(1)
            if len(potential_smiles) >= 2:
                return potential_smiles
        
        return None
    
    def _normalize_chemical_formula(self, formula: str) -> str:
        """Convert chemical formula with subscripts to normal text (C₄H₁₀O -> C4H10O)."""
        subscript_map = {
            '₀': '0', '₁': '1', '₂': '2', '₃': '3', '₄': '4', '₅': '5',
            '₆': '6', '₇': '7', '₈': '8', '₉': '9',
            '₊': '+', '₋': '-', '₌': '=', '₍': '(', '₎': ')'
        }
        normalized = formula
        for subscript, normal in subscript_map.items():
            normalized = normalized.replace(subscript, normal)
        return normalized
    
    def _search_web_for_smiles(self, molecule_name: str) -> Optional[str]:
        """Search the web for SMILES notation - use AI instead of web scraping."""
        # Web scraping is unreliable, so we'll skip it and go straight to AI
        # The AI provider should handle molecule name lookups better
        return None
    
    def _call_openai(self, system_prompt: str, user_prompt: str, temperature: float = 0.7) -> Optional[str]:
        """Call OpenAI API."""
        if not self.openai_client:
            return None
        
        try:
            response = self.openai_client.chat.completions.create(
                model="gpt-4",
                messages=[
                    {"role": "system", "content": system_prompt},
                    {"role": "user", "content": user_prompt}
                ],
                temperature=temperature
            )
            return response.choices[0].message.content
        except Exception as e:
            print(f"OpenAI error: {e}")
            return None
    
    def _call_gemini(self, system_prompt: str, user_prompt: str, temperature: float = 0.7) -> Optional[str]:
        """Call Google Gemini API (FREE TIER) with timeout."""
        if not self.gemini_client:
            return None
        
        try:
            # Combine system and user prompt for Gemini
            full_prompt = f"{system_prompt}\n\n{user_prompt}"
            
            # Check if using old API (google.generativeai)
            if hasattr(self, 'gemini_use_old') and self.gemini_use_old:
                # Old google.generativeai API
                try:
                    import google.generativeai as genai
                    response = self.gemini_client.generate_content(
                        full_prompt,
                        generation_config=genai.types.GenerationConfig(
                            temperature=temperature,
                            max_output_tokens=500
                        )
                    )
                    return response.text
                except Exception as e:
                    print(f"Gemini API error: {e}")
                    return None
            else:
                # New API (google.genai)
                try:
                    response = self.gemini_client.generate_content(
                        full_prompt,
                        config={'temperature': temperature, 'max_output_tokens': 500}
                    )
                    return response.text
                except Exception as e:
                    print(f"Gemini API error: {e}")
                    return None
        except Exception as e:
            print(f"Gemini error: {e}")
            return None
    
    def _call_huggingface(self, system_prompt: str, user_prompt: str, temperature: float = 0.7) -> Optional[str]:
        """Call Hugging Face Inference API (FREE TIER)."""
        if not self.hf_client:
            return None
        
        try:
            # Use a chemistry-friendly model
            # Try meta-llama/Llama-2-7b-chat-hf or microsoft/DialoGPT-medium
            # For chemistry, we'll use a general model
            model = "microsoft/DialoGPT-large"  # Free model
            
            # Combine prompts
            full_prompt = f"{system_prompt}\n\nUser: {user_prompt}\nAssistant:"
            
            response = self.hf_client.text_generation(
                full_prompt,
                max_new_tokens=200,
                temperature=temperature,
                return_full_text=False
            )
            return response
        except Exception as e:
            print(f"Hugging Face error: {e}")
            return None
    
    def _call_ai_provider(self, system_prompt: str, user_prompt: str, temperature: float = 0.7) -> Optional[str]:
        """Call the appropriate AI provider."""
        if self.provider == "openai":
            return self._call_openai(system_prompt, user_prompt, temperature)
        elif self.provider == "gemini":
            return self._call_gemini(system_prompt, user_prompt, temperature)
        elif self.provider == "huggingface":
            return self._call_huggingface(system_prompt, user_prompt, temperature)
        else:
            return None
    
    def process_chat_message(self, message: ChatMessage) -> ChatResponse:
        """Process a chat message and return AI response with optional structure modifications."""
        # Fast path: Check for common molecules first (instant response)
        message_lower = message.message.lower().strip()
        if any(word in message_lower for word in ['create', 'generate', 'make', 'show', 'build']):
            smiles = self._fallback_generate_smiles(message.message)
            if smiles:
                structure = MolecularService.smiles_to_structure(smiles)
                if structure:
                    return ChatResponse(
                        response=f"I've created {message.message}. The molecule is now displayed in the 3D viewer.",
                        structure=structure,
                        suggestions=[]
                    )
        
        system_prompt = self._get_chemistry_system_prompt()
        
        user_prompt = message.message
        if message.structure:
            smiles = MolecularService.structure_to_smiles(message.structure)
            if smiles:
                user_prompt += f"\n\nCurrent structure (SMILES): {smiles}"
        
        # Try AI provider first (only for complex queries)
        ai_response = None
        if self.provider != "fallback" and not any(mol in message_lower for mol in ['water', 'h2o', 'methane', 'benzene', 'aspirin', 'caffeine']):
            ai_response = self._call_ai_provider(system_prompt, user_prompt, temperature=0.7)
        
        # Fallback to simple responses if AI fails or for simple queries
        if not ai_response:
            ai_response = self._generate_fallback_response(message.message, message.structure)
        
        # Try to extract SMILES from response
        structure = None
        if message.structure:
            structure = message.structure
        
        smiles_from_response = self._extract_smiles_from_text(ai_response)
        if smiles_from_response:
            new_structure = MolecularService.smiles_to_structure(smiles_from_response)
            if new_structure:
                structure = new_structure
        
        return ChatResponse(
            response=ai_response,
            structure=structure,
            suggestions=self._extract_suggestions(ai_response)
        )
    
    def _generate_fallback_response(self, message: str, structure: Optional[MolecularStructure]) -> str:
        """Generate a fallback response when no AI provider is available."""
        message_lower = message.lower()
        
        # Check for common requests
        if any(word in message_lower for word in ['create', 'generate', 'make', 'show', 'build']):
            smiles = self._fallback_generate_smiles(message)
            if smiles:
                return f"I can help you create that molecule. Here's the SMILES notation: {smiles}\n\nNote: For more advanced AI features, you can set up a free API key for Google Gemini (GEMINI_API_KEY) or Hugging Face (HF_API_KEY) in your .env file."
            else:
                return "I can help with common molecules like water (H2O), benzene, aspirin, caffeine, etc. For more complex molecules, please provide a SMILES string or set up a free AI API key (Gemini or Hugging Face)."
        
        elif any(word in message_lower for word in ['what', 'explain', 'tell me', 'describe']):
            return "I can help explain molecular structures and properties. For detailed chemistry explanations, consider setting up a free Google Gemini API key (GEMINI_API_KEY) in your .env file."
        
        elif any(word in message_lower for word in ['modify', 'change', 'add', 'remove']):
            return "I can help modify structures. Use the structure editor panel to add or remove atoms and bonds. For AI-powered suggestions, set up a free API key."
        
        else:
            return "I'm here to help with molecular modeling! You can:\n- Create molecules using SMILES notation\n- Edit structures manually\n- Ask about common molecules\n\nFor advanced AI features, set up a free Google Gemini API key (GEMINI_API_KEY) in your .env file. Get one at: https://makersuite.google.com/app/apikey"
    
    def _extract_suggestions(self, response: str) -> list:
        """Extract actionable suggestions from AI response."""
        suggestions = []
        response_lower = response.lower()
        if "add" in response_lower or "modify" in response_lower:
            suggestions.append("Apply suggested modifications")
        if "optimize" in response_lower:
            suggestions.append("Optimize structure")
        return suggestions
    
    def generate_structure_from_text(self, description: str) -> Optional[MolecularStructure]:
        """Generate molecular structure from natural language description."""
        # Extract molecule name from description FIRST (before checking common molecules)
        # Remove common verbs like "create", "generate", "make", "show"
        molecule_name = re.sub(r'\b(create|generate|make|show|build|draw|display)\b', '', description, flags=re.IGNORECASE).strip()
        molecule_name = molecule_name.strip('.,!?').strip()
        
        # PRIORITY 1: Try PubChem FIRST (fastest, most reliable for known compounds)
        if molecule_name and len(molecule_name) > 2:
            print(f"Searching PubChem for: {molecule_name}")
            pubchem_result = self.pubchem_service.search_by_name(molecule_name)
            if pubchem_result and pubchem_result.get('smiles'):
                structure = MolecularService.smiles_to_structure(pubchem_result['smiles'])
                if structure:
                    print(f"Found in PubChem (CID: {pubchem_result.get('cid')}): {pubchem_result.get('smiles')}")
                    return structure
        
        # PRIORITY 2: Try fallback (instant for common molecules)
        smiles = self._fallback_generate_smiles(description)
        if smiles:
            structure = MolecularService.smiles_to_structure(smiles)
            if structure:
                return structure
        
        # Check if it's a chemical formula (contains subscripts or numbers)
        # Normalize chemical formulas (C₄H₁₀O -> C4H10O)
        normalized_formula = self._normalize_chemical_formula(molecule_name)
        is_formula = normalized_formula != molecule_name or (any(c.isdigit() for c in molecule_name) and any(c.isupper() for c in molecule_name) and len(molecule_name) <= 20)
        
        if is_formula:
            formula_to_use = normalized_formula if normalized_formula != molecule_name else molecule_name
            
            # Check common formulas first
            formula_lower = formula_to_use.lower().replace(' ', '')
            if formula_lower == 'c4h10o':
                # This is tert-Butanol or butanol isomers - use tert-Butanol
                smiles = 'CC(C)(C)O'
                structure = MolecularService.smiles_to_structure(smiles)
                if structure:
                    return structure
            
            # Try to use AI to convert formula to SMILES
            if self.provider != "fallback":
                prompt = f"Convert this chemical formula to SMILES notation: {formula_to_use}. Provide ONLY the SMILES string, no explanation, no labels."
                system_prompt = "You are a chemistry expert. Convert chemical formulas to SMILES. Respond with ONLY the SMILES string, nothing else."
                ai_response = self._call_ai_provider(system_prompt, prompt, temperature=0.1)
                if ai_response:
                    smiles = self._extract_smiles_from_text(ai_response)
                    if smiles:
                        structure = MolecularService.smiles_to_structure(smiles)
                        if structure:
                            return structure
        
        # If it looks like a molecule name (has capital letters, not just SMILES), search online
        if molecule_name and (any(c.isupper() for c in molecule_name) or any(c.isdigit() for c in molecule_name)) and len(molecule_name) > 2:
            print(f"Searching web for molecule: {molecule_name}")
            
            # Try web search first
            smiles = self._search_web_for_smiles(molecule_name)
            if smiles:
                structure = MolecularService.smiles_to_structure(smiles)
                if structure:
                    print(f"Found SMILES from web search: {smiles}")
                    return structure
            
            # Try AI provider for molecule lookup
            if self.provider != "fallback":
                # Better prompt for molecule name lookup
                prompt = f"What is the SMILES notation for {molecule_name}? Provide ONLY the SMILES string, no explanation."
                system_prompt = "You are a chemistry expert. When asked for a molecule's SMILES notation, respond with ONLY the SMILES string. Do not include any other text, explanations, or labels."
                
                ai_response = self._call_ai_provider(system_prompt, prompt, temperature=0.1)
                
                if ai_response:
                    # Clean the response - remove common prefixes
                    cleaned = ai_response.strip()
                    # Remove common prefixes like "SMILES:", "The SMILES is", etc.
                    cleaned = re.sub(r'^(SMILES|Canonical\s+SMILES|The\s+SMILES)[:\s=]+', '', cleaned, flags=re.IGNORECASE)
                    cleaned = cleaned.strip()
                    
                    # Extract SMILES from response
                    smiles = self._extract_smiles_from_text(cleaned)
                    if not smiles:
                        # Try to clean the response more aggressively
                        # Remove everything except valid SMILES characters
                        smiles = re.sub(r'[^A-Za-z0-9@\[\]()=+\-\\/]', '', cleaned)
                    
                    # Validate SMILES looks reasonable (at least 2 chars, has letters)
                    if smiles and len(smiles) >= 2 and any(c.isalpha() for c in smiles):
                        structure = MolecularService.smiles_to_structure(smiles)
                        if structure:
                            print(f"Found SMILES from AI: {smiles}")
                            return structure
                        else:
                            print(f"Invalid SMILES from AI: {smiles}")
        
        # Fallback: try AI provider with original description
        if self.provider != "fallback":
            prompt = f"SMILES for: {description}"
            system_prompt = "Respond with ONLY the SMILES string, nothing else."
            
            ai_response = self._call_ai_provider(system_prompt, prompt, temperature=0.3)
            
            if ai_response:
                # Extract SMILES from response
                smiles = self._extract_smiles_from_text(ai_response)
                if not smiles:
                    # Try to clean the response
                    smiles = ai_response.strip()
                    smiles = re.sub(r'[^A-Za-z0-9@\[\]()=+\-\\/]', '', smiles)
                
                if smiles:
                    structure = MolecularService.smiles_to_structure(smiles)
                    if structure:
                        return structure
        
        return None
