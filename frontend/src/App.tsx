import { useMemo, useCallback } from 'react';
import { useMolecularStore } from './store/molecularStore';
import MolecularViewer from './components/MolecularViewer';
import StructureEditor from './components/StructureEditor';
import AIChat from './components/AIChat';
import PropertyPanel from './components/PropertyPanel';
import ReactionSimulator from './components/ReactionSimulator';
import RepresentationSelector from './components/RepresentationSelector';
import { molecularAPI } from './services/api';
import toast from 'react-hot-toast';
import { Toaster } from 'react-hot-toast';

function App() {
  const { structure, selectedAtomId, setStructure, setSelectedAtomId } = useMolecularStore();

  const handleAtomClick = useCallback((atomId: number) => {
    if (selectedAtomId === atomId) {
      setSelectedAtomId(null);
    } else {
      setSelectedAtomId(atomId);
    }
  }, [selectedAtomId, setSelectedAtomId]);

  const handleLoadFromSmiles = useCallback(async (smiles: string) => {
    try {
      const response = await molecularAPI.structureFromSmiles(smiles, false);
      if ('structure' in response) {
        setStructure(response.structure);
        toast.success('Molecule loaded');
      } else {
        toast.error('Failed to load molecule');
      }
    } catch (error) {
      toast.error('Failed to load molecule');
    }
  }, [setStructure]);

  const handleCopySmiles = useCallback(async () => {
    if (structure) {
      try {
        const smilesData = await molecularAPI.structureToSmiles(structure);
        navigator.clipboard.writeText(smilesData.smiles);
        toast.success('Copied!');
      } catch (error) {
        toast.error('Failed to copy');
      }
    }
  }, [structure]);

  const handleAtomDelete = useCallback((atomId: number) => {
    if (structure) {
      const newAtoms = structure.atoms.filter(a => a.id !== atomId);
      const newBonds = structure.bonds.filter(
        b => b.atom1_id !== atomId && b.atom2_id !== atomId
      );
      setStructure({ ...structure, atoms: newAtoms, bonds: newBonds });
      setSelectedAtomId(null);
    }
  }, [structure, setStructure, setSelectedAtomId]);

  const hasStructure = useMemo(() => !!structure, [structure]);

  return (
    <div className="h-screen w-screen flex flex-col bg-gradient-to-br from-slate-50 to-slate-100 text-slate-900 overflow-hidden">
      <Toaster 
        position="top-right"
        toastOptions={{
          duration: 2000,
          style: {
            background: '#fff',
            color: '#1e293b',
            border: '1px solid #e2e8f0',
            boxShadow: '0 4px 6px -1px rgba(0, 0, 0, 0.1)',
          },
          success: {
            iconTheme: {
              primary: '#10b981',
              secondary: '#fff',
            },
          },
          error: {
            iconTheme: {
              primary: '#ef4444',
              secondary: '#fff',
            },
          },
        }}
      />
      
      {/* Minimal Header */}
      <header className="bg-white/90 backdrop-blur-md border-b border-slate-200/60 px-5 py-2.5 z-10">
        <div className="flex items-center justify-between">
          <div className="flex items-center gap-2.5">
            <div className="w-7 h-7 bg-gradient-to-br from-blue-500 to-indigo-600 rounded-md flex items-center justify-center shadow-sm">
              <span className="text-white text-xs font-semibold">M</span>
            </div>
            <h1 className="text-lg font-semibold text-slate-800">
              Molecular Draw
            </h1>
          </div>
          <div className="flex items-center gap-2">
            <input
              type="text"
              placeholder="SMILES..."
              onKeyPress={(e) => {
                if (e.key === 'Enter') {
                  const value = e.currentTarget.value.trim();
                  if (value) {
                    handleLoadFromSmiles(value);
                    e.currentTarget.value = '';
                  }
                }
              }}
              className="bg-slate-50/80 border border-slate-200 text-slate-700 rounded-md px-3 py-1.5 text-xs w-40 focus:outline-none focus:ring-1 focus:ring-blue-400 focus:border-blue-300 transition-all placeholder:text-slate-400"
            />
            {hasStructure && (
              <button
                onClick={handleCopySmiles}
                className="bg-slate-100 hover:bg-slate-200 text-slate-700 px-3 py-1.5 rounded-md text-xs font-medium transition-colors"
              >
                Copy
              </button>
            )}
          </div>
        </div>
      </header>

      {/* Main Content - Clean Layout */}
      <div className="flex-1 flex overflow-hidden bg-slate-50">
        {/* Left Sidebar - AI Chat */}
        <div className="w-80 bg-white/70 backdrop-blur-sm border-r border-slate-200/60 flex flex-col">
          <AIChat />
        </div>

        {/* Center - 3D Viewer */}
        <div className="flex-1 flex flex-col relative bg-gradient-to-br from-slate-50 to-slate-100">
          <div className="absolute top-3 left-3 z-10">
            <RepresentationSelector />
          </div>
          <div className="flex-1 relative">
            <MolecularViewer onAtomClick={handleAtomClick} />
          </div>
        </div>

        {/* Right Sidebar - Tools */}
        <div className="w-80 bg-white/70 backdrop-blur-sm border-l border-slate-200/60 flex flex-col overflow-hidden">
          <div className="flex-1 overflow-y-auto p-3 space-y-3">
            <ReactionSimulator />
            <StructureEditor onAtomDelete={handleAtomDelete} />
            <PropertyPanel />
          </div>
        </div>
      </div>
    </div>
  );
}

export default App;

