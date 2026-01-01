import { useMemo, useCallback } from 'react';
import { useMolecularStore } from './store/molecularStore';
import MolecularViewer from './components/MolecularViewer';
import StructureEditor from './components/StructureEditor';
import AIChat from './components/AIChat';
import PropertyPanel from './components/PropertyPanel';
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
      
      {/* Modern Header */}
      <header className="bg-white/80 backdrop-blur-md border-b border-slate-200/50 px-6 py-3 shadow-sm z-10">
        <div className="flex items-center justify-between">
          <div className="flex items-center gap-3">
            <div className="w-8 h-8 bg-gradient-to-br from-blue-500 to-purple-600 rounded-lg flex items-center justify-center">
              <span className="text-white text-sm font-bold">M</span>
            </div>
            <h1 className="text-xl font-semibold bg-gradient-to-r from-slate-900 to-slate-700 bg-clip-text text-transparent">
              Molecular Draw
            </h1>
          </div>
          <div className="flex items-center gap-3">
            <div className="relative">
              <input
                type="text"
                placeholder="SMILES (e.g., CCO)..."
                onKeyPress={(e) => {
                  if (e.key === 'Enter') {
                    const value = e.currentTarget.value.trim();
                    if (value) {
                      handleLoadFromSmiles(value);
                      e.currentTarget.value = '';
                    }
                  }
                }}
                className="bg-slate-50 border border-slate-200 text-slate-900 rounded-lg px-4 py-2 text-sm w-56 focus:outline-none focus:ring-2 focus:ring-blue-500 focus:border-transparent transition-all"
              />
            </div>
            {hasStructure && (
              <button
                onClick={handleCopySmiles}
                className="bg-gradient-to-r from-blue-500 to-blue-600 hover:from-blue-600 hover:to-blue-700 text-white px-4 py-2 rounded-lg text-sm font-medium shadow-sm hover:shadow-md transition-all active:scale-95"
              >
                Copy SMILES
              </button>
            )}
          </div>
        </div>
      </header>

      {/* Main Content - Modern Layout */}
      <div className="flex-1 flex overflow-hidden">
        {/* Left Sidebar - AI Chat */}
        <div className="w-96 bg-white/60 backdrop-blur-sm border-r border-slate-200/50 flex flex-col shadow-sm">
          <AIChat />
        </div>

        {/* Center - 3D Viewer */}
        <div className="flex-1 flex flex-col bg-slate-50 relative">
          <div className="flex-1 relative">
            <MolecularViewer onAtomClick={handleAtomClick} />
          </div>
        </div>

        {/* Right Sidebar - Editor & Properties */}
        <div className="w-96 bg-white/60 backdrop-blur-sm border-l border-slate-200/50 flex flex-col overflow-hidden shadow-sm">
          <div className="flex-1 overflow-y-auto p-4 space-y-4">
            <StructureEditor onAtomDelete={handleAtomDelete} />
            <PropertyPanel />
          </div>
        </div>
      </div>
    </div>
  );
}

export default App;

