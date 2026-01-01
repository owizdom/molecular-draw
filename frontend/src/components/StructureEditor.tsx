import { useState, useCallback, useMemo } from 'react';
import { MolecularStructure, Atom, Bond } from '../utils/molecular';
import { useMolecularStore } from '../store/molecularStore';

interface StructureEditorProps {
  onAtomDelete?: (atomId: number) => void;
}

const ELEMENTS = [
  // Period 1
  'H', 'He',
  // Period 2
  'Li', 'Be', 'B', 'C', 'N', 'O', 'F', 'Ne',
  // Period 3
  'Na', 'Mg', 'Al', 'Si', 'P', 'S', 'Cl', 'Ar',
  // Period 4
  'K', 'Ca', 'Sc', 'Ti', 'V', 'Cr', 'Mn', 'Fe', 'Co', 'Ni', 'Cu', 'Zn', 'Ga', 'Ge', 'As', 'Se', 'Br', 'Kr',
  // Period 5
  'Rb', 'Sr', 'Y', 'Zr', 'Nb', 'Mo', 'Tc', 'Ru', 'Rh', 'Pd', 'Ag', 'Cd', 'In', 'Sn', 'Sb', 'Te', 'I', 'Xe',
  // Period 6
  'Cs', 'Ba', 'La', 'Ce', 'Pr', 'Nd', 'Pm', 'Sm', 'Eu', 'Gd', 'Tb', 'Dy', 'Ho', 'Er', 'Tm', 'Yb', 'Lu',
  'Hf', 'Ta', 'W', 'Re', 'Os', 'Ir', 'Pt', 'Au', 'Hg', 'Tl', 'Pb', 'Bi', 'Po', 'At', 'Rn',
  // Period 7
  'Fr', 'Ra', 'Ac', 'Th', 'Pa', 'U', 'Np', 'Pu', 'Am', 'Cm', 'Bk', 'Cf', 'Es', 'Fm', 'Md', 'No', 'Lr',
  'Rf', 'Db', 'Sg', 'Bh', 'Hs', 'Mt', 'Ds', 'Rg', 'Cn', 'Nh', 'Fl', 'Mc', 'Lv', 'Ts', 'Og'
];

export default function StructureEditor({ onAtomDelete }: StructureEditorProps) {
  const { structure, selectedAtomId, updateStructure } = useMolecularStore();
  const [selectedElement, setSelectedElement] = useState<string>('C');
  const [editingMode, setEditingMode] = useState<'add' | 'delete'>('add');

  const handleAddAtom = useCallback(() => {
    if (!structure) return;
    
    // Get the next available ID
    const maxId = structure.atoms.length > 0 
      ? Math.max(...structure.atoms.map(a => a.id ?? 0))
      : -1;
    const newAtomId = maxId + 1;
    
    let newAtom: Atom;
    let newBonds: Bond[] = [...structure.bonds];
    
    // If an atom is selected, connect the new atom to it
    if (selectedAtomId !== null && selectedAtomId !== undefined) {
      const selectedAtom = structure.atoms.find(a => a.id === selectedAtomId);
      if (selectedAtom) {
        // Position new atom near the selected atom (1.5 units away, typical bond length)
        const offset = 1.5;
        const angle = Math.random() * Math.PI * 2;
        const elevation = Math.random() * Math.PI;
        
        newAtom = {
          element: selectedElement,
          x: selectedAtom.x + offset * Math.sin(elevation) * Math.cos(angle),
          y: selectedAtom.y + offset * Math.sin(elevation) * Math.sin(angle),
          z: selectedAtom.z + offset * Math.cos(elevation),
          id: newAtomId
        };
        
        // Create a bond between the selected atom and the new atom
        newBonds.push({
          atom1_id: selectedAtomId,
          atom2_id: newAtomId,
          bond_type: 1 // Single bond
        });
      } else {
        // Selected atom not found, place randomly
        newAtom = {
          element: selectedElement,
          x: Math.random() * 5 - 2.5,
          y: Math.random() * 5 - 2.5,
          z: Math.random() * 5 - 2.5,
          id: newAtomId
        };
      }
    } else {
      // No atom selected, place randomly
      newAtom = {
        element: selectedElement,
        x: Math.random() * 5 - 2.5,
        y: Math.random() * 5 - 2.5,
        z: Math.random() * 5 - 2.5,
        id: newAtomId
      };
    }
    
    updateStructure({ 
      ...structure, 
      atoms: [...structure.atoms, newAtom],
      bonds: newBonds
    });
  }, [structure, selectedElement, selectedAtomId, updateStructure]);

  const handleDeleteAtom = useCallback((atomId: number) => {
    if (!structure) return;
    const newAtoms = structure.atoms.filter(a => a.id !== atomId);
    const newBonds = structure.bonds.filter(
      b => b.atom1_id !== atomId && b.atom2_id !== atomId
    );
    updateStructure({ ...structure, atoms: newAtoms, bonds: newBonds });
  }, [structure, updateStructure]);

  const handleClear = useCallback(() => {
    updateStructure({ atoms: [], bonds: [] });
  }, [updateStructure]);

  const stats = useMemo(() => ({
    atoms: structure?.atoms.length || 0,
    bonds: structure?.bonds.length || 0,
  }), [structure]);

  if (!structure) {
    return (
      <div className="p-4 bg-white/60 border border-slate-200/50 rounded-lg">
        <p className="text-sm text-slate-500">No structure loaded. Create one using the AI chat!</p>
      </div>
    );
  }

  return (
    <div className="bg-white/80 border border-slate-200/60 rounded-lg p-3 space-y-3">
      <div className="flex items-center justify-between">
        <h3 className="text-xs font-semibold text-slate-800">Editor</h3>
        <div className="text-xs text-slate-400">
          {stats.atoms} · {stats.bonds}
        </div>
      </div>
      
      {/* Mode Toggle */}
      <div className="flex gap-1 p-0.5 bg-slate-100 rounded-md">
        <button
          onClick={() => setEditingMode('add')}
          className={`flex-1 px-2 py-1 rounded text-xs font-medium transition-colors ${
            editingMode === 'add'
              ? 'bg-white text-slate-800 shadow-sm'
              : 'text-slate-500 hover:text-slate-700'
          }`}
        >
          Add
        </button>
        <button
          onClick={() => setEditingMode('delete')}
          className={`flex-1 px-2 py-1 rounded text-xs font-medium transition-colors ${
            editingMode === 'delete'
              ? 'bg-white text-slate-800 shadow-sm'
              : 'text-slate-500 hover:text-slate-700'
          }`}
        >
          Delete
        </button>
      </div>

      {editingMode === 'add' && (
        <div className="space-y-2">
          <div>
            <label className="block text-xs font-medium text-slate-600 mb-1">Element</label>
            <select
              value={selectedElement}
              onChange={(e) => setSelectedElement(e.target.value)}
              className="w-full bg-white border border-slate-200/60 text-slate-700 rounded-md px-2.5 py-1.5 text-xs focus:outline-none focus:ring-1 focus:ring-blue-400 focus:border-blue-300"
            >
              {ELEMENTS.map(el => (
                <option key={el} value={el}>{el}</option>
              ))}
            </select>
          </div>

          {selectedAtomId !== null && selectedAtomId !== undefined ? (
            <div className="p-2 bg-blue-50 rounded border border-blue-200/60">
              <p className="text-xs text-blue-700 mb-1.5">
                Atom {selectedAtomId} - will connect
              </p>
              <button
                onClick={handleAddAtom}
                className="w-full bg-blue-500 hover:bg-blue-600 text-white px-3 py-1.5 rounded-md text-xs font-medium transition-colors"
              >
                Add Connected
              </button>
            </div>
          ) : (
            <div className="space-y-1.5">
              <p className="text-xs text-slate-400">
                Click atom to connect
              </p>
              <button
                onClick={handleAddAtom}
                className="w-full bg-blue-500 hover:bg-blue-600 text-white px-3 py-1.5 rounded-md text-xs font-medium transition-colors"
              >
                Add Atom
              </button>
            </div>
          )}
        </div>
      )}

      {editingMode === 'delete' && (
        <div className="space-y-2">
          <p className="text-xs text-slate-400">Click atoms to select</p>
          {selectedAtomId !== null && selectedAtomId !== undefined && (
            <div className="space-y-1.5 p-2 bg-slate-50 rounded border border-slate-200/60">
              <p className="text-xs font-medium text-slate-700">Atom {selectedAtomId}</p>
              <button
                onClick={() => onAtomDelete?.(selectedAtomId)}
                className="w-full bg-red-500 hover:bg-red-600 text-white px-3 py-1.5 rounded-md text-xs font-medium transition-colors"
              >
                Delete
              </button>
            </div>
          )}
        </div>
      )}

      <div className="pt-2 border-t border-slate-200/60">
        <button
          onClick={handleClear}
          className="w-full bg-slate-100 hover:bg-slate-200 text-slate-600 px-3 py-1.5 rounded-md text-xs font-medium transition-colors"
        >
          Clear
        </button>
      </div>
    </div>
  );
}
