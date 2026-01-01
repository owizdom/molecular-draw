import { useEffect, useState, useCallback } from 'react';
import { useMolecularStore } from '../store/molecularStore';
import { molecularAPI } from '../services/api';
import toast from 'react-hot-toast';

export default function PropertyPanel() {
  const { structure, updateStructure } = useMolecularStore();
  const [properties, setProperties] = useState<Record<string, any> | null>(null);
  const [loading, setLoading] = useState(false);
  const [smiles, setSmiles] = useState<string>('');

  useEffect(() => {
    if (!structure) {
      setProperties(null);
      setSmiles('');
      return;
    }

    const fetchProperties = async () => {
      setLoading(true);
      try {
        const smilesData = await molecularAPI.structureToSmiles(structure);
        setSmiles(smilesData.smiles);
        const props = await molecularAPI.getProperties(smilesData.smiles, undefined);
        setProperties(props);
      } catch (error) {
        console.error('Error fetching properties:', error);
        setProperties(null);
      } finally {
        setLoading(false);
      }
    };

    fetchProperties();
  }, [structure]);

  const handleOptimize = useCallback(async () => {
    if (!structure) return;
    try {
      const result = await molecularAPI.optimizeStructure(structure, false);
      if ('structure' in result && result.structure) {
        updateStructure(result.structure);
        toast.success('Optimized');
      } else {
        toast.error('Optimization failed');
      }
    } catch (error) {
      toast.error('Optimization failed');
    }
  }, [structure, updateStructure]);

  const handleExport = useCallback(async (format: 'mol' | 'sdf' | 'json') => {
    if (!structure) {
      toast.error('No structure to export');
      return;
    }
    
    try {
      const blob = await molecularAPI.exportStructure(structure, format);
      const url = window.URL.createObjectURL(blob);
      const a = document.createElement('a');
      a.href = url;
      a.download = `molecule.${format}`;
      document.body.appendChild(a);
      a.click();
      window.URL.revokeObjectURL(url);
      document.body.removeChild(a);
      toast.success(`Exported as ${format.toUpperCase()}`);
    } catch (error) {
      toast.error('Export failed');
    }
  }, [structure]);

  if (!structure) {
    return (
      <div className="bg-white/80 border border-slate-200/60 rounded-lg p-3">
        <h3 className="text-xs font-semibold text-slate-800 mb-1">Properties</h3>
        <p className="text-xs text-slate-400">No molecule loaded</p>
      </div>
    );
  }

  return (
    <div className="bg-white/80 border border-slate-200/60 rounded-lg p-3 space-y-3">
      <h3 className="text-xs font-semibold text-slate-800">Properties</h3>

      {loading ? (
        <div className="flex items-center gap-2 text-xs text-slate-400">
          <div className="w-3 h-3 border-2 border-slate-300 border-t-blue-500 rounded-full animate-spin"></div>
          Loading...
        </div>
      ) : (
        <>
          <div className="space-y-2">
            {smiles && (
              <div>
                <label className="text-xs font-medium text-slate-600 block mb-1">SMILES</label>
                <div className="bg-slate-50 border border-slate-200/60 p-1.5 rounded text-xs font-mono break-all text-slate-700">
                  {smiles}
                </div>
                <button
                  onClick={() => {
                    navigator.clipboard.writeText(smiles);
                    toast.success('Copied');
                  }}
                  className="mt-1 text-xs text-blue-500 hover:text-blue-600"
                >
                  Copy
                </button>
              </div>
            )}

            {properties && (
              <div className="grid grid-cols-2 gap-2">
                {properties.formula && (
                  <div>
                    <label className="text-xs font-medium text-slate-500 block mb-0.5">Formula</label>
                    <div className="text-xs font-semibold text-slate-800">{properties.formula}</div>
                  </div>
                )}

                {properties.molecular_weight && (
                  <div>
                    <label className="text-xs font-medium text-slate-500 block mb-0.5">Weight</label>
                    <div className="text-xs font-semibold text-slate-800">{properties.molecular_weight.toFixed(2)}</div>
                  </div>
                )}

                {properties.num_atoms !== undefined && (
                  <div>
                    <label className="text-xs font-medium text-slate-500 block mb-0.5">Atoms</label>
                    <div className="text-xs font-semibold text-slate-800">{properties.num_atoms}</div>
                  </div>
                )}

                {properties.num_bonds !== undefined && (
                  <div>
                    <label className="text-xs font-medium text-slate-500 block mb-0.5">Bonds</label>
                    <div className="text-xs font-semibold text-slate-800">{properties.num_bonds}</div>
                  </div>
                )}

                {properties.num_rings !== undefined && (
                  <div>
                    <label className="text-xs font-medium text-slate-500 block mb-0.5">Rings</label>
                    <div className="text-xs font-semibold text-slate-800">{properties.num_rings}</div>
                  </div>
                )}
              </div>
            )}
          </div>

          <div className="pt-2 border-t border-slate-200/60">
            <button
              onClick={handleOptimize}
              className="w-full bg-green-500 hover:bg-green-600 text-white px-3 py-1.5 rounded-md text-xs font-medium transition-colors"
            >
              Optimize
            </button>
          </div>

          <div className="pt-2 border-t border-slate-200/60">
            <label className="text-xs font-medium text-slate-600 block mb-1">Export</label>
            <div className="grid grid-cols-3 gap-1.5">
              <button
                onClick={() => handleExport('mol')}
                className="px-2 py-1 bg-slate-50 hover:bg-slate-100 border border-slate-200/60 text-slate-600 hover:text-slate-800 rounded text-xs font-medium transition-colors"
              >
                MOL
              </button>
              <button
                onClick={() => handleExport('sdf')}
                className="px-2 py-1 bg-slate-50 hover:bg-slate-100 border border-slate-200/60 text-slate-600 hover:text-slate-800 rounded text-xs font-medium transition-colors"
              >
                SDF
              </button>
              <button
                onClick={() => handleExport('json')}
                className="px-2 py-1 bg-slate-50 hover:bg-slate-100 border border-slate-200/60 text-slate-600 hover:text-slate-800 rounded text-xs font-medium transition-colors"
              >
                JSON
              </button>
            </div>
          </div>
        </>
      )}
    </div>
  );
}
