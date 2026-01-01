import { useEffect, useState, useCallback, useMemo } from 'react';
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
      <div className="p-4 bg-white/60 border border-slate-200/50 rounded-lg">
        <h3 className="text-sm font-semibold text-slate-900 mb-2">Properties</h3>
        <p className="text-xs text-slate-500">No molecule loaded</p>
      </div>
    );
  }

  return (
    <div className="p-4 bg-white/60 border border-slate-200/50 rounded-lg space-y-4 shadow-sm">
      <h3 className="text-sm font-semibold text-slate-900">Properties</h3>

      {loading ? (
        <div className="flex items-center gap-2 text-sm text-slate-500">
          <div className="w-4 h-4 border-2 border-slate-300 border-t-blue-500 rounded-full animate-spin"></div>
          Loading...
        </div>
      ) : (
        <>
          <div className="space-y-3">
            {smiles && (
              <div>
                <label className="text-xs font-medium text-slate-700 block mb-1.5">SMILES</label>
                <div className="bg-slate-50 border border-slate-200 p-2.5 rounded-lg text-xs font-mono break-all text-slate-700">
                  {smiles}
                </div>
                <button
                  onClick={() => {
                    navigator.clipboard.writeText(smiles);
                    toast.success('Copied');
                  }}
                  className="mt-1.5 text-xs text-blue-600 hover:text-blue-700 font-medium"
                >
                  Copy
                </button>
              </div>
            )}

            {properties && (
              <div className="grid grid-cols-2 gap-3">
                {properties.formula && (
                  <div>
                    <label className="text-xs font-medium text-slate-500 block mb-1">Formula</label>
                    <div className="text-sm font-semibold text-slate-900">{properties.formula}</div>
                  </div>
                )}

                {properties.molecular_weight && (
                  <div>
                    <label className="text-xs font-medium text-slate-500 block mb-1">Weight</label>
                    <div className="text-sm font-semibold text-slate-900">{properties.molecular_weight.toFixed(2)} g/mol</div>
                  </div>
                )}

                {properties.num_atoms !== undefined && (
                  <div>
                    <label className="text-xs font-medium text-slate-500 block mb-1">Atoms</label>
                    <div className="text-sm font-semibold text-slate-900">{properties.num_atoms}</div>
                  </div>
                )}

                {properties.num_bonds !== undefined && (
                  <div>
                    <label className="text-xs font-medium text-slate-500 block mb-1">Bonds</label>
                    <div className="text-sm font-semibold text-slate-900">{properties.num_bonds}</div>
                  </div>
                )}

                {properties.num_rings !== undefined && (
                  <div>
                    <label className="text-xs font-medium text-slate-500 block mb-1">Rings</label>
                    <div className="text-sm font-semibold text-slate-900">{properties.num_rings}</div>
                  </div>
                )}
              </div>
            )}
          </div>

          <div className="pt-3 border-t border-slate-200 space-y-2">
            <button
              onClick={handleOptimize}
              className="w-full bg-gradient-to-r from-green-500 to-green-600 hover:from-green-600 hover:to-green-700 text-white px-4 py-2 rounded-lg text-sm font-medium shadow-sm hover:shadow-md transition-all active:scale-95"
            >
              Optimize Structure
            </button>
          </div>

          <div className="pt-3 border-t border-slate-200 space-y-2">
            <label className="text-xs font-medium text-slate-700 block">Export</label>
            <div className="grid grid-cols-3 gap-2">
              <button
                onClick={() => handleExport('mol')}
                className="px-3 py-2 bg-white border border-slate-200 hover:border-blue-500 hover:bg-blue-50 text-slate-700 hover:text-blue-700 rounded-lg text-xs font-medium transition-all active:scale-95"
              >
                MOL
              </button>
              <button
                onClick={() => handleExport('sdf')}
                className="px-3 py-2 bg-white border border-slate-200 hover:border-blue-500 hover:bg-blue-50 text-slate-700 hover:text-blue-700 rounded-lg text-xs font-medium transition-all active:scale-95"
              >
                SDF
              </button>
              <button
                onClick={() => handleExport('json')}
                className="px-3 py-2 bg-white border border-slate-200 hover:border-blue-500 hover:bg-blue-50 text-slate-700 hover:text-blue-700 rounded-lg text-xs font-medium transition-all active:scale-95"
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
