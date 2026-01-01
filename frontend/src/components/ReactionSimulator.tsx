import { useState, useCallback } from 'react';
import { useMolecularStore } from '../store/molecularStore';
import { molecularAPI } from '../services/api';
import { MolecularStructure } from '../utils/molecular';
import toast from 'react-hot-toast';

export default function ReactionSimulator() {
  const { structure, setStructure, setReactionCollision, setIsReactionAnimating } = useMolecularStore();
  const [reactant1, setReactant1] = useState<{ smiles: string; structure: MolecularStructure | null } | null>(null);
  const [reactant2, setReactant2] = useState<{ smiles: string; structure: MolecularStructure | null } | null>(null);
  const [reactionResult, setReactionResult] = useState<any>(null);
  const [isAnimating, setIsAnimating] = useState(false);
  const [loading, setLoading] = useState(false);

  const handleSetReactant1 = useCallback(async () => {
    if (!structure) {
      toast.error('No molecule loaded. Create or load a molecule first.');
      return;
    }
    
    try {
      const smilesData = await molecularAPI.structureToSmiles(structure);
      setReactant1({
        smiles: smilesData.smiles,
        structure: structure
      });
      toast.success('Reactant 1 set');
    } catch (error) {
      toast.error('Failed to set reactant 1');
    }
  }, [structure]);

  const handleSetReactant2 = useCallback(async () => {
    if (!structure) {
      toast.error('No molecule loaded. Create or load a molecule first.');
      return;
    }
    
    try {
      const smilesData = await molecularAPI.structureToSmiles(structure);
      setReactant2({
        smiles: smilesData.smiles,
        structure: structure
      });
      toast.success('Reactant 2 set');
    } catch (error) {
      toast.error('Failed to set reactant 2');
    }
  }, [structure]);

  const handleCombine = useCallback(async () => {
    if (!reactant1 || !reactant2) {
      toast.error('Please set both reactants first');
      return;
    }

    setLoading(true);
    setIsAnimating(true);
    setReactionResult(null);

    try {
      const result = await molecularAPI.combineMolecules(
        {
          smiles: reactant1.smiles,
          structure: reactant1.structure || undefined
        },
        {
          smiles: reactant2.smiles,
          structure: reactant2.structure || undefined
        }
      );

      setReactionResult(result);
      
      // Set collision data for animation
      if (result.collision) {
        setReactionCollision(result.collision);
        setIsReactionAnimating(true);
      }
      
      // Show collision animation
      setTimeout(() => {
        setIsAnimating(false);
        setIsReactionAnimating(false);
        setReactionCollision(null);
        
        // Display first product if available
        if (result.products && result.products.length > 0 && result.products[0].structure) {
          setStructure(result.products[0].structure);
          toast.success(`Reaction complete! ${result.reaction_type} reaction detected.`);
        }
      }, 2000); // Animation duration
      
    } catch (error: any) {
      toast.error(error.response?.data?.detail || 'Reaction failed');
      setIsAnimating(false);
    } finally {
      setLoading(false);
    }
  }, [reactant1, reactant2, setStructure]);

  const handleClear = useCallback(() => {
    setReactant1(null);
    setReactant2(null);
    setReactionResult(null);
    setIsAnimating(false);
  }, []);

  return (
    <div className="bg-white/80 border border-slate-200/60 rounded-lg p-3 space-y-2.5">
      <div className="flex items-center justify-between">
        <h3 className="text-xs font-semibold text-slate-800">Reactions</h3>
        {(reactant1 || reactant2) && (
          <button
            onClick={handleClear}
            className="text-xs text-slate-400 hover:text-slate-600"
          >
            Clear
          </button>
        )}
      </div>

      <div className="space-y-2">
        {/* Reactant 1 */}
        <div className="space-y-1">
          <label className="block text-xs font-medium text-slate-600">Reactant 1</label>
          {reactant1 ? (
            <div className="p-1.5 bg-blue-50 rounded border border-blue-200/60">
              <p className="text-xs text-blue-700 font-mono truncate">{reactant1.smiles}</p>
              <button
                onClick={() => setReactant1(null)}
                className="text-xs text-blue-500 hover:text-blue-700 mt-0.5"
              >
                Remove
              </button>
            </div>
          ) : (
            <button
              onClick={handleSetReactant1}
              disabled={!structure}
              className="w-full bg-slate-50 hover:bg-slate-100 text-slate-600 px-2.5 py-1.5 rounded-md text-xs font-medium transition-colors disabled:opacity-40 disabled:cursor-not-allowed"
            >
              Set R1
            </button>
          )}
        </div>

        {/* Reactant 2 */}
        <div className="space-y-1">
          <label className="block text-xs font-medium text-slate-600">Reactant 2</label>
          {reactant2 ? (
            <div className="p-1.5 bg-green-50 rounded border border-green-200/60">
              <p className="text-xs text-green-700 font-mono truncate">{reactant2.smiles}</p>
              <button
                onClick={() => setReactant2(null)}
                className="text-xs text-green-500 hover:text-green-700 mt-0.5"
              >
                Remove
              </button>
            </div>
          ) : (
            <button
              onClick={handleSetReactant2}
              disabled={!structure}
              className="w-full bg-slate-50 hover:bg-slate-100 text-slate-600 px-2.5 py-1.5 rounded-md text-xs font-medium transition-colors disabled:opacity-40 disabled:cursor-not-allowed"
            >
              Set R2
            </button>
          )}
        </div>

        {/* Combine Button */}
        <button
          onClick={handleCombine}
          disabled={!reactant1 || !reactant2 || loading}
          className={`w-full px-3 py-1.5 rounded-md text-xs font-medium transition-colors ${
            isAnimating
              ? 'bg-yellow-500 text-white'
              : 'bg-purple-500 hover:bg-purple-600 text-white'
          } disabled:opacity-50 disabled:cursor-not-allowed`}
        >
          {isAnimating ? 'Colliding...' : loading ? 'Processing...' : 'Combine'}
        </button>

        {/* Reaction Result */}
        {reactionResult && (
          <div className="p-2 bg-slate-50 rounded border border-slate-200/60 space-y-1.5">
            <div className="text-xs font-medium text-slate-700">
              Type: <span className="text-purple-600">{reactionResult.reaction_type}</span>
            </div>
            {reactionResult.reaction_equation && (
              <div className="text-xs text-slate-600 font-mono bg-white p-1.5 rounded text-[10px]">
                {reactionResult.reaction_equation}
              </div>
            )}
            {reactionResult.products && reactionResult.products.length > 0 && (
              <div className="space-y-1">
                <p className="text-xs font-medium text-slate-600">Products:</p>
                {reactionResult.products.map((product: any, idx: number) => (
                  <button
                    key={idx}
                    onClick={() => {
                      if (product.structure) {
                        setStructure(product.structure);
                        toast.success(`Loaded ${product.name || 'Product'}`);
                      }
                    }}
                    className="w-full text-left text-xs text-slate-600 hover:text-slate-800 p-1.5 bg-white rounded border border-slate-200/60 hover:border-purple-300 transition-colors"
                  >
                    {product.name || `P${idx + 1}`}: {product.smiles}
                  </button>
                ))}
              </div>
            )}
          </div>
        )}
      </div>
    </div>
  );
}

