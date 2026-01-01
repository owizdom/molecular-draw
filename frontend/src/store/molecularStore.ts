import { create } from 'zustand';
import { MolecularStructure } from '../utils/molecular';
import { TaskProgress } from '../hooks/useSSE';

interface MolecularStore {
  // Structure state
  structure: MolecularStructure | null;
  selectedAtomId: number | null;
  
  // Task state
  currentTaskId: string | null;
  taskProgress: TaskProgress | null;
  isTaskRunning: boolean;
  
  // Reaction state
  reactionCollision: any | null;
  isReactionAnimating: boolean;
  
  // Visualization state
  representationMode: 'ball-and-stick' | 'wireframe' | 'stick' | 'van-der-waals' | 'line';
  
  // Actions
  setStructure: (structure: MolecularStructure | null) => void;
  setSelectedAtomId: (atomId: number | null) => void;
  setCurrentTaskId: (taskId: string | null) => void;
  setTaskProgress: (progress: TaskProgress | null) => void;
  setIsTaskRunning: (running: boolean) => void;
  setReactionCollision: (collision: any | null) => void;
  setIsReactionAnimating: (animating: boolean) => void;
  setRepresentationMode: (mode: 'ball-and-stick' | 'wireframe' | 'stick' | 'van-der-waals' | 'line') => void;
  
  // Helper actions
  updateStructure: (structure: MolecularStructure) => void;
  clearTask: () => void;
}

export const useMolecularStore = create<MolecularStore>((set) => ({
  // Initial state
  structure: null,
  selectedAtomId: null,
  currentTaskId: null,
  taskProgress: null,
  isTaskRunning: false,
  reactionCollision: null,
  isReactionAnimating: false,
  representationMode: 'ball-and-stick',
  
  // Actions
  setStructure: (structure) => set({ structure }),
  setSelectedAtomId: (atomId) => set({ selectedAtomId: atomId }),
  setCurrentTaskId: (taskId) => set({ currentTaskId: taskId }),
  setTaskProgress: (progress) => set({ 
    taskProgress: progress,
    isTaskRunning: progress?.status === 'processing' || progress?.status === 'pending'
  }),
  setIsTaskRunning: (running) => set({ isTaskRunning: running }),
  setReactionCollision: (collision) => set({ reactionCollision: collision }),
  setIsReactionAnimating: (animating) => set({ isReactionAnimating: animating }),
  setRepresentationMode: (mode) => set({ representationMode: mode }),
  
  // Helper actions
  updateStructure: (structure) => set({ structure }),
  clearTask: () => set({ 
    currentTaskId: null, 
    taskProgress: null, 
    isTaskRunning: false 
  }),
}));

