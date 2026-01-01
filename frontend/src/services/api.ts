import axios from 'axios';
import { MolecularStructure } from '../utils/molecular';
import { TaskProgress } from '../hooks/useSSE';

const API_BASE_URL = 'http://localhost:8000';

export interface StructureRequest {
  smiles?: string;
  structure?: MolecularStructure;
}

export interface StructureResponse {
  structure: MolecularStructure;
  molecular_weight?: number;
  formula?: string;
  properties?: Record<string, any>;
}

export interface ChatMessage {
  message: string;
  structure?: MolecularStructure;
}

export interface ChatResponse {
  response: string;
  structure?: MolecularStructure;
  suggestions?: string[];
}

export interface TaskResponse {
  task_id: string;
  status: string;
  message?: string;
}

const api = axios.create({
  baseURL: API_BASE_URL,
  headers: {
    'Content-Type': 'application/json',
  },
});

export const molecularAPI = {
  structureFromSmiles: async (smiles: string, asyncTask: boolean = false): Promise<StructureResponse | TaskResponse> => {
    const response = await api.post(`/structure/from-smiles?async_task=${asyncTask}`, { smiles });
    return response.data;
  },

  structureToSmiles: async (structure: MolecularStructure): Promise<{ smiles: string }> => {
    const response = await api.post('/structure/to-smiles', { structure });
    return response.data;
  },

  getProperties: async (smiles?: string, structure?: MolecularStructure, asyncTask: boolean = false): Promise<Record<string, any> | TaskResponse> => {
    const response = await api.post(`/structure/properties?async_task=${asyncTask}`, { smiles, structure });
    return response.data;
  },

  validateStructure: async (structure: MolecularStructure, asyncTask: boolean = false): Promise<{ is_valid: boolean; error?: string } | TaskResponse> => {
    const response = await api.post(`/structure/validate?async_task=${asyncTask}`, { structure });
    return response.data;
  },

  chat: async (message: string, structure?: MolecularStructure, asyncTask: boolean = true): Promise<ChatResponse | TaskResponse> => {
    const response = await api.post(`/chat?async_task=${asyncTask}`, { message, structure });
    return response.data;
  },

  generateStructure: async (description: string, asyncTask: boolean = true): Promise<{ structure: MolecularStructure } | TaskResponse> => {
    const response = await api.post(`/structure/generate?async_task=${asyncTask}`, { description });
    return response.data;
  },

  optimizeStructure: async (structure: MolecularStructure, asyncTask: boolean = true): Promise<{ structure: MolecularStructure } | TaskResponse> => {
    const response = await api.post(`/structure/optimize?async_task=${asyncTask}`, { structure });
    return response.data;
  },

  getTaskStatus: async (taskId: string): Promise<TaskProgress> => {
    const response = await api.get(`/tasks/${taskId}`);
    return response.data;
  },

  exportStructure: async (structure: MolecularStructure, format: 'mol' | 'sdf' | 'json'): Promise<Blob> => {
    const response = await api.post(`/structure/export/${format}`, { structure }, {
      responseType: 'blob',
    });
    return response.data;
  },

  combineMolecules: async (reactant1: { smiles: string; structure?: MolecularStructure }, 
                           reactant2: { smiles: string; structure?: MolecularStructure }): Promise<any> => {
    const response = await api.post('/reaction/combine', {
      reactant1: {
        smiles: reactant1.smiles,
        structure: reactant1.structure
      },
      reactant2: {
        smiles: reactant2.smiles,
        structure: reactant2.structure
      }
    });
    return response.data;
  },

  searchPubChem: async (name?: string, cid?: number): Promise<any> => {
    const params = new URLSearchParams();
    if (name) params.append('name', name);
    if (cid) params.append('cid', cid.toString());
    const response = await api.get(`/structure/search/pubchem?${params.toString()}`);
    return response.data;
  },
};

