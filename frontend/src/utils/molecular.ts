export interface Atom {
  element: string;
  x: number;
  y: number;
  z: number;
  id?: number;
}

export interface Bond {
  atom1_id: number;
  atom2_id: number;
  bond_type: number;
}

export interface MolecularStructure {
  atoms: Atom[];
  bonds: Bond[];
  smiles?: string;
}

export const ELEMENT_COLORS: Record<string, string> = {
  H: '#ffffff',
  C: '#909090',
  N: '#3050f8',
  O: '#ff0d0d',
  F: '#90e050',
  P: '#ff8000',
  S: '#ffff30',
  Cl: '#1ff01f',
  Br: '#a62929',
  I: '#940094',
  default: '#ff1493'
};

export const ELEMENT_RADII: Record<string, number> = {
  H: 0.31,
  C: 0.77,
  N: 0.71,
  O: 0.66,
  F: 0.57,
  P: 1.07,
  S: 1.05,
  Cl: 1.02,
  Br: 1.20,
  I: 1.39,
  default: 0.5
};

export function getElementColor(element: string): string {
  return ELEMENT_COLORS[element] || ELEMENT_COLORS.default;
}

export function getElementRadius(element: string): number {
  return ELEMENT_RADII[element] || ELEMENT_RADII.default;
}

