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

// Van der Waals radii in Angstroms (converted to relative units)
export const VAN_DER_WAALS_RADII: Record<string, number> = {
  H: 1.2,
  He: 1.4,
  Li: 1.82,
  Be: 1.53,
  B: 1.92,
  C: 1.7,
  N: 1.55,
  O: 1.52,
  F: 1.47,
  Ne: 1.54,
  Na: 2.27,
  Mg: 1.73,
  Al: 1.84,
  Si: 2.1,
  P: 1.8,
  S: 1.8,
  Cl: 1.75,
  Ar: 1.88,
  K: 2.75,
  Ca: 2.31,
  Sc: 2.11,
  Ti: 2.0,
  V: 2.0,
  Cr: 2.0,
  Mn: 2.0,
  Fe: 2.0,
  Co: 2.0,
  Ni: 1.63,
  Cu: 1.4,
  Zn: 1.39,
  Ga: 1.87,
  Ge: 2.11,
  As: 1.85,
  Se: 1.9,
  Br: 1.85,
  Kr: 2.02,
  Rb: 3.03,
  Sr: 2.49,
  Y: 2.0,
  Zr: 2.0,
  Nb: 2.0,
  Mo: 2.0,
  Tc: 2.0,
  Ru: 2.0,
  Rh: 2.0,
  Pd: 1.63,
  Ag: 1.72,
  Cd: 1.58,
  In: 1.93,
  Sn: 2.17,
  Sb: 2.06,
  Te: 2.06,
  I: 1.98,
  Xe: 2.16,
  Cs: 3.43,
  Ba: 2.68,
  La: 2.0,
  Ce: 2.0,
  Pr: 2.0,
  Nd: 2.0,
  Pm: 2.0,
  Sm: 2.0,
  Eu: 2.0,
  Gd: 2.0,
  Tb: 2.0,
  Dy: 2.0,
  Ho: 2.0,
  Er: 2.0,
  Tm: 2.0,
  Yb: 2.0,
  Lu: 2.0,
  Hf: 2.0,
  Ta: 2.0,
  W: 2.0,
  Re: 2.0,
  Os: 2.0,
  Ir: 2.0,
  Pt: 1.75,
  Au: 1.66,
  Hg: 1.55,
  Tl: 1.96,
  Pb: 2.02,
  Bi: 2.07,
  Po: 1.97,
  At: 2.02,
  Rn: 2.2,
  Fr: 3.48,
  Ra: 2.83,
  Ac: 2.0,
  Th: 2.0,
  Pa: 2.0,
  U: 1.86,
  Np: 2.0,
  Pu: 2.0,
  Am: 2.0,
  Cm: 2.0,
  Bk: 2.0,
  Cf: 2.0,
  Es: 2.0,
  Fm: 2.0,
  Md: 2.0,
  No: 2.0,
  Lr: 2.0,
  Rf: 2.0,
  Db: 2.0,
  Sg: 2.0,
  Bh: 2.0,
  Hs: 2.0,
  Mt: 2.0,
  Ds: 2.0,
  Rg: 2.0,
  Cn: 2.0,
  Nh: 2.0,
  Fl: 2.0,
  Mc: 2.0,
  Lv: 2.0,
  Ts: 2.0,
  Og: 2.0,
  default: 1.7
};

export function getVanDerWaalsRadius(element: string): number {
  // Convert Angstroms to relative units (divide by ~2 for scaling)
  return (VAN_DER_WAALS_RADII[element] || VAN_DER_WAALS_RADII.default) / 2;
}

