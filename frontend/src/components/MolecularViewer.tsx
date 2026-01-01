import { Canvas } from '@react-three/fiber';
import { OrbitControls, PerspectiveCamera } from '@react-three/drei';
import { MolecularStructure, getElementColor, getElementRadius } from '../utils/molecular';
import { useMemo, memo } from 'react';
import * as THREE from 'three';
import { useMolecularStore } from '../store/molecularStore';

interface MolecularViewerProps {
  onAtomClick?: (atomId: number) => void;
}

const AtomSphere = memo(({ atom, onClick, isSelected }: { atom: any; onClick?: () => void; isSelected: boolean }) => {
  const color = getElementColor(atom.element);
  const radius = getElementRadius(atom.element) * 0.5;

  return (
    <mesh
      position={[atom.x, atom.y, atom.z]}
      onClick={onClick}
      onPointerOver={(e) => {
        e.stopPropagation();
        document.body.style.cursor = 'pointer';
      }}
      onPointerOut={() => {
        document.body.style.cursor = 'default';
      }}
    >
      <sphereGeometry args={[radius, 32, 32]} />
      <meshStandardMaterial
        color={isSelected ? '#6366f1' : color}
        metalness={0.3}
        roughness={0.4}
        emissive={isSelected ? '#6366f1' : '#000000'}
        emissiveIntensity={isSelected ? 0.3 : 0}
      />
    </mesh>
  );
});
AtomSphere.displayName = 'AtomSphere';

const BondCylinder = memo(({ bond, atoms }: { bond: any; atoms: any[] }) => {
  const atom1 = atoms.find(a => a.id === bond.atom1_id);
  const atom2 = atoms.find(a => a.id === bond.atom2_id);
  
  if (!atom1 || !atom2) return null;

  const start = new THREE.Vector3(atom1.x, atom1.y, atom1.z);
  const end = new THREE.Vector3(atom2.x, atom2.y, atom2.z);
  const distance = start.distanceTo(end);
  const midpoint = new THREE.Vector3().addVectors(start, end).multiplyScalar(0.5);
  
  const direction = new THREE.Vector3().subVectors(end, start).normalize();
  const quaternion = new THREE.Quaternion().setFromUnitVectors(
    new THREE.Vector3(0, 1, 0),
    direction
  );

  const radius = 0.1;
  const color = bond.bond_type === 2 ? '#94a3b8' : bond.bond_type === 3 ? '#64748b' : '#cbd5e1';

  return (
    <mesh position={[midpoint.x, midpoint.y, midpoint.z]} quaternion={quaternion}>
      <cylinderGeometry args={[radius, radius, distance, 16]} />
      <meshStandardMaterial color={color} metalness={0.2} roughness={0.6} />
    </mesh>
  );
});
BondCylinder.displayName = 'BondCylinder';

function Molecule3D({ onAtomClick }: MolecularViewerProps) {
  const { structure, selectedAtomId } = useMolecularStore();
  
  const { atoms, bonds } = useMemo(() => {
    if (!structure) return { atoms: [], bonds: [] };
    
    // Ensure atoms have IDs
    const atomsWithIds = structure.atoms.map((atom, idx) => ({
      ...atom,
      id: atom.id ?? idx
    }));
    
    return { atoms: atomsWithIds, bonds: structure.bonds };
  }, [structure]);

  if (!structure || atoms.length === 0) {
    return (
      <div className="flex items-center justify-center h-full text-gray-500">
        No molecule loaded. Use the chat to create one!
      </div>
    );
  }

  return (
    <Canvas>
      <PerspectiveCamera makeDefault position={[10, 10, 10]} />
      <ambientLight intensity={0.5} />
      <directionalLight position={[10, 10, 5]} intensity={1} />
      <pointLight position={[-10, -10, -5]} intensity={0.5} />
      
      {atoms.map((atom) => (
        <AtomSphere
          key={atom.id}
          atom={atom}
          onClick={() => onAtomClick?.(atom.id!)}
          isSelected={selectedAtomId === atom.id}
        />
      ))}
      
      {bonds.map((bond, idx) => (
        <BondCylinder key={idx} bond={bond} atoms={atoms} />
      ))}
      
      <OrbitControls enableDamping dampingFactor={0.05} />
    </Canvas>
  );
}

export default function MolecularViewer(props: MolecularViewerProps) {
  return (
    <div className="w-full h-full bg-gradient-to-br from-slate-50 to-slate-100">
      <Molecule3D onAtomClick={props.onAtomClick} />
    </div>
  );
}

