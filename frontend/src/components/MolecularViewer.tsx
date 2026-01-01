import { Canvas, useFrame } from '@react-three/fiber';
import { OrbitControls, PerspectiveCamera } from '@react-three/drei';
import { MolecularStructure, getElementColor, getElementRadius, getVanDerWaalsRadius } from '../utils/molecular';
import { useMemo, memo, useRef, useState, useEffect } from 'react';
import * as THREE from 'three';
import { useMolecularStore } from '../store/molecularStore';

interface MolecularViewerProps {
  onAtomClick?: (atomId: number) => void;
}

const AtomSphere = memo(({ 
  atom, 
  onClick, 
  isSelected, 
  representationMode 
}: { 
  atom: any; 
  onClick?: () => void; 
  isSelected: boolean;
  representationMode: string;
}) => {
  const color = getElementColor(atom.element);
  
  // Calculate radius based on representation mode
  let radius: number;
  let opacity = 1;
  
  switch (representationMode) {
    case 'wireframe':
      radius = 0.15;
      break;
    case 'stick':
      // No atoms in stick mode
      return null;
    case 'van-der-waals':
      radius = getVanDerWaalsRadius(atom.element);
      opacity = 0.7;
      break;
    case 'line':
      radius = 0.1;
      break;
    case 'ball-and-stick':
    default:
      radius = getElementRadius(atom.element) * 0.5;
      break;
  }

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
        transparent={opacity < 1}
        opacity={opacity}
      />
    </mesh>
  );
});
AtomSphere.displayName = 'AtomSphere';

const BondCylinder = memo(({ 
  bond, 
  atoms, 
  representationMode 
}: { 
  bond: any; 
  atoms: any[]; 
  representationMode: string;
}) => {
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

  // Calculate radius and color based on representation mode
  let radius: number;
  let color: string;
  
  switch (representationMode) {
    case 'wireframe':
      radius = 0.05;
      color = bond.bond_type === 2 ? '#94a3b8' : bond.bond_type === 3 ? '#64748b' : '#cbd5e1';
      break;
    case 'stick':
      radius = 0.1;
      color = bond.bond_type === 2 ? '#94a3b8' : bond.bond_type === 3 ? '#64748b' : '#cbd5e1';
      break;
    case 'van-der-waals':
      // No bonds in van der Waals mode
      return null;
    case 'line':
      radius = 0.02;
      color = bond.bond_type === 2 ? '#94a3b8' : bond.bond_type === 3 ? '#64748b' : '#cbd5e1';
      break;
    case 'ball-and-stick':
    default:
      radius = 0.1;
      color = bond.bond_type === 2 ? '#94a3b8' : bond.bond_type === 3 ? '#64748b' : '#cbd5e1';
      break;
  }

  return (
    <mesh position={[midpoint.x, midpoint.y, midpoint.z]} quaternion={quaternion}>
      <cylinderGeometry args={[radius, radius, distance, 16]} />
      <meshStandardMaterial color={color} metalness={0.2} roughness={0.6} />
    </mesh>
  );
});
BondCylinder.displayName = 'BondCylinder';

// Animated molecule group for collision
const AnimatedMoleculeGroup = memo(({ 
  structure, 
  initialPos, 
  targetPos, 
  isAnimating,
  onAtomClick,
  selectedAtomId
}: { 
  structure: MolecularStructure; 
  initialPos: { x: number; y: number; z: number };
  targetPos: { x: number; y: number; z: number };
  isAnimating: boolean;
  onAtomClick?: (atomId: number) => void;
  selectedAtomId: number | null;
}) => {
  const groupRef = useRef<THREE.Group>(null);
  const [position, setPosition] = useState(initialPos);
  const [progress, setProgress] = useState(0);

  useEffect(() => {
    if (isAnimating) {
      setProgress(0);
      const startTime = Date.now();
      const duration = 2000; // 2 seconds

      const animate = () => {
        const elapsed = Date.now() - startTime;
        const newProgress = Math.min(elapsed / duration, 1);
        setProgress(newProgress);

        // Ease in-out animation
        const eased = newProgress < 0.5
          ? 2 * newProgress * newProgress
          : 1 - Math.pow(-2 * newProgress + 2, 2) / 2;

        const newPos = {
          x: initialPos.x + (targetPos.x - initialPos.x) * eased,
          y: initialPos.y + (targetPos.y - initialPos.y) * eased,
          z: initialPos.z + (targetPos.z - initialPos.z) * eased,
        };
        setPosition(newPos);

        if (newProgress < 1) {
          requestAnimationFrame(animate);
        }
      };

      animate();
    } else {
      setPosition(initialPos);
      setProgress(0);
    }
  }, [isAnimating, initialPos, targetPos]);

  useFrame(() => {
    if (groupRef.current) {
      groupRef.current.position.set(position.x, position.y, position.z);
    }
  });

  const atomsWithIds = useMemo(() => 
    structure.atoms.map((atom, idx) => ({
      ...atom,
      id: atom.id ?? idx
    })),
    [structure]
  );

  // Calculate center of mass for relative positioning
  const centerOfMass = useMemo(() => {
    if (structure.atoms.length === 0) return { x: 0, y: 0, z: 0 };
    const sum = structure.atoms.reduce((acc, atom) => ({
      x: acc.x + atom.x,
      y: acc.y + atom.y,
      z: acc.z + atom.z
    }), { x: 0, y: 0, z: 0 });
    return {
      x: sum.x / structure.atoms.length,
      y: sum.y / structure.atoms.length,
      z: sum.z / structure.atoms.length
    };
  }, [structure]);

  const relativeAtoms = useMemo(() => 
    atomsWithIds.map(atom => ({
      ...atom,
      x: atom.x - centerOfMass.x,
      y: atom.y - centerOfMass.y,
      z: atom.z - centerOfMass.z
    })),
    [atomsWithIds, centerOfMass]
  );

  return (
    <group ref={groupRef}>
      {relativeAtoms.map((atom) => (
        <AtomSphere
          key={atom.id}
          atom={atom}
          onClick={() => onAtomClick?.(atom.id!)}
          isSelected={selectedAtomId === atom.id}
          representationMode="ball-and-stick"
        />
      ))}
      {structure.bonds.map((bond, idx) => (
        <BondCylinder key={idx} bond={bond} atoms={relativeAtoms} representationMode="ball-and-stick" />
      ))}
      {isAnimating && progress > 0.8 && (
        <mesh position={[0, 0, 0]}>
          <sphereGeometry args={[0.5, 16, 16]} />
          <meshStandardMaterial color="#ffaa00" emissive="#ffaa00" emissiveIntensity={0.5} transparent opacity={0.6} />
        </mesh>
      )}
    </group>
  );
});
AnimatedMoleculeGroup.displayName = 'AnimatedMoleculeGroup';

function Molecule3D({ onAtomClick }: MolecularViewerProps) {
  const { structure, selectedAtomId, reactionCollision, isReactionAnimating, representationMode } = useMolecularStore();
  
  const { atoms, bonds } = useMemo(() => {
    if (!structure) return { atoms: [], bonds: [] };
    
    // Ensure atoms have IDs
    const atomsWithIds = structure.atoms.map((atom, idx) => ({
      ...atom,
      id: atom.id ?? idx
    }));
    
    return { atoms: atomsWithIds, bonds: structure.bonds };
  }, [structure]);

  // Show collision animation if reaction is happening
  if (isReactionAnimating && reactionCollision) {
    const { reactant1_structure, reactant2_structure, reactant1_initial, reactant2_initial, collision_point } = reactionCollision;
    
    return (
      <Canvas>
        <PerspectiveCamera makeDefault position={[10, 10, 10]} />
        <ambientLight intensity={0.5} />
        <directionalLight position={[10, 10, 5]} intensity={1} />
        <pointLight position={[-10, -10, -5]} intensity={0.5} />
        
        {reactant1_structure && (
          <AnimatedMoleculeGroup
            structure={reactant1_structure}
            initialPos={reactant1_initial}
            targetPos={collision_point}
            isAnimating={isReactionAnimating}
            onAtomClick={onAtomClick}
            selectedAtomId={selectedAtomId}
          />
        )}
        
        {reactant2_structure && (
          <AnimatedMoleculeGroup
            structure={reactant2_structure}
            initialPos={reactant2_initial}
            targetPos={collision_point}
            isAnimating={isReactionAnimating}
            onAtomClick={onAtomClick}
            selectedAtomId={selectedAtomId}
          />
        )}
        
        <OrbitControls enableDamping dampingFactor={0.05} />
      </Canvas>
    );
  }

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
          representationMode={representationMode}
        />
      ))}
      
      {bonds.map((bond, idx) => (
        <BondCylinder key={idx} bond={bond} atoms={atoms} representationMode={representationMode} />
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

