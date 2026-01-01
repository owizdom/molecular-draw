import { useState, useRef, useEffect } from 'react';
import { useMolecularStore } from '../store/molecularStore';

export default function RepresentationSelector() {
  const { representationMode, setRepresentationMode } = useMolecularStore();
  const [isOpen, setIsOpen] = useState(false);
  const dropdownRef = useRef<HTMLDivElement>(null);

  const modes = [
    { value: 'ball-and-stick', label: 'Ball & Stick', icon: '●' },
    { value: 'wireframe', label: 'Wireframe', icon: '□' },
    { value: 'stick', label: 'Stick', icon: '—' },
    { value: 'van-der-waals', label: 'Van der Waals', icon: '○' },
    { value: 'line', label: 'Line', icon: '━' },
  ];

  const currentMode = modes.find(m => m.value === representationMode) || modes[0];

  // Close dropdown when clicking outside
  useEffect(() => {
    const handleClickOutside = (event: MouseEvent) => {
      if (dropdownRef.current && !dropdownRef.current.contains(event.target as Node)) {
        setIsOpen(false);
      }
    };

    if (isOpen) {
      document.addEventListener('mousedown', handleClickOutside);
    }

    return () => {
      document.removeEventListener('mousedown', handleClickOutside);
    };
  }, [isOpen]);

  const handleSelect = (mode: typeof modes[0]) => {
    setRepresentationMode(mode.value as any);
    setIsOpen(false);
  };

  return (
    <div ref={dropdownRef} className="relative">
      <button
        onClick={() => setIsOpen(!isOpen)}
        className="bg-white/90 backdrop-blur-sm border border-slate-200/60 rounded-md shadow-sm px-2.5 py-1.5 flex items-center gap-1.5 text-xs font-medium text-slate-700 hover:bg-slate-50 transition-colors"
      >
        <span className="text-xs">{currentMode.icon}</span>
        <span>{currentMode.label}</span>
        <svg 
          className={`w-3 h-3 transition-transform ${isOpen ? 'rotate-180' : ''}`} 
          fill="none" 
          stroke="currentColor" 
          viewBox="0 0 24 24"
        >
          <path strokeLinecap="round" strokeLinejoin="round" strokeWidth={2} d="M19 9l-7 7-7-7" />
        </svg>
      </button>

      {isOpen && (
        <div className="absolute top-full left-0 mt-1 bg-white/95 backdrop-blur-sm border border-slate-200/60 rounded-md shadow-lg p-1 min-w-[140px] z-50">
          {modes.map((mode) => (
            <button
              key={mode.value}
              onClick={() => handleSelect(mode)}
              className={`w-full flex items-center gap-2 px-2.5 py-1.5 rounded text-xs font-medium transition-colors ${
                representationMode === mode.value
                  ? 'bg-blue-500 text-white'
                  : 'text-slate-600 hover:bg-slate-100'
              }`}
            >
              <span className="text-xs">{mode.icon}</span>
              <span>{mode.label}</span>
            </button>
          ))}
        </div>
      )}
    </div>
  );
}

