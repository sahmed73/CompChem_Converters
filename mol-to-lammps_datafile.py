#!/usr/bin/env python3
"""
MOL to LAMMPS Data File Converter

This script converts a MOL file to LAMMPS data file format with atom and bond sections.
"""

import sys
import os
from collections import defaultdict

class MOLToLAMMPS:
    def __init__(self):
        self.atoms = []
        self.bonds = []
        self.atom_types = {}
        self.bond_types = {}
        self.box_bounds = {'xlo': 0, 'xhi': 50, 'ylo': 0, 'yhi': 50, 'zlo': 0, 'zhi': 50}
        
    def read_mol_file(self, filename):
        """Read and parse MOL file (supports both V2000 and V3000 formats)"""
        try:
            with open(filename, 'r') as f:
                lines = f.readlines()
        except FileNotFoundError:
            print(f"Error: File '{filename}' not found.")
            return False
            
        if len(lines) < 4:
            print("Error: Invalid MOL file format.")
            return False
        
        # Check if it's V3000 format
        is_v3000 = any('V3000' in line for line in lines[:10])
        
        if is_v3000:
            return self._read_v3000_format(lines)
        else:
            return self._read_v2000_format(lines)
    
    def _read_v2000_format(self, lines):
        """Read V2000 format MOL file"""
        # Parse counts line (4th line)
        counts_line = lines[3].strip()
        if len(counts_line) < 6:
            print("Error: Invalid counts line in MOL file.")
            return False
            
        num_atoms = int(counts_line[:3])
        num_bonds = int(counts_line[3:6])
        
        print(f"Reading V2000 format: {num_atoms} atoms and {num_bonds} bonds...")
        
        # Parse atoms (starting from line 4, 0-indexed)
        atom_type_counter = 1
        for i in range(4, 4 + num_atoms):
            if i >= len(lines):
                print(f"Error: Not enough atom lines in MOL file.")
                return False
                
            line = lines[i].strip().split()
            if len(line) < 4:
                print(f"Error: Invalid atom line {i-3}: {lines[i].strip()}")
                return False
                
            x = float(line[0])
            y = float(line[1])
            z = float(line[2])
            element = line[3]
            
            # Assign atom type ID
            if element not in self.atom_types:
                self.atom_types[element] = atom_type_counter
                atom_type_counter += 1
                
            atom_id = len(self.atoms) + 1
            self.atoms.append({
                'id': atom_id,
                'type': self.atom_types[element],
                'element': element,
                'x': x,
                'y': y,
                'z': z,
                'charge': 0.0  # Default charge
            })
            
        # Parse bonds
        bond_type_counter = 1
        for i in range(4 + num_atoms, 4 + num_atoms + num_bonds):
            if i >= len(lines):
                print(f"Error: Not enough bond lines in MOL file.")
                return False
                
            line = lines[i].strip()
            if len(line) < 6:  # Minimum length for bond line
                print(f"Error: Invalid bond line {i-3-num_atoms}: {line}")
                return False
            
            # Parse bond line - MOL V2000 format can have concatenated atom numbers
            try:
                # Try normal split first
                parts = line.split()
                if len(parts) >= 3:
                    # Check if first part is a concatenated atom pair
                    first_part = parts[0]
                    if len(first_part) > 3:  # Likely concatenated atoms like "96100"
                        # For concatenated format, split in the middle
                        mid = len(first_part) // 2
                        atom1 = int(first_part[:mid])
                        atom2 = int(first_part[mid:])
                        bond_order = int(parts[1])
                    else:
                        # Normal format
                        atom1 = int(parts[0])
                        atom2 = int(parts[1])
                        bond_order = int(parts[2])
                else:
                    raise ValueError("Not enough parts in bond line")
                        
            except (ValueError, IndexError) as e:
                print(f"Error: Invalid bond line {i-3-num_atoms}: '{line}' - {e}")
                return False
            
            # Validate atom IDs
            if atom1 < 1 or atom1 > num_atoms or atom2 < 1 or atom2 > num_atoms:
                print(f"Error: Invalid atom IDs in bond line {i-3-num_atoms}: atoms {atom1}-{atom2} (valid range: 1-{num_atoms})")
                return False
            
            # Fix bond order 0 - treat as single bond
            if bond_order == 0:
                bond_order = 1
                print(f"Warning: Bond order 0 found, treating as single bond (atoms {atom1}-{atom2})")
            
            # Assign bond type (based on bond order)
            if bond_order not in self.bond_types:
                self.bond_types[bond_order] = bond_type_counter
                bond_type_counter += 1
                
            bond_id = len(self.bonds) + 1
            self.bonds.append({
                'id': bond_id,
                'type': self.bond_types[bond_order],
                'atom1': atom1,
                'atom2': atom2,
                'order': bond_order
            })
        
        self._calculate_box_bounds()
        return True
    
    def _read_v3000_format(self, lines):
        """Read V3000 format MOL file"""
        print("Reading V3000 format MOL file...")
        
        # Find the COUNTS line and sections
        counts_line = None
        atom_section_start = None
        bond_section_start = None
        atom_section_end = None
        bond_section_end = None
        
        for i, line in enumerate(lines):
            if 'COUNTS' in line:
                counts_line = line
            elif 'BEGIN ATOM' in line:
                atom_section_start = i + 1
            elif 'BEGIN BOND' in line:
                bond_section_start = i + 1
            elif 'END ATOM' in line and atom_section_start:
                atom_section_end = i
            elif 'END BOND' in line and bond_section_start:
                bond_section_end = i
        
        if not counts_line:
            print("Error: Could not find COUNTS line in V3000 format.")
            return False
        
        # Parse counts - format: M  V30 COUNTS 116 115 0 0 0
        parts = counts_line.split()
        counts_idx = -1
        for i, part in enumerate(parts):
            if part == 'COUNTS':
                counts_idx = i
                break
        
        if counts_idx == -1 or len(parts) < counts_idx + 3:
            print(f"Error: Invalid COUNTS line format: {counts_line.strip()}")
            return False
            
        num_atoms = int(parts[counts_idx + 1])
        num_bonds = int(parts[counts_idx + 2])
        
        print(f"Found {num_atoms} atoms and {num_bonds} bonds...")
        
        # Parse atoms
        atom_type_counter = 1
        if atom_section_start and atom_section_end:
            for i in range(atom_section_start, atom_section_end):
                line = lines[i].strip()
                if not line or not line.startswith('M  V30'):
                    continue
                    
                # Parse V3000 atom line: M  V30 1 C -4.825800 -8.138300 0.248500 0
                parts = line.split()
                if len(parts) < 7:  # Need at least: M V30 id element x y z
                    continue
                
                try:
                    atom_id = int(parts[2])
                    element = parts[3]
                    x = float(parts[4])
                    y = float(parts[5])
                    z = float(parts[6])
                except (ValueError, IndexError) as e:
                    print(f"Warning: Could not parse atom line: {line.strip()}")
                    continue
                
                # Assign atom type ID
                if element not in self.atom_types:
                    self.atom_types[element] = atom_type_counter
                    atom_type_counter += 1
                
                self.atoms.append({
                    'id': atom_id,
                    'type': self.atom_types[element],
                    'element': element,
                    'x': x,
                    'y': y,
                    'z': z,
                    'charge': 0.0
                })
        
        # Parse bonds
        bond_type_counter = 1
        if bond_section_start and bond_section_end:
            for i in range(bond_section_start, bond_section_end):
                line = lines[i].strip()
                if not line or not line.startswith('M  V30'):
                    continue
                
                # Parse V3000 bond line: M  V30 1 1 1 2
                parts = line.split()
                if len(parts) < 6:  # Need at least: M V30 id order atom1 atom2
                    continue
                
                try:
                    bond_id = int(parts[2])
                    bond_order = int(parts[3])
                    atom1 = int(parts[4])
                    atom2 = int(parts[5])
                except (ValueError, IndexError) as e:
                    print(f"Warning: Could not parse bond line: {line.strip()}")
                    continue
                
                # Fix bond order 0 - treat as single bond
                if bond_order == 0:
                    bond_order = 1
                    print(f"Warning: Bond order 0 found, treating as single bond (atoms {atom1}-{atom2})")
                
                # Assign bond type
                if bond_order not in self.bond_types:
                    self.bond_types[bond_order] = bond_type_counter
                    bond_type_counter += 1
                
                self.bonds.append({
                    'id': bond_id,
                    'type': self.bond_types[bond_order],
                    'atom1': atom1,
                    'atom2': atom2,
                    'order': bond_order
                })
        
        print(f"Successfully parsed {len(self.atoms)} atoms and {len(self.bonds)} bonds")
        self._calculate_box_bounds()
        return True
    
    def _calculate_box_bounds(self):
        """Calculate simulation box bounds"""
        if self.atoms:
            xs = [atom['x'] for atom in self.atoms]
            ys = [atom['y'] for atom in self.atoms]
            zs = [atom['z'] for atom in self.atoms]
            
            padding = 10.0
            self.box_bounds = {
                'xlo': min(xs) - padding,
                'xhi': max(xs) + padding,
                'ylo': min(ys) - padding,
                'yhi': max(ys) + padding,
                'zlo': min(zs) - padding,
                'zhi': max(zs) + padding
            }
        
    def write_lammps_data(self, filename):
        """Write LAMMPS data file with proper formatting for OVITO compatibility"""
        try:
            with open(filename, 'w') as f:
                # Header with clear description - OVITO compatible
                f.write("LAMMPS data file converted from MOL format\n\n")
                
                # Counts
                f.write(f"{len(self.atoms)} atoms\n")
                if self.bonds:
                    f.write(f"{len(self.bonds)} bonds\n")
                f.write("0 angles\n")
                f.write("0 dihedrals\n")
                f.write("0 impropers\n\n")
                
                # Types
                f.write(f"{len(self.atom_types)} atom types\n")
                if self.bond_types:
                    f.write(f"{len(self.bond_types)} bond types\n")
                f.write("0 angle types\n")
                f.write("0 dihedral types\n")
                f.write("0 improper types\n\n")
                
                # Box bounds - OVITO compatible format
                f.write(f"{self.box_bounds['xlo']:.6f} {self.box_bounds['xhi']:.6f} xlo xhi\n")
                f.write(f"{self.box_bounds['ylo']:.6f} {self.box_bounds['yhi']:.6f} ylo yhi\n")
                f.write(f"{self.box_bounds['zlo']:.6f} {self.box_bounds['zhi']:.6f} zlo zhi\n\n")
                
                # Masses section - required for OVITO
                f.write("Masses\n\n")
                mass_map = {
                    'H': 1.008, 'He': 4.003, 'Li': 6.941, 'Be': 9.012, 'B': 10.811,
                    'C': 12.011, 'N': 14.007, 'O': 15.999, 'F': 18.998, 'Ne': 20.180,
                    'Na': 22.990, 'Mg': 24.305, 'Al': 26.982, 'Si': 28.086, 'P': 30.974,
                    'S': 32.065, 'Cl': 35.453, 'Ar': 39.948, 'K': 39.098, 'Ca': 40.078,
                    'Sc': 44.956, 'Ti': 47.867, 'V': 50.942, 'Cr': 51.996, 'Mn': 54.938,
                    'Fe': 55.845, 'Co': 58.933, 'Ni': 58.693, 'Cu': 63.546, 'Zn': 65.38,
                    'Ga': 69.723, 'Ge': 72.64, 'As': 74.922, 'Se': 78.96, 'Br': 79.904,
                    'Kr': 83.798, 'Rb': 85.468, 'Sr': 87.62, 'Y': 88.906, 'Zr': 91.224,
                    'Nb': 92.906, 'Mo': 95.96, 'Tc': 98.0, 'Ru': 101.07, 'Rh': 102.906,
                    'Pd': 106.42, 'Ag': 107.868, 'Cd': 112.411, 'In': 114.818, 'Sn': 118.71,
                    'Sb': 121.76, 'Te': 127.6, 'I': 126.904, 'Xe': 131.293
                }
                
                for element, type_id in sorted(self.atom_types.items(), key=lambda x: x[1]):
                    mass = mass_map.get(element, 12.011)  # Default to carbon mass if unknown
                    f.write(f"{type_id} {mass:.3f}  # {element}\n")
                f.write("\n")
                
                # Atoms section - OVITO compatible atomic style
                f.write("Atoms  # atomic\n\n")
                for atom in self.atoms:
                    # Format: atom-ID atom-type x y z
                    f.write(f"{atom['id']} {atom['type']} ")
                    f.write(f"{atom['x']:.6f} {atom['y']:.6f} {atom['z']:.6f}\n")
                f.write("\n")
                
                # Bonds section
                if self.bonds:
                    f.write("Bonds\n\n")
                    for bond in self.bonds:
                        f.write(f"{bond['id']} {bond['type']} {bond['atom1']} {bond['atom2']}\n")
                    f.write("\n")
                    
            print(f"LAMMPS data file written to: {filename}")
            
            # Print summary
            print("\nConversion Summary:")
            print(f"Atoms: {len(self.atoms)}")
            print(f"Bonds: {len(self.bonds)}")
            print("Atom types:")
            for element, type_id in sorted(self.atom_types.items(), key=lambda x: x[1]):
                count = sum(1 for atom in self.atoms if atom['element'] == element)
                print(f"  Type {type_id}: {element} ({count} atoms)")
            if self.bond_types:
                print("Bond types:")
                for order, type_id in sorted(self.bond_types.items(), key=lambda x: x[1]):
                    count = sum(1 for bond in self.bonds if bond['order'] == order)
                    bond_name = {1: "single", 2: "double", 3: "triple"}.get(order, f"order-{order}")
                    print(f"  Type {type_id}: {bond_name} bonds ({count} bonds)")
            
            print(f"\nBox dimensions:")
            print(f"  X: {self.box_bounds['xlo']:.3f} to {self.box_bounds['xhi']:.3f}")
            print(f"  Y: {self.box_bounds['ylo']:.3f} to {self.box_bounds['yhi']:.3f}")
            print(f"  Z: {self.box_bounds['zlo']:.3f} to {self.box_bounds['zhi']:.3f}")
                    
        except Exception as e:
            print(f"Error writing LAMMPS file: {e}")
            return False
            
        return True

def convert_single_file(input_file, output_file):
    """Convert a single MOL file to LAMMPS data file"""
    converter = MOLToLAMMPS()
    
    print(f"Converting {input_file} to {output_file}...")
    
    if converter.read_mol_file(input_file):
        if converter.write_lammps_data(output_file):
            print("Conversion completed successfully!")
            return True
        else:
            print("Error: Failed to write LAMMPS data file.")
            return False
    else:
        print("Error: Failed to read MOL file.")
        return False

def convert_folder(folder_path):
    """Convert all MOL files in a folder to LAMMPS data files"""
    if not os.path.exists(folder_path):
        print(f"Error: Folder '{folder_path}' does not exist.")
        return
        
    if not os.path.isdir(folder_path):
        print(f"Error: '{folder_path}' is not a directory.")
        return
    
    mol_files = [f for f in os.listdir(folder_path) if f.lower().endswith('.mol')]
    
    if not mol_files:
        print(f"No .mol files found in folder '{folder_path}'")
        return
        
    print(f"Found {len(mol_files)} MOL files in '{folder_path}'")
    print("Converting files...")
    
    success_count = 0
    for mol_file in mol_files:
        input_path = os.path.join(folder_path, mol_file)
        output_file = mol_file[:-4] + '.data'  # Replace .mol with .data
        output_path = os.path.join(folder_path, output_file)
        
        print(f"\n--- Processing {mol_file} ---")
        if convert_single_file(input_path, output_path):
            success_count += 1
            
    print(f"\n=== Batch Conversion Complete ===")
    print(f"Successfully converted {success_count}/{len(mol_files)} files")

def main():
    """Main function"""
    if len(sys.argv) == 2:
        # Single argument - check if it's a file or folder
        input_path = sys.argv[1]
        
        if os.path.isfile(input_path):
            # Single file mode
            if not input_path.lower().endswith('.mol'):
                print("Error: Input file must have .mol extension")
                return
            output_file = input_path[:-4] + '.data'  # Replace .mol with .data
            convert_single_file(input_path, output_file)
            
        elif os.path.isdir(input_path):
            # Folder mode
            convert_folder(input_path)
            
        else:
            print(f"Error: '{input_path}' is not a valid file or directory.")
            
    elif len(sys.argv) == 3:
        # Two arguments - input and output file
        input_file = sys.argv[1]
        output_file = sys.argv[2]
        
        if not os.path.exists(input_file):
            print(f"Error: Input file '{input_file}' does not exist.")
            return
            
        convert_single_file(input_file, output_file)
        
    else:
        print("Usage:")
        print("  python mol_to_lammps.py input.mol output.data    # Convert single file")
        print("  python mol_to_lammps.py input.mol               # Convert single file (auto output name)")
        print("  python mol_to_lammps.py folder_path/            # Convert all .mol files in folder")
        print("\nExamples:")
        print("  python mol_to_lammps.py benzene.mol benzene.data")
        print("  python mol_to_lammps.py benzene.mol")
        print("  python mol_to_lammps.py ./molecules/")
        return

if __name__ == "__main__":
    main()