# CompChem Converters

A collection of computational chemistry file format converters for molecular dynamics and quantum chemistry applications.

## Contents

### MOL to LAMMPS Data File Converter

Converts MOL files (both V2000 and V3000 formats) to LAMMPS data files with proper atom and bond sections.

**Features:**
- Supports both V2000 and V3000 MOL file formats
- Generates LAMMPS data files compatible with OVITO visualization
- Automatic atom type assignment based on element
- Bond type assignment based on bond order
- Proper mass assignment for all elements
- Automatic simulation box calculation with padding
- Batch processing for multiple files

**Usage:**

```bash
# Convert single file with custom output name
python mol-to-lammps_datafile.py input.mol output.data

# Convert single file with automatic output name
python mol-to-lammps_datafile.py input.mol

# Convert all .mol files in a folder
python mol-to-lammps_datafile.py folder_path/
```

**Examples:**

```bash
python mol-to-lammps_datafile.py benzene.mol benzene.data
python mol-to-lammps_datafile.py benzene.mol
python mol-to-lammps_datafile.py ./molecules/
```

## Requirements

- Python 3.6+
- No external dependencies required

## Installation

1. Clone this repository:
   ```bash
   git clone https://github.com/sahmed73/CompChem_Converters.git
   cd CompChem_Converters
   ```

2. Make the scripts executable:
   ```bash
   chmod +x mol-to-lammps_datafile.py
   ```

## Output Format

The LAMMPS data files generated include:
- Header with atom/bond counts and types
- Box dimensions with automatic padding
- Masses section with accurate atomic masses
- Atoms section in "full" style format
- Bonds section (if present in MOL file)

## Supported Elements

The converter includes atomic masses for elements H through Xe. Unknown elements default to carbon mass (12.011 amu).

## File Formats

### Input: MOL Files
- V2000 format (traditional MOL format)
- V3000 format (extended MOL format for large molecules)

### Output: LAMMPS Data Files
- Compatible with LAMMPS molecular dynamics simulator
- OVITO-compatible for visualization
- Full atom style with charges (default: 0.0)

## Contributing

Feel free to contribute by:
- Adding support for new file formats
- Improving existing converters
- Adding unit tests
- Reporting bugs or feature requests

## License

[Add your license information here]

## Author

[Add your name and contact information]

## Changelog

### Version 1.0
- Initial release with MOL to LAMMPS converter
- Support for V2000 and V3000 MOL formats
- Batch processing capability
- OVITO compatibility
