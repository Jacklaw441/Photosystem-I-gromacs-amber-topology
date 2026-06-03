import numpy as np
from scipy.spatial import cKDTree
import argparse

DEFAULT_COFACTORS = [
    'ECH', '45D', 'EQ3', 'C7Z', 'CLA', 
    'PQN', 'BCR', 'QLA', 'LHG', 'LMG', 
    'SQD', 'LMT', 'SF4'
]

def parse_args():
    parser = argparse.ArgumentParser(
        description="""
        Filter solvated water molecules from a GROMACS .g96 file based on 
        their distance to cofactor atoms. Crystal waters (below --solvent_start 
        atom number) are always retained regardless of distance. Solvated waters 
        within --outer_threshold nm of any cofactor atom are removed.
        """,
        formatter_class=argparse.RawDescriptionHelpFormatter
    )
    parser.add_argument(
        '-i', '--input',
        required=True,
        help='Path to input .g96 file'
    )
    parser.add_argument(
        '-o', '--output',
        required=True,
        help='Path to output .g96 file'
    )
    parser.add_argument(
        '--solvent_start',
        type=int,
        required=True,
        help=('Atom number at which solvated water molecules begin. '
              'All SOL atoms below this number are treated as crystal '
              'waters and retained unconditionally. This value is '
              'system-specific and must be set by the user.')
    )
    parser.add_argument(
        '--outer_threshold',
        type=float,
        default=1.75,
        help='Maximum distance (nm) from cofactor atoms for water exclusion (default: 1.75)'
    )
    parser.add_argument(
        '--inner_threshold',
        type=float,
        default=0.0,
        help='Minimum distance (nm) from cofactor atoms for water exclusion (default: 0.0)'
    )
    parser.add_argument(
        '--cofactors',
        nargs='+',
        default=DEFAULT_COFACTORS,
        help=(f'List of cofactor residue names to use for distance calculations. '
              f'Default: {" ".join(DEFAULT_COFACTORS)}')
    )
    return parser.parse_args()


def parse_g96_file(file_path):
    atoms = []
    with open(file_path, 'r') as file:
        lines = file.readlines()
        for line in lines:
            parts = line.split()
            if len(parts) >= 7 and line[0] != 'G':
                atom_info = {
                    'residue_number': int(parts[0]),
                    'residue': parts[1],
                    'atom_name': parts[2],
                    'atom_number': int(parts[3]),
                    'x': float(parts[4]),
                    'y': float(parts[5]),
                    'z': float(parts[6])
                }
                atoms.append(atom_info)
    return atoms, lines


def filter_water_molecules(atoms, cofactors, solvent_start,
                           inner_threshold=0.0, outer_threshold=1.75):
    """
    Filter solvated water molecules based on distance to cofactor atoms.

    Crystal waters (SOL atoms with atom_number < solvent_start) are always
    retained. Solvated waters within outer_threshold nm of any cofactor atom
    are removed.

    Parameters
    ----------
    atoms : list of dict
        Parsed atom records from the .g96 file
    cofactors : list of str
        Residue names to treat as cofactors
    solvent_start : int
        Atom number at which solvated waters begin
    inner_threshold : float
        Minimum distance from cofactor for exclusion (nm)
    outer_threshold : float
        Maximum distance from cofactor for exclusion (nm)

    Returns
    -------
    kept_atoms : list of dict
        Atoms to write to the output file
    removed_waters : list of dict
        Water molecules that were removed
    """
    cofactors_sol = cofactors + ['SOL']

    protein_atoms = [atom for atom in atoms if atom['residue'] not in cofactors_sol]
    crystal_waters = [atom for atom in atoms if atom['residue'] == 'SOL' 
                      and atom['atom_number'] < solvent_start]
    solvated_waters = [atom for atom in atoms if atom['residue'] == 'SOL' 
                       and atom['atom_number'] >= solvent_start]
    cofactor_atoms = [atom for atom in atoms if atom['residue'] in cofactors]

    if not cofactor_atoms:
        print('Warning: no cofactor atoms found with the specified residue names. '
              'No waters will be removed.')
        return atoms, []

    cofactor_coords = np.array([
        [atom['x'], atom['y'], atom['z']] for atom in cofactor_atoms
    ])
    cofactor_tree = cKDTree(cofactor_coords)

    kept_waters = []
    removed_waters = []
    current_molecule = []
    exclude_molecule = False

    for water in solvated_waters:
        current_molecule.append(water)
        water_coord = np.array([water['x'], water['y'], water['z']])
        distance, _ = cofactor_tree.query(water_coord)

        if inner_threshold <= distance <= outer_threshold:
            exclude_molecule = True

        # HW2 is the last atom in a TIP3P water molecule
        # so we flush the current molecule when we see it
        if water['atom_name'] == 'HW2':
            if exclude_molecule:
                for atom in current_molecule:
                    removed_waters.append(atom)
            else:
                for atom in current_molecule:
                    kept_waters.append(atom)
            current_molecule = []
            exclude_molecule = False

    kept_atoms = protein_atoms + crystal_waters + kept_waters
    return kept_atoms, removed_waters


def write_g96_file(atoms, output_path):
    with open(output_path, 'w') as file:
        for atom_info in atoms:
            new_line = (
                f"{atom_info['residue_number']:>5} "
                f"{atom_info['residue']:<5} "
                f"{atom_info['atom_name']:<5} "
                f"{atom_info['atom_number']:>6} "
                f"{atom_info['x']:>14.9f} "
                f"{atom_info['y']:>14.9f} "
                f"{atom_info['z']:>14.9f}\n"
            )
            file.write(new_line)


def main():
    args = parse_args()

    print(f'Reading {args.input}...')
    atoms, original_lines = parse_g96_file(args.input)
    print(f'  {len(atoms)} atoms read')

    print(f'Filtering waters with cofactors: {", ".join(args.cofactors)}')
    print(f'  Solvent start atom number: {args.solvent_start}')
    print(f'  Exclusion distance: {args.inner_threshold}–{args.outer_threshold} nm')

    kept_atoms, removed_waters = filter_water_molecules(
        atoms,
        cofactors=args.cofactors,
        solvent_start=args.solvent_start,
        inner_threshold=args.inner_threshold,
        outer_threshold=args.outer_threshold
    )

    n_removed = len(removed_waters) // 3
    print(f'  {n_removed} solvated water molecules removed')
    print(f'  {len(kept_atoms)} atoms remaining')

    print(f'Writing {args.output}...')
    write_g96_file(kept_atoms, output_path=args.output)
    print('Done.')


if __name__ == '__main__':
    main()
