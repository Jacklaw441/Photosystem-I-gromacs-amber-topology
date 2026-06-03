import argparse
from collections import defaultdict

# =============================================================================
# COORDINATION MAPPING
# Edit this dictionary to match your system before running.
#
# Format:
#   cluster_index: { 'FE_name': acys_resid, ... }
#
# cluster_index  — integer, corresponds to the order SF4 residues appear
#                  in the topology (1 = first SF4 residue, etc.)
# FE_name        — iron atom name in the SF4 residue (FE1, FE2, FE3, FE4)
# acys_resid     — PDB residue ID of the coordinating ACYS residue
# =============================================================================

COORDINATION = {
    1: {
        'FE1': 565,
        'FE2': 583,
        'FE3': 574,
        'FE4': 556,
    },
    2: {
        'FE1': 51,
        'FE2': 21,
        'FE3': 54,
        'FE4': 48,
    },
    3: {
        'FE1': 58,
        'FE2': 17,
        'FE3': 14,
        'FE4': 11,
    },
}

# =============================================================================
# SF4 internal bonding
# Each Fe is NOT bonded to the inorganic S of the same index:
#   Fe1 bonds S2, S3, S4  (not S1)
#   Fe2 bonds S1, S3, S4  (not S2)
#   Fe3 bonds S1, S2, S4  (not S3)
#   Fe4 bonds S1, S2, S3  (not S4)
# Do not edit unless your SF4 topology differs from the standard cubane.
# =============================================================================

SF4_FE_TO_S = {
    'FE1': ['S2', 'S3', 'S4'],
    'FE2': ['S1', 'S3', 'S4'],
    'FE3': ['S1', 'S2', 'S4'],
    'FE4': ['S1', 'S2', 'S3'],
}

SF4_S_TO_FE = {
    'S1': ['FE2', 'FE3', 'FE4'],
    'S2': ['FE1', 'FE3', 'FE4'],
    'S3': ['FE1', 'FE2', 'FE4'],
    'S4': ['FE1', 'FE2', 'FE3'],
}

SF4_FE_NAMES = ['FE1', 'FE2', 'FE3', 'FE4']
SF4_S_NAMES  = ['S1',  'S2',  'S3',  'S4' ]


def parse_args():
    parser = argparse.ArgumentParser(
        description="""
        Generate cross-residue SF4-ACYS bonded terms (bonds, angles, dihedrals)
        and insert them into a GROMACS topology file. These terms span the
        Fe-SG coordination bond and cannot be generated automatically by pdb2gmx.

        Edit the COORDINATION dictionary at the top of this script to match
        your system before running.
        """,
        formatter_class=argparse.RawDescriptionHelpFormatter
    )
    parser.add_argument(
        '-p', '--topology',
        required=True,
        help='Path to input GROMACS topology file (.top or .itp)'
    )
    parser.add_argument(
        '-o', '--output',
        required=True,
        help='Path to output topology file'
    )
    return parser.parse_args()


def parse_topology(topology_path):
    """
    Parse a GROMACS topology file and extract atom records.

    Returns
    -------
    lines : list of str
        All lines in the file
    atoms : dict
        { atom_number (int): {
            'atom_number': int,
            'atom_name':   str,
            'resid':       int,
            'resname':     str,
          }
        }
    section_indices : dict
        { section_name (str): int } — line index of each [ section ] header
    """
    atoms = {}
    section_indices = {}
    lines = []

    with open(topology_path, 'r') as f:
        lines = f.readlines()

    in_atoms = False
    for i, line in enumerate(lines):
        stripped = line.strip()

        if stripped.startswith('[') and stripped.endswith(']'):
            section_name = stripped[1:-1].strip().lower()
            section_indices[section_name] = i
            in_atoms = (section_name == 'atoms')
            continue

        if in_atoms and stripped and not stripped.startswith(';'):
            parts = stripped.split()
            if len(parts) >= 6:
                try:
                    atom_number = int(parts[0])
                    atom_type   = parts[1]
                    resid       = int(parts[2])
                    resname     = parts[3]
                    atom_name   = parts[4]
                    atoms[atom_number] = {
                        'atom_number': atom_number,
                        'atom_type':   atom_type,
                        'resid':       resid,
                        'resname':     resname,
                        'atom_name':   atom_name,
                    }
                except ValueError:
                    continue

    return lines, atoms, section_indices


def build_lookup(atoms):
    """
    Build lookup dictionaries for fast atom retrieval.

    Returns
    -------
    by_resid_atomname : dict
        { (resid, atom_name): atom_number }
    sf4_clusters : dict
        { resid (int): { atom_name (str): atom_number (int) } }
    """
    by_resid_atomname = {}
    sf4_clusters      = defaultdict(dict)

    for atom_number, atom in atoms.items():
        by_resid_atomname[(atom['resid'], atom['atom_name'])] = atom_number
        if atom['resname'] == 'SF4':
            sf4_clusters[atom['resid']][atom['atom_name']] = atom_number

    return by_resid_atomname, sf4_clusters


def get_acys_atoms(resid, by_resid_atomname):
    """
    Retrieve atom numbers for key ACYS atoms given a residue ID.

    Returns dict with keys: SG, CB, HB1, HB2, CA
    Raises KeyError if any atom is not found.
    """
    return {
        'SG':  by_resid_atomname[(resid, 'SG')],
        'CB':  by_resid_atomname[(resid, 'CB')],
        'HB1': by_resid_atomname[(resid, 'HB1')],
        'HB2': by_resid_atomname[(resid, 'HB2')],
        'CA':  by_resid_atomname[(resid, 'CA')],
    }


def generate_bonds(fe_atom, sg_atom):
    """Generate Fe-SG bond term."""
    return [(fe_atom, sg_atom, 1)]


def generate_angles(fe_name, fe_atom, sg_atom, cb_atom, cluster_atoms):
    """
    Generate angle terms across the Fe-SG bond.

    Pattern 1 — Fe-SG-CB
    Pattern 2 — X-Fe-SG for every cluster atom X except:
                - Fe itself
                - the inorganic S with the same index as Fe
    """
    angles = []
    excluded_s = 'S' + fe_name[2]  # FE1 -> S1, FE2 -> S2, etc.

    # Pattern 1
    angles.append((fe_atom, sg_atom, cb_atom, 1))

    # Pattern 2
    for atom_name, atom_number in cluster_atoms.items():
        if atom_name == fe_name:
            continue
        if atom_name == excluded_s:
            continue
        angles.append((atom_number, fe_atom, sg_atom, 1))

    return angles


def get_cluster_bonds(atom_name, cluster_atoms):
    """
    Return atom numbers of all atoms bonded to atom_name within the SF4 cluster.
    """
    bonded = []
    if atom_name in SF4_FE_NAMES:
        for other_fe in SF4_FE_NAMES:
            if other_fe != atom_name:
                bonded.append(cluster_atoms[other_fe])
        for s_name in SF4_FE_TO_S[atom_name]:
            bonded.append(cluster_atoms[s_name])
    elif atom_name in SF4_S_NAMES:
        for fe_name in SF4_S_TO_FE[atom_name]:
            bonded.append(cluster_atoms[fe_name])
    return bonded


def generate_dihedrals(fe_name, fe_atom, sg_atom, cb_atom,
                       hb1_atom, hb2_atom, ca_atom, cluster_atoms):
    """
    Generate dihedral terms across the Fe-SG bond.

    Pattern 1 — Fe-SG-CB-X    where X is HB1, HB2, CA
    Pattern 2 — X-Fe-SG-CB    where X is every cluster atom except
                               Fe itself and S(same index as Fe)
    Pattern 3 — Y-X-Fe-SG     where X is from pattern 2's X list,
                               Y is every atom bonded to X in the
                               cluster, excluding Fe itself
    """
    dihedrals  = []
    excluded_s = 'S' + fe_name[2]

    # Build X list (shared between patterns 2 and 3)
    x_list = [
        (name, num)
        for name, num in cluster_atoms.items()
        if name != fe_name and name != excluded_s
    ]

    # Pattern 1
    for x_atom in [hb1_atom, hb2_atom, ca_atom]:
        dihedrals.append((fe_atom, sg_atom, cb_atom, x_atom, 9))

    # Pattern 2
    for x_name, x_atom in x_list:
        dihedrals.append((x_atom, fe_atom, sg_atom, cb_atom, 9))

    # Pattern 3
    for x_name, x_atom in x_list:
        for y_atom in get_cluster_bonds(x_name, cluster_atoms):
            if y_atom == fe_atom:
                continue
            dihedrals.append((y_atom, x_atom, fe_atom, sg_atom, 9))

    return dihedrals


def format_bonds(bonds):
    lines = [
        '; SF4-ACYS bridge bonds\n',
        '; {:>6} {:>6}   func\n'.format('Fe', 'SG')
    ]
    for a1, a2, func in bonds:
        lines.append(f'{a1:>6} {a2:>6}     {func}\n')
    return lines


def format_angles(angles, cluster_idx):
    lines = [f'; SF4-ACYS bridge angles — cluster {cluster_idx}\n']
    for a1, a2, a3, func in angles:
        lines.append(f'{a1:>6} {a2:>6} {a3:>6}     {func}\n')
    return lines


def format_dihedrals(dihedrals, cluster_idx):
    lines = [f'; SF4-ACYS bridge dihedrals — cluster {cluster_idx}\n']
    for a1, a2, a3, a4, func in dihedrals:
        lines.append(f'{a1:>6} {a2:>6} {a3:>6} {a4:>6}     {func}\n')
    return lines


def insert_terms(topology_lines, section_indices,
                 all_bonds, all_angles, all_dihedrals):
    """
    Insert generated terms into the topology after the existing
    [ bonds ], [ angles ], and [ dihedrals ] sections respectively.
    """
    lines  = list(topology_lines)
    offset = 0

    for section_name, new_lines in [
        ('bonds',     format_bonds(all_bonds)),
        ('angles',    all_angles),
        ('dihedrals', all_dihedrals),
    ]:
        if section_name not in section_indices:
            print(f'Warning: [ {section_name} ] section not found. Skipping.')
            continue

        section_start = section_indices[section_name] + offset

        # Find next section header
        insert_at = len(lines)
        for i in range(section_start + 1, len(lines)):
            s = lines[i].strip()
            if s.startswith('[') and s.endswith(']'):
                insert_at = i
                break

        block = ['\n'] + new_lines
        lines[insert_at:insert_at] = block
        offset += len(block)

    return lines


def main():
    args = parse_args()

    print(f'Reading topology: {args.topology}')
    topology_lines, atoms, section_indices = parse_topology(args.topology)
    print(f'  {len(atoms)} atoms parsed')

    _, sf4_clusters = build_lookup(atoms)
    by_resid_atomname, _ = build_lookup(atoms)

    sorted_sf4 = sorted(sf4_clusters.items())
    print(f'  {len(sorted_sf4)} SF4 clusters found in topology')

    if len(sorted_sf4) != len(COORDINATION):
        print(f'Warning: {len(sorted_sf4)} SF4 clusters in topology but '
              f'{len(COORDINATION)} defined in COORDINATION dictionary.')

    all_bonds      = []
    all_angle_lines    = []
    all_dihedral_lines = []

    for cluster_idx, fe_to_acys in sorted(COORDINATION.items()):
        try:
            sf4_resid, cluster_atoms = sorted_sf4[cluster_idx - 1]
        except IndexError:
            print(f'Error: cluster index {cluster_idx} not found in topology.')
            continue

        print(f'\nCluster {cluster_idx} (SF4 resid {sf4_resid}):')

        cluster_bonds     = []
        cluster_angles    = []
        cluster_dihedrals = []

        for fe_name, acys_resid in fe_to_acys.items():
            fe_atom = cluster_atoms.get(fe_name)
            if fe_atom is None:
                print(f'  Warning: {fe_name} not found in cluster {cluster_idx}')
                continue

            try:
                acys = get_acys_atoms(acys_resid, by_resid_atomname)
            except KeyError as e:
                print(f'  Warning: atom {e} not found in ACYS resid {acys_resid}')
                continue

            print(f'  {fe_name} (atom {fe_atom}) -> '
                  f'ACYS {acys_resid} SG (atom {acys["SG"]})')

            cluster_bonds.extend(
                generate_bonds(fe_atom, acys['SG'])
            )
            cluster_angles.extend(
                generate_angles(
                    fe_name, fe_atom, acys['SG'],
                    acys['CB'], cluster_atoms
                )
            )
            cluster_dihedrals.extend(
                generate_dihedrals(
                    fe_name, fe_atom, acys['SG'], acys['CB'],
                    acys['HB1'], acys['HB2'], acys['CA'],
                    cluster_atoms
                )
            )

        all_bonds.extend(cluster_bonds)
        all_angle_lines.extend(format_angles(cluster_angles, cluster_idx))
        all_dihedral_lines.extend(format_dihedrals(cluster_dihedrals, cluster_idx))

    print(f'\nGenerated:')
    print(f'  {len(all_bonds)} bonds')
    print(f'  {len(all_angle_lines) - len(COORDINATION)} angles')
    print(f'  {len(all_dihedral_lines) - len(COORDINATION)} dihedrals')

    print(f'\nInserting terms into topology...')
    new_lines = insert_terms(
        topology_lines, section_indices,
        all_bonds, all_angle_lines, all_dihedral_lines
    )

    print(f'Writing output: {args.output}')
    with open(args.output, 'w') as f:
        f.writelines(new_lines)
    print('Done.')


if __name__ == '__main__':
    main()
