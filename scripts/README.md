# Scripts

This directory contains two utility scripts developed to assist with 
the preparation of a GROMACS topology for Photosystem I (PSI). Both 
scripts address limitations in the standard GROMACS preprocessing 
workflow when dealing with non-standard cofactors and cross-residue 
covalent bonds.

---

## filter_waters.py

Filters solvated water molecules from a GROMACS `.g96` file based on 
their distance to cofactor atoms. This script was developed to address 
a limitation in `gmx solvate`, which does not always respect cofactor atom 
positions when placing solvent molecules, resulting in waters 
being placed unphysically close to cofactors.

Crystal waters present in the original PDB structure should be 
retained regardless of their distance to any cofactor. Only waters 
added by `gmx solvate` are subject to filtering.

### Required Packages

- Python 3.x
- NumPy
- SciPy

# Scripts

This directory contains two utility scripts developed to assist with 
the preparation of a GROMACS topology for Photosystem I (PSI). Both 
scripts address limitations in the standard GROMACS preprocessing 
workflow when dealing with non-standard cofactors and cross-residue 
covalent bonds.

---

## filter_waters.py

Filters solvated water molecules from a GROMACS `.g96` file based on 
their distance to cofactor atoms. This script was developed to address 
a limitation in `gmx solvate`, which does not respect cofactor atom 
positions when placing solvent molecules, often resulting in waters 
being placed unphysically close to cofactors.

Crystal waters present in the original PDB structure are always 
retained regardless of their distance to any cofactor. Only waters 
added by `gmx solvate` are subject to filtering.

### Dependencies

- Python 3.x
- NumPy
- SciPy

Install with:
```bash
pip install numpy scipy
```

### Usage

```bash
python filter_waters.py -i <input.g96> -o <output.g96> --solvent_start <atom_number> [options]
```

### Arguments

| Argument | Required | Default | Description |
|----------|----------|---------|-------------|
| `-i`, `--input` | Yes | — | Path to input `.g96` file |
| `-o`, `--output` | Yes | — | Path to output `.g96` file |
| `--solvent_start` | Yes | — | Atom number at which solvated waters begin. All SOL atoms below this number are treated as crystal waters and retained unconditionally. This value is system-specific — see note below. |
| `--outer_threshold` | No | 1.75 | Maximum distance (nm) from any cofactor atom within which solvated waters are removed |
| `--inner_threshold` | No | 0.0 | Minimum distance (nm) from cofactor atoms for exclusion |
| `--cofactors` | No | ECH 45D EQ3 C7Z CLA PQN BCR QLA LHG LMG SQD LMT SF4 | Space-separated list of cofactor residue names to use for distance calculations |

### Finding your --solvent_start value

Open your `.g96` file and locate the first SOL atom that was added 
by `gmx solvate` that was NOT present in the original crystal 
structure. The atom number on that line is your `--solvent_start` 
value. This is most easily found by checking the atom number of the 
last crystal water in the original PDB and adding 1 
(validate with .g96 or .gro!).

### Notes

- Crystal waters (SOL atoms with atom number below `--solvent_start`) 
  are always retained regardless of distance to cofactors.
- The script removes whole water molecules only. If any atom in a 
  TIP3P water molecule (OW, HW1, HW2) falls within the exclusion 
  distance of a cofactor, the entire molecule is removed.
- This script assumes TIP3P water topology where HW2 is the last 
  atom of each water molecule. If using a different water model the 
  molecule flushing logic may need to be adjusted.
- Distances are in nanometers, typical unit for a GROMACS structure file.

### Example

```bash
python filter_waters.py \
    -i solv_wt.g96 \
    -o solv_wt_filtered.g96 \
    --solvent_start 53901 \
    --outer_threshold 1.75 \
    --cofactors CLA SF4 BCR PQN
```

Expected output:
```
Reading solv_wt.g96...
  142563 atoms read
Filtering waters with cofactors: CLA, SF4, BCR, PQN
  Solvent start atom number: 53901
  Exclusion distance: 0.0–1.75 nm
  47 solvated water molecules removed
  142422 atoms remaining
Writing solv_wt_filtered.g96...
Done.
```

---

## add_SF4_ACYS_bonds.py

Generates and inserts the cross-residue bonded terms required to 
describe the covalent coordination between the SF4 iron-sulfur 
clusters and their coordinating ACYS cysteine residues. These terms 
span the Fe–SG bond and include bonds, angles, and dihedrals that 
`pdb2gmx` cannot generate automatically because they cross residue 
boundaries.

### Dependencies

- Python 3.x

### Before running

Open the script and edit the `COORDINATION` dictionary at the top 
of the file to match your system. This dictionary defines which ACYS 
residue (identified by PDB residue ID) coordinates which iron atom 
in each SF4 cluster. 

```python
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
```

The cluster index (1, 2, 3...) corresponds to the order in which 
SF4 residues appear in the topology file. ACYS residue IDs are 
PDB residue numbers as they appear in the `[ atoms ]` section of 
the topology. PDB residues are relabelled in the structures found on 
this page as such: A3001 -> 1 SF4; C3002 -> 1 SF4; C3003 -> 1 SF4 
(this might be different for your system).

### Usage

```bash
python add_SF4_ACYS_bonds.py -p <input.top> -o <output.top>
```

### Arguments

| Argument | Required | Description |
|----------|----------|-------------|
| `-p`, `--topology` | Yes | Path to input GROMACS topology file |
| `-o`, `--output` | Yes | Path to output topology file |

### What terms are generated

For each Fe–SG coordination bond the following terms are added:

**Bonds**
- One Fe–SG bond per coordination pair (12 total for 3 clusters)

**Angles**
- Fe–SG–CB: the angle from the iron through the sulfur to the 
  beta carbon of the cysteine
- S–Fe–SG: for three S atoms in the SF4 cluster bound to the iron.
- Fe–Fe–SG: for the other 3 iron atoms in the cluster other than the iron itself.
**Dihedrals** (function type 9)
- Fe–SG–CB–X: where X is HB1, HB2, and CA of the ACYS residue
- X–Fe–SG–CB: where X is every cluster atom in the angle list 
- Y–X–Fe–SG: where X is every cluster atom in the angle list and 
  Y is every atom bonded to X within the cluster, excluding the 
  central Fe

### SF4 internal bonding assumed

The script assumes the following SF4 cubane bonding topology, 
where each Fe is NOT bonded to the inorganic S of the same index:

| Iron | Bonded inorganic sulfurs |
|------|--------------------------|
| FE1  | S2, S3, S4 (not S1)      |
| FE2  | S1, S3, S4 (not S2)      |
| FE3  | S1, S2, S4 (not S3)      |
| FE4  | S1, S2, S3 (not S4)      |

If your SF4 topology differs from this, edit the `SF4_FE_TO_S` and 
`SF4_S_TO_FE` dictionaries in the script accordingly.

### Example

```bash
python add_SF4_ACYS_bonds.py \
    -p topol.top \
    -o topol_SF4.top
```

Expected output:
```
Reading topology: topol.top
  34946 atoms parsed
  3 SF4 clusters found in topology

Cluster 1 (SF4 resid 1):
  FE1 (atom 34923) -> ACYS 565 SG (atom 20260)
  FE2 (atom 34924) -> ACYS 583 SG (atom 8999)
  FE3 (atom 34925) -> ACYS 574 SG (atom 8897)
  FE4 (atom 34926) -> ACYS 556 SG (atom 20158)

Cluster 2 (SF4 resid 2):
  FE1 (atom 34931) -> ACYS 51 SG (atom 23686)
  FE2 (atom 34932) -> ACYS 21 SG (atom 23258)
  FE3 (atom 34933) -> ACYS 54 SG (atom 23742)
  FE4 (atom 34934) -> ACYS 48 SG (atom 23653)

Cluster 3 (SF4 resid 3):
  FE1 (atom 34939) -> ACYS 58 SG (atom 23791)
  FE2 (atom 34940) -> ACYS 17 SG (atom 23198)
  FE3 (atom 34941) -> ACYS 14 SG (atom 23157)
  FE4 (atom 34942) -> ACYS 11 SG (atom 23121)

Generated:
  12 bonds
  84 angles
  360 dihedrals

Inserting terms into topology...
Writing output: topol_SF4.top
Done.
```

### Notes

- The script always writes a new output file and never modifies 
  the input topology in place.
- Generated terms are inserted before the next section header 
  after the existing `[ bonds ]`, `[ angles ]`, and `[ dihedrals ]` 
  sections respectively, and are annotated with comments identifying 
  which cluster they belong to.
- If a warning is printed for a missing atom or residue, check that 
  the ACYS residue IDs in the `COORDINATION` dictionary match the 
  residue numbering in your topology file. Note that GROMACS may 
  renumber residues during preprocessing — the residue IDs in the 
  topology `[ atoms ]` section are the ones that matter, not the 
  original PDB residue numbers.
- The expected bond, angle, and dihedral counts for a system with 
  N SF4 clusters are: 4N bonds, 28N angles, and 120N dihedrals. 
  For 3 clusters: 12 bonds, 84 angles, 360 dihedrals. If your 
  counts differ, check the `COORDINATION` dictionary for missing 
  or duplicate entries.

---

## Recommended workflow

These two scripts fit into the overall system preparation workflow 
as follows:

```
gmx pdb2gmx       →  generate initial protein topology
gmx editconf      →  define simulation box
gmx solvate       →  add solvent
filter_waters.py  →  remove waters clashing with cofactors
gmx grompp        →  prepare for ion addition
gmx genion        →  add ions
add_SF4_ACYS_bonds.py  →  add cross-residue SF4-ACYS bonded terms
gmx grompp        →  prepare for energy minimisation
gmx mdrun         →  energy minimisation → NVT → NPT → production MD
```
