# PSI-GROMACS (AMBER03) topology

A GROMACS-compatible topology and modified AMBER03 force field for 
molecular dynamics simulation of Photosystem I (PSI), based on the 
crystal structure 5OY0. This repository provides all files necessary 
to reproduce the PSI simulation system.

---

## Contents

This repository provides:
- A modified AMBER03 force field accommodating the non-standard 
  residues and cofactors present in PSI
- Individual cofactor topology files for chlorophyll a (CLA), 
  carotenoids, quinones, and lipids
- Protein topology excerpts for the SF4 iron-sulfur clusters and 
  coordinating cysteine residues (ACYS)
- A script for adding cross-residue SF4–ACYS bonded terms that 
  cannot be generated automatically by pdb2gmx
- GROMACS parameter files (.mdp) for energy minimisation, 
  NVT equilibration, NPT equilibration, and NPT production MD

---

## Background

Photosystem I is a trimeric membrane protein complex responsible for 
light-driven electron transfer during photosynthesis. The system 
presents several non-standard topology challenges for molecular 
dynamics simulation:

- Several chlorophyll a residues have incomplete phytol tails in the 
  deposited crystal structure, requiring manual structural repair
- Chlorophyll a and the [4Fe-4S] iron-sulfur clusters lack 
  well-established AMBER-compatible parameters
- Each [4Fe-4S] cluster is covalently coordinated by four cysteine 
  residues (ACYS) across formally separate protein chains, 
  requiring cross-residue bonded terms that pdb2gmx cannot generate 
  automatically

---

## Requirements

- GROMACS (made with 2021, but should be generally compatible with most versions)
- AMBER03 force field (distributed with GROMACS)
- Python 3.x (for the SF4-ACYS bonding script)
- The 5OY0 crystal structure, available from the RCSB Protein Data 
  Bank (https://www.rcsb.org/structure/5OY0); or the OPM reoriented structure.

---

## Usage

### 1. Obtain the crystal structure
Download 5OY0 from the RCSB PDB. Chains A-M are used for simulation; 
chains a-m and 1-0 should be removed. Chain identifier information must 
be stripped from the PDB file prior to topology generation. Additionally,
you will need to remove any non-amino acid residues temporarily (then 
add them back after protein.itp is generated).

### 2. Install the modified force field
Copy the `amber03_PSI.ff` directory into your GROMACS force field 
directory, or into your working directory.

### 3. Generate the protein topology
Run pdb2gmx with the modified force field:

```bash
gmx pdb2gmx -f 5OY0_monomer_nochain.pdb -ff amber03_PSI -water tip3p
```
When prompted, you must specify where the protein chain is not continuous 
(most prompts are not continuous).

### 4. Add cofactor topologies
The cofactor .itp files in the `cofactors/` directory should be 
included in your system topology file (topol.top). Add an include 
statement for each cofactor (CLA, BCR, PQN, etc...); topol.top is an 
example of the overall topology file.

A structure file with protonated cofactors is provided. If you are 
using a different structure, you will need to protonate each cofactor
(IQmol was used for these structures). Additionally, CLA molecules are 
incomplete in the PDB structure. They were rebuilt in PyMol, then WinCoot
was used in an attempt to orient the tails as close to the electron density
(.cif from the PDB website) as possible. This is an approximate starting 
structure, and equillibration will help produce a more proper starting position.

All protonated cofactors should be added back to the full structure. Note that
this structure file has some residue numbers rearranged.

### 5. Add SF4–ACYS cross-residue bonded terms
Run the bonding script to add the inter-residue Fe–S bonds, angles, 
and dihedrals that pdb2gmx cannot generate automatically:

```bash
python scripts/add_SF4_ACYS_bonds.py topol.top
```
This will create protein.itp, which has the needed SF4-ACYS terms.

### 6. Define the simulation box
Define a periodic simulation box around the protein. The `editconf` 
command centres the protein and sets the box dimensions. Typically, 
a dodecahedral box is recommended, but this work used a cubic box for 
familiarity:

```bash
gmx editconf -f system.gro -o system_box.gro -c -d 1.0 -bt cubic
```

The `-d 1.0` flag sets a minimum distance of 1.0 nm between the protein 
and the box edge. Adjust this value if your system requires a larger 
buffer.

### 7. Solvate the system
Add explicit TIP3P water molecules to fill the simulation box:

```bash
gmx solvate -cp system_box.gro -cs spc216.gro -o system_solv.gro -p topol_SF4.top
```

Note that `gmx solvate` will automatically update the `[ molecules ]` 
directive in `protein.top` to reflect the number of water molecules added.

### 8. Add ions
First, generate a .tpr file for the ion addition step using a minimal 
.mdp file:

```bash
gmx grompp -f mdp/ions.mdp -c system_solv.gro -p topol_SF4.top -o ions.tpr
```

Then add ions to neutralise the system charge.

```bash
gmx genion -s ions.tpr -o system_ions.gro -p topol_SF4.top -pname NA -nname CL -neutral
```

See `scripts/README.md` for full usage instructions and a description 
of the terms added.

### 6. Run MD

The `.mdp` files in the `mdp/` directory. Are provided. The order 
they should be run is the following:

```bash
gmx grompp -f mdp/minim.mdp -c system.gro -p topol_SF4.top -o em.tpr
gmx mdrun -v -deffnm em

gmx grompp -f mdp/nvt.mdp -c em.gro -p topol_SF4.top -o nvt.tpr
gmx mdrun -deffnm nvt

gmx grompp -f mdp/npt.mdp -c nvt.gro -p topol_SF4.top -o npt.tpr
gmx mdrun -deffnm npt

gmx grompp -f mdp/md.mdp -c npt.gro -p topol_SF4.top -o md.tpr
gmx mdrun -deffnm md
```

---

## Force Field Modifications

A complete change log of all modifications made to the stock AMBER03 
force field is provided in `forcefield/CHANGES.md`. In summary:

| File | Modification |
|------|-------------|
| ffnonbonded.itp | Added non-bonded parameters for fe1, sf1, and cofactor atom types |
| ffbonded.itp | Added bonded parameters for SF4, ACYS, and cofactor interactions |
| aminoacids.rtp | Added SF4 and ACYS residue definitions |
| posre_ffbonded.itp | Updated to reflect new bonded parameters |
| posre_ffnonbonded.itp | Updated to reflect new non-bonded parameters |

---

## Contact
Slipchenko Research Lab
Purdue University
