# MDP Files

The following parameter files are provided for each stage of the 
simulation workflow. They should be run in the order listed below.

| File | Stage | Use |
|------|-------|-----|
| ions.mdp | Pre-processing | Minimal parameter file used as input for gmx grompp before running gmx genion. Not used for simulation. |
| em.mdp | Energy minimisation | Steepest descent minimisation until maximum force converges below 1000 kJ mol⁻¹ nm⁻¹. Run on the solvated, ionised system before equilibration. |
| nvt.mdp | Equilibration stage 1 | NVT ensemble at 300 K for 500 ps with position restraints on protein heavy atoms and cofactors. |
| npt.mdp | Equilibration stage 2 | NPT ensemble at 300 K and 1 bar for 500 ps with position restraints on protein heavy atoms and cofactors. |
| md.mdp | Production MD | NPT ensemble at 300 K and 1 bar, 1 ns per run. Run consecutively to achieve desired total simulation time. |

---
