# GMXFlow

**Protein-ligand molecular dynamics with GROMACS, from a PDB file to publication-ready plots.**

GMXFlow chains the three stages of a typical protein-ligand MD study:

| Stage | Tool | Output |
|---|---|---|
| **1. Set up and run** | `runner/gromacs-runner.sh` | CHARMM36 protein topology, CGenFF ligand parameters, solvated and neutralised box, minimisation, NVT and NPT equilibration, production trajectory |
| **2. Analyse** | `runner/gromacs-analyzer.sh` | Backbone and ligand RMSD, RMSF, radius of gyration, SASA, PCA and free energy landscape (`.xvg` files) |
| **3. Plot** | `GMXPlotter.ipynb` [![Open in Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/yboulaamane/GMXFlow/blob/main/GMXPlotter.ipynb) | Figures for each analysis, plus hydrogen bonds and MM/PBSA per-residue decomposition if you provide them, exported as one ZIP |

The scripts work for any protein (one chain or several) and any ligand: nothing depends on your system's GROMACS group numbers.

## Quick start

Requirements: GROMACS, Open Babel, Perl and Python 3, plus the [third-party files](#third-party-files) below.

```bash
# 1. In a folder with your complex PDB, the runner/ files and the third-party files:
PROTEIN=complex.pdb LIG=LIG MD_NS=100 bash gromacs-runner.sh
#    It stops once: upload lig_fix.mol2 to the CGenFF server, save the result
#    as lig_fix.str in the folder, and run the same command again.

# 2. When the production run has finished:
LIG=LIG bash gromacs-analyzer.sh

# 3. Open the notebook in Colab and point it at your .xvg files
#    (the example data folders show the layout it expects).
```

`lig` stands for your ligand's residue name in lower case. To run the production step on a GPU cluster, adapt `runner/submit_gromacs_gpu.sh` (an example SLURM job).

## Settings

Both scripts take their settings from environment variables:

| Variable | Default | Meaning |
|---|---|---|
| `PROTEIN` | `protein.pdb` | Input complex. `ATOM` records are the protein; other `HETATM` records (crystal waters, ions) are left out |
| `LIG` | `UNK` | Residue name of the ligand in the PDB |
| `FF` | `charmm36-jul2022` | Force field folder, without `.ff` |
| `WATER` | `tip3p` | Water model |
| `BOX`, `DIST` | `cubic`, `1.0` | Box type and solute-to-edge distance (nm) |
| `CONC` | `0.15` | NaCl concentration (mol/L) on top of neutralisation |
| `MD_NS` | *(md.mdp)* | Production length in ns; `md.mdp` as shipped runs 10 ns |
| `DEFFNM` | `md_0_100` | Name of the production run files |
| `GMX` | `gmx` | GROMACS binary, e.g. `gmx_mpi` |
| `GMX_MODULE` | *(none)* | Environment module to load on a cluster |

Groups are chosen by name (`Protein`, `Backbone`, `C-alpha`, the ligand, `Protein_<LIG>`, `SOL`), so the scripts work for any protein, one chain or several, and any ligand name. The analysis reads the trajectory's real end time for the final frame. `submit_gromacs_gpu.sh` is an example SLURM job: change its `#SBATCH` lines to your cluster.

## Example data

The folders at the top level (`rmsd/`, `rmsf/`, `pca/`, `fel/`, `mmpbsa/` and so on) hold example output, so you can open the notebook and see every plot before running anything.

## Third-party files

These are needed but not included in this repository; download them from their authors:

| File | Source |
|---|---|
| `charmm36-jul2022.ff` (unpacked folder) | CHARMM36 force field for GROMACS, MacKerell lab: https://mackerell.umaryland.edu/charmm_ff.shtml#gromacs |
| `cgenff_charmm2gmx_py3_nx2.py` | CGenFF to GROMACS conversion script, same MacKerell lab page (GNU AGPL-3.0) |
| `sort_mol2_bonds.pl` | Justin Lemkul's GROMACS protein-ligand tutorial: http://www.mdtutorials.com/gmx/complex/ |
| `xpm2txt.py` (optional) | Converts the free energy landscape (`FES.xpm`) to text; the analysis script skips that step without it |

The `.mdp` parameter files are adapted from Justin Lemkul's protein-ligand complex tutorial.

