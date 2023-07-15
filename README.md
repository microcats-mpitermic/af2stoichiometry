## "Co-evolution at protein-protein interfaces may guide inference of stoichiometry of oligomeric protein complexes by de novo structure prediction"
### Code and source data
Kilian and Bischofs, 2023 - Max-Planck-Institute for Terrestrial Microbiology

### Requirements
All code was written in Jupyter Notebooks using Python 3.9.13 as part of the [Anaconda](https://www.anaconda.com/) distribution.
To execute the code, the following packages were used:
* matplotlib 3.5.2 - plotting and visualization
* pandas 1.4.4 - handling the tabular data
* numpy 1.21.5 - array processing for plotting
* [evcouplings V0.1.1](https://github.com/debbiemarkslab/EVcouplings) - functions for co-evolution and structure analysis

### Raw data generation
Structure predictions were generated with [Alphafold V2.2.4](https://github.com/deepmind/alphafold) on the HPC "Raven" at the MPCDF, Garching. MSA assembly and co-evolution analysis was carried out by plmDCA with [evcouplings V0.1.1](https://github.com/debbiemarkslab/EVcouplings) as described in the manuscript.

### Data structure
Folders in this repository are sorted by target for structure and co-evolution analysis.
Each folder contains multiple subfolders:
* ChimeraX scripts - .cxc files that can be used to plot co-evolving residue pairs from plmDCA onto PDB structures in [UCSF ChimeraX](https://www.cgl.ucsf.edu/chimerax/)
* Confidence scores - Alphafold2 confidence scores (pTM/ipTM/pLDDT/PAE) in .json file format, derived from the .pkl output of the Alphafold2 pipeline directly
  * ranking_debug.json - weighted total confidence scores for all predictions that were made
  * result_model_X_multimer_v2_pred_Y.json - confidence scores for top-ranking model
* Contact data - filtered data of significantly co-evolving residues mapped per-chain to identify which correspond to native contacts in a PDB structure
  * intra: within a chain (tertiary structure)
  * inter: from a specific chain to any other chain (quaternary structure)
* EVcouplings - raw output from plmDCA on input MSAs, list of all co-evolving residues with associated cn-score and probability (long-range contacts residues at least AAs apart only)
* MSAs - input MSAs for plmDCA built by EVCouplings using jackhmmer
* Notebooks - Jupyter Notebooks for co-evolution and structure analysis, and mapping of co-evolving residues on structures
  * metrics_plotting: plot AF2 confidence scores and scoring of co-evolving residue pairs
  * complex_contact_maps:  co-evolving residue mapping on each individual complex with specific stoichiometry (see below)
* PDB - structure predictions from Alphafold2 or crystal structures with PDB ID used for analysis

### Running the notebooks
1. Import the libraries in the top cell and provide an absolute path under %run to *funcdef.py* (in the top level of this repo).
2. Provide the absolute paths to the PDB folder and then a list of strings of the names to each PDB file to be loaded. Then, provide a list of chain pair tuples for distance map computation. Then, provide a path to the folder containing the co-evolution analysis data and the filenames of the plmDCA *couplingscores_longrange.csv* tables.
3. Execute the following cells to load the data and create the distance maps. This may take some time.
4. For plotting, *plot_contact_map* from the EVCouplings package requires a list of co-evolving residues, which you provide by indexing the dictionary of co-evolving residues (intra_A, intra_B, inter, etc.) with the filename of the *couplingscores_longrange.csv* table. Multimer data is provided similarly. Distance maps are provided by from either a) the monomer_distancemaps dictionary with a key in the format *'filename.pdb/chainID'* or b) the intra_distancemaps dictionary for homomeric distances with a key in the *'filename.pdb/chainID1chainID2'*. If multiple chains are to be analyzed they are pooled by the DistanceMap.aggregate function by providing multiple input distance maps.
5. Distance mapping and writing of the .cxc scripts for visualization in ChimeraX is subsequently carried out on chains provided by the user in the *coupling_scores_compared* function. The mapped co-evolving residues (in .csv format) and .cxc script is finally written to the disk at a specific output path.
