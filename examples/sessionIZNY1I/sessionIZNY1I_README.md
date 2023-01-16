# ProteinLens README file
This README lists the contents of the results folder downloaded for session XXXXXX

## Input files
- `sessionXXXXXX_config.json`: contains all input options you chose for the session

- `sessionXXXXXX_graph_atoms.csv`: lists the atoms that were used as nodes in the constructed graph
>The data in this file gives information on the __ID__ of each atom (its internal numbering), __PDB number__ (in modified PDB file), the __atom element__, the __PDB chain__ it is located on, information on its __respective residue name__ and __residue number__ (in the original PDB file), and the __xyz coordinates__ of the atom itself.

- `sessionXXXXXX_graph_bonds.csv`: lists the bonds which were detected between all atoms, i.e. the edges in the constructed graph
>The data in this file gives information on the __ID__ of each bond (its internal numbering), its __bond type__, the __bond weight__ it has been assigned, and the __bond distance__ associated with it (i.e. the length of the bond). Additionally, this file provides information on the two atoms connected by each bond, giving the __atom ID__ (internal numbering) for each atom and the __residue name__, __number__ and __chain__ (in the original PDB file) each atom is a part of.

- `sessionXXXXXX_modified.pdb`: the modified PDB file containing the chosen model, chains and only the atoms which have not been stripped off

- `sessionXXXXXX_original.pdb`: the original PDB file which has been uploaded by the user or downloaded from the PDB

## Results files

### Bond-to-bond propensity
- `sessionXXXXXX_propensity_bonds.csv`: results of the bond-to-bond propensity calculation for all weak bonds present in the graph
>The data in this file organises results across lines corresponding to each bond, identified by __bond ID__ (internal numbering) and __bond name__ (in the format: _atom1 chain atomName : atom2 chain atomName_). The next columns contain the calculated __raw__, __normalised__ and __adjusted propensity__ values. Finally, the last columns give information on the two atoms connected by the bond, the __distance__ between the bond and the source, the __bond weight__ and the calculated __quantile score__ and the quantile score compared to a test set.

- `sessionXXXXXX_propensity_residues.csv`: results of the bond-to-bond propensity calculations for every residue of the biomolecule
>The data in this file organises the propensity results by residue, summing all the propensities of the bonds contained in each residue. Each row contains the __residue number__, __chain__ and __name__. Then, the minimum __distance__ between the source and the residue and the __raw__ and __normalised propensity__ values are provided. The last two columns contain the calculated __quantile score__ and the quantile core compared to a test set.

- `sessionXXXXXX_propensity_scores.csv`: results of the scoring your sites of interest analysis
>The data in this file gives you entry for each scoring query. Each row contains the __query__ number and  __residues__ that were selected as a site of interest. Then, the __target_score__ of the site of interest, the __random_score__ of 1000 randomly sampled sites of the same size and the 95% confidence interval as __random_score_ci__ of this value are given.

### Markov transient times
- `sessionXXXXXX_transient_atoms.csv`: results of the Markov transient calculations for each atom present in the graph
>The data in this file gives information on each __atom ID__ (internal numbering) and __atom name__ (in the format: _residueID chain atomName_). The next columns contain the __transient half time__, the __distance__ of each atom from the source and the calculated __quantile score__.

- `sessionXXXXXX_transient_residues.csv`: results of the Markov transient calculations for every residue of the biomolecule
>The data in this file organises the Markov transients results by residue, averaging the transient half times of all atoms contained in each residue. Each row contains the, __residue number__, __chain__, and __distance__ between the source and the give residue. Then, more information on the residue is given in the column __residue name__, in the format _AminoAcid/residue number/chain_. The next columns contain the __transient half times__ and the calculated __quantile score__.

### Plot folders
- `sessionXXXXXX_propensity_plots`: contains plots similar to the ones in the __Relevant Residues__ section of the bond-bond propensity results (1B) at different cutoffs

 - `sessionXXXXXX_transient_plots`: contains plots similar to the ones in the __Relevant Residues__ section of the Markov transients results (2B) at different cutoffs
