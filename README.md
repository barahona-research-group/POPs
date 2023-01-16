# POPs: Propensity Optimised Paths

This is the official repository of POPs, a method to compute and score paths of optimised propensity that link the orthosteric site with the identified allosteric sites, and identifies crucial residues that contribute to those paths based on results from [bond-to-bond propensity analysis](https://doi.org/10.1038/ncomms12477).

<p align = 'center'>
  <img src = 'doc/graphical_abstract.png' width = '800' />
</p>

Allostery commonly refers to the mechanism that regulates protein activity through the binding of a molecule at a different, usually distal, site from the orthosteric site. The omnipresence of allosteric regulation in nature and its potential for drug design and screening render the study of allostery invaluable. Nevertheless, challenges remain as few computational methods are available to effectively predict allosteric sites, identify signalling pathways involved in allostery, or to aid with the design of suitable molecules targeting such sites. Recently, bond-to-bond propensity analysis has been shown successful at identifying allosteric sites for a large and diverse group of proteins from knowledge of the orthosteric sites and its ligands alone by using network analysis applied to energy-weighted atomistic protein graphs. To address the identification of signalling pathways, we propose here a method to compute and score paths of optimised propensity that link the orthosteric site with the identified allosteric sites, and identifies crucial residues that contribute to those paths.

## Main Work Flow

This code is only available in python. We apologise that it is not available in other languages for users who are not familiar with python and would like to use this code.

### 1. Run bond-to-bond propensity analysis

Bond-to-bond propensity analysis quantifies the non-local effect of instantaneous bond fluctuations propagating through the protein. It is available as a webserver - [ProteinLens](https://proteinlens.io/webserver/). Please follow the [Tutorial](https://proteinlens.io/webserver/tutorial) provided in [ProteinLens](https://proteinlens.io/webserver/) to complete bond-to-bond propensity analysis on your chosen protein and download the results. The results are required for POPs computation and an example of the result folder, [sessionIZN1I](https://github.com/nw97nan/POPs/tree/main/examples/sessionIZNY1I) the bond-to-bond propensity analysis result of h-Ras (PDB ID: 3K8Y), is provided here for illustration.

### 2. Run POPs computation and analysis

After completing bond-to-bond propensity analysis, [POP.py](https://github.com/nw97nan/POPs/blob/main/POP.py) is the only script required for POPs computation and analysis. Please follow the step by step instruction in [POP_example.ipynb](https://github.com/nw97nan/POPs/blob/main/examples/POP_example.ipynb) in the `examples/` directory to complete POPs computation and analysis. Running everything in [POP_example.ipynb](https://github.com/nw97nan/POPs/blob/main/examples/POP_example.ipynb) would give you two folders as the results and they should be the same as the folders given in [examples_completed](https://github.com/nw97nan/POPs/tree/main/examples/examples_completed) in the `examples/` directory. This will make sure everything is set up correctly and you can then proceed with your own protein(s) of interest.

## Contributor

- Nan Wu, GitHub: `nw97nan <https://github.com/nw97nan>`

## Cite

Please cite our paper :

Wu, N., Yaliraki, S. N., & Barahona, M. (2022). Prediction of Protein Allosteric Signalling Pathways and Functional Residues Through Paths of Optimised Propensity. Journal of Molecular Biology, 167749. https://doi.org/https://doi.org/10.1016/j.jmb.2022.167749

Originally appeared as a arXiv preprint:
arXiv:2207.07202; doi: https://doi.org/10.48550/arXiv.2207.07202

and the paper on bond-to-bond propensity analysis:

Amor, B. R. C., Schaub, M. T., Yaliraki, S. N., & Barahona, M. (2016). Prediction of allosteric sites and mediating interactions through bond-to-bond propensities. Nature Communications, 7, 1–13. https://doi.org/10.1038/ncomms12477

if you use this code in your own work.