# Pairwise Sequence Alignment and Evolutionary Clustering of Mammalian Hemoglobin α Chains

## Team Members:

- Abhigna Nimmagadda (UFID: 31864878)
- Rohith Kumar Ballem (UFID: 30969136)
- Ashfaq Ahmed Mohammed (UFID: 86835927)

## Abstract

Sequence alignment is fundamental to bioinformatics for comparing proteins and understanding evolutionary relationships. This project implements three alignment algorithms in Python: Needleman-Wunsch (global), Smith-Waterman (local), and Semi-global alignment. We apply these algorithms to perform pairwise sequence alignment on hemoglobin α chain sequences from mammalian species, generating alignment scores and similarity measures that enable basic evolutionary clustering analysis. The project compares algorithm performance under different scoring matrices (BLOSUM62, BLOSUM80, PAM250) and gap penalties, validates results against established tools, and uses alignment similarity scores to create simple distance matrices for basic clustering visualization.

## Plan of Action

1. Study and understand the various prediction algorithms by reading relevant papers.

2. Retrieve mammalian hemoglobin α sequences from NCBI Protein database using automated scripts (Biopython/Entrez).

3. Implement three algorithms: Needleman-Wunsch (global), Smith-Waterman (local), and Semi-global alignment with pluggable scoring matrices.

4. Test algorithm accuracy against ClustalW/MUSCLE on a small subset of sequences and analyze computational performance.

5. Generate pairwise similarity matrices and create simple distance-based clustering visualization.


## Workload distribution:

Workload will be distributed evenly amongst members. We plan to meet regularly and work together on all aspects of the project. Different parts of the code will be implemented separately, and we will combine our results. Documentation tasks will be split up, verified and proofread by each team member.


## References:

- Marco-Sola S, et al. (2021). Fast gap-affine pairwise alignment using the wavefront algorithm. Bioinformatics 37(4):456-463

- Günther M, Nettelblad C. (2023). The Long Road to Accurate Pairwise Sequence Alignment. PLoS Comput Biol 19(4):e1010999

- Zhang Z, et al. (2025). QuickEd: high-performance exact sequence alignment based on bound-and-align strategy. Bioinformatics 41(3):btaf112

- Cock PJA, et al. (2022). Biopython: freely available Python tools for computational molecular biology. Bioinformatics 38(23):5813-5815


## Datasets:

- NCBI Protein Database - Mammalian hemoglobin α chains 
https://www.ncbi.nlm.nih.gov/protein/?term=hemoglobin+alpha+AND+Mammalia[Organism]
