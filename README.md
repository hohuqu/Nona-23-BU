# NonaWorks 2023 - Software x Biology Hackathon

September 2023

## Team
The Team of Nearly Finite Possibilities

## Collaborators
* Quan Huu Ho – MS Biomedical Engineering, Boston University
* Phuong Khanh Tran (Jade) – MS Electrical and Computer Engineering, Boston University
* Payton J Thomas – Ph.D. Bioengineering, UC San Diego
* Taos Transue – Ph.D. Mathematics, University of Utah
* Yehuda Binik – BS Computer Science, Rhode Island College

## Project Title
BioGeneius

## Project Description
Developing a tool for enumerating sequences from a set of oligos submitted to a DNA printing core is a crucial project at the intersection of bioengineering and biosecurity. BioGenius leverages common cloning techniques to predict the potential DNA sequences that could be constructed from provided oligos, offering a valuable resource for ensuring the safety and security of genetic material synthesis. By enabling the rapid identification and evaluation of sequences, this project not only aids in responsible research but also contributes to safeguarding against the synthesis of potentially hazardous genetic material.

## Project Presentation
https://docs.google.com/presentation/d/1wh1UimAxOS5B4HWA8eUSIt-bICXVbQhdXA2GiGQE_eo/edit?usp=sharing

# Documentation

## Restriction-Ligation

Our restriction-ligation algorithm is included in "Twist Problem"/restriction_site.py.

The algorithm simulates the complete digestion of any number of input oligos using BioPython for all known TypeIIP endonucleases. This creates a large set of minimal DNA fragments, each annotated with its 5' and 3' end compatability. We also perform all alternative digests in the case when restriction sites interfere with one another. To reduce time complexity, our algorithm divides endonucleases into equivalence classes of isoschizomers and considers only one representative isoschizomer from each class. 

Our algorithm then traverses the tree of possible assemblies that may be made by ligating these fragments (given their end compatibilities). To make the tree of possible molecules finite, we stop ligating assemblies once they reach 25-30 bp in length. These 25-30bp assemblies may be compared to a database of hazardous DNA molecules (BLAST advises that search sequences be at least 22bp in length).

## Polymerase Cycling Assembly (PCA)

Our Polymerase Cycling Assembly (PCA) algorithm is included in "Twist Problem"/pca_gibson.py.

The algorithm identifies regions of complementarity of one oligo and the reverse complement of another oligo, as in PCA. The algorithm only considers regions of complementarity of at least 20bp, in accordance with conventional cloning wisdom. If sufficient overlaps are identified, the PCA reaction is simulated, and the resulting assembly is saved. All possible PCA assemblies produced from two input oligos are calculated. 

Note that this algorithm may only be applied pair-wise. 

## Gibson Assembly

Our Gibson assembly algorithm is included in "Twist Problem"/pca_gibson.py.

The algorithm finds homologous sequences at the opposite-polarity ends of two oligos of at least 15bp in length (in accordance with the GeneArt Seamless Cloning protocol for Gibson Assembly). Gibson assembly is simulated between these two oligos to form an assembly that is saved. All possible Gibson assemblies are calculated.

Note that this algorithm may only be applied pair-wise. 
