# sorl1_4way_RNASeq_6m
Analysis of a zebrafish model of an early-onset familial Alzheimer's disease (EOfAD) mutation in *sortilin-related receptor 1* (*sorl1*). Analysis was performed on a family of heterozygous EOfAD-like mutants (V1482Afs/+), heterozygous (putatively) null mutants (R122Pfs/+), mutants trans-heterozgous for the the EOfAD-like mutation and the null mutations (V1482Afs/R122Pfs), and their wild type (+/+) siblings in a 'four-way' analysis. Analysing the brains from siblings raised in three tanks, side by side in a recirculating water system with shared water supply reduces environmental and genetic variation, which improves our chances of seeing subtle effects due to *sorl1* genotype. 

- Differential gene expression analysis was performed using `edgeR`
- RUVseq was performed (RUVg) to remove 1 factor of unwanted variation. 
- Enrichment analysis was performed using `fry`, `camera` and `fgsea`. Then the raw p-values were combined to give an overall signifcance value by calculating the harmonic mean p-value. 
- random datasets were also subjected to enrichment testing, to confirm that the observed results due to *sorl1* genotype indeed exceeded chance. 
