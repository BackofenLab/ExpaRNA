[![Build Status](https://travis-ci.org/BackofenLab/ExpaRNA.svg?branch=master)](https://travis-ci.org/BackofenLab/ExpaRNA)

# ExpaRNA
Pairwise comparison of RNAs based on exact sequence-structure matches

The program finds the longest common subsequence of exact pattern matches (LCS-EPM).
This is the best co-linear arrangement of substructures common to two RNAs.

The Complexity of the algorithm is O( n^2 m^2 ) time and O (nm) space for two RNAs 
of lengths n and m.


**Motivation**: Specific functions of ribonucleic acid (RNA) molecules are often 
associated with different motifs in the RNA structure. The key feature that forms 
such an RNA motif is the combination of sequence and structure properties. 
ExpaRNA is a new RNA sequence€-structure comparison method 
which maintains exact matching substructures. Existing common substructures are 
treated as whole unit while variability is allowed between such structural motifs.

Based on a fast detectable set of overlapping and crossing substructure matches 
for two nested RNA secondary structures, our method ExpaRNA (exact pattern of 
alignment of RNA) computes the longest co-linear sequence of substructures 
common to two RNAs. Applied to different RNAs, our method correctly identifies 
sequence-€“structure similarities between two RNAs.

**Results**: We have compared ExpaRNA with two other alignment methods that 
work with given RNA structures, namely RNAforester and RNA_align. The results 
are in good agreement, but can be obtained in a fraction of running time, in 
particular for larger RNAs. We have also used ExpaRNA to speed up state-of-the-art 
Sankoff-style alignment tools like LocARNA, and observe a tradeoff between 
quality and speed. However, we get a speedup of 4.25 even in the highest quality 
setting, where the quality of the produced alignment is comparable to that of 
LocARNA alone. 

## Dependencies

- Vienna RNA package (developed with v1.8.5, compilation tested with v2.2.10)


## Contribution

Feel free to contribute to this project by writing [Issues](https://github.com/BackofenLab/ExpaRNA/issues) 
with feature requests or bug reports.

## Cite
If you use IntaRNA, please cite our [article](http://bioinformatics.oxfordjournals.org/content/25/16/2095):
```
doi: 10.1093/bioinformatics/btp065
```
