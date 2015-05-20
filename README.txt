These functions are free software; you can redistribute it and/or modify them under the terms of the GNU General Public License.

These functions are distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

Maintained by:  Stephen Douglass <smdougla@ucla.edu>
Summary:        A set of functions for discovering novel microRNAs from small RNA sequence data using a Naive Bayes Classifier approach.
Dependencies:   Bioinformatics toolbox
		Bowtie short read aligner


The functions should be executed in the order:
1. readInGenome
2. readInAlignments
3. generateStartSiteCoverage
4. computeShannonEntropy
5. computeStarStrand
6. NBayes
With the exception of computeStarStrand, which can be done any time after step 2 and before step 6.  Output variables from all functions must be saved for use by subsequent functions.

The user must provide the following files:
1. Small RNA sequence read alignments in .bow format
2. Known miRNA sequence read alignments in .bow format
3. Reference genome sequence in FASTA format
The reference sequence (3) must be the same as that used by Bowtie to align the sequences.  The user may wish to have two small RNA sequence read alignments, one for 'readInAlignments' and another for 'generateStartSiteCoverage'.  This is because all possible alignments should be used to generate the genome coverage for determining entropy, but extremely highly repetitive sequences may want to filtered in the alignment step by the user as they are unlikely to be true microRNAs.

Each function is documented more extensively at the top of the corresponding function code. Briefly, 'readInGenome' is needed to convert the FASTA formatted genome file to a MATLAB variable for use by 'readInAlignments', 'generateStartSiteCoverage', and 'computeStarStrand'.  'readInAlignments' converts both the small RNA sequence reads and known miRNA sequence reads to MATLAB variables which store the counts, length, and multiplicity of each sequence.  'generateStartSiteCoverage' outputs a MATLAB variable that is needed to compute the entropy of each sequence read.  'computeShannonEntropy' computes the spread of aligned sequence reads in
the form of Shannon Entropy.  'computeStarStrand' attempts to find a miRNA* sequence corresponding to each sequence read and records whether none or at least one putative miRNA* is found.  'NBayes' takes the data structures defined by the previous functions and computes the Naive Bayesian probability that each mapped location of each sequence (both known microRNAs and small RNA sequence reads) is a microRNA.
