# TAPER-Results
Scripts and data used for the evaulation of TAPER as well as comparisions to some other methods.

## CSV Files ##
Contains finalized data in .csv and .txt form.  
**res_NewErrRates.csv**: contains our main results on TAPER, evaluating varying error lengths, varying percentage of erroneous sequences, and more generalized testing with erroneous lengths/sequences drawn from normal distributions. Alignments used are 16S.B, Hackett, and small-10-aa.  
**res_other_methods.csv**: contains results of comparing TAPER to other methods (DivA and Divvier). DivA is run on small-10-aa and Divvier is run on both 16S.B and small-10-aa.  
**res_treecmp.csv**: contains results of inferred trees from running TAPER as well as the other methods.  
**res_16SK.csv**: contains our preliminary tests running TAPER with various values of the *k* parameter.  
**variedUnion.csv**: contains our preliminary tests running TAPER with the upper bounded unions of various *k* values.  

## Scripts - Generate FASTA ##  
Contains scripts used to generate .fasta files to get 16S.B alignments of varying diameters.  
**chooseMultipleAlignmentsByDistance.py**: traverses through a tree to find subtrees in a diameter range to produce desired .fasta files.  
**chooseAllAlignments.sh**: script to get .fasta files in all desired diameter ranges.  
**getAlignmentData.sh**: gets the tree diameter of a tree.


## Scripts - Error Simulation ##
Contains scripts for the error simulation process.  
**Divvier Scripts**: additional scripts needed to reformat Divvier outputs to count the error rates.  
**Tree Accuracy Scripts**: scripts to run FastTree and TreeCmp to get tree results.  
**generateErrorModel.py**: adds specified errors to an alignment file.  
**chooseSequences.py**: chooses sequences from an alignment file, resulting in a smaller alignment file output.  
**getErrorRates-filtered.py**: computes FP, FN, TP, TN counts for correction algorithm outputs.  
**generateData Scripts**: pipeline to generate error models, runs TAPER, and gets the error rates.  
**unionK.py**: initial script to compute the upper bounded union of running TAPER on various values of *k*.  


## Test Data ##
Contains compressed alignment files (.fasta.gz) and tree files (.tree) used for evaluation.  
**16S.B Orig**: original .fasta alignments for 16S.B dataset.  
**16S.B Orig Trees**: original .tree files corresponding to the original .fasta alignments for 16S.B.  
**16S.B Error**: error .fasta alignments with varying error lengths, varying percentage of erroneous sequences, and more generalized testing for 16S.B.  
**16S.B Error Trees**: error .tree files corresponding to the error .fasta alignments for 16S.B generalized testing.  
**16S.B Divver**: alignment results of running Divvier on 16S.B generalized testing.  
**16S.B Divvier Trees**: .tree files of Divvier results corresponding to the alignment files.  
**Hackett Orig**: original .fasta alignments for the Hackett dataset.  
**Hackett Orig Trees**: original .tree files corresponding to the original .fasta alignments for Hackett.  
**Hackett Error**: error .fasta alignments with varying error lengths and varying percentage of erroneous sequences for Hackett.  
**Hackett Error Trees**: error .tree files corresponding to the Hackett error .fasta files.  
**small-10-aa Orig**: original .fasta alignment for small-10-aa-RV100-BBA0039.   
**small-10-aa Orig Trees**: original .tree file corresponding to the .fasta alignment for small-10-aa-RV100-BBA0039.  
**small-10-aa Error**: error .fasta alignments with varying error lengths and varying percentage of erroneous sequences for small-10-aa-RV100-BBA0039.  
**small-10-aa Error Trees**: error .tree files corresponding to the error .fasta files for small-10-aa-RV100-BBA0039.  
**small-10-aa Divvier**: alignment results of running Divvier on small-10-aa error alignments with error length 88 and 5% of erroneous sequences.  
**small-10-aa Divvier Trees**: .tree files of Divvier results corresponding to the alignment files.  
**small-10-aa DivA**: alignment results of running Divvier on small-10-aa error alignments with error length 88 and 5% of erroneous sequences.  
**small-10-aa DivA Trees**: .tree files of DivA results corresponding to the alignment files.  


## multi-parameters ##
Data for correction-multi testing on a subset of Hackett and small-10-aa alignments.
