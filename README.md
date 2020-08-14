# AlignmentErrorRemoval-Results

## CSV Files ##
Contains the CSV data used in R plots.  
*res.csv* - contains results for general, varying error lengths, and varying percentage of erroneous sequences for 16S.B, Hackett, and small-10-aa. All of the results use 3 k values, and the general data draws parameters from normal distributions. The specific parameters are specified in the Figures section.  
*save_res.csv* - contains same data as res.csv, just used to ensure backup copy  
*variedUnion.csv* - contains results for the union of running the correction algorithm multiple times with varying k values  

## FASTA Data ##
Contains the .fasta files of alignments used in the error simulation process, organized by diameter ranges

## Figures ##
*General Figures* - contains figures of error parameters being drawn from a normal distribution
- Error lengths are drawn from normal distribution centered at 50 with standard deviation 10
- For number of sequences with error, x is drawn from normal distribution centered at 20 with standard deviation 2. The alignment will have N/x erroneous sequences, where N is the total number of sequences in the alignment.  

*Union Figures* - contains figures for testing the union of running the correction algorithm multiple times with varying k values.
- 2 k values: 
  - k = 5 to cover errors of length [0, 30]
  - k = 9 to cover all errors
- 3 k values:
  - k = 5 to cover errors of length [0, 30]
  - k = 9 to cover errors of length [0, 54]
  - k = 17 to cover all errors
- 4 k values:
  - k = 5 to cover errors of length [0, 20]
  - k = 7 to cover errors of length [0, 35]
  - k = 11 to cover errors of length [0, 66]
  - k = 17 to cover all errors  
  
*ErrLen Figures* - contains figures of varying error lengths, with the percentage of erroneous sequences being set to 5%  

*NumErrSeq Figures* - contains figures of varying percentage of erroneous sequences, with the error length being set to 88 characters.

## Script Results ##
Contains output of scripts, which have been processed into the CSV files

## Scripts To Generate FASTA Data Organized By Distance ##
Contains scripts used to generate the .fasta files located in the FASTA Data folder. 

## Scripts To Simulate Errors ##
Contains scripts used to run the error simulation process. The general idea behind most of the files are the same, which is given a directory of .fasta files, run the error simulation process. They differ by changing variables (e.g. error length or percentage of erroneous sequences). There are also util scripts used for calculating error rates, drawing from a normal distribution, etc.

## draw.R ##
The R program that generates the figures.
  
