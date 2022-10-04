Hello! Here are some useful things to know before using FALCON (Fast ALgorithm for Codon OptimizatioN).

FALCON is a multi-objective codon optimization algorithm. It optimizes nucleic acid (NA) sequences to be expressed in HEK or B-cells (Epstein-Barr virally immortalized B cells). Additionally, it supports the optimization for no tissue in particular (in this case, the codon usage tables are according to the entire human transcriptome and not HEK or B-cells). 

---------------------------------------------------------------------
                           FEATURES
---------------------------------------------------------------------

FALCON backtranslates amino acid sequences to NAs, taking into account the following features:
--Single-codon usage in highly expressed genes, also corrected for tRNA levels. 
--Bi-codon usage in highly expressed genes (also known as codon context), corrected for tRNA levels.
--Codon Autocorrelation Bias.
--Correction towards desired GC-content.
--Motif avoidance:
	--11 Restriction sites (BamHI, Bsp1407I, EcoRI, HindIII, KpnI, MscI, NcoI, NheI, SpeI, XbaI, XhoI)
	--Homopolymer stretches (A, T, C or G) >= 6
	--A, T, or AT runs >= 8
	--Pyrimidine stretches >= 10
--Minimum Free Energy (MFE) optimization at the 5' start.
	It creates 10 sub-strings of the first 20 amino acids and calculates their MFE. The sequence with the
	highest MFE is used as a starting string for the creation of the 10 full candidate strings.
--Tournament selection. Among 10 candidate sequences, it selects the one with the highest score based on:
	--Higher Codon Adaptation Index.
	--GC-content closer to the one specified by the user.
	--Lower number of CG dinucleotides.

FALCON is fast! It uses multiprocessing to backtransate your sequences as fast as your machine allows. 

-----------------------------------------------------------------------
                            INPUT
-----------------------------------------------------------------------
FALCON supports input in FASTA format:

>GeneName1
yourAAseq1
>GeneName2
yourAAseq2

Or the input in which every line is an entry. Each entry is composed of Name and AAsequence separated by a tab (note the absence of '>'):

GeneName1 \t yourAAseq1
GeneName2 \t yourAAseq2

Recommendations:
--Make sure that your input sequences do not have other letters than the 20 common amino acids + U (selenocystein). Two examples: 
1) if there are sequences that contain ‘X’ when the amino acid is ambiguous or undetermined, FALCON will crash; 
2) especially with files from BioMart, there are entries that can contain “Sequence unavailable” instead of an actual amino acid sequence. If FALCON does not crash, chances are good the output for those sequences will be of length 3 (i.e. the S will be back translated).
FALCON accepts the use of ‘*’ as a stop codon signal.

--If the input file contains two sequences with the exact same name, FALCON will only back translate the last entry. 
This has to do with how FALCON deals with the input. It saves the sequences with their names in a dictionary as Name:Seq (Key:Value) before processing them. In Python dictionaries, the keys must be unique. Thus, if another sequence with the same name is saved afterwards, it will overwrite the one that was previously stored.

--If your input file is not in FASTA format, make sure that the GeneNames do not start with a ‘>’. 
This is because FALCON uses the ‘>’ to know if it is dealing with a FASTA file. 

-----------------------------------------------------------------------------
                              OUTPUT
-----------------------------------------------------------------------------
FALCON creates an output file in which every line is an entry. Each entry is composed of Name and NAsequence, separated by a tab:

GeneName1 \t yourNAseq1
GeneName2 \t yourNAseq2

-----------------------------------------------------------------------------
                          BEFORE RUNNING FALCON
-----------------------------------------------------------------------------

--Make sure you have installed Python 3.7 or above. 
If you don't, I'd recommend installing it from the Anaconda environment, as it already includes many additional packages that you will not have to worry about installing. Python is automatically installed with the Anaconda distribution. 

link for Anaconda. It's free: https://www.anaconda.com/products/individual

FALCON uses these additional packages that are not included in the standard library: "scipy.optimize" and "seqfold" (JJTimmons (https://pypi.org/project/seqfold/)). 

To install any package, just open the terminal prompt and type (example for seqfold):
pip install seqfold

If the seqfold package is not installed, FALCON can still run. It will automatically disable the MFE optimization at the 5' start.
While it can be somewhat inconvenient having to install an additional package to run FALCON ("scipy.optimize" included in Anaconda btw), the advantage is that you can i) input your desired GC aim and ii) optimize the MFE of the start of your sequences. i) FALCON uses a 4-parameter logistic function to constantly correct the probabilities of the codons to choose, partially depending on the GC-content of the growing NA sequence. The scipy.optimize package allows FALCON to fit the values of the four parameters (A, B, C, D) based on the GC% you want the optimized sequences to have. The more a growing sequence deviates from your desired GC, the stronger the correction. ii) For every AAseq to be backtranslated: 10 candidate sub-strings of the first 20 AAs are generated. The minimum free energy is calculated and the one with the highest is chosen. This 60-nucleotide long sequence is used as a starting point to make the 10 candidate full-strings.
Rationale: The sequence with the highest minimum free energy should be the one forming the least thermodynamically stable secondary structure. This in turn should favor translation initiation.     

-----------------------------------------------------------------------------
                            RUNNING FALCON
-----------------------------------------------------------------------------

--Open the terminal. It is easiest if both FALCON and your input files are in the same folder you opened your terminal in.
--Type the following, separated by a single space: 

python3 FALCON_v1_1.py

--If FALCON is not in the same directory you opened the terminal in, either go to that folder or type the path. Example:

python3 Desktop/my_folder/FALCON_v1_1.py

--FALCON will guide you through the optimization options and will ask you for the name of the input/output files.
--You can monitor the progress in real time.
--After optimizing, FALCON will ask if you want to see the results printed on the screen. If you just have optimized a huge file, I don't recommend printing everything out haha.  

-----------------------------------------------------------------------------
                                CONTACT
-----------------------------------------------------------------------------
If you have comments or issues, feel free to let me know:

Miguel A. Hernandez
contact mail 1: miguel.hernandez@stud.uni-heidelberg.de
contact mail 2: miguel13hh@gmail.com 

