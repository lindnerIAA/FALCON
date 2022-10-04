#!/anaconda3/bin/python

# This is a script for multi-objective codon optimization of protein sequences.
# Contact:
# Miguel A. Hernandez
# mail 1: miguel.hernandez@stud.uni-heidelberg.de
# mail 2: miguel13hh@gmail.com

import random, re, math

#The MFE is calculated with seqfold package, developed by JJTimmons (https://pypi.org/project/seqfold/).
#Import sub-module of seqfold if available.
try:
    from seqfold import dg
    seq_fold_exists = True
except:
    seq_fold_exists = False

a_line = '-' ; a_space = ' ' #for output aesthetics.

#
#This chunk of code ONLY runs in the MAIN script (not in child parallel processes).
#
if __name__ == '__main__':
    import pathlib, time
    from scipy.optimize import leastsq
    import concurrent.futures
    #
    #-----------------------Dialogue with the user to set desired options------------
    #
    print(f"\n{a_space*30}Hello there! This is FALCON!\n{a_space*30}Let's get right into it!!\n\n")

    #Input file and check if exists in working directory
    InFilename = input('Please, write the name of your INPUT file (e.g. "sequences.txt", "AAseqs.fasta")\n>>> ')
    InFilename2 = pathlib.Path(InFilename)
    while not InFilename2.exists():
        InFilename = input("\nI can't see your file :( Please check spelling or if it's in the working directory.\nName of your input file:\n>>> ")
        InFilename2 = pathlib.Path(InFilename)

    #Enter output file and check if already exists.
    OutFilename = input('\nHow do you want to call your OUTPUT file? (e.g. "out_sequences.txt")\n>>> ')
    OutFilename2 = pathlib.Path(OutFilename)
    if OutFilename2.exists():
        answer = input("\nFile already exists, overwrite? (y/n)\n>>> "); options = ['y', 'n']
        while answer not in options:
            answer = input("\nI know you can do it!! ;) Please choose a valid option:\n(y/n)\n>>> ")
        if answer == 'n':
            while OutFilename2.exists():
                OutFilename = input("\nPlease enter a different name:\n>>> ")
                OutFilename2 = pathlib.Path(OutFilename)
        else:
            #delete the contents of the existing file
            with open(OutFilename, 'w') as f:
                pass

    #Optimization options
    #Expression System
    print(f"\n\n{a_line*10}'Expression System'{a_line*10}")
    ex_sys = input('Four optimization options:\n1) Homo sapiens (no tissue in particular)\n2) Homo sapiens (no tissue in particular, tRNA-corrected)\n3) B-cells (Epstein-Barr virally immortalized)\n4) HEK293T-cells\n\nPlease, enter the number of your option (e.g. "2" )\n>>> ')
    while ex_sys not in ['1', '2', '3', '4']:
        ex_sys = input('\nAnswer not in options: 1, 2, 3 or 4. Please enter a valid number:\n>>> ')
    if ex_sys == '1':
        str_ex_sys = "Homo sapiens (no tissue in particular)"
    elif ex_sys == '2':
        str_ex_sys = "Homo sapiens (no tissue in particular, tRNA-corrected)"
    elif ex_sys == '3':
        str_ex_sys = "B-cells"
    else:
        str_ex_sys = "HEK293T-cells"
    #Desired GC%
    print(f"\n\n{a_line*10}'Desired GC%'{a_line*10}")
    answer = input("\nDo you want to set a different GC% 'aim' for your sequences? Default is 55% (y/n)\n>>> ")
    while answer not in ['y', 'n']:
        answer = input("Answer not in options: y or n. Please enter a valid option.\n>>> ")
    if answer == 'y':
        des_GC = float(input('Alright. Please input desired GC% (e.g. 50, 57.5)\n>>> '))
    else:
        des_GC = 55

    #If seqfold module available, ask if it should be used. If not available, inform user.
    if seq_fold_exists:
        print(f"\n\n{a_line*10}'Minimum Free Energy (MFE)'{a_line*10}")
        answer = input("\nWith the aim to improve translation initiation, you can choose\nto increase the MFE of the first 60 nucleotides of each seq.\nThis favors the formation of less thermodinamically stable secondary structures.\nHowever, it will take around x46 longer.\n\nDo you want to optimize the start of your sequences? (y/n)\n>>> ")
        while answer not in ['y', 'n']:
            answer = input("Answer not in options: y or n. Please enter a valid option.\n>>> ")
        if answer == 'y':
            seq_fold = True
        else:
            seq_fold = False
    else:
        print(f"\n\n{a_line*30}\nThe seqfold module that performs MFE calculations has not\nbeen installed in your computer.\nThe MFE otpimization option is thus disabled\n{a_line*30}\n")
        seq_fold = False
    #Set the MaxThreshold according to desired GC.
    #In this script, there is both a MaxThreshold (60%) and a MinThreshold (48%).
    #They are the tolerable limits of GC% a backtranslated sequence can have. If a sequence is finished and its GC% falls outside
    #of this range, it is deleted and backtranslated again. For every 10 times that the seq is deleted, the respective threshold
    #is relaxed by 0.5% (i.e 60.5% for MaxThreshold).
    #Because the MaxThreshold is 60 by default, if the user sets the desired GC% to 60, the MaxThreshold has to be increased to avoid
    #conflicts (i.e waisting time deleting the backtranslated seq 10 times before relaxing the Threshold to 60.5%).
    if des_GC >= 60:
        MaxThreshold = des_GC+2
    else:
        MaxThreshold = 60

    #To avoid errors when ex_sys 1 or 2 are selected, the variable "CC_evaluation_dict" will be pre-defined (to None).
    #If ex_sys 3 or 4 are selected, the "CC_evaluation_dict" will anyways be updated below to its dictionary.
    CC_evaluation_dict = None

    input(f"\nGreat! Your options were:\nInput file: {InFilename}\nOutput file: {OutFilename}\nOptimize Sequences for: {str_ex_sys}\nDesired GC%: {des_GC}\nMFE optimization: {seq_fold}\n\nPress Enter to start optimizing")

    #---------------------------Single Codon Usage Dictionaries------
    #The weights of the stop codons '*' reflect the average usage in humans (no tissue in particular). Taken from: https://www.genscript.com/tools/codon-frequency-table
    #Single codon usage dictionary for expression systems 1 and 2 (Homo sapiens).
    #It reflects average codon usage in the transcriptome.
    if ex_sys == '1' or ex_sys == '2':
        codons_dict = {
        'L' : [[8, 13, 13, 20, 7, 40], ['TTA', 'TTG', 'CTT', 'CTC', 'CTA', 'CTG']],
        'R' : [[8, 18, 11, 20, 21, 21], ['CGT', 'CGC', 'CGA', 'CGG', 'AGA', 'AGG']],
        'S' : [[19, 22, 15, 5, 15, 24], ['TCT', 'TCC', 'TCA', 'TCG', 'AGT', 'AGC']],
        'A' : [[27, 40, 23, 11], ['GCT', 'GCC', 'GCA', 'GCG']],
        'G' : [[16, 34, 25, 25], ['GGT', 'GGC', 'GGA', 'GGG']],
        'P' : [[29, 32, 28, 11], ['CCT', 'CCC', 'CCA', 'CCG']],
        'T' : [[25, 36, 28, 11], ['ACT', 'ACC', 'ACA', 'ACG']],
        'V' : [[18, 24, 12, 46], ['GTT', 'GTC', 'GTA', 'GTG']],
        'I' : [[36, 47, 17], ['ATT', 'ATC', 'ATA']],
        'C' : [[46, 54], ['TGT', 'TGC']],
        'U' : [[46, 54], ['TGT', 'TGC']],
        'D' : [[46, 54], ['GAT', 'GAC']],
        'E' : [[42, 58], ['GAA', 'GAG']],
        'F' : [[46, 54], ['TTT', 'TTC']],
        'H' : [[42, 58], ['CAT', 'CAC']],
        'N' : [[47, 53], ['AAT', 'AAC']],
        'K' : [[43, 57], ['AAA', 'AAG']],
        'Q' : [[27, 73], ['CAA', 'CAG']],
        'Y' : [[44, 56], ['TAT', 'TAC']],
        'M' : [[100], ['ATG']],
        'W' : [[100], ['TGG']],
        '*' : [[20, 28, 52], ['TAG', 'TAA', 'TGA']]}

    #Single codon usage dictionary for B-cells. It reflects codon usage in highly expressed genes.
    #It is also corrected for tRNA abundance.
    elif ex_sys == '3':
        codons_dict = {
        'F':[ [33,67],["TTT","TTC"] ],
        'L': [ [14,15,8,39,10,14],["CTT","CTC","CTA","CTG","TTA","TTG"] ],
        'I': [ [48,39,13],["ATT","ATC","ATA"] ],
        'M': [ [100],["ATG"] ],
        'V': [ [18,18,16,48],["GTT","GTC","GTA","GTG"] ],
        'S': [ [22,18,15,6,11,28],["TCT","TCC","TCA","TCG","AGT","AGC"] ],
        'P': [ [25,24,40,11],["CCT","CCC","CCA","CCG"] ],
        'T': [ [34,28,26,12],["ACT","ACC","ACA","ACG"] ],
        'A': [ [34,30,22,14],["GCT","GCC","GCA","GCG"] ],
        'Y': [ [32,68],["TAT","TAC"] ],
        'H': [ [30,70],["CAT","CAC"] ],
        'Q': [ [23,77],["CAA","CAG"] ],
        'N': [ [34,66],["AAT","AAC"] ],
        'K': [ [41,59],["AAA","AAG"] ],
        'D': [ [36,64],["GAT","GAC"] ],
        'E': [ [42,58],["GAA","GAG"] ],
        'C': [ [34,66],["TGT","TGC"] ],
        'U': [ [34,66],["TGT","TGC"] ],
        'W': [ [100],["TGG"] ],
        'R': [ [11,15,12,19,23,20],["CGT","CGC","CGA","CGG","AGA","AGG"] ],
        'G': [ [15,40,21,24],["GGT","GGC","GGA","GGG"] ],
        '*': [ [20, 28, 52], ['TAG', 'TAA', 'TGA'] ]}

    #Single codon usage dictionary for HEK293 cells. It reflects codon usage in highly expressed genes.
    #It is also corrected for tRNA abundance.
    else:
        codons_dict = {
        'L' : [[13, 18, 10, 10, 6, 43], ['TTA', 'TTG', 'CTT', 'CTC', 'CTA', 'CTG']],
        'R' : [[16, 11, 19, 22, 15, 17], ['CGT', 'CGC', 'CGA', 'CGG', 'AGA', 'AGG']],
        'S' : [[23, 12, 19, 8, 7, 31], ['TCT', 'TCC', 'TCA', 'TCG', 'AGT', 'AGC']],
        'A' : [[40, 21, 27, 13], ['GCT', 'GCC', 'GCA', 'GCG']],
        'G' : [[10, 46, 22, 22], ['GGT', 'GGC', 'GGA', 'GGG']],
        'P' : [[16, 17, 56, 11], ['CCT', 'CCC', 'CCA', 'CCG']],
        'T' : [[38, 20, 27, 15], ['ACT', 'ACC', 'ACA', 'ACG']],
        'V' : [[13, 12, 10, 65], ['GTT', 'GTC', 'GTA', 'GTG']],
        'I' : [[57, 27, 16], ['ATT', 'ATC', 'ATA']],
        'C' : [[22, 78], ['TGT', 'TGC']],
        'U' : [[22, 78], ['TGT', 'TGC']],
        'D' : [[24, 76], ['GAT', 'GAC']],
        'E' : [[45, 55], ['GAA', 'GAG']],
        'F' : [[21, 79], ['TTT', 'TTC']],
        'H' : [[20, 80], ['CAT', 'CAC']],
        'N' : [[20, 80], ['AAT', 'AAC']],
        'K' : [[49, 51], ['AAA', 'AAG']],
        'Q' : [[24, 76], ['CAA', 'CAG']],
        'Y' : [[21, 79], ['TAT', 'TAC']],
        'M' : [[100], ['ATG']],
        'W' : [[100], ['TGG']],
        '*' : [[20, 28, 52], ['TAG', 'TAA', 'TGA']]}

    #-----Bicodon Usage Dictionaries (Codon Context)-------
    #Global Bicodon Usage for Homo sapiens. It reflects average bicodon usage in the transcriptome.
    if ex_sys == '1':
        CC_dict = {
        'AA' : { 'GCAGCA' : 25, 'GCAGCC' : 37, 'GCAGCG' : 9, 'GCAGCT' : 29, 'GCCGCA' : 11, 'GCCGCC' : 53, 'GCCGCG' : 22, 'GCCGCT' : 14, 'GCGGCA' : 12, 'GCGGCC' : 42, 'GCGGCG' : 31, 'GCGGCT' : 15, 'GCTGCA' : 23, 'GCTGCC' : 38, 'GCTGCG' : 9, 'GCTGCT' : 30} ,
        'AC' : { 'GCATGC' : 41, 'GCATGT' : 59, 'GCCTGC' : 62, 'GCCTGT' : 38, 'GCGTGC' : 62, 'GCGTGT' : 38, 'GCTTGC' : 45, 'GCTTGT' : 55} ,
        'AU' : { 'GCATGC' : 41, 'GCATGT' : 59, 'GCCTGC' : 62, 'GCCTGT' : 38, 'GCGTGC' : 62, 'GCGTGT' : 38, 'GCTTGC' : 45, 'GCTTGT' : 55} ,
        'AD' : { 'GCAGAC' : 45, 'GCAGAT' : 55, 'GCCGAC' : 75, 'GCCGAT' : 25, 'GCGGAC' : 65, 'GCGGAT' : 35, 'GCTGAC' : 52, 'GCTGAT' : 48} ,
        'AE' : { 'GCAGAA' : 45, 'GCAGAG' : 55, 'GCCGAA' : 18, 'GCCGAG' : 82, 'GCGGAA' : 23, 'GCGGAG' : 77, 'GCTGAA' : 42, 'GCTGAG' : 58} ,
        'AF' : { 'GCATTC' : 35, 'GCATTT' : 65, 'GCCTTC' : 65, 'GCCTTT' : 35, 'GCGTTC' : 55, 'GCGTTT' : 45, 'GCTTTC' : 46, 'GCTTTT' : 54} ,
        'AG' : { 'GCAGGA' : 28, 'GCAGGC' : 29, 'GCAGGG' : 27, 'GCAGGT' : 16, 'GCCGGA' : 12, 'GCCGGC' : 45, 'GCCGGG' : 35, 'GCCGGT' : 9, 'GCGGGA' : 14, 'GCGGGC' : 46, 'GCGGGG' : 30, 'GCGGGT' : 10, 'GCTGGA' : 26, 'GCTGGC' : 28, 'GCTGGG' : 28, 'GCTGGT' : 18} ,
        'AH' : { 'GCACAC' : 54, 'GCACAT' : 46, 'GCCCAC' : 66, 'GCCCAT' : 34, 'GCGCAC' : 80, 'GCGCAT' : 20, 'GCTCAC' : 51, 'GCTCAT' : 49} ,
        'AI' : { 'GCAATA' : 26, 'GCAATC' : 28, 'GCAATT' : 46, 'GCCATA' : 11, 'GCCATC' : 60, 'GCCATT' : 29, 'GCGATA' : 17, 'GCGATC' : 45, 'GCGATT' : 37, 'GCTATA' : 21, 'GCTATC' : 31, 'GCTATT' : 48} ,
        'AK' : { 'GCAAAA' : 50, 'GCAAAG' : 50, 'GCCAAA' : 32, 'GCCAAG' : 68, 'GCGAAA' : 34, 'GCGAAG' : 66, 'GCTAAA' : 53, 'GCTAAG' : 47} ,
        'AL' : { 'GCACTA' : 8, 'GCACTC' : 15, 'GCACTG' : 36, 'GCACTT' : 15, 'GCATTA' : 11, 'GCATTG' : 15, 'GCCCTA' : 5, 'GCCCTC' : 18, 'GCCCTG' : 47, 'GCCCTT' : 9, 'GCCTTA' : 5, 'GCCTTG' : 15, 'GCGCTA' : 4, 'GCGCTC' : 22, 'GCGCTG' : 59, 'GCGCTT' : 6, 'GCGTTA' : 2, 'GCGTTG' : 6, 'GCTCTA' : 7, 'GCTCTC' : 15, 'GCTCTG' : 34, 'GCTCTT' : 15, 'GCTTTA' : 10, 'GCTTTG' : 19} ,
        'AM' : { 'GCAATG' : 100, 'GCCATG' : 100, 'GCGATG' : 100, 'GCTATG' : 100} ,
        'AN' : { 'GCAAAC' : 41, 'GCAAAT' : 59, 'GCCAAC' : 64, 'GCCAAT' : 36, 'GCGAAC' : 55, 'GCGAAT' : 45, 'GCTAAC' : 41, 'GCTAAT' : 59} ,
        'AP' : { 'GCACCA' : 27, 'GCACCC' : 33, 'GCACCG' : 9, 'GCACCT' : 31, 'GCCCCA' : 26, 'GCCCCC' : 30, 'GCCCCG' : 19, 'GCCCCT' : 25, 'GCGCCA' : 12, 'GCGCCC' : 47, 'GCGCCG' : 25, 'GCGCCT' : 15, 'GCTCCA' : 30, 'GCTCCC' : 30, 'GCTCCG' : 9, 'GCTCCT' : 31} ,
        'AQ' : { 'GCACAA' : 30, 'GCACAG' : 70, 'GCCCAA' : 17, 'GCCCAG' : 83, 'GCGCAA' : 11, 'GCGCAG' : 89, 'GCTCAA' : 27, 'GCTCAG' : 73} ,
        'AR' : { 'GCAAGA' : 31, 'GCAAGG' : 22, 'GCACGA' : 10, 'GCACGC' : 12, 'GCACGG' : 18, 'GCACGT' : 8, 'GCCAGA' : 15, 'GCCAGG' : 23, 'GCCCGA' : 9, 'GCCCGC' : 21, 'GCCCGG' : 27, 'GCCCGT' : 6, 'GCGAGA' : 7, 'GCGAGG' : 13, 'GCGCGA' : 6, 'GCGCGC' : 35, 'GCGCGG' : 34, 'GCGCGT' : 6, 'GCTAGA' : 15, 'GCTAGG' : 10, 'GCTCGA' : 19, 'GCTCGC' : 19, 'GCTCGG' : 26, 'GCTCGT' : 12} ,
        'AS' : { 'GCAAGC' : 14, 'GCAAGT' : 16, 'GCATCA' : 18, 'GCATCC' : 22, 'GCATCG' : 5, 'GCATCT' : 26, 'GCCAGC' : 28, 'GCCAGT' : 13, 'GCCTCA' : 11, 'GCCTCC' : 24, 'GCCTCG' : 9, 'GCCTCT' : 15, 'GCGAGC' : 19, 'GCGAGT' : 8, 'GCGTCA' : 10, 'GCGTCC' : 33, 'GCGTCG' : 14, 'GCGTCT' : 15, 'GCTAGC' : 8, 'GCTAGT' : 8, 'GCTTCA' : 22, 'GCTTCC' : 28, 'GCTTCG' : 5, 'GCTTCT' : 29} ,
        'AT' : { 'GCAACA' : 36, 'GCAACC' : 25, 'GCAACG' : 8, 'GCAACT' : 31, 'GCCACA' : 25, 'GCCACC' : 41, 'GCCACG' : 15, 'GCCACT' : 20, 'GCGACA' : 22, 'GCGACC' : 37, 'GCGACG' : 21, 'GCGACT' : 20, 'GCTACA' : 33, 'GCTACC' : 30, 'GCTACG' : 7, 'GCTACT' : 30} ,
        'AV' : { 'GCAGTA' : 15, 'GCAGTC' : 19, 'GCAGTG' : 45, 'GCAGTT' : 22, 'GCCGTA' : 5, 'GCCGTC' : 26, 'GCCGTG' : 61, 'GCCGTT' : 8, 'GCGGTA' : 7, 'GCGGTC' : 21, 'GCGGTG' : 62, 'GCGGTT' : 10, 'GCTGTA' : 12, 'GCTGTC' : 21, 'GCTGTG' : 47, 'GCTGTT' : 20} ,
        'AW' : { 'GCATGG' : 100, 'GCCTGG' : 100, 'GCGTGG' : 100, 'GCTTGG' : 100} ,
        'AY' : { 'GCATAC' : 41, 'GCATAT' : 59, 'GCCTAC' : 62, 'GCCTAT' : 38, 'GCGTAC' : 58, 'GCGTAT' : 42, 'GCTTAC' : 43, 'GCTTAT' : 57} ,
        'CA' : { 'TGCGCA' : 12, 'TGCGCC' : 53, 'TGCGCG' : 21, 'TGCGCT' : 14, 'TGTGCA' : 24, 'TGTGCC' : 41, 'TGTGCG' : 7, 'TGTGCT' : 29} ,
        'CC' : { 'TGCTGC' : 63, 'TGCTGT' : 37, 'TGTTGC' : 48, 'TGTTGT' : 52} ,
        'CU' : { 'TGCTGC' : 63, 'TGCTGT' : 37, 'TGTTGC' : 48, 'TGTTGT' : 52} ,
        'CD' : { 'TGCGAC' : 71, 'TGCGAT' : 29, 'TGTGAC' : 54, 'TGTGAT' : 46} ,
        'CE' : { 'TGCGAA' : 20, 'TGCGAG' : 80, 'TGTGAA' : 50, 'TGTGAG' : 50} ,
        'CF' : { 'TGCTTC' : 62, 'TGCTTT' : 38, 'TGTTTC' : 48, 'TGTTTT' : 52} ,
        'CG' : { 'TGCGGA' : 13, 'TGCGGC' : 44, 'TGCGGG' : 33, 'TGCGGT' : 10, 'TGTGGA' : 24, 'TGTGGC' : 31, 'TGTGGG' : 30, 'TGTGGT' : 14} ,
        'CH' : { 'TGCCAC' : 64, 'TGCCAT' : 36, 'TGTCAC' : 52, 'TGTCAT' : 48} ,
        'CI' : { 'TGCATA' : 12, 'TGCATC' : 56, 'TGCATT' : 32, 'TGTATA' : 22, 'TGTATC' : 36, 'TGTATT' : 42} ,
        'CK' : { 'TGCAAA' : 42, 'TGCAAG' : 58, 'TGTAAA' : 54, 'TGTAAG' : 46} ,
        'CL' : { 'TGCCTA' : 5, 'TGCCTC' : 23, 'TGCCTG' : 42, 'TGCCTT' : 13, 'TGCTTA' : 5, 'TGCTTG' : 11, 'TGTCTA' : 8, 'TGTCTC' : 18, 'TGTCTG' : 28, 'TGTCTT' : 17, 'TGTTTA' : 12, 'TGTTTG' : 18} ,
        'CM' : { 'TGCATG' : 100, 'TGTATG' : 100} ,
        'CN' : { 'TGCAAC' : 59, 'TGCAAT' : 41, 'TGTAAC' : 42, 'TGTAAT' : 58} ,
        'CP' : { 'TGCCCA' : 25, 'TGCCCC' : 36, 'TGCCCG' : 15, 'TGCCCT' : 24, 'TGTCCA' : 29, 'TGTCCC' : 32, 'TGTCCG' : 7, 'TGTCCT' : 32} ,
        'CQ' : { 'TGCCAA' : 22, 'TGCCAG' : 78, 'TGTCAA' : 34, 'TGTCAG' : 66} ,
        'CR' : { 'TGCAGA' : 17, 'TGCAGG' : 21, 'TGCCGA' : 9, 'TGCCGC' : 25, 'TGCCGG' : 20, 'TGCCGT' : 7, 'TGTAGA' : 22, 'TGTAGG' : 15, 'TGTCGA' : 15, 'TGTCGC' : 18, 'TGTCGG' : 18, 'TGTCGT' : 11} ,
        'CS' : { 'TGCAGC' : 30, 'TGCAGT' : 16, 'TGCTCA' : 11, 'TGCTCC' : 24, 'TGCTCG' : 5, 'TGCTCT' : 14, 'TGTAGC' : 15, 'TGTAGT' : 15, 'TGTTCA' : 20, 'TGTTCC' : 22, 'TGTTCG' : 4, 'TGTTCT' : 24} ,
        'CT' : { 'TGCACA' : 26, 'TGCACC' : 40, 'TGCACG' : 14, 'TGCACT' : 21, 'TGTACA' : 30, 'TGTACC' : 34, 'TGTACG' : 6, 'TGTACT' : 30} ,
        'CV' : { 'TGCGTA' : 5, 'TGCGTC' : 26, 'TGCGTG' : 60, 'TGCGTT' : 9, 'TGTGTA' : 12, 'TGTGTC' : 24, 'TGTGTG' : 43, 'TGTGTT' : 20} ,
        'CW' : { 'TGCTGG' : 100, 'TGTTGG' : 100} ,
        'CY' : { 'TGCTAC' : 63, 'TGCTAT' : 37, 'TGTTAC' : 50, 'TGTTAT' : 50} ,
        'UA' : { 'TGCGCA' : 12, 'TGCGCC' : 53, 'TGCGCG' : 21, 'TGCGCT' : 14, 'TGTGCA' : 24, 'TGTGCC' : 41, 'TGTGCG' : 7, 'TGTGCT' : 29} ,
        'UC' : { 'TGCTGC' : 63, 'TGCTGT' : 37, 'TGTTGC' : 48, 'TGTTGT' : 52} ,
        'UU' : { 'TGCTGC' : 63, 'TGCTGT' : 37, 'TGTTGC' : 48, 'TGTTGT' : 52} ,
        'UD' : { 'TGCGAC' : 71, 'TGCGAT' : 29, 'TGTGAC' : 54, 'TGTGAT' : 46} ,
        'UE' : { 'TGCGAA' : 20, 'TGCGAG' : 80, 'TGTGAA' : 50, 'TGTGAG' : 50} ,
        'UF' : { 'TGCTTC' : 62, 'TGCTTT' : 38, 'TGTTTC' : 48, 'TGTTTT' : 52} ,
        'UG' : { 'TGCGGA' : 13, 'TGCGGC' : 44, 'TGCGGG' : 33, 'TGCGGT' : 10, 'TGTGGA' : 24, 'TGTGGC' : 31, 'TGTGGG' : 30, 'TGTGGT' : 14} ,
        'UH' : { 'TGCCAC' : 64, 'TGCCAT' : 36, 'TGTCAC' : 52, 'TGTCAT' : 48} ,
        'UI' : { 'TGCATA' : 12, 'TGCATC' : 56, 'TGCATT' : 32, 'TGTATA' : 22, 'TGTATC' : 36, 'TGTATT' : 42} ,
        'UK' : { 'TGCAAA' : 42, 'TGCAAG' : 58, 'TGTAAA' : 54, 'TGTAAG' : 46} ,
        'UL' : { 'TGCCTA' : 5, 'TGCCTC' : 23, 'TGCCTG' : 42, 'TGCCTT' : 13, 'TGCTTA' : 5, 'TGCTTG' : 11, 'TGTCTA' : 8, 'TGTCTC' : 18, 'TGTCTG' : 28, 'TGTCTT' : 17, 'TGTTTA' : 12, 'TGTTTG' : 18} ,
        'UM' : { 'TGCATG' : 100, 'TGTATG' : 100} ,
        'UN' : { 'TGCAAC' : 59, 'TGCAAT' : 41, 'TGTAAC' : 42, 'TGTAAT' : 58} ,
        'UP' : { 'TGCCCA' : 25, 'TGCCCC' : 36, 'TGCCCG' : 15, 'TGCCCT' : 24, 'TGTCCA' : 29, 'TGTCCC' : 32, 'TGTCCG' : 7, 'TGTCCT' : 32} ,
        'UQ' : { 'TGCCAA' : 22, 'TGCCAG' : 78, 'TGTCAA' : 34, 'TGTCAG' : 66} ,
        'UR' : { 'TGCAGA' : 17, 'TGCAGG' : 21, 'TGCCGA' : 9, 'TGCCGC' : 25, 'TGCCGG' : 20, 'TGCCGT' : 7, 'TGTAGA' : 22, 'TGTAGG' : 15, 'TGTCGA' : 15, 'TGTCGC' : 18, 'TGTCGG' : 18, 'TGTCGT' : 11} ,
        'US' : { 'TGCAGC' : 30, 'TGCAGT' : 16, 'TGCTCA' : 11, 'TGCTCC' : 24, 'TGCTCG' : 5, 'TGCTCT' : 14, 'TGTAGC' : 15, 'TGTAGT' : 15, 'TGTTCA' : 20, 'TGTTCC' : 22, 'TGTTCG' : 4, 'TGTTCT' : 24} ,
        'UT' : { 'TGCACA' : 26, 'TGCACC' : 40, 'TGCACG' : 14, 'TGCACT' : 21, 'TGTACA' : 30, 'TGTACC' : 34, 'TGTACG' : 6, 'TGTACT' : 30} ,
        'UV' : { 'TGCGTA' : 5, 'TGCGTC' : 26, 'TGCGTG' : 60, 'TGCGTT' : 9, 'TGTGTA' : 12, 'TGTGTC' : 24, 'TGTGTG' : 43, 'TGTGTT' : 20} ,
        'UW' : { 'TGCTGG' : 100, 'TGTTGG' : 100} ,
        'UY' : { 'TGCTAC' : 63, 'TGCTAT' : 37, 'TGTTAC' : 50, 'TGTTAT' : 50} ,
        'DA' : { 'GACGCA' : 15, 'GACGCC' : 47, 'GACGCG' : 20, 'GACGCT' : 17, 'GATGCA' : 25, 'GATGCC' : 37, 'GATGCG' : 6, 'GATGCT' : 32} ,
        'DC' : { 'GACTGC' : 55, 'GACTGT' : 45, 'GATTGC' : 42, 'GATTGT' : 58} ,
        'DU' : { 'GACTGC' : 55, 'GACTGT' : 45, 'GATTGC' : 42, 'GATTGT' : 58} ,
        'DD' : { 'GACGAC' : 68, 'GACGAT' : 32, 'GATGAC' : 48, 'GATGAT' : 52} ,
        'DE' : { 'GACGAA' : 23, 'GACGAG' : 77, 'GATGAA' : 52, 'GATGAG' : 48} ,
        'DF' : { 'GACTTC' : 59, 'GACTTT' : 41, 'GATTTC' : 44, 'GATTTT' : 56} ,
        'DG' : { 'GACGGA' : 13, 'GACGGC' : 46, 'GACGGG' : 30, 'GACGGT' : 11, 'GATGGA' : 28, 'GATGGC' : 29, 'GATGGG' : 25, 'GATGGT' : 18} ,
        'DH' : { 'GACCAC' : 64, 'GACCAT' : 36, 'GATCAC' : 47, 'GATCAT' : 53} ,
        'DI' : { 'GACATA' : 13, 'GACATC' : 53, 'GACATT' : 33, 'GATATA' : 21, 'GATATC' : 31, 'GATATT' : 47} ,
        'DK' : { 'GACAAA' : 40, 'GACAAG' : 60, 'GATAAA' : 60, 'GATAAG' : 40} ,
        'DL' : { 'GACCTA' : 6, 'GACCTC' : 21, 'GACCTG' : 41, 'GACCTT' : 11, 'GACTTA' : 7, 'GACTTG' : 14, 'GATCTA' : 8, 'GATCTC' : 16, 'GATCTG' : 25, 'GATCTT' : 18, 'GATTTA' : 13, 'GATTTG' : 20} ,
        'DM' : { 'GACATG' : 100, 'GATATG' : 100} ,
        'DN' : { 'GACAAC' : 59, 'GACAAT' : 41, 'GATAAC' : 41, 'GATAAT' : 59} ,
        'DP' : { 'GACCCA' : 24, 'GACCCC' : 37, 'GACCCG' : 12, 'GACCCT' : 27, 'GATCCA' : 29, 'GATCCC' : 28, 'GATCCG' : 6, 'GATCCT' : 37} ,
        'DQ' : { 'GACCAA' : 27, 'GACCAG' : 73, 'GATCAA' : 40, 'GATCAG' : 60} ,
        'DR' : { 'GACAGA' : 23, 'GACAGG' : 21, 'GACCGA' : 10, 'GACCGC' : 20, 'GACCGG' : 19, 'GACCGT' : 7, 'GATAGA' : 26, 'GATAGG' : 13, 'GATCGA' : 17, 'GATCGC' : 15, 'GATCGG' : 16, 'GATCGT' : 13} ,
        'DS' : { 'GACAGC' : 30, 'GACAGT' : 18, 'GACTCA' : 13, 'GACTCC' : 19, 'GACTCG' : 7, 'GACTCT' : 14, 'GATAGC' : 12, 'GATAGT' : 13, 'GATTCA' : 21, 'GATTCC' : 21, 'GATTCG' : 4, 'GATTCT' : 29} ,
        'DT' : { 'GACACA' : 28, 'GACACC' : 36, 'GACACG' : 14, 'GACACT' : 22, 'GATACA' : 33, 'GATACC' : 29, 'GATACG' : 7, 'GATACT' : 31} ,
        'DV' : { 'GACGTA' : 6, 'GACGTC' : 26, 'GACGTG' : 58, 'GACGTT' : 10, 'GATGTA' : 14, 'GATGTC' : 23, 'GATGTG' : 40, 'GATGTT' : 23} ,
        'DW' : { 'GACTGG' : 100, 'GATTGG' : 100} ,
        'DY' : { 'GACTAC' : 59, 'GACTAT' : 41, 'GATTAC' : 42, 'GATTAT' : 58} ,
        'EA' : { 'GAAGCA' : 29, 'GAAGCC' : 34, 'GAAGCG' : 6, 'GAAGCT' : 31, 'GAGGCA' : 20, 'GAGGCC' : 43, 'GAGGCG' : 14, 'GAGGCT' : 23} ,
        'EC' : { 'GAATGC' : 30, 'GAATGT' : 70, 'GAGTGC' : 50, 'GAGTGT' : 50} ,
        'EU' : { 'GAATGC' : 30, 'GAATGT' : 70, 'GAGTGC' : 50, 'GAGTGT' : 50} ,
        'ED' : { 'GAAGAC' : 41, 'GAAGAT' : 59, 'GAGGAC' : 59, 'GAGGAT' : 41} ,
        'EE' : { 'GAAGAA' : 55, 'GAAGAG' : 45, 'GAGGAA' : 33, 'GAGGAG' : 67} ,
        'EF' : { 'GAATTC' : 38, 'GAATTT' : 62, 'GAGTTC' : 53, 'GAGTTT' : 47} ,
        'EG' : { 'GAAGGA' : 32, 'GAAGGC' : 27, 'GAAGGG' : 23, 'GAAGGT' : 19, 'GAGGGA' : 18, 'GAGGGC' : 41, 'GAGGGG' : 25, 'GAGGGT' : 16} ,
        'EH' : { 'GAACAC' : 44, 'GAACAT' : 56, 'GAGCAC' : 64, 'GAGCAT' : 36} ,
        'EI' : { 'GAAATA' : 28, 'GAAATC' : 31, 'GAAATT' : 42, 'GAGATA' : 15, 'GAGATC' : 52, 'GAGATT' : 33} ,
        'EK' : { 'GAAAAA' : 51, 'GAAAAG' : 49, 'GAGAAA' : 40, 'GAGAAG' : 60} ,
        'EL' : { 'GAACTA' : 11, 'GAACTC' : 15, 'GAACTG' : 28, 'GAACTT' : 19, 'GAATTA' : 13, 'GAATTG' : 14, 'GAGCTA' : 7, 'GAGCTC' : 19, 'GAGCTG' : 51, 'GAGCTT' : 10, 'GAGTTA' : 5, 'GAGTTG' : 9} ,
        'EM' : { 'GAAATG' : 100, 'GAGATG' : 100} ,
        'EN' : { 'GAAAAC' : 41, 'GAAAAT' : 59, 'GAGAAC' : 57, 'GAGAAT' : 43} ,
        'EP' : { 'GAACCA' : 34, 'GAACCC' : 26, 'GAACCG' : 7, 'GAACCT' : 32, 'GAGCCA' : 24, 'GAGCCC' : 37, 'GAGCCG' : 14, 'GAGCCT' : 25} ,
        'EQ' : { 'GAACAA' : 41, 'GAACAG' : 59, 'GAGCAA' : 23, 'GAGCAG' : 77} ,
        'ER' : { 'GAAAGA' : 36, 'GAAAGG' : 25, 'GAACGA' : 10, 'GAACGC' : 9, 'GAACGG' : 12, 'GAACGT' : 8, 'GAGAGA' : 20, 'GAGAGG' : 23, 'GAGCGA' : 8, 'GAGCGC' : 19, 'GAGCGG' : 23, 'GAGCGT' : 7} ,
        'ES' : { 'GAAAGC' : 21, 'GAAAGT' : 24, 'GAATCA' : 16, 'GAATCC' : 15, 'GAATCG' : 3, 'GAATCT' : 21, 'GAGAGC' : 31, 'GAGAGT' : 17, 'GAGTCA' : 12, 'GAGTCC' : 19, 'GAGTCG' : 6, 'GAGTCT' : 15} ,
        'ET' : { 'GAAACA' : 36, 'GAAACC' : 26, 'GAAACG' : 7, 'GAAACT' : 30, 'GAGACA' : 26, 'GAGACC' : 38, 'GAGACG' : 15, 'GAGACT' : 21} ,
        'EV' : { 'GAAGTA' : 19, 'GAAGTC' : 21, 'GAAGTG' : 35, 'GAAGTT' : 25, 'GAGGTA' : 10, 'GAGGTC' : 22, 'GAGGTG' : 54, 'GAGGTT' : 14} ,
        'EW' : { 'GAATGG' : 100, 'GAGTGG' : 100} ,
        'EY' : { 'GAATAC' : 39, 'GAATAT' : 61, 'GAGTAC' : 59, 'GAGTAT' : 41} ,
        'FA' : { 'TTCGCA' : 13, 'TTCGCC' : 53, 'TTCGCG' : 18, 'TTCGCT' : 16, 'TTTGCA' : 27, 'TTTGCC' : 36, 'TTTGCG' : 5, 'TTTGCT' : 32} ,
        'FC' : { 'TTCTGC' : 57, 'TTCTGT' : 43, 'TTTTGC' : 42, 'TTTTGT' : 58} ,
        'FU' : { 'TTCTGC' : 57, 'TTCTGT' : 43, 'TTTTGC' : 42, 'TTTTGT' : 58} ,
        'FD' : { 'TTCGAC' : 68, 'TTCGAT' : 32, 'TTTGAC' : 47, 'TTTGAT' : 53} ,
        'FE' : { 'TTCGAA' : 23, 'TTCGAG' : 77, 'TTTGAA' : 52, 'TTTGAG' : 48} ,
        'FF' : { 'TTCTTC' : 63, 'TTCTTT' : 37, 'TTTTTC' : 50, 'TTTTTT' : 50} ,
        'FG' : { 'TTCGGA' : 14, 'TTCGGC' : 42, 'TTCGGG' : 32, 'TTCGGT' : 11, 'TTTGGA' : 31, 'TTTGGC' : 25, 'TTTGGG' : 25, 'TTTGGT' : 19} ,
        'FH' : { 'TTCCAC' : 63, 'TTCCAT' : 37, 'TTTCAC' : 43, 'TTTCAT' : 57} ,
        'FI' : { 'TTCATA' : 11, 'TTCATC' : 56, 'TTCATT' : 33, 'TTTATA' : 21, 'TTTATC' : 31, 'TTTATT' : 48} ,
        'FK' : { 'TTCAAA' : 39, 'TTCAAG' : 61, 'TTTAAA' : 54, 'TTTAAG' : 46} ,
        'FL' : { 'TTCCTA' : 6, 'TTCCTC' : 22, 'TTCCTG' : 43, 'TTCCTT' : 11, 'TTCTTA' : 6, 'TTCTTG' : 11, 'TTTCTA' : 10, 'TTTCTC' : 17, 'TTTCTG' : 27, 'TTTCTT' : 19, 'TTTTTA' : 12, 'TTTTTG' : 15} ,
        'FM' : { 'TTCATG' : 100, 'TTTATG' : 100} ,
        'FN' : { 'TTCAAC' : 62, 'TTCAAT' : 38, 'TTTAAC' : 44, 'TTTAAT' : 56} ,
        'FP' : { 'TTCCCA' : 26, 'TTCCCC' : 34, 'TTCCCG' : 13, 'TTCCCT' : 26, 'TTTCCA' : 33, 'TTTCCC' : 24, 'TTTCCG' : 5, 'TTTCCT' : 38} ,
        'FQ' : { 'TTCCAA' : 22, 'TTCCAG' : 78, 'TTTCAA' : 37, 'TTTCAG' : 63} ,
        'FR' : { 'TTCAGA' : 17, 'TTCAGG' : 17, 'TTCCGA' : 12, 'TTCCGC' : 23, 'TTCCGG' : 22, 'TTCCGT' : 8, 'TTTAGA' : 24, 'TTTAGG' : 17, 'TTTCGA' : 20, 'TTTCGC' : 12, 'TTTCGG' : 16, 'TTTCGT' : 11} ,
        'FS' : { 'TTCAGC' : 26, 'TTCAGT' : 16, 'TTCTCA' : 12, 'TTCTCC' : 25, 'TTCTCG' : 5, 'TTCTCT' : 16, 'TTTAGC' : 12, 'TTTAGT' : 14, 'TTTTCA' : 22, 'TTTTCC' : 20, 'TTTTCG' : 3, 'TTTTCT' : 29} ,
        'FT' : { 'TTCACA' : 23, 'TTCACC' : 42, 'TTCACG' : 12, 'TTCACT' : 23, 'TTTACA' : 32, 'TTTACC' : 30, 'TTTACG' : 5, 'TTTACT' : 33} ,
        'FV' : { 'TTCGTA' : 6, 'TTCGTC' : 26, 'TTCGTG' : 59, 'TTCGTT' : 9, 'TTTGTA' : 14, 'TTTGTC' : 22, 'TTTGTG' : 42, 'TTTGTT' : 22} ,
        'FW' : { 'TTCTGG' : 100, 'TTTTGG' : 100} ,
        'FY' : { 'TTCTAC' : 62, 'TTCTAT' : 38, 'TTTTAC' : 40, 'TTTTAT' : 60} ,
        'GA' : { 'GGAGCA' : 27, 'GGAGCC' : 37, 'GGAGCG' : 7, 'GGAGCT' : 29, 'GGCGCA' : 11, 'GGCGCC' : 48, 'GGCGCG' : 26, 'GGCGCT' : 14, 'GGGGCA' : 18, 'GGGGCC' : 44, 'GGGGCG' : 13, 'GGGGCT' : 25, 'GGTGCA' : 23, 'GGTGCC' : 39, 'GGTGCG' : 7, 'GGTGCT' : 31} ,
        'GC' : { 'GGATGC' : 45, 'GGATGT' : 55, 'GGCTGC' : 63, 'GGCTGT' : 37, 'GGGTGC' : 55, 'GGGTGT' : 45, 'GGTTGC' : 42, 'GGTTGT' : 58} ,
        'GU' : { 'GGATGC' : 45, 'GGATGT' : 55, 'GGCTGC' : 63, 'GGCTGT' : 37, 'GGGTGC' : 55, 'GGGTGT' : 45, 'GGTTGC' : 42, 'GGTTGT' : 58} ,
        'GD' : { 'GGAGAC' : 52, 'GGAGAT' : 48, 'GGCGAC' : 69, 'GGCGAT' : 31, 'GGGGAC' : 65, 'GGGGAT' : 35, 'GGTGAC' : 53, 'GGTGAT' : 47} ,
        'GE' : { 'GGAGAA' : 48, 'GGAGAG' : 52, 'GGCGAA' : 17, 'GGCGAG' : 83, 'GGGGAA' : 32, 'GGGGAG' : 68, 'GGTGAA' : 50, 'GGTGAG' : 50} ,
        'GF' : { 'GGATTC' : 37, 'GGATTT' : 63, 'GGCTTC' : 65, 'GGCTTT' : 35, 'GGGTTC' : 47, 'GGGTTT' : 53, 'GGTTTC' : 48, 'GGTTTT' : 52} ,
        'GG' : { 'GGAGGA' : 31, 'GGAGGC' : 31, 'GGAGGG' : 20, 'GGAGGT' : 18, 'GGCGGA' : 9, 'GGCGGC' : 56, 'GGCGGG' : 24, 'GGCGGT' : 10, 'GGGGGA' : 17, 'GGGGGC' : 53, 'GGGGGG' : 12, 'GGGGGT' : 18, 'GGTGGA' : 24, 'GGTGGC' : 34, 'GGTGGG' : 21, 'GGTGGT' : 21} ,
        'GH' : { 'GGACAC' : 52, 'GGACAT' : 48, 'GGCCAC' : 69, 'GGCCAT' : 31, 'GGGCAC' : 64, 'GGGCAT' : 36, 'GGTCAC' : 53, 'GGTCAT' : 47} ,
        'GI' : { 'GGAATA' : 24, 'GGAATC' : 33, 'GGAATT' : 42, 'GGCATA' : 10, 'GGCATC' : 60, 'GGCATT' : 30, 'GGGATA' : 17, 'GGGATC' : 48, 'GGGATT' : 35, 'GGTATA' : 19, 'GGTATC' : 38, 'GGTATT' : 43} ,
        'GK' : { 'GGAAAA' : 55, 'GGAAAG' : 45, 'GGCAAA' : 40, 'GGCAAG' : 60, 'GGGAAA' : 43, 'GGGAAG' : 57, 'GGTAAA' : 70, 'GGTAAG' : 30} ,
        'GL' : { 'GGACTA' : 9, 'GGACTC' : 19, 'GGACTG' : 27, 'GGACTT' : 19, 'GGATTA' : 12, 'GGATTG' : 13, 'GGCCTA' : 5, 'GGCCTC' : 25, 'GGCCTG' : 42, 'GGCCTT' : 11, 'GGCTTA' : 5, 'GGCTTG' : 12, 'GGGCTA' : 6, 'GGGCTC' : 22, 'GGGCTG' : 49, 'GGGCTT' : 12, 'GGGTTA' : 4, 'GGGTTG' : 7, 'GGTCTA' : 7, 'GGTCTC' : 19, 'GGTCTG' : 27, 'GGTCTT' : 17, 'GGTTTA' : 12, 'GGTTTG' : 18} ,
        'GM' : { 'GGAATG' : 100, 'GGCATG' : 100, 'GGGATG' : 100, 'GGTATG' : 100} ,
        'GN' : { 'GGAAAC' : 45, 'GGAAAT' : 55, 'GGCAAC' : 63, 'GGCAAT' : 37, 'GGGAAC' : 60, 'GGGAAT' : 40, 'GGTAAC' : 50, 'GGTAAT' : 50} ,
        'GP' : { 'GGACCA' : 30, 'GGACCC' : 34, 'GGACCG' : 6, 'GGACCT' : 29, 'GGCCCA' : 22, 'GGCCCC' : 39, 'GGCCCG' : 14, 'GGCCCT' : 25, 'GGGCCA' : 21, 'GGGCCC' : 40, 'GGGCCG' : 14, 'GGGCCT' : 25, 'GGTCCA' : 25, 'GGTCCC' : 35, 'GGTCCG' : 7, 'GGTCCT' : 34} ,
        'GQ' : { 'GGACAA' : 35, 'GGACAG' : 65, 'GGCCAA' : 21, 'GGCCAG' : 79, 'GGGCAA' : 20, 'GGGCAG' : 80, 'GGTCAA' : 33, 'GGTCAG' : 67} ,
        'GR' : { 'GGAAGA' : 38, 'GGAAGG' : 22, 'GGACGA' : 10, 'GGACGC' : 11, 'GGACGG' : 12, 'GGACGT' : 8, 'GGCAGA' : 14, 'GGCAGG' : 17, 'GGCCGA' : 10, 'GGCCGC' : 29, 'GGCCGG' : 22, 'GGCCGT' : 8, 'GGGAGA' : 20, 'GGGAGG' : 22, 'GGGCGA' : 9, 'GGGCGC' : 21, 'GGGCGG' : 21, 'GGGCGT' : 7, 'GGTAGA' : 18, 'GGTAGG' : 10, 'GGTCGA' : 17, 'GGTCGC' : 20, 'GGTCGG' : 19, 'GGTCGT' : 15} ,
        'GS' : { 'GGAAGC' : 25, 'GGAAGT' : 24, 'GGATCA' : 14, 'GGATCC' : 16, 'GGATCG' : 3, 'GGATCT' : 17, 'GGCAGC' : 30, 'GGCAGT' : 13, 'GGCTCA' : 10, 'GGCTCC' : 26, 'GGCTCG' : 6, 'GGCTCT' : 14, 'GGGAGC' : 31, 'GGGAGT' : 15, 'GGGTCA' : 12, 'GGGTCC' : 20, 'GGGTCG' : 5, 'GGGTCT' : 16, 'GGTAGC' : 13, 'GGTAGT' : 11, 'GGTTCA' : 18, 'GGTTCC' : 25, 'GGTTCG' : 4, 'GGTTCT' : 29} ,
        'GT' : { 'GGAACA' : 35, 'GGAACC' : 29, 'GGAACG' : 7, 'GGAACT' : 29, 'GGCACA' : 23, 'GGCACC' : 44, 'GGCACG' : 13, 'GGCACT' : 20, 'GGGACA' : 27, 'GGGACC' : 39, 'GGGACG' : 12, 'GGGACT' : 22, 'GGTACA' : 29, 'GGTACC' : 36, 'GGTACG' : 5, 'GGTACT' : 30} ,
        'GV' : { 'GGAGTA' : 17, 'GGAGTC' : 25, 'GGAGTG' : 34, 'GGAGTT' : 24, 'GGCGTA' : 5, 'GGCGTC' : 32, 'GGCGTG' : 53, 'GGCGTT' : 10, 'GGGGTA' : 9, 'GGGGTC' : 32, 'GGGGTG' : 42, 'GGGGTT' : 16, 'GGTGTA' : 12, 'GGTGTC' : 26, 'GGTGTG' : 40, 'GGTGTT' : 23} ,
        'GW' : { 'GGATGG' : 100, 'GGCTGG' : 100, 'GGGTGG' : 100, 'GGTTGG' : 100} ,
        'GY' : { 'GGATAC' : 42, 'GGATAT' : 58, 'GGCTAC' : 64, 'GGCTAT' : 36, 'GGGTAC' : 60, 'GGGTAT' : 40, 'GGTTAC' : 44, 'GGTTAT' : 56} ,
        'HA' : { 'CACGCA' : 17, 'CACGCC' : 47, 'CACGCG' : 19, 'CACGCT' : 17, 'CATGCA' : 26, 'CATGCC' : 36, 'CATGCG' : 7, 'CATGCT' : 32} ,
        'HC' : { 'CACTGC' : 61, 'CACTGT' : 39, 'CATTGC' : 47, 'CATTGT' : 53} ,
        'HU' : { 'CACTGC' : 61, 'CACTGT' : 39, 'CATTGC' : 47, 'CATTGT' : 53} ,
        'HD' : { 'CACGAC' : 66, 'CACGAT' : 34, 'CATGAC' : 49, 'CATGAT' : 51} ,
        'HE' : { 'CACGAA' : 26, 'CACGAG' : 74, 'CATGAA' : 52, 'CATGAG' : 48} ,
        'HF' : { 'CACTTC' : 61, 'CACTTT' : 39, 'CATTTC' : 46, 'CATTTT' : 54} ,
        'HG' : { 'CACGGA' : 17, 'CACGGC' : 42, 'CACGGG' : 29, 'CACGGT' : 12, 'CATGGA' : 29, 'CATGGC' : 29, 'CATGGG' : 23, 'CATGGT' : 19} ,
        'HH' : { 'CACCAC' : 66, 'CACCAT' : 34, 'CATCAC' : 49, 'CATCAT' : 51} ,
        'HI' : { 'CACATA' : 14, 'CACATC' : 55, 'CACATT' : 31, 'CATATA' : 22, 'CATATC' : 32, 'CATATT' : 46} ,
        'HK' : { 'CACAAA' : 39, 'CACAAG' : 61, 'CATAAA' : 51, 'CATAAG' : 49} ,
        'HL' : { 'CACCTA' : 6, 'CACCTC' : 21, 'CACCTG' : 44, 'CACCTT' : 12, 'CACTTA' : 6, 'CACTTG' : 12, 'CATCTA' : 8, 'CATCTC' : 17, 'CATCTG' : 26, 'CATCTT' : 18, 'CATTTA' : 13, 'CATTTG' : 18} ,
        'HM' : { 'CACATG' : 100, 'CATATG' : 100} ,
        'HN' : { 'CACAAC' : 60, 'CACAAT' : 40, 'CATAAC' : 41, 'CATAAT' : 59} ,
        'HP' : { 'CACCCA' : 25, 'CACCCC' : 33, 'CACCCG' : 17, 'CACCCT' : 25, 'CATCCA' : 31, 'CATCCC' : 26, 'CATCCG' : 7, 'CATCCT' : 36} ,
        'HQ' : { 'CACCAA' : 21, 'CACCAG' : 79, 'CATCAA' : 30, 'CATCAG' : 70} ,
        'HR' : { 'CACAGA' : 21, 'CACAGG' : 23, 'CACCGA' : 9, 'CACCGC' : 21, 'CACCGG' : 20, 'CACCGT' : 7, 'CATAGA' : 19, 'CATAGG' : 14, 'CATCGA' : 16, 'CATCGC' : 16, 'CATCGG' : 20, 'CATCGT' : 15} ,
        'HS' : { 'CACAGC' : 34, 'CACAGT' : 18, 'CACTCA' : 12, 'CACTCC' : 17, 'CACTCG' : 8, 'CACTCT' : 12, 'CATAGC' : 12, 'CATAGT' : 13, 'CATTCA' : 22, 'CATTCC' : 21, 'CATTCG' : 4, 'CATTCT' : 27} ,
        'HT' : { 'CACACA' : 28, 'CACACC' : 30, 'CACACG' : 17, 'CACACT' : 25, 'CATACA' : 26, 'CATACC' : 21, 'CATACG' : 6, 'CATACT' : 46} ,
        'HV' : { 'CACGTA' : 6, 'CACGTC' : 26, 'CACGTG' : 57, 'CACGTT' : 11, 'CATGTA' : 14, 'CATGTC' : 23, 'CATGTG' : 40, 'CATGTT' : 23} ,
        'HW' : { 'CACTGG' : 100, 'CATTGG' : 100} ,
        'HY' : { 'CACTAC' : 63, 'CACTAT' : 37, 'CATTAC' : 47, 'CATTAT' : 53} ,
        'IA' : { 'ATAGCA' : 36, 'ATAGCC' : 29, 'ATAGCG' : 5, 'ATAGCT' : 31, 'ATCGCA' : 16, 'ATCGCC' : 51, 'ATCGCG' : 14, 'ATCGCT' : 19, 'ATTGCA' : 27, 'ATTGCC' : 33, 'ATTGCG' : 5, 'ATTGCT' : 35} ,
        'IC' : { 'ATATGC' : 37, 'ATATGT' : 63, 'ATCTGC' : 60, 'ATCTGT' : 40, 'ATTTGC' : 43, 'ATTTGT' : 57} ,
        'IU' : { 'ATATGC' : 37, 'ATATGT' : 63, 'ATCTGC' : 60, 'ATCTGT' : 40, 'ATTTGC' : 43, 'ATTTGT' : 57} ,
        'ID' : { 'ATAGAC' : 43, 'ATAGAT' : 57, 'ATCGAC' : 66, 'ATCGAT' : 34, 'ATTGAC' : 47, 'ATTGAT' : 53} ,
        'IE' : { 'ATAGAA' : 61, 'ATAGAG' : 39, 'ATCGAA' : 25, 'ATCGAG' : 75, 'ATTGAA' : 55, 'ATTGAG' : 45} ,
        'IF' : { 'ATATTC' : 33, 'ATATTT' : 67, 'ATCTTC' : 63, 'ATCTTT' : 37, 'ATTTTC' : 45, 'ATTTTT' : 55} ,
        'IG' : { 'ATAGGA' : 37, 'ATAGGC' : 23, 'ATAGGG' : 19, 'ATAGGT' : 20, 'ATCGGA' : 15, 'ATCGGC' : 41, 'ATCGGG' : 31, 'ATCGGT' : 12, 'ATTGGA' : 31, 'ATTGGC' : 27, 'ATTGGG' : 22, 'ATTGGT' : 20} ,
        'IH' : { 'ATACAC' : 40, 'ATACAT' : 60, 'ATCCAC' : 64, 'ATCCAT' : 36, 'ATTCAC' : 40, 'ATTCAT' : 60} ,
        'II' : { 'ATAATA' : 26, 'ATAATC' : 25, 'ATAATT' : 49, 'ATCATA' : 12, 'ATCATC' : 57, 'ATCATT' : 31, 'ATTATA' : 22, 'ATTATC' : 30, 'ATTATT' : 48} ,
        'IK' : { 'ATAAAA' : 55, 'ATAAAG' : 45, 'ATCAAA' : 39, 'ATCAAG' : 61, 'ATTAAA' : 59, 'ATTAAG' : 41} ,
        'IL' : { 'ATACTA' : 10, 'ATACTC' : 11, 'ATACTG' : 23, 'ATACTT' : 19, 'ATATTA' : 18, 'ATATTG' : 18, 'ATCCTA' : 6, 'ATCCTC' : 22, 'ATCCTG' : 41, 'ATCCTT' : 11, 'ATCTTA' : 6, 'ATCTTG' : 12, 'ATTCTA' : 8, 'ATTCTC' : 15, 'ATTCTG' : 24, 'ATTCTT' : 18, 'ATTTTA' : 15, 'ATTTTG' : 19} ,
        'IM' : { 'ATAATG' : 100, 'ATCATG' : 100, 'ATTATG' : 100} ,
        'IN' : { 'ATAAAC' : 39, 'ATAAAT' : 61, 'ATCAAC' : 60, 'ATCAAT' : 40, 'ATTAAC' : 40, 'ATTAAT' : 60} ,
        'IP' : { 'ATACCA' : 36, 'ATACCC' : 23, 'ATACCG' : 5, 'ATACCT' : 36, 'ATCCCA' : 26, 'ATCCCC' : 35, 'ATCCCG' : 12, 'ATCCCT' : 27, 'ATTCCA' : 34, 'ATTCCC' : 23, 'ATTCCG' : 5, 'ATTCCT' : 38} ,
        'IQ' : { 'ATACAA' : 39, 'ATACAG' : 61, 'ATCCAA' : 21, 'ATCCAG' : 79, 'ATTCAA' : 35, 'ATTCAG' : 65} ,
        'IR' : { 'ATAAGA' : 37, 'ATAAGG' : 24, 'ATACGA' : 12, 'ATACGC' : 7, 'ATACGG' : 12, 'ATACGT' : 8, 'ATCAGA' : 17, 'ATCAGG' : 15, 'ATCCGA' : 12, 'ATCCGC' : 25, 'ATCCGG' : 22, 'ATCCGT' : 9, 'ATTAGA' : 19, 'ATTAGG' : 11, 'ATTCGA' : 23, 'ATTCGC' : 16, 'ATTCGG' : 17, 'ATTCGT' : 14} ,
        'IS' : { 'ATAAGC' : 15, 'ATAAGT' : 18, 'ATATCA' : 23, 'ATATCC' : 17, 'ATATCG' : 3, 'ATATCT' : 24, 'ATCAGC' : 27, 'ATCAGT' : 15, 'ATCTCA' : 12, 'ATCTCC' : 26, 'ATCTCG' : 5, 'ATCTCT' : 16, 'ATTAGC' : 10, 'ATTAGT' : 11, 'ATTTCA' : 24, 'ATTTCC' : 22, 'ATTTCG' : 3, 'ATTTCT' : 31} ,
        'IT' : { 'ATAACA' : 36, 'ATAACC' : 24, 'ATAACG' : 6, 'ATAACT' : 34, 'ATCACA' : 26, 'ATCACC' : 41, 'ATCACG' : 11, 'ATCACT' : 23, 'ATTACA' : 36, 'ATTACC' : 25, 'ATTACG' : 6, 'ATTACT' : 33} ,
        'IV' : { 'ATAGTA' : 20, 'ATAGTC' : 19, 'ATAGTG' : 34, 'ATAGTT' : 27, 'ATCGTA' : 6, 'ATCGTC' : 30, 'ATCGTG' : 53, 'ATCGTT' : 10, 'ATTGTA' : 15, 'ATTGTC' : 23, 'ATTGTG' : 38, 'ATTGTT' : 24} ,
        'IW' : { 'ATATGG' : 100, 'ATCTGG' : 100, 'ATTTGG' : 100} ,
        'IY' : { 'ATATAC' : 37, 'ATATAT' : 63, 'ATCTAC' : 61, 'ATCTAT' : 39, 'ATTTAC' : 41, 'ATTTAT' : 59} ,
        'KA' : { 'AAAGCA' : 30, 'AAAGCC' : 35, 'AAAGCG' : 5, 'AAAGCT' : 30, 'AAGGCA' : 23, 'AAGGCC' : 41, 'AAGGCG' : 11, 'AAGGCT' : 25} ,
        'KC' : { 'AAATGC' : 35, 'AAATGT' : 65, 'AAGTGC' : 52, 'AAGTGT' : 48} ,
        'KU' : { 'AAATGC' : 35, 'AAATGT' : 65, 'AAGTGC' : 52, 'AAGTGT' : 48} ,
        'KD' : { 'AAAGAC' : 43, 'AAAGAT' : 57, 'AAGGAC' : 55, 'AAGGAT' : 45} ,
        'KE' : { 'AAAGAA' : 60, 'AAAGAG' : 40, 'AAGGAA' : 41, 'AAGGAG' : 59} ,
        'KF' : { 'AAATTC' : 42, 'AAATTT' : 58, 'AAGTTC' : 52, 'AAGTTT' : 48} ,
        'KG' : { 'AAAGGA' : 35, 'AAAGGC' : 26, 'AAAGGG' : 20, 'AAAGGT' : 18, 'AAGGGA' : 22, 'AAGGGC' : 37, 'AAGGGG' : 24, 'AAGGGT' : 17} ,
        'KH' : { 'AAACAC' : 46, 'AAACAT' : 54, 'AAGCAC' : 59, 'AAGCAT' : 41} ,
        'KI' : { 'AAAATA' : 28, 'AAAATC' : 33, 'AAAATT' : 40, 'AAGATA' : 19, 'AAGATC' : 47, 'AAGATT' : 34} ,
        'KK' : { 'AAAAAA' : 45, 'AAAAAG' : 55, 'AAGAAA' : 44, 'AAGAAG' : 56} ,
        'KL' : { 'AAACTA' : 11, 'AAACTC' : 17, 'AAACTG' : 27, 'AAACTT' : 19, 'AAATTA' : 13, 'AAATTG' : 13, 'AAGCTA' : 8, 'AAGCTC' : 19, 'AAGCTG' : 43, 'AAGCTT' : 13, 'AAGTTA' : 7, 'AAGTTG' : 11} ,
        'KM' : { 'AAAATG' : 100, 'AAGATG' : 100} ,
        'KN' : { 'AAAAAC' : 41, 'AAAAAT' : 59, 'AAGAAC' : 54, 'AAGAAT' : 46} ,
        'KP' : { 'AAACCA' : 33, 'AAACCC' : 30, 'AAACCG' : 5, 'AAACCT' : 32, 'AAGCCA' : 29, 'AAGCCC' : 35, 'AAGCCG' : 11, 'AAGCCT' : 25} ,
        'KQ' : { 'AAACAA' : 44, 'AAACAG' : 56, 'AAGCAA' : 30, 'AAGCAG' : 70} ,
        'KR' : { 'AAAAGA' : 34, 'AAAAGG' : 24, 'AAACGA' : 11, 'AAACGC' : 10, 'AAACGG' : 12, 'AAACGT' : 9, 'AAGAGA' : 24, 'AAGAGG' : 27, 'AAGCGA' : 9, 'AAGCGC' : 16, 'AAGCGG' : 18, 'AAGCGT' : 6} ,
        'KS' : { 'AAAAGC' : 18, 'AAAAGT' : 19, 'AAATCA' : 20, 'AAATCC' : 17, 'AAATCG' : 4, 'AAATCT' : 23, 'AAGAGC' : 25, 'AAGAGT' : 16, 'AAGTCA' : 17, 'AAGTCC' : 20, 'AAGTCG' : 5, 'AAGTCT' : 17} ,
        'KT' : { 'AAAACA' : 36, 'AAAACC' : 28, 'AAAACG' : 7, 'AAAACT' : 29, 'AAGACA' : 29, 'AAGACC' : 35, 'AAGACG' : 13, 'AAGACT' : 23} ,
        'KV' : { 'AAAGTA' : 18, 'AAAGTC' : 22, 'AAAGTG' : 36, 'AAAGTT' : 24, 'AAGGTA' : 10, 'AAGGTC' : 24, 'AAGGTG' : 49, 'AAGGTT' : 17} ,
        'KW' : { 'AAATGG' : 100, 'AAGTGG' : 100} ,
        'KY' : { 'AAATAC' : 47, 'AAATAT' : 53, 'AAGTAC' : 59, 'AAGTAT' : 41} ,
        'LA' : { 'CTAGCA' : 32, 'CTAGCC' : 31, 'CTAGCG' : 7, 'CTAGCT' : 30, 'CTCGCA' : 15, 'CTCGCC' : 51, 'CTCGCG' : 13, 'CTCGCT' : 21, 'CTGGCA' : 18, 'CTGGCC' : 46, 'CTGGCG' : 13, 'CTGGCT' : 23, 'CTTGCA' : 28, 'CTTGCC' : 32, 'CTTGCG' : 5, 'CTTGCT' : 35, 'TTAGCA' : 35, 'TTAGCC' : 26, 'TTAGCG' : 5, 'TTAGCT' : 34, 'TTGGCA' : 27, 'TTGGCC' : 33, 'TTGGCG' : 8, 'TTGGCT' : 32} ,
        'LC' : { 'CTATGC' : 44, 'CTATGT' : 56, 'CTCTGC' : 61, 'CTCTGT' : 39, 'CTGTGC' : 60, 'CTGTGT' : 40, 'CTTTGC' : 45, 'CTTTGT' : 55, 'TTATGC' : 38, 'TTATGT' : 62, 'TTGTGC' : 45, 'TTGTGT' : 55} ,
        'LU' : { 'CTATGC' : 44, 'CTATGT' : 56, 'CTCTGC' : 61, 'CTCTGT' : 39, 'CTGTGC' : 60, 'CTGTGT' : 40, 'CTTTGC' : 45, 'CTTTGT' : 55, 'TTATGC' : 38, 'TTATGT' : 62, 'TTGTGC' : 45, 'TTGTGT' : 55} ,
        'LD' : { 'CTAGAC' : 47, 'CTAGAT' : 53, 'CTCGAC' : 62, 'CTCGAT' : 38, 'CTGGAC' : 60, 'CTGGAT' : 40, 'CTTGAC' : 43, 'CTTGAT' : 57, 'TTAGAC' : 35, 'TTAGAT' : 65, 'TTGGAC' : 42, 'TTGGAT' : 58} ,
        'LE' : { 'CTAGAA' : 51, 'CTAGAG' : 49, 'CTCGAA' : 28, 'CTCGAG' : 72, 'CTGGAA' : 29, 'CTGGAG' : 71, 'CTTGAA' : 58, 'CTTGAG' : 42, 'TTAGAA' : 62, 'TTAGAG' : 38, 'TTGGAA' : 49, 'TTGGAG' : 51} ,
        'LF' : { 'CTATTC' : 36, 'CTATTT' : 64, 'CTCTTC' : 65, 'CTCTTT' : 35, 'CTGTTC' : 52, 'CTGTTT' : 48, 'CTTTTC' : 50, 'CTTTTT' : 50, 'TTATTC' : 31, 'TTATTT' : 69, 'TTGTTC' : 38, 'TTGTTT' : 62} ,
        'LG' : { 'CTAGGA' : 32, 'CTAGGC' : 27, 'CTAGGG' : 23, 'CTAGGT' : 18, 'CTCGGA' : 15, 'CTCGGC' : 41, 'CTCGGG' : 33, 'CTCGGT' : 12, 'CTGGGA' : 18, 'CTGGGC' : 43, 'CTGGGG' : 26, 'CTGGGT' : 13, 'CTTGGA' : 32, 'CTTGGC' : 26, 'CTTGGG' : 23, 'CTTGGT' : 19, 'TTAGGA' : 39, 'TTAGGC' : 20, 'TTAGGG' : 20, 'TTAGGT' : 21, 'TTGGGA' : 30, 'TTGGGC' : 28, 'TTGGGG' : 24, 'TTGGGT' : 19} ,
        'LH' : { 'CTACAC' : 49, 'CTACAT' : 51, 'CTCCAC' : 61, 'CTCCAT' : 39, 'CTGCAC' : 70, 'CTGCAT' : 30, 'CTTCAC' : 42, 'CTTCAT' : 58, 'TTACAC' : 40, 'TTACAT' : 60, 'TTGCAC' : 51, 'TTGCAT' : 49} ,
        'LI' : { 'CTAATA' : 26, 'CTAATC' : 28, 'CTAATT' : 46, 'CTCATA' : 10, 'CTCATC' : 60, 'CTCATT' : 30, 'CTGATA' : 17, 'CTGATC' : 47, 'CTGATT' : 36, 'CTTATA' : 20, 'CTTATC' : 32, 'CTTATT' : 48, 'TTAATA' : 29, 'TTAATC' : 24, 'TTAATT' : 46, 'TTGATA' : 25, 'TTGATC' : 33, 'TTGATT' : 42} ,
        'LK' : { 'CTAAAA' : 48, 'CTAAAG' : 52, 'CTCAAA' : 32, 'CTCAAG' : 68, 'CTGAAA' : 34, 'CTGAAG' : 66, 'CTTAAA' : 56, 'CTTAAG' : 44, 'TTAAAA' : 56, 'TTAAAG' : 44, 'TTGAAA' : 50, 'TTGAAG' : 50} ,
        'LL' : { 'CTACTA' : 11, 'CTACTC' : 17, 'CTACTG' : 33, 'CTACTT' : 18, 'CTATTA' : 8, 'CTATTG' : 12, 'CTCCTA' : 6, 'CTCCTC' : 22, 'CTCCTG' : 44, 'CTCCTT' : 11, 'CTCTTA' : 6, 'CTCTTG' : 12, 'CTGCTA' : 5, 'CTGCTC' : 22, 'CTGCTG' : 53, 'CTGCTT' : 10, 'CTGTTA' : 3, 'CTGTTG' : 8, 'CTTCTA' : 9, 'CTTCTC' : 18, 'CTTCTG' : 30, 'CTTCTT' : 17, 'CTTTTA' : 10, 'CTTTTG' : 15, 'TTACTA' : 12, 'TTACTC' : 12, 'TTACTG' : 25, 'TTACTT' : 19, 'TTATTA' : 15, 'TTATTG' : 16, 'TTGCTA' : 10, 'TTGCTC' : 15, 'TTGCTG' : 39, 'TTGCTT' : 16, 'TTGTTA' : 8, 'TTGTTG' : 13} ,
        'LM' : { 'CTAATG' : 100, 'CTCATG' : 100, 'CTGATG' : 100, 'CTTATG' : 100, 'TTAATG' : 100, 'TTGATG' : 100} ,
        'LN' : { 'CTAAAC' : 39, 'CTAAAT' : 61, 'CTCAAC' : 64, 'CTCAAT' : 36, 'CTGAAC' : 56, 'CTGAAT' : 44, 'CTTAAC' : 39, 'CTTAAT' : 61, 'TTAAAC' : 33, 'TTAAAT' : 67, 'TTGAAC' : 41, 'TTGAAT' : 59} ,
        'LP' : { 'CTACCA' : 30, 'CTACCC' : 30, 'CTACCG' : 7, 'CTACCT' : 34, 'CTCCCA' : 27, 'CTCCCC' : 30, 'CTCCCG' : 15, 'CTCCCT' : 28, 'CTGCCA' : 18, 'CTGCCC' : 44, 'CTGCCG' : 13, 'CTGCCT' : 24, 'CTTCCA' : 32, 'CTTCCC' : 25, 'CTTCCG' : 6, 'CTTCCT' : 38, 'TTACCA' : 36, 'TTACCC' : 20, 'TTACCG' : 5, 'TTACCT' : 39, 'TTGCCA' : 33, 'TTGCCC' : 27, 'TTGCCG' : 8, 'TTGCCT' : 32} ,
        'LQ' : { 'CTACAA' : 32, 'CTACAG' : 68, 'CTCCAA' : 22, 'CTCCAG' : 78, 'CTGCAA' : 16, 'CTGCAG' : 84, 'CTTCAA' : 34, 'CTTCAG' : 66, 'TTACAA' : 40, 'TTACAG' : 60, 'TTGCAA' : 32, 'TTGCAG' : 68} ,
        'LR' : { 'CTAAGA' : 27, 'CTAAGG' : 22, 'CTACGA' : 11, 'CTACGC' : 13, 'CTACGG' : 18, 'CTACGT' : 9, 'CTCAGA' : 16, 'CTCAGG' : 23, 'CTCCGA' : 13, 'CTCCGC' : 18, 'CTCCGG' : 23, 'CTCCGT' : 8, 'CTGAGA' : 14, 'CTGAGG' : 20, 'CTGCGA' : 7, 'CTGCGC' : 24, 'CTGCGG' : 28, 'CTGCGT' : 7, 'CTTAGA' : 14, 'CTTAGG' : 10, 'CTTCGA' : 23, 'CTTCGC' : 15, 'CTTCGG' : 24, 'CTTCGT' : 14, 'TTAAGA' : 35, 'TTAAGG' : 23, 'TTACGA' : 14, 'TTACGC' : 7, 'TTACGG' : 13, 'TTACGT' : 8, 'TTGAGA' : 29, 'TTGAGG' : 25, 'TTGCGA' : 10, 'TTGCGC' : 12, 'TTGCGG' : 17, 'TTGCGT' : 7} ,
        'LS' : { 'CTAAGC' : 23, 'CTAAGT' : 23, 'CTATCA' : 16, 'CTATCC' : 17, 'CTATCG' : 3, 'CTATCT' : 18, 'CTCAGC' : 30, 'CTCAGT' : 15, 'CTCTCA' : 12, 'CTCTCC' : 23, 'CTCTCG' : 6, 'CTCTCT' : 15, 'CTGAGC' : 29, 'CTGAGT' : 14, 'CTGTCA' : 11, 'CTGTCC' : 24, 'CTGTCG' : 7, 'CTGTCT' : 15, 'CTTAGC' : 10, 'CTTAGT' : 11, 'CTTTCA' : 21, 'CTTTCC' : 23, 'CTTTCG' : 4, 'CTTTCT' : 31, 'TTAAGC' : 13, 'TTAAGT' : 20, 'TTATCA' : 22, 'TTATCC' : 16, 'TTATCG' : 3, 'TTATCT' : 26, 'TTGAGC' : 17, 'TTGAGT' : 16, 'TTGTCA' : 19, 'TTGTCC' : 20, 'TTGTCG' : 4, 'TTGTCT' : 22} ,
        'LT' : { 'CTAACA' : 35, 'CTAACC' : 25, 'CTAACG' : 8, 'CTAACT' : 32, 'CTCACA' : 21, 'CTCACC' : 43, 'CTCACG' : 13, 'CTCACT' : 23, 'CTGACA' : 24, 'CTGACC' : 40, 'CTGACG' : 14, 'CTGACT' : 23, 'CTTACA' : 28, 'CTTACC' : 27, 'CTTACG' : 6, 'CTTACT' : 38, 'TTAACA' : 39, 'TTAACC' : 22, 'TTAACG' : 5, 'TTAACT' : 34, 'TTGACA' : 33, 'TTGACC' : 28, 'TTGACG' : 9, 'TTGACT' : 31} ,
        'LV' : { 'CTAGTA' : 17, 'CTAGTC' : 18, 'CTAGTG' : 44, 'CTAGTT' : 21, 'CTCGTA' : 6, 'CTCGTC' : 30, 'CTCGTG' : 52, 'CTCGTT' : 12, 'CTGGTA' : 8, 'CTGGTC' : 22, 'CTGGTG' : 57, 'CTGGTT' : 13, 'CTTGTA' : 15, 'CTTGTC' : 23, 'CTTGTG' : 37, 'CTTGTT' : 25, 'TTAGTA' : 23, 'TTAGTC' : 16, 'TTAGTG' : 33, 'TTAGTT' : 28, 'TTGGTA' : 15, 'TTGGTC' : 20, 'TTGGTG' : 43, 'TTGGTT' : 22} ,
        'LW' : { 'CTATGG' : 100, 'CTCTGG' : 100, 'CTGTGG' : 100, 'CTTTGG' : 100, 'TTATGG' : 100, 'TTGTGG' : 100} ,
        'LY' : { 'CTATAC' : 45, 'CTATAT' : 55, 'CTCTAC' : 63, 'CTCTAT' : 37, 'CTGTAC' : 61, 'CTGTAT' : 39, 'CTTTAC' : 40, 'CTTTAT' : 60, 'TTATAC' : 36, 'TTATAT' : 64, 'TTGTAC' : 47, 'TTGTAT' : 53} ,
        'MA' : { 'ATGGCA' : 23, 'ATGGCC' : 36, 'ATGGCG' : 15, 'ATGGCT' : 26} ,
        'MC' : { 'ATGTGC' : 51, 'ATGTGT' : 49} ,
        'MU' : { 'ATGTGC' : 51, 'ATGTGT' : 49} ,
        'MD' : { 'ATGGAC' : 51, 'ATGGAT' : 49} ,
        'ME' : { 'ATGGAA' : 42, 'ATGGAG' : 58} ,
        'MF' : { 'ATGTTC' : 49, 'ATGTTT' : 51} ,
        'MG' : { 'ATGGGA' : 25, 'ATGGGC' : 34, 'ATGGGG' : 24, 'ATGGGT' : 17} ,
        'MH' : { 'ATGCAC' : 58, 'ATGCAT' : 42} ,
        'MI' : { 'ATGATA' : 19, 'ATGATC' : 43, 'ATGATT' : 38} ,
        'MK' : { 'ATGAAA' : 41, 'ATGAAG' : 59} ,
        'ML' : { 'ATGCTA' : 8, 'ATGCTC' : 15, 'ATGCTG' : 43, 'ATGCTT' : 13, 'ATGTTA' : 8, 'ATGTTG' : 13} ,
        'MM' : { 'ATGATG' : 100} ,
        'MN' : { 'ATGAAC' : 52, 'ATGAAT' : 48} ,
        'MP' : { 'ATGCCA' : 27, 'ATGCCC' : 33, 'ATGCCG' : 11, 'ATGCCT' : 29} ,
        'MQ' : { 'ATGCAA' : 27, 'ATGCAG' : 73} ,
        'MR' : { 'ATGAGA' : 26, 'ATGAGG' : 27, 'ATGCGA' : 8, 'ATGCGC' : 14, 'ATGCGG' : 19, 'ATGCGT' : 7} ,
        'MS' : { 'ATGAGC' : 22, 'ATGAGT' : 15, 'ATGTCA' : 15, 'ATGTCC' : 21, 'ATGTCG' : 7, 'ATGTCT' : 20} ,
        'MT' : { 'ATGACA' : 28, 'ATGACC' : 35, 'ATGACG' : 10, 'ATGACT' : 26} ,
        'MV' : { 'ATGGTA' : 12, 'ATGGTC' : 21, 'ATGGTG' : 49, 'ATGGTT' : 18} ,
        'MW' : { 'ATGTGG' : 100} ,
        'MY' : { 'ATGTAC' : 56, 'ATGTAT' : 44} ,
        'NA' : { 'AACGCA' : 18, 'AACGCC' : 48, 'AACGCG' : 17, 'AACGCT' : 17, 'AATGCA' : 28, 'AATGCC' : 36, 'AATGCG' : 5, 'AATGCT' : 31} ,
        'NC' : { 'AACTGC' : 58, 'AACTGT' : 42, 'AATTGC' : 44, 'AATTGT' : 56} ,
        'NU' : { 'AACTGC' : 58, 'AACTGT' : 42, 'AATTGC' : 44, 'AATTGT' : 56} ,
        'ND' : { 'AACGAC' : 66, 'AACGAT' : 34, 'AATGAC' : 49, 'AATGAT' : 51} ,
        'NE' : { 'AACGAA' : 29, 'AACGAG' : 71, 'AATGAA' : 56, 'AATGAG' : 44} ,
        'NF' : { 'AACTTC' : 60, 'AACTTT' : 40, 'AATTTC' : 43, 'AATTTT' : 57} ,
        'NG' : { 'AACGGA' : 18, 'AACGGC' : 38, 'AACGGG' : 33, 'AACGGT' : 12, 'AATGGA' : 32, 'AATGGC' : 28, 'AATGGG' : 23, 'AATGGT' : 17} ,
        'NH' : { 'AACCAC' : 63, 'AACCAT' : 37, 'AATCAC' : 46, 'AATCAT' : 54} ,
        'NI' : { 'AACATA' : 15, 'AACATC' : 54, 'AACATT' : 32, 'AATATA' : 23, 'AATATC' : 31, 'AATATT' : 46} ,
        'NK' : { 'AACAAA' : 43, 'AACAAG' : 57, 'AATAAA' : 59, 'AATAAG' : 41} ,
        'NL' : { 'AACCTA' : 7, 'AACCTC' : 23, 'AACCTG' : 38, 'AACCTT' : 12, 'AACTTA' : 7, 'AACTTG' : 14, 'AATCTA' : 9, 'AATCTC' : 16, 'AATCTG' : 23, 'AATCTT' : 19, 'AATTTA' : 15, 'AATTTG' : 18} ,
        'NM' : { 'AACATG' : 100, 'AATATG' : 100} ,
        'NN' : { 'AACAAC' : 60, 'AACAAT' : 40, 'AATAAC' : 42, 'AATAAT' : 58} ,
        'NP' : { 'AACCCA' : 27, 'AACCCC' : 36, 'AACCCG' : 11, 'AACCCT' : 27, 'AATCCA' : 35, 'AATCCC' : 24, 'AATCCG' : 5, 'AATCCT' : 35} ,
        'NQ' : { 'AACCAA' : 27, 'AACCAG' : 73, 'AATCAA' : 41, 'AATCAG' : 59} ,
        'NR' : { 'AACAGA' : 25, 'AACAGG' : 21, 'AACCGA' : 10, 'AACCGC' : 19, 'AACCGG' : 18, 'AACCGT' : 7, 'AATAGA' : 28, 'AATAGG' : 16, 'AATCGA' : 16, 'AATCGC' : 14, 'AATCGG' : 14, 'AATCGT' : 12} ,
        'NS' : { 'AACAGC' : 32, 'AACAGT' : 17, 'AACTCA' : 13, 'AACTCC' : 19, 'AACTCG' : 5, 'AACTCT' : 14, 'AATAGC' : 14, 'AATAGT' : 14, 'AATTCA' : 24, 'AATTCC' : 19, 'AATTCG' : 3, 'AATTCT' : 26} ,
        'NT' : { 'AACACA' : 30, 'AACACC' : 34, 'AACACG' : 13, 'AACACT' : 23, 'AATACA' : 36, 'AATACC' : 26, 'AATACG' : 7, 'AATACT' : 31} ,
        'NV' : { 'AACGTA' : 6, 'AACGTC' : 27, 'AACGTG' : 55, 'AACGTT' : 11, 'AATGTA' : 15, 'AATGTC' : 24, 'AATGTG' : 37, 'AATGTT' : 23} ,
        'NW' : { 'AACTGG' : 100, 'AATTGG' : 100} ,
        'NY' : { 'AACTAC' : 61, 'AACTAT' : 39, 'AATTAC' : 45, 'AATTAT' : 55} ,
        'PA' : { 'CCAGCA' : 26, 'CCAGCC' : 38, 'CCAGCG' : 8, 'CCAGCT' : 27, 'CCCGCA' : 14, 'CCCGCC' : 46, 'CCCGCG' : 27, 'CCCGCT' : 13, 'CCGGCA' : 14, 'CCGGCC' : 48, 'CCGGCG' : 21, 'CCGGCT' : 18, 'CCTGCA' : 25, 'CCTGCC' : 38, 'CCTGCG' : 9, 'CCTGCT' : 28} ,
        'PC' : { 'CCATGC' : 50, 'CCATGT' : 50, 'CCCTGC' : 66, 'CCCTGT' : 34, 'CCGTGC' : 64, 'CCGTGT' : 36, 'CCTTGC' : 49, 'CCTTGT' : 51} ,
        'PU' : { 'CCATGC' : 50, 'CCATGT' : 50, 'CCCTGC' : 66, 'CCCTGT' : 34, 'CCGTGC' : 64, 'CCGTGT' : 36, 'CCTTGC' : 49, 'CCTTGT' : 51} ,
        'PD' : { 'CCAGAC' : 48, 'CCAGAT' : 52, 'CCCGAC' : 73, 'CCCGAT' : 27, 'CCGGAC' : 66, 'CCGGAT' : 34, 'CCTGAC' : 50, 'CCTGAT' : 50} ,
        'PE' : { 'CCAGAA' : 46, 'CCAGAG' : 54, 'CCCGAA' : 20, 'CCCGAG' : 80, 'CCGGAA' : 26, 'CCGGAG' : 74, 'CCTGAA' : 42, 'CCTGAG' : 58} ,
        'PF' : { 'CCATTC' : 39, 'CCATTT' : 61, 'CCCTTC' : 66, 'CCCTTT' : 34, 'CCGTTC' : 60, 'CCGTTT' : 40, 'CCTTTC' : 52, 'CCTTTT' : 48} ,
        'PG' : { 'CCAGGA' : 29, 'CCAGGC' : 30, 'CCAGGG' : 26, 'CCAGGT' : 16, 'CCCGGA' : 15, 'CCCGGC' : 42, 'CCCGGG' : 35, 'CCCGGT' : 9, 'CCGGGA' : 18, 'CCGGGC' : 41, 'CCGGGG' : 30, 'CCGGGT' : 11, 'CCTGGA' : 27, 'CCTGGC' : 29, 'CCTGGG' : 28, 'CCTGGT' : 16} ,
        'PH' : { 'CCACAC' : 52, 'CCACAT' : 48, 'CCCCAC' : 66, 'CCCCAT' : 34, 'CCGCAC' : 77, 'CCGCAT' : 23, 'CCTCAC' : 50, 'CCTCAT' : 50} ,
        'PI' : { 'CCAATA' : 26, 'CCAATC' : 30, 'CCAATT' : 45, 'CCCATA' : 11, 'CCCATC' : 62, 'CCCATT' : 28, 'CCGATA' : 18, 'CCGATC' : 45, 'CCGATT' : 37, 'CCTATA' : 22, 'CCTATC' : 34, 'CCTATT' : 44} ,
        'PK' : { 'CCAAAA' : 44, 'CCAAAG' : 56, 'CCCAAA' : 37, 'CCCAAG' : 63, 'CCGAAA' : 34, 'CCGAAG' : 66, 'CCTAAA' : 53, 'CCTAAG' : 47} ,
        'PL' : { 'CCACTA' : 8, 'CCACTC' : 18, 'CCACTG' : 36, 'CCACTT' : 15, 'CCATTA' : 9, 'CCATTG' : 13, 'CCCCTA' : 5, 'CCCCTC' : 20, 'CCCCTG' : 43, 'CCCCTT' : 10, 'CCCTTA' : 5, 'CCCTTG' : 16, 'CCGCTA' : 4, 'CCGCTC' : 22, 'CCGCTG' : 57, 'CCGCTT' : 7, 'CCGTTA' : 3, 'CCGTTG' : 7, 'CCTCTA' : 8, 'CCTCTC' : 17, 'CCTCTG' : 32, 'CCTCTT' : 16, 'CCTTTA' : 10, 'CCTTTG' : 17} ,
        'PM' : { 'CCAATG' : 100, 'CCCATG' : 100, 'CCGATG' : 100, 'CCTATG' : 100} ,
        'PN' : { 'CCAAAC' : 41, 'CCAAAT' : 59, 'CCCAAC' : 65, 'CCCAAT' : 35, 'CCGAAC' : 53, 'CCGAAT' : 47, 'CCTAAC' : 42, 'CCTAAT' : 58} ,
        'PP' : { 'CCACCA' : 30, 'CCACCC' : 30, 'CCACCG' : 10, 'CCACCT' : 30, 'CCCCCA' : 34, 'CCCCCC' : 19, 'CCCCCG' : 21, 'CCCCCT' : 26, 'CCGCCA' : 16, 'CCGCCC' : 36, 'CCGCCG' : 31, 'CCGCCT' : 16, 'CCTCCA' : 32, 'CCTCCC' : 26, 'CCTCCG' : 9, 'CCTCCT' : 33} ,
        'PQ' : { 'CCACAA' : 28, 'CCACAG' : 72, 'CCCCAA' : 19, 'CCCCAG' : 81, 'CCGCAA' : 13, 'CCGCAG' : 87, 'CCTCAA' : 29, 'CCTCAG' : 71} ,
        'PR' : { 'CCAAGA' : 28, 'CCAAGG' : 24, 'CCACGA' : 10, 'CCACGC' : 12, 'CCACGG' : 18, 'CCACGT' : 7, 'CCCAGA' : 17, 'CCCAGG' : 28, 'CCCCGA' : 10, 'CCCCGC' : 17, 'CCCCGG' : 23, 'CCCCGT' : 5, 'CCGAGA' : 9, 'CCGAGG' : 18, 'CCGCGA' : 7, 'CCGCGC' : 30, 'CCGCGG' : 31, 'CCGCGT' : 5, 'CCTAGA' : 14, 'CCTAGG' : 11, 'CCTCGA' : 19, 'CCTCGC' : 16, 'CCTCGG' : 28, 'CCTCGT' : 11} ,
        'PS' : { 'CCAAGC' : 17, 'CCAAGT' : 16, 'CCATCA' : 19, 'CCATCC' : 21, 'CCATCG' : 5, 'CCATCT' : 22, 'CCCAGC' : 35, 'CCCAGT' : 15, 'CCCTCA' : 12, 'CCCTCC' : 18, 'CCCTCG' : 8, 'CCCTCT' : 12, 'CCGAGC' : 24, 'CCGAGT' : 11, 'CCGTCA' : 11, 'CCGTCC' : 27, 'CCGTCG' : 13, 'CCGTCT' : 13, 'CCTAGC' : 9, 'CCTAGT' : 9, 'CCTTCA' : 23, 'CCTTCC' : 26, 'CCTTCG' : 5, 'CCTTCT' : 28} ,
        'PT' : { 'CCAACA' : 35, 'CCAACC' : 29, 'CCAACG' : 8, 'CCAACT' : 28, 'CCCACA' : 24, 'CCCACC' : 41, 'CCCACG' : 16, 'CCCACT' : 19, 'CCGACA' : 22, 'CCGACC' : 39, 'CCGACG' : 18, 'CCGACT' : 21, 'CCTACA' : 33, 'CCTACC' : 30, 'CCTACG' : 8, 'CCTACT' : 29} ,
        'PV' : { 'CCAGTA' : 15, 'CCAGTC' : 21, 'CCAGTG' : 42, 'CCAGTT' : 23, 'CCCGTA' : 5, 'CCCGTC' : 26, 'CCCGTG' : 60, 'CCCGTT' : 9, 'CCGGTA' : 7, 'CCGGTC' : 24, 'CCGGTG' : 59, 'CCGGTT' : 11, 'CCTGTA' : 12, 'CCTGTC' : 24, 'CCTGTG' : 44, 'CCTGTT' : 20} ,
        'PW' : { 'CCATGG' : 100, 'CCCTGG' : 100, 'CCGTGG' : 100, 'CCTTGG' : 100} ,
        'PY' : { 'CCATAC' : 43, 'CCATAT' : 57, 'CCCTAC' : 62, 'CCCTAT' : 38, 'CCGTAC' : 56, 'CCGTAT' : 44, 'CCTTAC' : 47, 'CCTTAT' : 53} ,
        'QA' : { 'CAAGCA' : 30, 'CAAGCC' : 34, 'CAAGCG' : 6, 'CAAGCT' : 30, 'CAGGCA' : 21, 'CAGGCC' : 42, 'CAGGCG' : 13, 'CAGGCT' : 24} ,
        'QC' : { 'CAATGC' : 37, 'CAATGT' : 63, 'CAGTGC' : 53, 'CAGTGT' : 47} ,
        'QU' : { 'CAATGC' : 37, 'CAATGT' : 63, 'CAGTGC' : 53, 'CAGTGT' : 47} ,
        'QD' : { 'CAAGAC' : 43, 'CAAGAT' : 57, 'CAGGAC' : 59, 'CAGGAT' : 41} ,
        'QE' : { 'CAAGAA' : 57, 'CAAGAG' : 43, 'CAGGAA' : 35, 'CAGGAG' : 65} ,
        'QF' : { 'CAATTC' : 40, 'CAATTT' : 60, 'CAGTTC' : 51, 'CAGTTT' : 49} ,
        'QG' : { 'CAAGGA' : 34, 'CAAGGC' : 26, 'CAAGGG' : 21, 'CAAGGT' : 18, 'CAGGGA' : 21, 'CAGGGC' : 41, 'CAGGGG' : 21, 'CAGGGT' : 17} ,
        'QH' : { 'CAACAC' : 40, 'CAACAT' : 60, 'CAGCAC' : 61, 'CAGCAT' : 39} ,
        'QI' : { 'CAAATA' : 28, 'CAAATC' : 31, 'CAAATT' : 41, 'CAGATA' : 18, 'CAGATC' : 48, 'CAGATT' : 35} ,
        'QK' : { 'CAAAAA' : 49, 'CAAAAG' : 51, 'CAGAAA' : 41, 'CAGAAG' : 59} ,
        'QL' : { 'CAACTA' : 12, 'CAACTC' : 17, 'CAACTG' : 29, 'CAACTT' : 22, 'CAATTA' : 10, 'CAATTG' : 9, 'CAGCTA' : 7, 'CAGCTC' : 21, 'CAGCTG' : 44, 'CAGCTT' : 12, 'CAGTTA' : 6, 'CAGTTG' : 9} ,
        'QM' : { 'CAAATG' : 100, 'CAGATG' : 100} ,
        'QN' : { 'CAAAAC' : 41, 'CAAAAT' : 59, 'CAGAAC' : 54, 'CAGAAT' : 46} ,
        'QP' : { 'CAACCA' : 34, 'CAACCC' : 27, 'CAACCG' : 6, 'CAACCT' : 33, 'CAGCCA' : 24, 'CAGCCC' : 36, 'CAGCCG' : 13, 'CAGCCT' : 26} ,
        'QQ' : { 'CAACAA' : 38, 'CAACAG' : 62, 'CAGCAA' : 23, 'CAGCAG' : 77} ,
        'QR' : { 'CAAAGA' : 40, 'CAAAGG' : 28, 'CAACGA' : 9, 'CAACGC' : 8, 'CAACGG' : 9, 'CAACGT' : 7, 'CAGAGA' : 22, 'CAGAGG' : 25, 'CAGCGA' : 8, 'CAGCGC' : 18, 'CAGCGG' : 20, 'CAGCGT' : 7} ,
        'QS' : { 'CAAAGC' : 28, 'CAAAGT' : 27, 'CAATCA' : 16, 'CAATCC' : 12, 'CAATCG' : 2, 'CAATCT' : 15, 'CAGAGC' : 30, 'CAGAGT' : 17, 'CAGTCA' : 14, 'CAGTCC' : 18, 'CAGTCG' : 5, 'CAGTCT' : 16} ,
        'QT' : { 'CAAACA' : 36, 'CAAACC' : 26, 'CAAACG' : 7, 'CAAACT' : 31, 'CAGACA' : 28, 'CAGACC' : 36, 'CAGACG' : 14, 'CAGACT' : 23} ,
        'QV' : { 'CAAGTA' : 17, 'CAAGTC' : 21, 'CAAGTG' : 38, 'CAAGTT' : 24, 'CAGGTA' : 9, 'CAGGTC' : 23, 'CAGGTG' : 51, 'CAGGTT' : 16} ,
        'QW' : { 'CAATGG' : 100, 'CAGTGG' : 100} ,
        'QY' : { 'CAATAC' : 41, 'CAATAT' : 59, 'CAGTAC' : 56, 'CAGTAT' : 44} ,
        'RA' : { 'AGAGCA' : 29, 'AGAGCC' : 34, 'AGAGCG' : 6, 'AGAGCT' : 31, 'AGGGCA' : 24, 'AGGGCC' : 41, 'AGGGCG' : 10, 'AGGGCT' : 25, 'CGAGCA' : 25, 'CGAGCC' : 40, 'CGAGCG' : 9, 'CGAGCT' : 27, 'CGCGCA' : 9, 'CGCGCC' : 53, 'CGCGCG' : 27, 'CGCGCT' : 11, 'CGGGCA' : 21, 'CGGGCC' : 45, 'CGGGCG' : 16, 'CGGGCT' : 19, 'CGTGCA' : 20, 'CGTGCC' : 45, 'CGTGCG' : 9, 'CGTGCT' : 25} ,
        'RC' : { 'AGATGC' : 41, 'AGATGT' : 59, 'AGGTGC' : 54, 'AGGTGT' : 46, 'CGATGC' : 47, 'CGATGT' : 53, 'CGCTGC' : 68, 'CGCTGT' : 32, 'CGGTGC' : 62, 'CGGTGT' : 38, 'CGTTGC' : 51, 'CGTTGT' : 49} ,
        'RU' : { 'AGATGC' : 41, 'AGATGT' : 59, 'AGGTGC' : 54, 'AGGTGT' : 46, 'CGATGC' : 47, 'CGATGT' : 53, 'CGCTGC' : 68, 'CGCTGT' : 32, 'CGGTGC' : 62, 'CGGTGT' : 38, 'CGTTGC' : 51, 'CGTTGT' : 49} ,
        'RD' : { 'AGAGAC' : 46, 'AGAGAT' : 54, 'AGGGAC' : 58, 'AGGGAT' : 42, 'CGAGAC' : 53, 'CGAGAT' : 47, 'CGCGAC' : 77, 'CGCGAT' : 23, 'CGGGAC' : 67, 'CGGGAT' : 33, 'CGTGAC' : 56, 'CGTGAT' : 44} ,
        'RE' : { 'AGAGAA' : 55, 'AGAGAG' : 45, 'AGGGAA' : 40, 'AGGGAG' : 60, 'CGAGAA' : 44, 'CGAGAG' : 56, 'CGCGAA' : 12, 'CGCGAG' : 88, 'CGGGAA' : 28, 'CGGGAG' : 72, 'CGTGAA' : 40, 'CGTGAG' : 60} ,
        'RF' : { 'AGATTC' : 38, 'AGATTT' : 62, 'AGGTTC' : 49, 'AGGTTT' : 51, 'CGATTC' : 45, 'CGATTT' : 55, 'CGCTTC' : 72, 'CGCTTT' : 28, 'CGGTTC' : 56, 'CGGTTT' : 44, 'CGTTTC' : 53, 'CGTTTT' : 47} ,
        'RG' : { 'AGAGGA' : 36, 'AGAGGC' : 27, 'AGAGGG' : 19, 'AGAGGT' : 18, 'AGGGGA' : 25, 'AGGGGC' : 40, 'AGGGGG' : 17, 'AGGGGT' : 18, 'CGAGGA' : 26, 'CGAGGC' : 33, 'CGAGGG' : 25, 'CGAGGT' : 16, 'CGCGGA' : 10, 'CGCGGC' : 51, 'CGCGGG' : 31, 'CGCGGT' : 8, 'CGGGGA' : 17, 'CGGGGC' : 48, 'CGGGGG' : 19, 'CGGGGT' : 15, 'CGTGGA' : 21, 'CGTGGC' : 35, 'CGTGGG' : 26, 'CGTGGT' : 18} ,
        'RH' : { 'AGACAC' : 49, 'AGACAT' : 51, 'AGGCAC' : 58, 'AGGCAT' : 42, 'CGACAC' : 50, 'CGACAT' : 50, 'CGCCAC' : 70, 'CGCCAT' : 30, 'CGGCAC' : 70, 'CGGCAT' : 30, 'CGTCAC' : 50, 'CGTCAT' : 50} ,
        'RI' : { 'AGAATA' : 24, 'AGAATC' : 27, 'AGAATT' : 49, 'AGGATA' : 21, 'AGGATC' : 43, 'AGGATT' : 36, 'CGAATA' : 22, 'CGAATC' : 37, 'CGAATT' : 41, 'CGCATA' : 7, 'CGCATC' : 71, 'CGCATT' : 22, 'CGGATA' : 16, 'CGGATC' : 54, 'CGGATT' : 30, 'CGTATA' : 17, 'CGTATC' : 43, 'CGTATT' : 40} ,
        'RK' : { 'AGAAAA' : 56, 'AGAAAG' : 44, 'AGGAAA' : 44, 'AGGAAG' : 56, 'CGAAAA' : 49, 'CGAAAG' : 51, 'CGCAAA' : 29, 'CGCAAG' : 71, 'CGGAAA' : 36, 'CGGAAG' : 64, 'CGTAAA' : 51, 'CGTAAG' : 49} ,
        'RL' : { 'AGACTA' : 10, 'AGACTC' : 18, 'AGACTG' : 27, 'AGACTT' : 19, 'AGATTA' : 13, 'AGATTG' : 13, 'AGGCTA' : 8, 'AGGCTC' : 20, 'AGGCTG' : 42, 'AGGCTT' : 14, 'AGGTTA' : 7, 'AGGTTG' : 9, 'CGACTA' : 10, 'CGACTC' : 20, 'CGACTG' : 33, 'CGACTT' : 17, 'CGATTA' : 10, 'CGATTG' : 11, 'CGCCTA' : 5, 'CGCCTC' : 25, 'CGCCTG' : 49, 'CGCCTT' : 9, 'CGCTTA' : 3, 'CGCTTG' : 8, 'CGGCTA' : 6, 'CGGCTC' : 23, 'CGGCTG' : 53, 'CGGCTT' : 9, 'CGGTTA' : 3, 'CGGTTG' : 6, 'CGTCTA' : 8, 'CGTCTC' : 19, 'CGTCTG' : 34, 'CGTCTT' : 15, 'CGTTTA' : 10, 'CGTTTG' : 15} ,
        'RM' : { 'AGAATG' : 100, 'AGGATG' : 100, 'CGAATG' : 100, 'CGCATG' : 100, 'CGGATG' : 100, 'CGTATG' : 100} ,
        'RN' : { 'AGAAAC' : 43, 'AGAAAT' : 57, 'AGGAAC' : 54, 'AGGAAT' : 46, 'CGAAAC' : 45, 'CGAAAT' : 55, 'CGCAAC' : 70, 'CGCAAT' : 30, 'CGGAAC' : 62, 'CGGAAT' : 38, 'CGTAAC' : 49, 'CGTAAT' : 51} ,
        'RP' : { 'AGACCA' : 31, 'AGACCC' : 29, 'AGACCG' : 6, 'AGACCT' : 34, 'AGGCCA' : 26, 'AGGCCC' : 34, 'AGGCCG' : 12, 'AGGCCT' : 28, 'CGACCA' : 29, 'CGACCC' : 32, 'CGACCG' : 9, 'CGACCT' : 30, 'CGCCCA' : 23, 'CGCCCC' : 37, 'CGCCCG' : 19, 'CGCCCT' : 21, 'CGGCCA' : 18, 'CGGCCC' : 44, 'CGGCCG' : 18, 'CGGCCT' : 20, 'CGTCCA' : 28, 'CGTCCC' : 31, 'CGTCCG' : 9, 'CGTCCT' : 31} ,
        'RQ' : { 'AGACAA' : 42, 'AGACAG' : 58, 'AGGCAA' : 29, 'AGGCAG' : 71, 'CGACAA' : 31, 'CGACAG' : 69, 'CGCCAA' : 16, 'CGCCAG' : 84, 'CGGCAA' : 16, 'CGGCAG' : 84, 'CGTCAA' : 27, 'CGTCAG' : 73} ,
        'RR' : { 'AGAAGA' : 41, 'AGAAGG' : 25, 'AGACGA' : 9, 'AGACGC' : 9, 'AGACGG' : 10, 'AGACGT' : 6, 'AGGAGA' : 27, 'AGGAGG' : 29, 'AGGCGA' : 8, 'AGGCGC' : 13, 'AGGCGG' : 17, 'AGGCGT' : 5, 'CGAAGA' : 31, 'CGAAGG' : 25, 'CGACGA' : 9, 'CGACGC' : 13, 'CGACGG' : 15, 'CGACGT' : 7, 'CGCAGA' : 11, 'CGCAGG' : 21, 'CGCCGA' : 8, 'CGCCGC' : 31, 'CGCCGG' : 23, 'CGCCGT' : 7, 'CGGAGA' : 13, 'CGGAGG' : 21, 'CGGCGA' : 8, 'CGGCGC' : 25, 'CGGCGG' : 28, 'CGGCGT' : 6, 'CGTAGA' : 15, 'CGTAGG' : 13, 'CGTCGA' : 14, 'CGTCGC' : 22, 'CGTCGG' : 24, 'CGTCGT' : 12} ,
        'RS' : { 'AGAAGC' : 25, 'AGAAGT' : 25, 'AGATCA' : 17, 'AGATCC' : 13, 'AGATCG' : 3, 'AGATCT' : 17, 'AGGAGC' : 32, 'AGGAGT' : 17, 'AGGTCA' : 14, 'AGGTCC' : 17, 'AGGTCG' : 4, 'AGGTCT' : 15, 'CGAAGC' : 28, 'CGAAGT' : 24, 'CGATCA' : 14, 'CGATCC' : 15, 'CGATCG' : 3, 'CGATCT' : 16, 'CGCAGC' : 34, 'CGCAGT' : 10, 'CGCTCA' : 10, 'CGCTCC' : 26, 'CGCTCG' : 10, 'CGCTCT' : 11, 'CGGAGC' : 39, 'CGGAGT' : 16, 'CGGTCA' : 11, 'CGGTCC' : 18, 'CGGTCG' : 5, 'CGGTCT' : 11, 'CGTAGC' : 13, 'CGTAGT' : 10, 'CGTTCA' : 18, 'CGTTCC' : 27, 'CGTTCG' : 5, 'CGTTCT' : 27} ,
        'RT' : { 'AGAACA' : 36, 'AGAACC' : 24, 'AGAACG' : 7, 'AGAACT' : 33, 'AGGACA' : 29, 'AGGACC' : 33, 'AGGACG' : 12, 'AGGACT' : 26, 'CGAACA' : 32, 'CGAACC' : 30, 'CGAACG' : 8, 'CGAACT' : 30, 'CGCACA' : 22, 'CGCACC' : 42, 'CGCACG' : 21, 'CGCACT' : 14, 'CGGACA' : 25, 'CGGACC' : 39, 'CGGACG' : 16, 'CGGACT' : 20, 'CGTACA' : 26, 'CGTACC' : 36, 'CGTACG' : 8, 'CGTACT' : 30} ,
        'RV' : { 'AGAGTA' : 17, 'AGAGTC' : 23, 'AGAGTG' : 33, 'AGAGTT' : 27, 'AGGGTA' : 12, 'AGGGTC' : 32, 'AGGGTG' : 40, 'AGGGTT' : 17, 'CGAGTA' : 15, 'CGAGTC' : 24, 'CGAGTG' : 42, 'CGAGTT' : 19, 'CGCGTA' : 4, 'CGCGTC' : 27, 'CGCGTG' : 62, 'CGCGTT' : 7, 'CGGGTA' : 10, 'CGGGTC' : 28, 'CGGGTG' : 51, 'CGGGTT' : 12, 'CGTGTA' : 11, 'CGTGTC' : 24, 'CGTGTG' : 48, 'CGTGTT' : 17} ,
        'RW' : { 'AGATGG' : 100, 'AGGTGG' : 100, 'CGATGG' : 100, 'CGCTGG' : 100, 'CGGTGG' : 100, 'CGTTGG' : 100} ,
        'RY' : { 'AGATAC' : 54, 'AGATAT' : 46, 'AGGTAC' : 57, 'AGGTAT' : 43, 'CGATAC' : 53, 'CGATAT' : 47, 'CGCTAC' : 68, 'CGCTAT' : 32, 'CGGTAC' : 65, 'CGGTAT' : 35, 'CGTTAC' : 53, 'CGTTAT' : 47} ,
        'SA' : { 'AGCGCA' : 13, 'AGCGCC' : 53, 'AGCGCG' : 18, 'AGCGCT' : 16, 'AGTGCA' : 24, 'AGTGCC' : 40, 'AGTGCG' : 6, 'AGTGCT' : 31, 'TCAGCA' : 29, 'TCAGCC' : 35, 'TCAGCG' : 7, 'TCAGCT' : 29, 'TCCGCA' : 17, 'TCCGCC' : 45, 'TCCGCG' : 21, 'TCCGCT' : 17, 'TCGGCA' : 14, 'TCGGCC' : 49, 'TCGGCG' : 18, 'TCGGCT' : 19, 'TCTGCA' : 29, 'TCTGCC' : 34, 'TCTGCG' : 7, 'TCTGCT' : 31} ,
        'SC' : { 'AGCTGC' : 62, 'AGCTGT' : 38, 'AGTTGC' : 40, 'AGTTGT' : 60, 'TCATGC' : 50, 'TCATGT' : 50, 'TCCTGC' : 62, 'TCCTGT' : 38, 'TCGTGC' : 60, 'TCGTGT' : 40, 'TCTTGC' : 45, 'TCTTGT' : 55} ,
        'SU' : { 'AGCTGC' : 62, 'AGCTGT' : 38, 'AGTTGC' : 40, 'AGTTGT' : 60, 'TCATGC' : 50, 'TCATGT' : 50, 'TCCTGC' : 62, 'TCCTGT' : 38, 'TCGTGC' : 60, 'TCGTGT' : 40, 'TCTTGC' : 45, 'TCTTGT' : 55} ,
        'SD' : { 'AGCGAC' : 69, 'AGCGAT' : 31, 'AGTGAC' : 52, 'AGTGAT' : 48, 'TCAGAC' : 45, 'TCAGAT' : 55, 'TCCGAC' : 67, 'TCCGAT' : 33, 'TCGGAC' : 63, 'TCGGAT' : 37, 'TCTGAC' : 46, 'TCTGAT' : 54} ,
        'SE' : { 'AGCGAA' : 22, 'AGCGAG' : 78, 'AGTGAA' : 52, 'AGTGAG' : 48, 'TCAGAA' : 51, 'TCAGAG' : 49, 'TCCGAA' : 26, 'TCCGAG' : 74, 'TCGGAA' : 29, 'TCGGAG' : 71, 'TCTGAA' : 49, 'TCTGAG' : 51} ,
        'SF' : { 'AGCTTC' : 63, 'AGCTTT' : 37, 'AGTTTC' : 45, 'AGTTTT' : 55, 'TCATTC' : 36, 'TCATTT' : 64, 'TCCTTC' : 63, 'TCCTTT' : 37, 'TCGTTC' : 50, 'TCGTTT' : 50, 'TCTTTC' : 48, 'TCTTTT' : 52} ,
        'SG' : { 'AGCGGA' : 13, 'AGCGGC' : 47, 'AGCGGG' : 29, 'AGCGGT' : 11, 'AGTGGA' : 27, 'AGTGGC' : 30, 'AGTGGG' : 25, 'AGTGGT' : 18, 'TCAGGA' : 32, 'TCAGGC' : 25, 'TCAGGG' : 25, 'TCAGGT' : 17, 'TCCGGA' : 18, 'TCCGGC' : 36, 'TCCGGG' : 35, 'TCCGGT' : 12, 'TCGGGA' : 18, 'TCGGGC' : 41, 'TCGGGG' : 30, 'TCGGGT' : 11, 'TCTGGA' : 31, 'TCTGGC' : 25, 'TCTGGG' : 27, 'TCTGGT' : 17} ,
        'SH' : { 'AGCCAC' : 66, 'AGCCAT' : 34, 'AGTCAC' : 50, 'AGTCAT' : 50, 'TCACAC' : 51, 'TCACAT' : 49, 'TCCCAC' : 62, 'TCCCAT' : 38, 'TCGCAC' : 74, 'TCGCAT' : 26, 'TCTCAC' : 46, 'TCTCAT' : 54} ,
        'SI' : { 'AGCATA' : 13, 'AGCATC' : 56, 'AGCATT' : 31, 'AGTATA' : 22, 'AGTATC' : 33, 'AGTATT' : 45, 'TCAATA' : 30, 'TCAATC' : 27, 'TCAATT' : 43, 'TCCATA' : 13, 'TCCATC' : 57, 'TCCATT' : 30, 'TCGATA' : 20, 'TCGATC' : 40, 'TCGATT' : 40, 'TCTATA' : 24, 'TCTATC' : 33, 'TCTATT' : 43} ,
        'SK' : { 'AGCAAA' : 43, 'AGCAAG' : 57, 'AGTAAA' : 62, 'AGTAAG' : 38, 'TCAAAA' : 52, 'TCAAAG' : 48, 'TCCAAA' : 41, 'TCCAAG' : 59, 'TCGAAA' : 40, 'TCGAAG' : 60, 'TCTAAA' : 53, 'TCTAAG' : 47} ,
        'SL' : { 'AGCCTA' : 5, 'AGCCTC' : 25, 'AGCCTG' : 42, 'AGCCTT' : 12, 'AGCTTA' : 5, 'AGCTTG' : 11, 'AGTCTA' : 8, 'AGTCTC' : 16, 'AGTCTG' : 25, 'AGTCTT' : 17, 'AGTTTA' : 15, 'AGTTTG' : 19, 'TCACTA' : 8, 'TCACTC' : 15, 'TCACTG' : 32, 'TCACTT' : 16, 'TCATTA' : 13, 'TCATTG' : 16, 'TCCCTA' : 6, 'TCCCTC' : 19, 'TCCCTG' : 41, 'TCCCTT' : 10, 'TCCTTA' : 7, 'TCCTTG' : 16, 'TCGCTA' : 4, 'TCGCTC' : 20, 'TCGCTG' : 56, 'TCGCTT' : 8, 'TCGTTA' : 4, 'TCGTTG' : 8, 'TCTCTA' : 8, 'TCTCTC' : 16, 'TCTCTG' : 30, 'TCTCTT' : 16, 'TCTTTA' : 12, 'TCTTTG' : 18} ,
        'SM' : { 'AGCATG' : 100, 'AGTATG' : 100, 'TCAATG' : 100, 'TCCATG' : 100, 'TCGATG' : 100, 'TCTATG' : 100} ,
        'SN' : { 'AGCAAC' : 60, 'AGCAAT' : 40, 'AGTAAC' : 47, 'AGTAAT' : 53, 'TCAAAC' : 38, 'TCAAAT' : 62, 'TCCAAC' : 62, 'TCCAAT' : 38, 'TCGAAC' : 50, 'TCGAAT' : 50, 'TCTAAC' : 42, 'TCTAAT' : 58} ,
        'SP' : { 'AGCCCA' : 22, 'AGCCCC' : 39, 'AGCCCG' : 13, 'AGCCCT' : 26, 'AGTCCA' : 28, 'AGTCCC' : 30, 'AGTCCG' : 6, 'AGTCCT' : 37, 'TCACCA' : 34, 'TCACCC' : 28, 'TCACCG' : 7, 'TCACCT' : 31, 'TCCCCA' : 33, 'TCCCCC' : 23, 'TCCCCG' : 18, 'TCCCCT' : 25, 'TCGCCA' : 19, 'TCGCCC' : 42, 'TCGCCG' : 20, 'TCGCCT' : 19, 'TCTCCA' : 35, 'TCTCCC' : 24, 'TCTCCG' : 7, 'TCTCCT' : 34} ,
        'SQ' : { 'AGCCAA' : 24, 'AGCCAG' : 76, 'AGTCAA' : 37, 'AGTCAG' : 63, 'TCACAA' : 36, 'TCACAG' : 64, 'TCCCAA' : 22, 'TCCCAG' : 78, 'TCGCAA' : 16, 'TCGCAG' : 84, 'TCTCAA' : 33, 'TCTCAG' : 67} ,
        'SR' : { 'AGCAGA' : 18, 'AGCAGG' : 21, 'AGCCGA' : 9, 'AGCCGC' : 24, 'AGCCGG' : 20, 'AGCCGT' : 7, 'AGTAGA' : 21, 'AGTAGG' : 13, 'AGTCGA' : 16, 'AGTCGC' : 17, 'AGTCGG' : 19, 'AGTCGT' : 13, 'TCAAGA' : 36, 'TCAAGG' : 24, 'TCACGA' : 9, 'TCACGC' : 10, 'TCACGG' : 13, 'TCACGT' : 8, 'TCCAGA' : 21, 'TCCAGG' : 26, 'TCCCGA' : 10, 'TCCCGC' : 17, 'TCCCGG' : 21, 'TCCCGT' : 6, 'TCGAGA' : 13, 'TCGAGG' : 21, 'TCGCGA' : 6, 'TCGCGC' : 26, 'TCGCGG' : 29, 'TCGCGT' : 5, 'TCTAGA' : 19, 'TCTAGG' : 15, 'TCTCGA' : 18, 'TCTCGC' : 14, 'TCTCGG' : 24, 'TCTCGT' : 11} ,
        'SS' : { 'AGCAGC' : 35, 'AGCAGT' : 17, 'AGCTCA' : 10, 'AGCTCC' : 20, 'AGCTCG' : 5, 'AGCTCT' : 12, 'AGTAGC' : 16, 'AGTAGT' : 14, 'AGTTCA' : 21, 'AGTTCC' : 20, 'AGTTCG' : 3, 'AGTTCT' : 26, 'TCAAGC' : 13, 'TCAAGT' : 16, 'TCATCA' : 21, 'TCATCC' : 20, 'TCATCG' : 4, 'TCATCT' : 26, 'TCCAGC' : 25, 'TCCAGT' : 14, 'TCCTCA' : 15, 'TCCTCC' : 23, 'TCCTCG' : 7, 'TCCTCT' : 15, 'TCGAGC' : 16, 'TCGAGT' : 10, 'TCGTCA' : 12, 'TCGTCC' : 31, 'TCGTCG' : 14, 'TCGTCT' : 16, 'TCTAGC' : 9, 'TCTAGT' : 9, 'TCTTCA' : 25, 'TCTTCC' : 24, 'TCTTCG' : 5, 'TCTTCT' : 28} ,
        'ST' : { 'AGCACA' : 28, 'AGCACC' : 40, 'AGCACG' : 11, 'AGCACT' : 22, 'AGTACA' : 34, 'AGTACC' : 29, 'AGTACG' : 6, 'AGTACT' : 32, 'TCAACA' : 36, 'TCAACC' : 24, 'TCAACG' : 6, 'TCAACT' : 34, 'TCCACA' : 27, 'TCCACC' : 37, 'TCCACG' : 15, 'TCCACT' : 21, 'TCGACA' : 24, 'TCGACC' : 34, 'TCGACG' : 15, 'TCGACT' : 28, 'TCTACA' : 34, 'TCTACC' : 28, 'TCTACG' : 6, 'TCTACT' : 32} ,
        'SV' : { 'AGCGTA' : 5, 'AGCGTC' : 32, 'AGCGTG' : 52, 'AGCGTT' : 10, 'AGTGTA' : 13, 'AGTGTC' : 26, 'AGTGTG' : 38, 'AGTGTT' : 23, 'TCAGTA' : 17, 'TCAGTC' : 19, 'TCAGTG' : 42, 'TCAGTT' : 22, 'TCCGTA' : 7, 'TCCGTC' : 25, 'TCCGTG' : 58, 'TCCGTT' : 10, 'TCGGTA' : 7, 'TCGGTC' : 22, 'TCGGTG' : 61, 'TCGGTT' : 10, 'TCTGTA' : 13, 'TCTGTC' : 22, 'TCTGTG' : 44, 'TCTGTT' : 21} ,
        'SW' : { 'AGCTGG' : 100, 'AGTTGG' : 100, 'TCATGG' : 100, 'TCCTGG' : 100, 'TCGTGG' : 100, 'TCTTGG' : 100} ,
        'SY' : { 'AGCTAC' : 63, 'AGCTAT' : 37, 'AGTTAC' : 52, 'AGTTAT' : 48, 'TCATAC' : 43, 'TCATAT' : 57, 'TCCTAC' : 61, 'TCCTAT' : 39, 'TCGTAC' : 60, 'TCGTAT' : 40, 'TCTTAC' : 51, 'TCTTAT' : 49} ,
        'TA' : { 'ACAGCA' : 28, 'ACAGCC' : 36, 'ACAGCG' : 7, 'ACAGCT' : 28, 'ACCGCA' : 16, 'ACCGCC' : 50, 'ACCGCG' : 17, 'ACCGCT' : 17, 'ACGGCA' : 18, 'ACGGCC' : 48, 'ACGGCG' : 16, 'ACGGCT' : 18, 'ACTGCA' : 29, 'ACTGCC' : 33, 'ACTGCG' : 7, 'ACTGCT' : 31} ,
        'TC' : { 'ACATGC' : 44, 'ACATGT' : 56, 'ACCTGC' : 62, 'ACCTGT' : 38, 'ACGTGC' : 56, 'ACGTGT' : 44, 'ACTTGC' : 44, 'ACTTGT' : 56} ,
        'TU' : { 'ACATGC' : 44, 'ACATGT' : 56, 'ACCTGC' : 62, 'ACCTGT' : 38, 'ACGTGC' : 56, 'ACGTGT' : 44, 'ACTTGC' : 44, 'ACTTGT' : 56} ,
        'TD' : { 'ACAGAC' : 47, 'ACAGAT' : 53, 'ACCGAC' : 70, 'ACCGAT' : 30, 'ACGGAC' : 62, 'ACGGAT' : 38, 'ACTGAC' : 48, 'ACTGAT' : 52} ,
        'TE' : { 'ACAGAA' : 50, 'ACAGAG' : 50, 'ACCGAA' : 25, 'ACCGAG' : 75, 'ACGGAA' : 30, 'ACGGAG' : 70, 'ACTGAA' : 51, 'ACTGAG' : 49} ,
        'TF' : { 'ACATTC' : 38, 'ACATTT' : 62, 'ACCTTC' : 66, 'ACCTTT' : 34, 'ACGTTC' : 50, 'ACGTTT' : 50, 'ACTTTC' : 48, 'ACTTTT' : 52} ,
        'TG' : { 'ACAGGA' : 31, 'ACAGGC' : 27, 'ACAGGG' : 22, 'ACAGGT' : 20, 'ACCGGA' : 17, 'ACCGGC' : 42, 'ACCGGG' : 30, 'ACCGGT' : 11, 'ACGGGA' : 18, 'ACGGGC' : 42, 'ACGGGG' : 27, 'ACGGGT' : 13, 'ACTGGA' : 38, 'ACTGGC' : 24, 'ACTGGG' : 22, 'ACTGGT' : 17} ,
        'TH' : { 'ACACAC' : 52, 'ACACAT' : 48, 'ACCCAC' : 65, 'ACCCAT' : 35, 'ACGCAC' : 73, 'ACGCAT' : 27, 'ACTCAC' : 53, 'ACTCAT' : 47} ,
        'TI' : { 'ACAATA' : 28, 'ACAATC' : 29, 'ACAATT' : 43, 'ACCATA' : 12, 'ACCATC' : 60, 'ACCATT' : 28, 'ACGATA' : 20, 'ACGATC' : 43, 'ACGATT' : 36, 'ACTATA' : 22, 'ACTATC' : 32, 'ACTATT' : 46} ,
        'TK' : { 'ACAAAA' : 51, 'ACAAAG' : 49, 'ACCAAA' : 37, 'ACCAAG' : 63, 'ACGAAA' : 40, 'ACGAAG' : 60, 'ACTAAA' : 59, 'ACTAAG' : 41} ,
        'TL' : { 'ACACTA' : 9, 'ACACTC' : 15, 'ACACTG' : 34, 'ACACTT' : 15, 'ACATTA' : 13, 'ACATTG' : 14, 'ACCCTA' : 5, 'ACCCTC' : 21, 'ACCCTG' : 41, 'ACCCTT' : 11, 'ACCTTA' : 6, 'ACCTTG' : 15, 'ACGCTA' : 5, 'ACGCTC' : 18, 'ACGCTG' : 56, 'ACGCTT' : 7, 'ACGTTA' : 5, 'ACGTTG' : 9, 'ACTCTA' : 9, 'ACTCTC' : 15, 'ACTCTG' : 30, 'ACTCTT' : 15, 'ACTTTA' : 13, 'ACTTTG' : 19} ,
        'TM' : { 'ACAATG' : 100, 'ACCATG' : 100, 'ACGATG' : 100, 'ACTATG' : 100} ,
        'TN' : { 'ACAAAC' : 40, 'ACAAAT' : 60, 'ACCAAC' : 63, 'ACCAAT' : 37, 'ACGAAC' : 51, 'ACGAAT' : 49, 'ACTAAC' : 43, 'ACTAAT' : 57} ,
        'TP' : { 'ACACCA' : 32, 'ACACCC' : 29, 'ACACCG' : 8, 'ACACCT' : 32, 'ACCCCA' : 33, 'ACCCCC' : 26, 'ACCCCG' : 14, 'ACCCCT' : 27, 'ACGCCA' : 20, 'ACGCCC' : 42, 'ACGCCG' : 19, 'ACGCCT' : 20, 'ACTCCA' : 35, 'ACTCCC' : 25, 'ACTCCG' : 7, 'ACTCCT' : 33} ,
        'TQ' : { 'ACACAA' : 33, 'ACACAG' : 67, 'ACCCAA' : 22, 'ACCCAG' : 78, 'ACGCAA' : 16, 'ACGCAG' : 84, 'ACTCAA' : 34, 'ACTCAG' : 66} ,
        'TR' : { 'ACAAGA' : 32, 'ACAAGG' : 24, 'ACACGA' : 10, 'ACACGC' : 11, 'ACACGG' : 16, 'ACACGT' : 8, 'ACCAGA' : 20, 'ACCAGG' : 26, 'ACCCGA' : 9, 'ACCCGC' : 18, 'ACCCGG' : 20, 'ACCCGT' : 7, 'ACGAGA' : 14, 'ACGAGG' : 19, 'ACGCGA' : 7, 'ACGCGC' : 26, 'ACGCGG' : 29, 'ACGCGT' : 5, 'ACTAGA' : 20, 'ACTAGG' : 12, 'ACTCGA' : 18, 'ACTCGC' : 16, 'ACTCGG' : 22, 'ACTCGT' : 12} ,
        'TS' : { 'ACAAGC' : 15, 'ACAAGT' : 17, 'ACATCA' : 21, 'ACATCC' : 18, 'ACATCG' : 4, 'ACATCT' : 24, 'ACCAGC' : 30, 'ACCAGT' : 14, 'ACCTCA' : 14, 'ACCTCC' : 20, 'ACCTCG' : 7, 'ACCTCT' : 14, 'ACGAGC' : 17, 'ACGAGT' : 10, 'ACGTCA' : 16, 'ACGTCC' : 28, 'ACGTCG' : 12, 'ACGTCT' : 17, 'ACTAGC' : 9, 'ACTAGT' : 10, 'ACTTCA' : 25, 'ACTTCC' : 23, 'ACTTCG' : 5, 'ACTTCT' : 29} ,
        'TT' : { 'ACAACA' : 35, 'ACAACC' : 28, 'ACAACG' : 7, 'ACAACT' : 31, 'ACCACA' : 25, 'ACCACC' : 42, 'ACCACG' : 13, 'ACCACT' : 20, 'ACGACA' : 26, 'ACGACC' : 35, 'ACGACG' : 16, 'ACGACT' : 23, 'ACTACA' : 36, 'ACTACC' : 27, 'ACTACG' : 8, 'ACTACT' : 29} ,
        'TV' : { 'ACAGTA' : 15, 'ACAGTC' : 21, 'ACAGTG' : 43, 'ACAGTT' : 21, 'ACCGTA' : 6, 'ACCGTC' : 28, 'ACCGTG' : 56, 'ACCGTT' : 10, 'ACGGTA' : 7, 'ACGGTC' : 23, 'ACGGTG' : 59, 'ACGGTT' : 11, 'ACTGTA' : 14, 'ACTGTC' : 23, 'ACTGTG' : 43, 'ACTGTT' : 20} ,
        'TW' : { 'ACATGG' : 100, 'ACCTGG' : 100, 'ACGTGG' : 100, 'ACTTGG' : 100} ,
        'TY' : { 'ACATAC' : 49, 'ACATAT' : 51, 'ACCTAC' : 67, 'ACCTAT' : 33, 'ACGTAC' : 60, 'ACGTAT' : 40, 'ACTTAC' : 57, 'ACTTAT' : 43} ,
        'VA' : { 'GTAGCA' : 30, 'GTAGCC' : 32, 'GTAGCG' : 5, 'GTAGCT' : 33, 'GTCGCA' : 16, 'GTCGCC' : 48, 'GTCGCG' : 15, 'GTCGCT' : 21, 'GTGGCA' : 19, 'GTGGCC' : 46, 'GTGGCG' : 10, 'GTGGCT' : 26, 'GTTGCA' : 28, 'GTTGCC' : 30, 'GTTGCG' : 4, 'GTTGCT' : 37} ,
        'VC' : { 'GTATGC' : 39, 'GTATGT' : 61, 'GTCTGC' : 59, 'GTCTGT' : 41, 'GTGTGC' : 57, 'GTGTGT' : 43, 'GTTTGC' : 44, 'GTTTGT' : 56} ,
        'VU' : { 'GTATGC' : 39, 'GTATGT' : 61, 'GTCTGC' : 59, 'GTCTGT' : 41, 'GTGTGC' : 57, 'GTGTGT' : 43, 'GTTTGC' : 44, 'GTTTGT' : 56} ,
        'VD' : { 'GTAGAC' : 41, 'GTAGAT' : 59, 'GTCGAC' : 59, 'GTCGAT' : 41, 'GTGGAC' : 58, 'GTGGAT' : 42, 'GTTGAC' : 43, 'GTTGAT' : 57} ,
        'VE' : { 'GTAGAA' : 58, 'GTAGAG' : 42, 'GTCGAA' : 33, 'GTCGAG' : 67, 'GTGGAA' : 35, 'GTGGAG' : 65, 'GTTGAA' : 60, 'GTTGAG' : 40} ,
        'VF' : { 'GTATTC' : 32, 'GTATTT' : 68, 'GTCTTC' : 64, 'GTCTTT' : 36, 'GTGTTC' : 50, 'GTGTTT' : 50, 'GTTTTC' : 49, 'GTTTTT' : 51} ,
        'VG' : { 'GTAGGA' : 36, 'GTAGGC' : 22, 'GTAGGG' : 22, 'GTAGGT' : 19, 'GTCGGA' : 17, 'GTCGGC' : 34, 'GTCGGG' : 35, 'GTCGGT' : 14, 'GTGGGA' : 19, 'GTGGGC' : 42, 'GTGGGG' : 24, 'GTGGGT' : 15, 'GTTGGA' : 32, 'GTTGGC' : 25, 'GTTGGG' : 21, 'GTTGGT' : 22} ,
        'VH' : { 'GTACAC' : 43, 'GTACAT' : 57, 'GTCCAC' : 62, 'GTCCAT' : 38, 'GTGCAC' : 67, 'GTGCAT' : 33, 'GTTCAC' : 46, 'GTTCAT' : 54} ,
        'VI' : { 'GTAATA' : 27, 'GTAATC' : 24, 'GTAATT' : 49, 'GTCATA' : 11, 'GTCATC' : 58, 'GTCATT' : 31, 'GTGATA' : 16, 'GTGATC' : 45, 'GTGATT' : 38, 'GTTATA' : 21, 'GTTATC' : 31, 'GTTATT' : 48} ,
        'VK' : { 'GTAAAA' : 54, 'GTAAAG' : 46, 'GTCAAA' : 38, 'GTCAAG' : 62, 'GTGAAA' : 38, 'GTGAAG' : 62, 'GTTAAA' : 61, 'GTTAAG' : 39} ,
        'VL' : { 'GTACTA' : 10, 'GTACTC' : 12, 'GTACTG' : 26, 'GTACTT' : 18, 'GTATTA' : 16, 'GTATTG' : 18, 'GTCCTA' : 7, 'GTCCTC' : 22, 'GTCCTG' : 40, 'GTCCTT' : 12, 'GTCTTA' : 6, 'GTCTTG' : 13, 'GTGCTA' : 6, 'GTGCTC' : 19, 'GTGCTG' : 50, 'GTGCTT' : 10, 'GTGTTA' : 5, 'GTGTTG' : 10, 'GTTCTA' : 9, 'GTTCTC' : 16, 'GTTCTG' : 25, 'GTTCTT' : 19, 'GTTTTA' : 14, 'GTTTTG' : 18} ,
        'VM' : { 'GTAATG' : 100, 'GTCATG' : 100, 'GTGATG' : 100, 'GTTATG' : 100} ,
        'VN' : { 'GTAAAC' : 36, 'GTAAAT' : 64, 'GTCAAC' : 61, 'GTCAAT' : 39, 'GTGAAC' : 54, 'GTGAAT' : 46, 'GTTAAC' : 40, 'GTTAAT' : 60} ,
        'VP' : { 'GTACCA' : 32, 'GTACCC' : 27, 'GTACCG' : 6, 'GTACCT' : 35, 'GTCCCA' : 25, 'GTCCCC' : 35, 'GTCCCG' : 11, 'GTCCCT' : 29, 'GTGCCA' : 21, 'GTGCCC' : 43, 'GTGCCG' : 11, 'GTGCCT' : 26, 'GTTCCA' : 32, 'GTTCCC' : 26, 'GTTCCG' : 5, 'GTTCCT' : 38} ,
        'VQ' : { 'GTACAA' : 38, 'GTACAG' : 62, 'GTCCAA' : 22, 'GTCCAG' : 78, 'GTGCAA' : 18, 'GTGCAG' : 82, 'GTTCAA' : 37, 'GTTCAG' : 63} ,
        'VR' : { 'GTAAGA' : 30, 'GTAAGG' : 19, 'GTACGA' : 14, 'GTACGC' : 12, 'GTACGG' : 16, 'GTACGT' : 10, 'GTCAGA' : 20, 'GTCAGG' : 20, 'GTCCGA' : 12, 'GTCCGC' : 20, 'GTCCGG' : 20, 'GTCCGT' : 9, 'GTGAGA' : 16, 'GTGAGG' : 20, 'GTGCGA' : 8, 'GTGCGC' : 24, 'GTGCGG' : 25, 'GTGCGT' : 7, 'GTTAGA' : 16, 'GTTAGG' : 10, 'GTTCGA' : 23, 'GTTCGC' : 16, 'GTTCGG' : 18, 'GTTCGT' : 17} ,
        'VS' : { 'GTAAGC' : 14, 'GTAAGT' : 12, 'GTATCA' : 21, 'GTATCC' : 20, 'GTATCG' : 3, 'GTATCT' : 29, 'GTCAGC' : 28, 'GTCAGT' : 15, 'GTCTCA' : 11, 'GTCTCC' : 25, 'GTCTCG' : 4, 'GTCTCT' : 16, 'GTGAGC' : 22, 'GTGAGT' : 10, 'GTGTCA' : 13, 'GTGTCC' : 27, 'GTGTCG' : 6, 'GTGTCT' : 21, 'GTTAGC' : 8, 'GTTAGT' : 9, 'GTTTCA' : 22, 'GTTTCC' : 24, 'GTTTCG' : 3, 'GTTTCT' : 33} ,
        'VT' : { 'GTAACA' : 36, 'GTAACC' : 24, 'GTAACG' : 6, 'GTAACT' : 34, 'GTCACA' : 22, 'GTCACC' : 45, 'GTCACG' : 9, 'GTCACT' : 24, 'GTGACA' : 26, 'GTGACC' : 38, 'GTGACG' : 11, 'GTGACT' : 25, 'GTTACA' : 32, 'GTTACC' : 28, 'GTTACG' : 5, 'GTTACT' : 35} ,
        'VV' : { 'GTAGTA' : 22, 'GTAGTC' : 17, 'GTAGTG' : 36, 'GTAGTT' : 26, 'GTCGTA' : 8, 'GTCGTC' : 33, 'GTCGTG' : 47, 'GTCGTT' : 12, 'GTGGTA' : 10, 'GTGGTC' : 25, 'GTGGTG' : 51, 'GTGGTT' : 15, 'GTTGTA' : 17, 'GTTGTC' : 23, 'GTTGTG' : 34, 'GTTGTT' : 26} ,
        'VW' : { 'GTATGG' : 100, 'GTCTGG' : 100, 'GTGTGG' : 100, 'GTTTGG' : 100} ,
        'VY' : { 'GTATAC' : 38, 'GTATAT' : 62, 'GTCTAC' : 62, 'GTCTAT' : 38, 'GTGTAC' : 59, 'GTGTAT' : 41, 'GTTTAC' : 42, 'GTTTAT' : 58} ,
        'WA' : { 'TGGGCA' : 24, 'TGGGCC' : 40, 'TGGGCG' : 11, 'TGGGCT' : 26} ,
        'WC' : { 'TGGTGC' : 56, 'TGGTGT' : 44} ,
        'WU' : { 'TGGTGC' : 56, 'TGGTGT' : 44} ,
        'WD' : { 'TGGGAC' : 55, 'TGGGAT' : 45} ,
        'WE' : { 'TGGGAA' : 42, 'TGGGAG' : 58} ,
        'WF' : { 'TGGTTC' : 52, 'TGGTTT' : 48} ,
        'WG' : { 'TGGGGA' : 25, 'TGGGGC' : 37, 'TGGGGG' : 20, 'TGGGGT' : 18} ,
        'WH' : { 'TGGCAC' : 57, 'TGGCAT' : 43} ,
        'WI' : { 'TGGATA' : 18, 'TGGATC' : 45, 'TGGATT' : 38} ,
        'WK' : { 'TGGAAA' : 42, 'TGGAAG' : 58} ,
        'WL' : { 'TGGCTA' : 8, 'TGGCTC' : 19, 'TGGCTG' : 44, 'TGGCTT' : 14, 'TGGTTA' : 6, 'TGGTTG' : 9} ,
        'WM' : { 'TGGATG' : 100} ,
        'WN' : { 'TGGAAC' : 51, 'TGGAAT' : 49} ,
        'WP' : { 'TGGCCA' : 28, 'TGGCCC' : 34, 'TGGCCG' : 12, 'TGGCCT' : 26} ,
        'WQ' : { 'TGGCAA' : 27, 'TGGCAG' : 73} ,
        'WR' : { 'TGGAGA' : 25, 'TGGAGG' : 28, 'TGGCGA' : 8, 'TGGCGC' : 15, 'TGGCGG' : 17, 'TGGCGT' : 7} ,
        'WS' : { 'TGGAGC' : 29, 'TGGAGT' : 18, 'TGGTCA' : 14, 'TGGTCC' : 19, 'TGGTCG' : 5, 'TGGTCT' : 16} ,
        'WT' : { 'TGGACA' : 30, 'TGGACC' : 33, 'TGGACG' : 13, 'TGGACT' : 24} ,
        'WV' : { 'TGGGTA' : 12, 'TGGGTC' : 25, 'TGGGTG' : 46, 'TGGGTT' : 17} ,
        'WW' : { 'TGGTGG' : 100} ,
        'WY' : { 'TGGTAC' : 57, 'TGGTAT' : 43} ,
        'YA' : { 'TACGCA' : 17, 'TACGCC' : 46, 'TACGCG' : 21, 'TACGCT' : 15, 'TATGCA' : 29, 'TATGCC' : 35, 'TATGCG' : 6, 'TATGCT' : 30} ,
        'YC' : { 'TACTGC' : 59, 'TACTGT' : 41, 'TATTGC' : 49, 'TATTGT' : 51} ,
        'YU' : { 'TACTGC' : 59, 'TACTGT' : 41, 'TATTGC' : 49, 'TATTGT' : 51} ,
        'YD' : { 'TACGAC' : 66, 'TACGAT' : 34, 'TATGAC' : 50, 'TATGAT' : 50} ,
        'YE' : { 'TACGAA' : 27, 'TACGAG' : 73, 'TATGAA' : 54, 'TATGAG' : 46} ,
        'YF' : { 'TACTTC' : 62, 'TACTTT' : 38, 'TATTTC' : 45, 'TATTTT' : 55} ,
        'YG' : { 'TACGGA' : 18, 'TACGGC' : 40, 'TACGGG' : 30, 'TACGGT' : 12, 'TATGGA' : 31, 'TATGGC' : 29, 'TATGGG' : 22, 'TATGGT' : 18} ,
        'YH' : { 'TACCAC' : 63, 'TACCAT' : 37, 'TATCAC' : 50, 'TATCAT' : 50} ,
        'YI' : { 'TACATA' : 13, 'TACATC' : 56, 'TACATT' : 31, 'TATATA' : 20, 'TATATC' : 36, 'TATATT' : 45} ,
        'YK' : { 'TACAAA' : 44, 'TACAAG' : 56, 'TATAAA' : 57, 'TATAAG' : 43} ,
        'YL' : { 'TACCTA' : 7, 'TACCTC' : 19, 'TACCTG' : 45, 'TACCTT' : 11, 'TACTTA' : 6, 'TACTTG' : 12, 'TATCTA' : 8, 'TATCTC' : 16, 'TATCTG' : 25, 'TATCTT' : 17, 'TATTTA' : 15, 'TATTTG' : 20} ,
        'YM' : { 'TACATG' : 100, 'TATATG' : 100} ,
        'YN' : { 'TACAAC' : 61, 'TACAAT' : 39, 'TATAAC' : 43, 'TATAAT' : 57} ,
        'YP' : { 'TACCCA' : 30, 'TACCCC' : 30, 'TACCCG' : 15, 'TACCCT' : 24, 'TATCCA' : 36, 'TATCCC' : 24, 'TATCCG' : 6, 'TATCCT' : 35} ,
        'YQ' : { 'TACCAA' : 24, 'TACCAG' : 76, 'TATCAA' : 35, 'TATCAG' : 65} ,
        'YR' : { 'TACAGA' : 22, 'TACAGG' : 19, 'TACCGA' : 11, 'TACCGC' : 23, 'TACCGG' : 19, 'TACCGT' : 7, 'TATAGA' : 25, 'TATAGG' : 17, 'TATCGA' : 17, 'TATCGC' : 14, 'TATCGG' : 16, 'TATCGT' : 11} ,
        'YS' : { 'TACAGC' : 33, 'TACAGT' : 17, 'TACTCA' : 12, 'TACTCC' : 19, 'TACTCG' : 7, 'TACTCT' : 11, 'TATAGC' : 15, 'TATAGT' : 13, 'TATTCA' : 21, 'TATTCC' : 21, 'TATTCG' : 4, 'TATTCT' : 25} ,
        'YT' : { 'TACACA' : 28, 'TACACC' : 36, 'TACACG' : 17, 'TACACT' : 19, 'TATACA' : 32, 'TATACC' : 30, 'TATACG' : 7, 'TATACT' : 31} ,
        'YV' : { 'TACGTA' : 7, 'TACGTC' : 24, 'TACGTG' : 59, 'TACGTT' : 10, 'TATGTA' : 14, 'TATGTC' : 23, 'TATGTG' : 42, 'TATGTT' : 21} ,
        'YW' : { 'TACTGG' : 100, 'TATTGG' : 100} ,
        'YY' : { 'TACTAC' : 62, 'TACTAT' : 38, 'TATTAC' : 50, 'TATTAT' : 50}}

    #Global Bicodon Usage for Homo sapiens. It reflects average bicodon usage in the transcriptome and it is
    #corrected for tRNA abundance (average between HEK and brain cells).
    elif ex_sys == '2':
        CC_dict = {
        'AA' : { 'GCAGCA' : 28, 'GCAGCC' : 19, 'GCAGCG' : 19, 'GCAGCT' : 35, 'GCCGCA' : 11, 'GCCGCC' : 53, 'GCCGCG' : 22, 'GCCGCT' : 14, 'GCGGCA' : 21, 'GCGGCC' : 21, 'GCGGCG' : 30, 'GCGGCT' : 28, 'GCTGCA' : 27, 'GCTGCC' : 19, 'GCTGCG' : 19, 'GCTGCT' : 35 } ,
        'AC' : { 'GCATGC' : 55, 'GCATGT' : 45, 'GCCTGC' : 62, 'GCCTGT' : 38, 'GCGTGC' : 66, 'GCGTGT' : 34, 'GCTTGC' : 54, 'GCTTGT' : 46 } ,
        'AU' : { 'GCATGC' : 55, 'GCATGT' : 45, 'GCCTGC' : 62, 'GCCTGT' : 38, 'GCGTGC' : 66, 'GCGTGT' : 34, 'GCTTGC' : 54, 'GCTTGT' : 46 } ,
        'AD' : { 'GCAGAC' : 73, 'GCAGAT' : 28, 'GCCGAC' : 75, 'GCCGAT' : 25, 'GCGGAC' : 83, 'GCGGAT' : 18, 'GCTGAC' : 76, 'GCTGAT' : 24 } ,
        'AE' : { 'GCAGAA' : 45, 'GCAGAG' : 55, 'GCCGAA' : 18, 'GCCGAG' : 82, 'GCGGAA' : 34, 'GCGGAG' : 66, 'GCTGAA' : 44, 'GCTGAG' : 56 } ,
        'AF' : { 'GCATTC' : 68, 'GCATTT' : 33, 'GCCTTC' : 65, 'GCCTTT' : 35, 'GCGTTC' : 78, 'GCGTTT' : 23, 'GCTTTC' : 73, 'GCTTTT' : 27 } ,
        'AG' : { 'GCAGGA' : 19, 'GCAGGC' : 42, 'GCAGGG' : 31, 'GCAGGT' : 8, 'GCCGGA' : 12, 'GCCGGC' : 45, 'GCCGGG' : 35, 'GCCGGT' : 9, 'GCGGGA' : 12, 'GCGGGC' : 50, 'GCGGGG' : 32, 'GCGGGT' : 5, 'GCTGGA' : 20, 'GCTGGC' : 40, 'GCTGGG' : 31, 'GCTGGT' : 9 } ,
        'AH' : { 'GCACAC' : 77, 'GCACAT' : 23, 'GCCCAC' : 66, 'GCCCAT' : 34, 'GCGCAC' : 90, 'GCGCAT' : 10, 'GCTCAC' : 76, 'GCTCAT' : 25 } ,
        'AI' : { 'GCAATA' : 26, 'GCAATC' : 24, 'GCAATT' : 50, 'GCCATA' : 11, 'GCCATC' : 60, 'GCCATT' : 29, 'GCGATA' : 21, 'GCGATC' : 32, 'GCGATT' : 47, 'GCTATA' : 24, 'GCTATC' : 27, 'GCTATT' : 48 } ,
        'AK' : { 'GCAAAA' : 49, 'GCAAAG' : 51, 'GCCAAA' : 32, 'GCCAAG' : 68, 'GCGAAA' : 41, 'GCGAAG' : 59, 'GCTAAA' : 51, 'GCTAAG' : 49 } ,
        'AL' : { 'GCACTA' : 12, 'GCACTC' : 8, 'GCACTG' : 33, 'GCACTT' : 16, 'GCATTA' : 16, 'GCATTG' : 17, 'GCCCTA' : 5, 'GCCCTC' : 18, 'GCCCTG' : 47, 'GCCCTT' : 9, 'GCCTTA' : 5, 'GCCTTG' : 15, 'GCGCTA' : 9, 'GCGCTC' : 11, 'GCGCTG' : 45, 'GCGCTT' : 11, 'GCGTTA' : 11, 'GCGTTG' : 13, 'GCTCTA' : 12, 'GCTCTC' : 7, 'GCTCTG' : 31, 'GCTCTT' : 16, 'GCTTTA' : 15, 'GCTTTG' : 19 } ,
        'AM' : { 'GCAATG' : 100, 'GCCATG' : 100, 'GCGATG' : 100, 'GCTATG' : 100 } ,
        'AN' : { 'GCAAAC' : 58, 'GCAAAT' : 42, 'GCCAAC' : 64, 'GCCAAT' : 36, 'GCGAAC' : 66, 'GCGAAT' : 34, 'GCTAAC' : 55, 'GCTAAT' : 45 } ,
        'AP' : { 'GCACCA' : 35, 'GCACCC' : 25, 'GCACCG' : 15, 'GCACCT' : 25, 'GCCCCA' : 26, 'GCCCCC' : 30, 'GCCCCG' : 19, 'GCCCCT' : 25, 'GCGCCA' : 28, 'GCGCCC' : 31, 'GCGCCG' : 23, 'GCGCCT' : 17, 'GCTCCA' : 34, 'GCTCCC' : 25, 'GCTCCG' : 16, 'GCTCCT' : 26 } ,
        'AQ' : { 'GCACAA' : 31, 'GCACAG' : 69, 'GCCCAA' : 17, 'GCCCAG' : 83, 'GCGCAA' : 21, 'GCGCAG' : 79, 'GCTCAA' : 31, 'GCTCAG' : 69 } ,
        'AR' : { 'GCAAGA' : 33, 'GCAAGG' : 22, 'GCACGA' : 12, 'GCACGC' : 6, 'GCACGG' : 16, 'GCACGT' : 10, 'GCCAGA' : 15, 'GCCAGG' : 23, 'GCCCGA' : 9, 'GCCCGC' : 21, 'GCCCGG' : 27, 'GCCCGT' : 6, 'GCGAGA' : 22, 'GCGAGG' : 18, 'GCGCGA' : 10, 'GCGCGC' : 17, 'GCGCGG' : 24, 'GCGCGT' : 9, 'GCTAGA' : 24, 'GCTAGG' : 16, 'GCTCGA' : 17, 'GCTCGC' : 9, 'GCTCGG' : 21, 'GCTCGT' : 13 } ,
        'AS' : { 'GCAAGC' : 22, 'GCAAGT' : 13, 'GCATCA' : 19, 'GCATCC' : 11, 'GCATCG' : 10, 'GCATCT' : 24, 'GCCAGC' : 28, 'GCCAGT' : 13, 'GCCTCA' : 11, 'GCCTCC' : 24, 'GCCTCG' : 9, 'GCCTCT' : 15, 'GCGAGC' : 25, 'GCGAGT' : 9, 'GCGTCA' : 16, 'GCGTCC' : 17, 'GCGTCG' : 14, 'GCGTCT' : 19, 'GCTAGC' : 18, 'GCTAGT' : 11, 'GCTTCA' : 21, 'GCTTCC' : 14, 'GCTTCG' : 11, 'GCTTCT' : 26 } ,
        'AT' : { 'GCAACA' : 32, 'GCAACC' : 13, 'GCAACG' : 16, 'GCAACT' : 40, 'GCCACA' : 25, 'GCCACC' : 41, 'GCCACG' : 15, 'GCCACT' : 20, 'GCGACA' : 25, 'GCGACC' : 19, 'GCGACG' : 22, 'GCGACT' : 35, 'GCTACA' : 31, 'GCTACC' : 15, 'GCTACG' : 16, 'GCTACT' : 38 } ,
        'AV' : { 'GCAGTA' : 19, 'GCAGTC' : 9, 'GCAGTG' : 52, 'GCAGTT' : 19, 'GCCGTA' : 5, 'GCCGTC' : 26, 'GCCGTG' : 61, 'GCCGTT' : 8, 'GCGGTA' : 15, 'GCGGTC' : 11, 'GCGGTG' : 61, 'GCGGTT' : 13, 'GCTGTA' : 18, 'GCTGTC' : 11, 'GCTGTG' : 52, 'GCTGTT' : 20 } ,
        'AW' : { 'GCATGG' : 100, 'GCCTGG' : 100, 'GCGTGG' : 100, 'GCTTGG' : 100 } ,
        'AY' : { 'GCATAC' : 54, 'GCATAT' : 46, 'GCCTAC' : 62, 'GCCTAT' : 38, 'GCGTAC' : 64, 'GCGTAT' : 36, 'GCTTAC' : 52, 'GCTTAT' : 48 } ,
        'CA' : { 'TGCGCA' : 21, 'TGCGCC' : 26, 'TGCGCG' : 25, 'TGCGCT' : 27, 'TGTGCA' : 26, 'TGTGCC' : 20, 'TGTGCG' : 16, 'TGTGCT' : 38 } ,
        'CC' : { 'TGCTGC' : 65, 'TGCTGT' : 35, 'TGTTGC' : 74, 'TGTTGT' : 26 } ,
        'CU' : { 'TGCTGC' : 65, 'TGCTGT' : 35, 'TGTTGC' : 74, 'TGTTGT' : 26 } ,
        'CD' : { 'TGCGAC' : 86, 'TGCGAT' : 14, 'TGTGAC' : 77, 'TGTGAT' : 23 } ,
        'CE' : { 'TGCGAA' : 33, 'TGCGAG' : 67, 'TGTGAA' : 47, 'TGTGAG' : 53 } ,
        'CF' : { 'TGCTTC' : 81, 'TGCTTT' : 19, 'TGTTTC' : 74, 'TGTTTT' : 26 } ,
        'CG' : { 'TGCGGA' : 12, 'TGCGGC' : 49, 'TGCGGG' : 34, 'TGCGGT' : 5, 'TGTGGA' : 15, 'TGTGGC' : 45, 'TGTGGG' : 33, 'TGTGGT' : 7 } ,
        'CH' : { 'TGCCAC' : 82, 'TGCCAT' : 18, 'TGTCAC' : 76, 'TGTCAT' : 24 } ,
        'CI' : { 'TGCATA' : 19, 'TGCATC' : 38, 'TGCATT' : 42, 'TGTATA' : 19, 'TGTATC' : 19, 'TGTATT' : 61 } ,
        'CK' : { 'TGCAAA' : 45, 'TGCAAG' : 55, 'TGTAAA' : 51, 'TGTAAG' : 49 } ,
        'CL' : { 'TGCCTA' : 10, 'TGCCTC' : 12, 'TGCCTG' : 36, 'TGCCTT' : 15, 'TGCTTA' : 13, 'TGCTTG' : 15, 'TGTCTA' : 9, 'TGTCTC' : 9, 'TGTCTG' : 33, 'TGTCTT' : 15, 'TGTTTA' : 16, 'TGTTTG' : 18 } ,
        'CM' : { 'TGCATG' : 100, 'TGTATG' : 100 } ,
        'CN' : { 'TGCAAC' : 66, 'TGCAAT' : 34, 'TGTAAC' : 71, 'TGTAAT' : 29 } ,
        'CP' : { 'TGCCCA' : 33, 'TGCCCC' : 27, 'TGCCCG' : 18, 'TGCCCT' : 22, 'TGTCCA' : 53, 'TGTCCC' : 16, 'TGTCCG' : 10, 'TGTCCT' : 21 } ,
        'CQ' : { 'TGCCAA' : 27, 'TGCCAG' : 73, 'TGTCAA' : 29, 'TGTCAG' : 71 } ,
        'CR' : { 'TGCAGA' : 26, 'TGCAGG' : 22, 'TGCCGA' : 12, 'TGCCGC' : 13, 'TGCCGG' : 18, 'TGCCGT' : 10, 'TGTAGA' : 32, 'TGTAGG' : 20, 'TGTCGA' : 13, 'TGTCGC' : 9, 'TGTCGG' : 16, 'TGTCGT' : 10 } ,
        'CS' : { 'TGCAGC' : 30, 'TGCAGT' : 14, 'TGCTCA' : 16, 'TGCTCC' : 12, 'TGCTCG' : 10, 'TGCTCT' : 18, 'TGTAGC' : 29, 'TGTAGT' : 8, 'TGTTCA' : 21, 'TGTTCC' : 11, 'TGTTCG' : 6, 'TGTTCT' : 25 } ,
        'CT' : { 'TGCACA' : 27, 'TGCACC' : 20, 'TGCACG' : 19, 'TGCACT' : 34, 'TGTACA' : 27, 'TGTACC' : 17, 'TGTACG' : 12, 'TGTACT' : 44 } ,
        'CV' : { 'TGCGTA' : 14, 'TGCGTC' : 13, 'TGCGTG' : 59, 'TGCGTT' : 13, 'TGTGTA' : 16, 'TGTGTC' : 12, 'TGTGTG' : 56, 'TGTGTT' : 16 } ,
        'CW' : { 'TGCTGG' : 100, 'TGTTGG' : 100 } ,
        'CY' : { 'TGCTAC' : 65, 'TGCTAT' : 35, 'TGTTAC' : 75, 'TGTTAT' : 25 } ,
        'UA' : { 'TGCGCA' : 21, 'TGCGCC' : 26, 'TGCGCG' : 25, 'TGCGCT' : 27, 'TGTGCA' : 26, 'TGTGCC' : 20, 'TGTGCG' : 16, 'TGTGCT' : 38 } ,
        'UC' : { 'TGCTGC' : 65, 'TGCTGT' : 35, 'TGTTGC' : 74, 'TGTTGT' : 26 } ,
        'UU' : { 'TGCTGC' : 65, 'TGCTGT' : 35, 'TGTTGC' : 74, 'TGTTGT' : 26 } ,
        'UD' : { 'TGCGAC' : 86, 'TGCGAT' : 14, 'TGTGAC' : 77, 'TGTGAT' : 23 } ,
        'UE' : { 'TGCGAA' : 33, 'TGCGAG' : 67, 'TGTGAA' : 47, 'TGTGAG' : 53 } ,
        'UF' : { 'TGCTTC' : 81, 'TGCTTT' : 19, 'TGTTTC' : 74, 'TGTTTT' : 26 } ,
        'UG' : { 'TGCGGA' : 12, 'TGCGGC' : 49, 'TGCGGG' : 34, 'TGCGGT' : 5, 'TGTGGA' : 15, 'TGTGGC' : 45, 'TGTGGG' : 33, 'TGTGGT' : 7 } ,
        'UH' : { 'TGCCAC' : 82, 'TGCCAT' : 18, 'TGTCAC' : 76, 'TGTCAT' : 24 } ,
        'UI' : { 'TGCATA' : 19, 'TGCATC' : 38, 'TGCATT' : 42, 'TGTATA' : 19, 'TGTATC' : 19, 'TGTATT' : 61 } ,
        'UK' : { 'TGCAAA' : 45, 'TGCAAG' : 55, 'TGTAAA' : 51, 'TGTAAG' : 49 } ,
        'UL' : { 'TGCCTA' : 10, 'TGCCTC' : 12, 'TGCCTG' : 36, 'TGCCTT' : 15, 'TGCTTA' : 13, 'TGCTTG' : 15, 'TGTCTA' : 9, 'TGTCTC' : 9, 'TGTCTG' : 33, 'TGTCTT' : 15, 'TGTTTA' : 16, 'TGTTTG' : 18 } ,
        'UM' : { 'TGCATG' : 100, 'TGTATG' : 100 } ,
        'UN' : { 'TGCAAC' : 66, 'TGCAAT' : 34, 'TGTAAC' : 71, 'TGTAAT' : 29 } ,
        'UP' : { 'TGCCCA' : 33, 'TGCCCC' : 27, 'TGCCCG' : 18, 'TGCCCT' : 22, 'TGTCCA' : 53, 'TGTCCC' : 16, 'TGTCCG' : 10, 'TGTCCT' : 21 } ,
        'UQ' : { 'TGCCAA' : 27, 'TGCCAG' : 73, 'TGTCAA' : 29, 'TGTCAG' : 71 } ,
        'UR' : { 'TGCAGA' : 26, 'TGCAGG' : 22, 'TGCCGA' : 12, 'TGCCGC' : 13, 'TGCCGG' : 18, 'TGCCGT' : 10, 'TGTAGA' : 32, 'TGTAGG' : 20, 'TGTCGA' : 13, 'TGTCGC' : 9, 'TGTCGG' : 16, 'TGTCGT' : 10 } ,
        'US' : { 'TGCAGC' : 30, 'TGCAGT' : 14, 'TGCTCA' : 16, 'TGCTCC' : 12, 'TGCTCG' : 10, 'TGCTCT' : 18, 'TGTAGC' : 29, 'TGTAGT' : 8, 'TGTTCA' : 21, 'TGTTCC' : 11, 'TGTTCG' : 6, 'TGTTCT' : 25 } ,
        'UT' : { 'TGCACA' : 27, 'TGCACC' : 20, 'TGCACG' : 19, 'TGCACT' : 34, 'TGTACA' : 27, 'TGTACC' : 17, 'TGTACG' : 12, 'TGTACT' : 44 } ,
        'UV' : { 'TGCGTA' : 14, 'TGCGTC' : 13, 'TGCGTG' : 59, 'TGCGTT' : 13, 'TGTGTA' : 16, 'TGTGTC' : 12, 'TGTGTG' : 56, 'TGTGTT' : 16 } ,
        'UW' : { 'TGCTGG' : 100, 'TGTTGG' : 100 } ,
        'UY' : { 'TGCTAC' : 65, 'TGCTAT' : 35, 'TGTTAC' : 75, 'TGTTAT' : 25 } ,
        'DA' : { 'GACGCA' : 24, 'GACGCC' : 24, 'GACGCG' : 26, 'GACGCT' : 27, 'GATGCA' : 25, 'GATGCC' : 37, 'GATGCG' : 6, 'GATGCT' : 32 } ,
        'DC' : { 'GACTGC' : 55, 'GACTGT' : 45, 'GATTGC' : 42, 'GATTGT' : 58 } ,
        'DU' : { 'GACTGC' : 55, 'GACTGT' : 45, 'GATTGC' : 42, 'GATTGT' : 58 } ,
        'DD' : { 'GACGAC' : 84, 'GACGAT' : 16, 'GATGAC' : 48, 'GATGAT' : 52 } ,
        'DE' : { 'GACGAA' : 35, 'GACGAG' : 65, 'GATGAA' : 52, 'GATGAG' : 48 } ,
        'DF' : { 'GACTTC' : 80, 'GACTTT' : 21, 'GATTTC' : 44, 'GATTTT' : 56 } ,
        'DG' : { 'GACGGA' : 16, 'GACGGC' : 47, 'GACGGG' : 32, 'GACGGT' : 6, 'GATGGA' : 28, 'GATGGC' : 29, 'GATGGG' : 25, 'GATGGT' : 18 } ,
        'DH' : { 'GACCAC' : 82, 'GACCAT' : 18, 'GATCAC' : 47, 'GATCAT' : 53 } ,
        'DI' : { 'GACATA' : 22, 'GACATC' : 41, 'GACATT' : 37, 'GATATA' : 21, 'GATATC' : 31, 'GATATT' : 47 } ,
        'DK' : { 'GACAAA' : 45, 'GACAAG' : 55, 'GATAAA' : 60, 'GATAAG' : 40 } ,
        'DL' : { 'GACCTA' : 12, 'GACCTC' : 11, 'GACCTG' : 32, 'GACCTT' : 15, 'GACTTA' : 14, 'GACTTG' : 17, 'GATCTA' : 8, 'GATCTC' : 16, 'GATCTG' : 25, 'GATCTT' : 18, 'GATTTA' : 13, 'GATTTG' : 20 } ,
        'DM' : { 'GACATG' : 100, 'GATATG' : 100 } ,
        'DN' : { 'GACAAC' : 59, 'GACAAT' : 41, 'GATAAC' : 41, 'GATAAT' : 59 } ,
        'DP' : { 'GACCCA' : 27, 'GACCCC' : 30, 'GACCCG' : 18, 'GACCCT' : 25, 'GATCCA' : 29, 'GATCCC' : 28, 'GATCCG' : 6, 'GATCCT' : 37 } ,
        'DQ' : { 'GACCAA' : 34, 'GACCAG' : 66, 'GATCAA' : 40, 'GATCAG' : 60 } ,
        'DR' : { 'GACAGA' : 26, 'GACAGG' : 21, 'GACCGA' : 13, 'GACCGC' : 10, 'GACCGG' : 18, 'GACCGT' : 12, 'GATAGA' : 26, 'GATAGG' : 13, 'GATCGA' : 17, 'GATCGC' : 15, 'GATCGG' : 16, 'GATCGT' : 13 } ,
        'DS' : { 'GACAGC' : 27, 'GACAGT' : 17, 'GACTCA' : 17, 'GACTCC' : 9, 'GACTCG' : 12, 'GACTCT' : 17, 'GATAGC' : 12, 'GATAGT' : 13, 'GATTCA' : 21, 'GATTCC' : 21, 'GATTCG' : 4, 'GATTCT' : 29 } ,
        'DT' : { 'GACACA' : 29, 'GACACC' : 18, 'GACACG' : 21, 'GACACT' : 31, 'GATACA' : 33, 'GATACC' : 29, 'GATACG' : 7, 'GATACT' : 31 } ,
        'DV' : { 'GACGTA' : 17, 'GACGTC' : 13, 'GACGTG' : 53, 'GACGTT' : 17, 'GATGTA' : 14, 'GATGTC' : 23, 'GATGTG' : 40, 'GATGTT' : 23 } ,
        'DW' : { 'GACTGG' : 100, 'GATTGG' : 100 } ,
        'DY' : { 'GACTAC' : 57, 'GACTAT' : 43, 'GATTAC' : 42, 'GATTAT' : 58 } ,
        'EA' : { 'GAAGCA' : 31, 'GAAGCC' : 17, 'GAAGCG' : 19, 'GAAGCT' : 33, 'GAGGCA' : 26, 'GAGGCC' : 22, 'GAGGCG' : 23, 'GAGGCT' : 29 } ,
        'EC' : { 'GAATGC' : 42, 'GAATGT' : 58, 'GAGTGC' : 52, 'GAGTGT' : 48 } ,
        'EU' : { 'GAATGC' : 42, 'GAATGT' : 58, 'GAGTGC' : 52, 'GAGTGT' : 48 } ,
        'ED' : { 'GAAGAC' : 71, 'GAAGAT' : 30, 'GAGGAC' : 80, 'GAGGAT' : 21 } ,
        'EE' : { 'GAAGAA' : 51, 'GAAGAG' : 49, 'GAGGAA' : 40, 'GAGGAG' : 60 } ,
        'EF' : { 'GAATTC' : 69, 'GAATTT' : 31, 'GAGTTC' : 77, 'GAGTTT' : 24 } ,
        'EG' : { 'GAAGGA' : 26, 'GAAGGC' : 36, 'GAAGGG' : 28, 'GAAGGT' : 9, 'GAGGGA' : 20, 'GAGGGC' : 42, 'GAGGGG' : 30, 'GAGGGT' : 8 } ,
        'EH' : { 'GAACAC' : 72, 'GAACAT' : 28, 'GAGCAC' : 82, 'GAGCAT' : 18 } ,
        'EI' : { 'GAAATA' : 29, 'GAAATC' : 30, 'GAAATT' : 40, 'GAGATA' : 23, 'GAGATC' : 41, 'GAGATT' : 36 } ,
        'EK' : { 'GAAAAA' : 50, 'GAAAAG' : 50, 'GAGAAA' : 45, 'GAGAAG' : 55 } ,
        'EL' : { 'GAACTA' : 15, 'GAACTC' : 8, 'GAACTG' : 25, 'GAACTT' : 19, 'GAATTA' : 17, 'GAATTG' : 17, 'GAGCTA' : 13, 'GAGCTC' : 9, 'GAGCTG' : 37, 'GAGCTT' : 14, 'GAGTTA' : 12, 'GAGTTG' : 14 } ,
        'EM' : { 'GAAATG' : 100, 'GAGATG' : 100 } ,
        'EN' : { 'GAAAAC' : 49, 'GAAAAT' : 51, 'GAGAAC' : 57, 'GAGAAT' : 43 } ,
        'EP' : { 'GAACCA' : 32, 'GAACCC' : 25, 'GAACCG' : 16, 'GAACCT' : 28, 'GAGCCA' : 26, 'GAGCCC' : 30, 'GAGCCG' : 19, 'GAGCCT' : 24 } ,
        'EQ' : { 'GAACAA' : 42, 'GAACAG' : 58, 'GAGCAA' : 33, 'GAGCAG' : 67 } ,
        'ER' : { 'GAAAGA' : 32, 'GAAAGG' : 23, 'GAACGA' : 14, 'GAACGC' : 5, 'GAACGG' : 15, 'GAACGT' : 12, 'GAGAGA' : 23, 'GAGAGG' : 22, 'GAGCGA' : 13, 'GAGCGC' : 10, 'GAGCGG' : 21, 'GAGCGT' : 12 } ,
        'ES' : { 'GAAAGC' : 22, 'GAAAGT' : 21, 'GAATCA' : 18, 'GAATCC' : 8, 'GAATCG' : 11, 'GAATCT' : 21, 'GAGAGC' : 27, 'GAGAGT' : 17, 'GAGTCA' : 16, 'GAGTCC' : 10, 'GAGTCG' : 12, 'GAGTCT' : 18 } ,
        'ET' : { 'GAAACA' : 34, 'GAAACC' : 13, 'GAAACG' : 18, 'GAAACT' : 35, 'GAGACA' : 29, 'GAGACC' : 19, 'GAGACG' : 23, 'GAGACT' : 30 } ,
        'EV' : { 'GAAGTA' : 24, 'GAAGTC' : 10, 'GAAGTG' : 40, 'GAAGTT' : 25, 'GAGGTA' : 20, 'GAGGTC' : 11, 'GAGGTG' : 49, 'GAGGTT' : 20 } ,
        'EW' : { 'GAATGG' : 100, 'GAGTGG' : 100 } ,
        'EY' : { 'GAATAC' : 47, 'GAATAT' : 53, 'GAGTAC' : 56, 'GAGTAT' : 44 } ,
        'FA' : { 'TTCGCA' : 22, 'TTCGCC' : 27, 'TTCGCG' : 24, 'TTCGCT' : 28, 'TTTGCA' : 27, 'TTTGCC' : 36, 'TTTGCG' : 5, 'TTTGCT' : 32 } ,
        'FC' : { 'TTCTGC' : 61, 'TTCTGT' : 39, 'TTTTGC' : 42, 'TTTTGT' : 58 } ,
        'FU' : { 'TTCTGC' : 61, 'TTCTGT' : 39, 'TTTTGC' : 42, 'TTTTGT' : 58 } ,
        'FD' : { 'TTCGAC' : 84, 'TTCGAT' : 16, 'TTTGAC' : 47, 'TTTGAT' : 53 } ,
        'FE' : { 'TTCGAA' : 34, 'TTCGAG' : 66, 'TTTGAA' : 52, 'TTTGAG' : 48 } ,
        'FF' : { 'TTCTTC' : 82, 'TTCTTT' : 19, 'TTTTTC' : 50, 'TTTTTT' : 50 } ,
        'FG' : { 'TTCGGA' : 13, 'TTCGGC' : 48, 'TTCGGG' : 34, 'TTCGGT' : 6, 'TTTGGA' : 31, 'TTTGGC' : 25, 'TTTGGG' : 25, 'TTTGGT' : 19 } ,
        'FH' : { 'TTCCAC' : 82, 'TTCCAT' : 19, 'TTTCAC' : 43, 'TTTCAT' : 57 } ,
        'FI' : { 'TTCATA' : 19, 'TTCATC' : 39, 'TTCATT' : 42, 'TTTATA' : 21, 'TTTATC' : 31, 'TTTATT' : 48 } ,
        'FK' : { 'TTCAAA' : 44, 'TTCAAG' : 56, 'TTTAAA' : 54, 'TTTAAG' : 46 } ,
        'FL' : { 'TTCCTA' : 11, 'TTCCTC' : 11, 'TTCCTG' : 36, 'TTCCTT' : 14, 'TTCTTA' : 13, 'TTCTTG' : 15, 'TTTCTA' : 10, 'TTTCTC' : 17, 'TTTCTG' : 27, 'TTTCTT' : 19, 'TTTTTA' : 12, 'TTTTTG' : 15 } ,
        'FM' : { 'TTCATG' : 100, 'TTTATG' : 100 } ,
        'FN' : { 'TTCAAC' : 66, 'TTCAAT' : 34, 'TTTAAC' : 44, 'TTTAAT' : 56 } ,
        'FP' : { 'TTCCCA' : 33, 'TTCCCC' : 26, 'TTCCCG' : 17, 'TTCCCT' : 24, 'TTTCCA' : 33, 'TTTCCC' : 24, 'TTTCCG' : 5, 'TTTCCT' : 38 } ,
        'FQ' : { 'TTCCAA' : 28, 'TTCCAG' : 72, 'TTTCAA' : 37, 'TTTCAG' : 63 } ,
        'FR' : { 'TTCAGA' : 26, 'TTCAGG' : 20, 'TTCCGA' : 13, 'TTCCGC' : 12, 'TTCCGG' : 19, 'TTCCGT' : 11, 'TTTAGA' : 24, 'TTTAGG' : 17, 'TTTCGA' : 20, 'TTTCGC' : 12, 'TTTCGG' : 16, 'TTTCGT' : 11 } ,
        'FS' : { 'TTCAGC' : 27, 'TTCAGT' : 14, 'TTCTCA' : 16, 'TTCTCC' : 13, 'TTCTCG' : 10, 'TTCTCT' : 19, 'TTTAGC' : 12, 'TTTAGT' : 14, 'TTTTCA' : 22, 'TTTTCC' : 20, 'TTTTCG' : 3, 'TTTTCT' : 29 } ,
        'FT' : { 'TTCACA' : 26, 'TTCACC' : 21, 'TTCACG' : 18, 'TTCACT' : 35, 'TTTACA' : 32, 'TTTACC' : 30, 'TTTACG' : 5, 'TTTACT' : 33 } ,
        'FV' : { 'TTCGTA' : 15, 'TTCGTC' : 13, 'TTCGTG' : 58, 'TTCGTT' : 14, 'TTTGTA' : 14, 'TTTGTC' : 22, 'TTTGTG' : 42, 'TTTGTT' : 22 } ,
        'FW' : { 'TTCTGG' : 100, 'TTTTGG' : 100 } ,
        'FY' : { 'TTCTAC' : 63, 'TTCTAT' : 37, 'TTTTAC' : 40, 'TTTTAT' : 60 } ,
        'GA' : { 'GGAGCA' : 29, 'GGAGCC' : 19, 'GGAGCG' : 18, 'GGAGCT' : 35, 'GGCGCA' : 22, 'GGCGCC' : 24, 'GGCGCG' : 29, 'GGCGCT' : 25, 'GGGGCA' : 25, 'GGGGCC' : 22, 'GGGGCG' : 22, 'GGGGCT' : 30, 'GGTGCA' : 23, 'GGTGCC' : 39, 'GGTGCG' : 7, 'GGTGCT' : 31 } ,
        'GC' : { 'GGATGC' : 56, 'GGATGT' : 44, 'GGCTGC' : 58, 'GGCTGT' : 42, 'GGGTGC' : 55, 'GGGTGT' : 45, 'GGTTGC' : 42, 'GGTTGT' : 58 } ,
        'GU' : { 'GGATGC' : 56, 'GGATGT' : 44, 'GGCTGC' : 58, 'GGCTGT' : 42, 'GGGTGC' : 55, 'GGGTGT' : 45, 'GGTTGC' : 42, 'GGTTGT' : 58 } ,
        'GD' : { 'GGAGAC' : 76, 'GGAGAT' : 24, 'GGCGAC' : 85, 'GGCGAT' : 16, 'GGGGAC' : 83, 'GGGGAT' : 18, 'GGTGAC' : 53, 'GGTGAT' : 47 } ,
        'GE' : { 'GGAGAA' : 47, 'GGAGAG' : 53, 'GGCGAA' : 32, 'GGCGAG' : 68, 'GGGGAA' : 40, 'GGGGAG' : 60, 'GGTGAA' : 50, 'GGTGAG' : 50 } ,
        'GF' : { 'GGATTC' : 69, 'GGATTT' : 32, 'GGCTTC' : 83, 'GGCTTT' : 18, 'GGGTTC' : 74, 'GGGTTT' : 27, 'GGTTTC' : 48, 'GGTTTT' : 52 } ,
        'GG' : { 'GGAGGA' : 21, 'GGAGGC' : 43, 'GGAGGG' : 27, 'GGAGGT' : 9, 'GGCGGA' : 16, 'GGCGGC' : 49, 'GGCGGG' : 29, 'GGCGGT' : 5, 'GGGGGA' : 19, 'GGGGGC' : 49, 'GGGGGG' : 23, 'GGGGGT' : 9, 'GGTGGA' : 24, 'GGTGGC' : 34, 'GGTGGG' : 21, 'GGTGGT' : 21 } ,
        'GH' : { 'GGACAC' : 76, 'GGACAT' : 24, 'GGCCAC' : 85, 'GGCCAT' : 16, 'GGGCAC' : 82, 'GGGCAT' : 18, 'GGTCAC' : 53, 'GGTCAT' : 47 } ,
        'GI' : { 'GGAATA' : 25, 'GGAATC' : 27, 'GGAATT' : 48, 'GGCATA' : 21, 'GGCATC' : 45, 'GGCATT' : 34, 'GGGATA' : 24, 'GGGATC' : 39, 'GGGATT' : 37, 'GGTATA' : 19, 'GGTATC' : 38, 'GGTATT' : 43 } ,
        'GK' : { 'GGAAAA' : 52, 'GGAAAG' : 48, 'GGCAAA' : 45, 'GGCAAG' : 55, 'GGGAAA' : 46, 'GGGAAG' : 54, 'GGTAAA' : 70, 'GGTAAG' : 30 } ,
        'GL' : { 'GGACTA' : 12, 'GGACTC' : 10, 'GGACTG' : 28, 'GGACTT' : 18, 'GGATTA' : 16, 'GGATTG' : 16, 'GGCCTA' : 12, 'GGCCTC' : 12, 'GGCCTG' : 32, 'GGCCTT' : 15, 'GGCTTA' : 13, 'GGCTTG' : 16, 'GGGCTA' : 12, 'GGGCTC' : 11, 'GGGCTG' : 36, 'GGGCTT' : 15, 'GGGTTA' : 12, 'GGGTTG' : 13, 'GGTCTA' : 7, 'GGTCTC' : 19, 'GGTCTG' : 27, 'GGTCTT' : 17, 'GGTTTA' : 12, 'GGTTTG' : 18 } ,
        'GM' : { 'GGAATG' : 100, 'GGCATG' : 100, 'GGGATG' : 100, 'GGTATG' : 100 } ,
        'GN' : { 'GGAAAC' : 60, 'GGAAAT' : 40, 'GGCAAC' : 59, 'GGCAAT' : 41, 'GGGAAC' : 59, 'GGGAAT' : 41, 'GGTAAC' : 50, 'GGTAAT' : 50 } ,
        'GP' : { 'GGACCA' : 37, 'GGACCC' : 25, 'GGACCG' : 14, 'GGACCT' : 25, 'GGCCCA' : 25, 'GGCCCC' : 31, 'GGCCCG' : 19, 'GGCCCT' : 25, 'GGGCCA' : 25, 'GGGCCC' : 31, 'GGGCCG' : 19, 'GGGCCT' : 24, 'GGTCCA' : 25, 'GGTCCC' : 35, 'GGTCCG' : 7, 'GGTCCT' : 34 } ,
        'GQ' : { 'GGACAA' : 33, 'GGACAG' : 67, 'GGCCAA' : 33, 'GGCCAG' : 67, 'GGGCAA' : 31, 'GGGCAG' : 69, 'GGTCAA' : 33, 'GGTCAG' : 67 } ,
        'GR' : { 'GGAAGA' : 37, 'GGAAGG' : 22, 'GGACGA' : 12, 'GGACGC' : 5, 'GGACGG' : 13, 'GGACGT' : 10, 'GGCAGA' : 20, 'GGCAGG' : 19, 'GGCCGA' : 14, 'GGCCGC' : 14, 'GGCCGG' : 20, 'GGCCGT' : 13, 'GGGAGA' : 24, 'GGGAGG' : 22, 'GGGCGA' : 13, 'GGGCGC' : 11, 'GGGCGG' : 19, 'GGGCGT' : 12, 'GGTAGA' : 18, 'GGTAGG' : 10, 'GGTCGA' : 17, 'GGTCGC' : 20, 'GGTCGG' : 19, 'GGTCGT' : 15 } ,
        'GS' : { 'GGAAGC' : 28, 'GGAAGT' : 17, 'GGATCA' : 18, 'GGATCC' : 8, 'GGATCG' : 9, 'GGATCT' : 20, 'GGCAGC' : 26, 'GGCAGT' : 16, 'GGCTCA' : 15, 'GGCTCC' : 13, 'GGCTCG' : 13, 'GGCTCT' : 17, 'GGGAGC' : 27, 'GGGAGT' : 16, 'GGGTCA' : 16, 'GGGTCC' : 10, 'GGGTCG' : 12, 'GGGTCT' : 19, 'GGTAGC' : 13, 'GGTAGT' : 11, 'GGTTCA' : 18, 'GGTTCC' : 25, 'GGTTCG' : 4, 'GGTTCT' : 29 } ,
        'GT' : { 'GGAACA' : 32, 'GGAACC' : 14, 'GGAACG' : 15, 'GGAACT' : 39, 'GGCACA' : 28, 'GGCACC' : 22, 'GGCACG' : 22, 'GGCACT' : 29, 'GGGACA' : 29, 'GGGACC' : 20, 'GGGACG' : 21, 'GGGACT' : 31, 'GGTACA' : 29, 'GGTACC' : 36, 'GGTACG' : 5, 'GGTACT' : 30 } ,
        'GV' : { 'GGAGTA' : 20, 'GGAGTC' : 13, 'GGAGTG' : 47, 'GGAGTT' : 20, 'GGCGTA' : 17, 'GGCGTC' : 16, 'GGCGTG' : 48, 'GGCGTT' : 19, 'GGGGTA' : 19, 'GGGGTC' : 16, 'GGGGTG' : 44, 'GGGGTT' : 21, 'GGTGTA' : 12, 'GGTGTC' : 26, 'GGTGTG' : 40, 'GGTGTT' : 23 } ,
        'GW' : { 'GGATGG' : 100, 'GGCTGG' : 100, 'GGGTGG' : 100, 'GGTTGG' : 100 } ,
        'GY' : { 'GGATAC' : 55, 'GGATAT' : 45, 'GGCTAC' : 58, 'GGCTAT' : 42, 'GGGTAC' : 57, 'GGGTAT' : 43, 'GGTTAC' : 44, 'GGTTAT' : 56 } ,
        'HA' : { 'CACGCA' : 24, 'CACGCC' : 24, 'CACGCG' : 25, 'CACGCT' : 27, 'CATGCA' : 26, 'CATGCC' : 36, 'CATGCG' : 7, 'CATGCT' : 32 } ,
        'HC' : { 'CACTGC' : 59, 'CACTGT' : 41, 'CATTGC' : 47, 'CATTGT' : 53 } ,
        'HU' : { 'CACTGC' : 59, 'CACTGT' : 41, 'CATTGC' : 47, 'CATTGT' : 53 } ,
        'HD' : { 'CACGAC' : 83, 'CACGAT' : 17, 'CATGAC' : 49, 'CATGAT' : 51 } ,
        'HE' : { 'CACGAA' : 36, 'CACGAG' : 64, 'CATGAA' : 52, 'CATGAG' : 48 } ,
        'HF' : { 'CACTTC' : 81, 'CACTTT' : 20, 'CATTTC' : 46, 'CATTTT' : 54 } ,
        'HG' : { 'CACGGA' : 17, 'CACGGC' : 45, 'CACGGG' : 32, 'CACGGT' : 6, 'CATGGA' : 29, 'CATGGC' : 29, 'CATGGG' : 23, 'CATGGT' : 19 } ,
        'HH' : { 'CACCAC' : 83, 'CACCAT' : 17, 'CATCAC' : 49, 'CATCAT' : 51 } ,
        'HI' : { 'CACATA' : 22, 'CACATC' : 41, 'CACATT' : 37, 'CATATA' : 22, 'CATATC' : 32, 'CATATT' : 46 } ,
        'HK' : { 'CACAAA' : 44, 'CACAAG' : 56, 'CATAAA' : 51, 'CATAAG' : 49 } ,
        'HL' : { 'CACCTA' : 12, 'CACCTC' : 10, 'CACCTG' : 34, 'CACCTT' : 15, 'CACTTA' : 13, 'CACTTG' : 16, 'CATCTA' : 8, 'CATCTC' : 17, 'CATCTG' : 26, 'CATCTT' : 18, 'CATTTA' : 13, 'CATTTG' : 18 } ,
        'HM' : { 'CACATG' : 100, 'CATATG' : 100 } ,
        'HN' : { 'CACAAC' : 60, 'CACAAT' : 40, 'CATAAC' : 41, 'CATAAT' : 59 } ,
        'HP' : { 'CACCCA' : 28, 'CACCCC' : 27, 'CACCCG' : 20, 'CACCCT' : 24, 'CATCCA' : 31, 'CATCCC' : 26, 'CATCCG' : 7, 'CATCCT' : 36 } ,
        'HQ' : { 'CACCAA' : 30, 'CACCAG' : 70, 'CATCAA' : 30, 'CATCAG' : 70 } ,
        'HR' : { 'CACAGA' : 25, 'CACAGG' : 22, 'CACCGA' : 13, 'CACCGC' : 10, 'CACCGG' : 18, 'CACCGT' : 11, 'CATAGA' : 19, 'CATAGG' : 14, 'CATCGA' : 16, 'CATCGC' : 16, 'CATCGG' : 20, 'CATCGT' : 15 } ,
        'HS' : { 'CACAGC' : 29, 'CACAGT' : 17, 'CACTCA' : 16, 'CACTCC' : 8, 'CACTCG' : 13, 'CACTCT' : 17, 'CATAGC' : 12, 'CATAGT' : 13, 'CATTCA' : 22, 'CATTCC' : 21, 'CATTCG' : 4, 'CATTCT' : 27 } ,
        'HT' : { 'CACACA' : 29, 'CACACC' : 15, 'CACACG' : 23, 'CACACT' : 33, 'CATACA' : 26, 'CATACC' : 21, 'CATACG' : 6, 'CATACT' : 46 } ,
        'HV' : { 'CACGTA' : 17, 'CACGTC' : 13, 'CACGTG' : 53, 'CACGTT' : 17, 'CATGTA' : 14, 'CATGTC' : 23, 'CATGTG' : 40, 'CATGTT' : 23 } ,
        'HW' : { 'CACTGG' : 100, 'CATTGG' : 100 } ,
        'HY' : { 'CACTAC' : 60, 'CACTAT' : 40, 'CATTAC' : 47, 'CATTAT' : 53 } ,
        'IA' : { 'ATAGCA' : 32, 'ATAGCC' : 14, 'ATAGCG' : 16, 'ATAGCT' : 38, 'ATCGCA' : 22, 'ATCGCC' : 26, 'ATCGCG' : 19, 'ATCGCT' : 33, 'ATTGCA' : 29, 'ATTGCC' : 17, 'ATTGCG' : 17, 'ATTGCT' : 37 } ,
        'IC' : { 'ATATGC' : 58, 'ATATGT' : 42, 'ATCTGC' : 77, 'ATCTGT' : 23, 'ATTTGC' : 52, 'ATTTGT' : 48 } ,
        'IU' : { 'ATATGC' : 58, 'ATATGT' : 42, 'ATCTGC' : 77, 'ATCTGT' : 23, 'ATTTGC' : 52, 'ATTTGT' : 48 } ,
        'ID' : { 'ATAGAC' : 72, 'ATAGAT' : 28, 'ATCGAC' : 83, 'ATCGAT' : 17, 'ATTGAC' : 74, 'ATTGAT' : 27 } ,
        'IE' : { 'ATAGAA' : 53, 'ATAGAG' : 47, 'ATCGAA' : 35, 'ATCGAG' : 65, 'ATTGAA' : 50, 'ATTGAG' : 50 } ,
        'IF' : { 'ATATTC' : 67, 'ATATTT' : 34, 'ATCTTC' : 82, 'ATCTTT' : 19, 'ATTTTC' : 73, 'ATTTTT' : 28 } ,
        'IG' : { 'ATAGGA' : 23, 'ATAGGC' : 40, 'ATAGGG' : 27, 'ATAGGT' : 10, 'ATCGGA' : 11, 'ATCGGC' : 50, 'ATCGGG' : 33, 'ATCGGT' : 6, 'ATTGGA' : 22, 'ATTGGC' : 39, 'ATTGGG' : 28, 'ATTGGT' : 10 } ,
        'IH' : { 'ATACAC' : 70, 'ATACAT' : 30, 'ATCCAC' : 82, 'ATCCAT' : 18, 'ATTCAC' : 70, 'ATTCAT' : 30 } ,
        'II' : { 'ATAATA' : 24, 'ATAATC' : 19, 'ATAATT' : 57, 'ATCATA' : 15, 'ATCATC' : 31, 'ATCATT' : 54, 'ATTATA' : 25, 'ATTATC' : 27, 'ATTATT' : 48 } ,
        'IK' : { 'ATAAAA' : 52, 'ATAAAG' : 48, 'ATCAAA' : 44, 'ATCAAG' : 56, 'ATTAAA' : 54, 'ATTAAG' : 46 } ,
        'IL' : { 'ATACTA' : 12, 'ATACTC' : 6, 'ATACTG' : 28, 'ATACTT' : 17, 'ATATTA' : 19, 'ATATTG' : 18, 'ATCCTA' : 9, 'ATCCTC' : 11, 'ATCCTG' : 39, 'ATCCTT' : 12, 'ATCTTA' : 13, 'ATCTTG' : 15, 'ATTCTA' : 12, 'ATTCTC' : 8, 'ATTCTG' : 25, 'ATTCTT' : 18, 'ATTTTA' : 18, 'ATTTTG' : 19 } ,
        'IM' : { 'ATAATG' : 100, 'ATCATG' : 100, 'ATTATG' : 100 } ,
        'IN' : { 'ATAAAC' : 62, 'ATAAAT' : 38, 'ATCAAC' : 78, 'ATCAAT' : 22, 'ATTAAC' : 53, 'ATTAAT' : 47 } ,
        'IP' : { 'ATACCA' : 45, 'ATACCC' : 17, 'ATACCG' : 12, 'ATACCT' : 26, 'ATCCCA' : 48, 'ATCCCC' : 19, 'ATCCCG' : 13, 'ATCCCT' : 20, 'ATTCCA' : 35, 'ATTCCC' : 21, 'ATTCCG' : 14, 'ATTCCT' : 30 } ,
        'IQ' : { 'ATACAA' : 33, 'ATACAG' : 67, 'ATCCAA' : 23, 'ATCCAG' : 77, 'ATTCAA' : 35, 'ATTCAG' : 65 } ,
        'IR' : { 'ATAAGA' : 38, 'ATAAGG' : 24, 'ATACGA' : 12, 'ATACGC' : 4, 'ATACGG' : 13, 'ATACGT' : 10, 'ATCAGA' : 29, 'ATCAGG' : 19, 'ATCCGA' : 12, 'ATCCGC' : 13, 'ATCCGG' : 18, 'ATCCGT' : 10, 'ATTAGA' : 26, 'ATTAGG' : 17, 'ATTCGA' : 19, 'ATTCGC' : 8, 'ATTCGG' : 16, 'ATTCGT' : 14 } ,
        'IS' : { 'ATAAGC' : 25, 'ATAAGT' : 12, 'ATATCA' : 22, 'ATATCC' : 9, 'ATATCG' : 8, 'ATATCT' : 24, 'ATCAGC' : 34, 'ATCAGT' : 8, 'ATCTCA' : 17, 'ATCTCC' : 13, 'ATCTCG' : 7, 'ATCTCT' : 21, 'ATTAGC' : 18, 'ATTAGT' : 12, 'ATTTCA' : 22, 'ATTTCC' : 11, 'ATTTCG' : 10, 'ATTTCT' : 26 } ,
        'IT' : { 'ATAACA' : 31, 'ATAACC' : 12, 'ATAACG' : 13, 'ATAACT' : 44, 'ATCACA' : 25, 'ATCACC' : 20, 'ATCACG' : 14, 'ATCACT' : 40, 'ATTACA' : 33, 'ATTACC' : 13, 'ATTACG' : 16, 'ATTACT' : 39 } ,
        'IV' : { 'ATAGTA' : 21, 'ATAGTC' : 10, 'ATAGTG' : 49, 'ATAGTT' : 21, 'ATCGTA' : 13, 'ATCGTC' : 15, 'ATCGTG' : 60, 'ATCGTT' : 11, 'ATTGTA' : 20, 'ATTGTC' : 12, 'ATTGTG' : 47, 'ATTGTT' : 22 } ,
        'IW' : { 'ATATGG' : 100, 'ATCTGG' : 100, 'ATTTGG' : 100 } ,
        'IY' : { 'ATATAC' : 58, 'ATATAT' : 42, 'ATCTAC' : 78, 'ATCTAT' : 22, 'ATTTAC' : 51, 'ATTTAT' : 49 } ,
        'KA' : { 'AAAGCA' : 31, 'AAAGCC' : 18, 'AAAGCG' : 18, 'AAAGCT' : 34, 'AAGGCA' : 27, 'AAGGCC' : 20, 'AAGGCG' : 21, 'AAGGCT' : 31 } ,
        'KC' : { 'AAATGC' : 46, 'AAATGT' : 54, 'AAGTGC' : 55, 'AAGTGT' : 45 } ,
        'KU' : { 'AAATGC' : 46, 'AAATGT' : 54, 'AAGTGC' : 55, 'AAGTGT' : 45 } ,
        'KD' : { 'AAAGAC' : 72, 'AAAGAT' : 28, 'AAGGAC' : 78, 'AAGGAT' : 23 } ,
        'KE' : { 'AAAGAA' : 53, 'AAAGAG' : 47, 'AAGGAA' : 44, 'AAGGAG' : 56 } ,
        'KF' : { 'AAATTC' : 71, 'AAATTT' : 29, 'AAGTTC' : 76, 'AAGTTT' : 24 } ,
        'KG' : { 'AAAGGA' : 26, 'AAAGGC' : 38, 'AAAGGG' : 27, 'AAAGGT' : 9, 'AAGGGA' : 19, 'AAGGGC' : 43, 'AAGGGG' : 29, 'AAGGGT' : 9 } ,
        'KH' : { 'AAACAC' : 73, 'AAACAT' : 27, 'AAGCAC' : 80, 'AAGCAT' : 21 } ,
        'KI' : { 'AAAATA' : 29, 'AAAATC' : 30, 'AAAATT' : 41, 'AAGATA' : 25, 'AAGATC' : 37, 'AAGATT' : 38 } ,
        'KK' : { 'AAAAAA' : 47, 'AAAAAG' : 53, 'AAGAAA' : 47, 'AAGAAG' : 53 } ,
        'KL' : { 'AAACTA' : 14, 'AAACTC' : 8, 'AAACTG' : 26, 'AAACTT' : 19, 'AAATTA' : 17, 'AAATTG' : 16, 'AAGCTA' : 13, 'AAGCTC' : 9, 'AAGCTG' : 34, 'AAGCTT' : 15, 'AAGTTA' : 13, 'AAGTTG' : 15 } ,
        'KM' : { 'AAAATG' : 100, 'AAGATG' : 100 } ,
        'KN' : { 'AAAAAC' : 51, 'AAAAAT' : 49, 'AAGAAC' : 58, 'AAGAAT' : 42 } ,
        'KP' : { 'AAACCA' : 33, 'AAACCC' : 26, 'AAACCG' : 14, 'AAACCT' : 27, 'AAGCCA' : 30, 'AAGCCC' : 28, 'AAGCCG' : 17, 'AAGCCT' : 24 } ,
        'KQ' : { 'AAACAA' : 41, 'AAACAG' : 59, 'AAGCAA' : 35, 'AAGCAG' : 65 } ,
        'KR' : { 'AAAAGA' : 32, 'AAAAGG' : 23, 'AAACGA' : 14, 'AAACGC' : 5, 'AAACGG' : 14, 'AAACGT' : 12, 'AAGAGA' : 27, 'AAGAGG' : 24, 'AAGCGA' : 13, 'AAGCGC' : 8, 'AAGCGG' : 17, 'AAGCGT' : 11 } ,
        'KS' : { 'AAAAGC' : 21, 'AAAAGT' : 17, 'AAATCA' : 20, 'AAATCC' : 8, 'AAATCG' : 11, 'AAATCT' : 22, 'AAGAGC' : 25, 'AAGAGT' : 16, 'AAGTCA' : 19, 'AAGTCC' : 10, 'AAGTCG' : 11, 'AAGTCT' : 19 } ,
        'KT' : { 'AAAACA' : 33, 'AAAACC' : 14, 'AAAACG' : 17, 'AAAACT' : 35, 'AAGACA' : 30, 'AAGACC' : 18, 'AAGACG' : 20, 'AAGACT' : 32 } ,
        'KV' : { 'AAAGTA' : 22, 'AAAGTC' : 11, 'AAAGTG' : 43, 'AAAGTT' : 23, 'AAGGTA' : 18, 'AAGGTC' : 12, 'AAGGTG' : 50, 'AAGGTT' : 20 } ,
        'KW' : { 'AAATGG' : 100, 'AAGTGG' : 100 } ,
        'KY' : { 'AAATAC' : 52, 'AAATAT' : 48, 'AAGTAC' : 58, 'AAGTAT' : 42 } ,
        'LA' : { 'CTAGCA' : 31, 'CTAGCC' : 15, 'CTAGCG' : 17, 'CTAGCT' : 37, 'CTCGCA' : 15, 'CTCGCC' : 51, 'CTCGCG' : 13, 'CTCGCT' : 21, 'CTGGCA' : 25, 'CTGGCC' : 23, 'CTGGCG' : 22, 'CTGGCT' : 31, 'CTTGCA' : 29, 'CTTGCC' : 16, 'CTTGCG' : 16, 'CTTGCT' : 39, 'TTAGCA' : 33, 'TTAGCC' : 13, 'TTAGCG' : 17, 'TTAGCT' : 37, 'TTGGCA' : 29, 'TTGGCC' : 17, 'TTGGCG' : 18, 'TTGGCT' : 37 } ,
        'LC' : { 'CTATGC' : 59, 'CTATGT' : 41, 'CTCTGC' : 61, 'CTCTGT' : 39, 'CTGTGC' : 60, 'CTGTGT' : 40, 'CTTTGC' : 58, 'CTTTGT' : 42, 'TTATGC' : 52, 'TTATGT' : 48, 'TTGTGC' : 56, 'TTGTGT' : 44 } ,
        'LU' : { 'CTATGC' : 59, 'CTATGT' : 41, 'CTCTGC' : 61, 'CTCTGT' : 39, 'CTGTGC' : 60, 'CTGTGT' : 40, 'CTTTGC' : 58, 'CTTTGT' : 42, 'TTATGC' : 52, 'TTATGT' : 48, 'TTGTGC' : 56, 'TTGTGT' : 44 } ,
        'LD' : { 'CTAGAC' : 74, 'CTAGAT' : 27, 'CTCGAC' : 62, 'CTCGAT' : 38, 'CTGGAC' : 80, 'CTGGAT' : 20, 'CTTGAC' : 72, 'CTTGAT' : 28, 'TTAGAC' : 68, 'TTAGAT' : 33, 'TTGGAC' : 71, 'TTGGAT' : 29 } ,
        'LE' : { 'CTAGAA' : 48, 'CTAGAG' : 52, 'CTCGAA' : 28, 'CTCGAG' : 72, 'CTGGAA' : 38, 'CTGGAG' : 62, 'CTTGAA' : 52, 'CTTGAG' : 48, 'TTAGAA' : 54, 'TTAGAG' : 46, 'TTGGAA' : 47, 'TTGGAG' : 53 } ,
        'LF' : { 'CTATTC' : 68, 'CTATTT' : 32, 'CTCTTC' : 65, 'CTCTTT' : 35, 'CTGTTC' : 76, 'CTGTTT' : 24, 'CTTTTC' : 75, 'CTTTTT' : 25, 'TTATTC' : 66, 'TTATTT' : 35, 'TTGTTC' : 69, 'TTGTTT' : 31 } ,
        'LG' : { 'CTAGGA' : 21, 'CTAGGC' : 41, 'CTAGGG' : 29, 'CTAGGT' : 9, 'CTCGGA' : 15, 'CTCGGC' : 41, 'CTCGGG' : 33, 'CTCGGT' : 12, 'CTGGGA' : 16, 'CTGGGC' : 47, 'CTGGGG' : 30, 'CTGGGT' : 7, 'CTTGGA' : 21, 'CTTGGC' : 41, 'CTTGGG' : 29, 'CTTGGT' : 10, 'TTAGGA' : 25, 'TTAGGC' : 37, 'TTAGGG' : 27, 'TTAGGT' : 11, 'TTGGGA' : 20, 'TTGGGC' : 41, 'TTGGGG' : 29, 'TTGGGT' : 9 } ,
        'LH' : { 'CTACAC' : 75, 'CTACAT' : 26, 'CTCCAC' : 61, 'CTCCAT' : 39, 'CTGCAC' : 85, 'CTGCAT' : 15, 'CTTCAC' : 71, 'CTTCAT' : 29, 'TTACAC' : 70, 'TTACAT' : 30, 'TTGCAC' : 76, 'TTGCAT' : 25 } ,
        'LI' : { 'CTAATA' : 25, 'CTAATC' : 22, 'CTAATT' : 53, 'CTCATA' : 10, 'CTCATC' : 60, 'CTCATT' : 30, 'CTGATA' : 23, 'CTGATC' : 36, 'CTGATT' : 41, 'CTTATA' : 22, 'CTTATC' : 25, 'CTTATT' : 53, 'TTAATA' : 28, 'TTAATC' : 23, 'TTAATT' : 50, 'TTGATA' : 25, 'TTGATC' : 27, 'TTGATT' : 48 } ,
        'LK' : { 'CTAAAA' : 48, 'CTAAAG' : 52, 'CTCAAA' : 32, 'CTCAAG' : 68, 'CTGAAA' : 42, 'CTGAAG' : 58, 'CTTAAA' : 52, 'CTTAAG' : 48, 'TTAAAA' : 52, 'TTAAAG' : 48, 'TTGAAA' : 49, 'TTGAAG' : 51 } ,
        'LL' : { 'CTACTA' : 13, 'CTACTC' : 9, 'CTACTG' : 33, 'CTACTT' : 17, 'CTATTA' : 14, 'CTATTG' : 15, 'CTCCTA' : 6, 'CTCCTC' : 22, 'CTCCTG' : 44, 'CTCCTT' : 11, 'CTCTTA' : 6, 'CTCTTG' : 12, 'CTGCTA' : 11, 'CTGCTC' : 11, 'CTGCTG' : 39, 'CTGCTT' : 14, 'CTGTTA' : 12, 'CTGTTG' : 14, 'CTTCTA' : 12, 'CTTCTC' : 9, 'CTTCTG' : 31, 'CTTCTT' : 16, 'CTTTTA' : 15, 'CTTTTG' : 17, 'TTACTA' : 14, 'TTACTC' : 6, 'TTACTG' : 27, 'TTACTT' : 18, 'TTATTA' : 18, 'TTATTG' : 18, 'TTGCTA' : 13, 'TTGCTC' : 7, 'TTGCTG' : 34, 'TTGCTT' : 16, 'TTGTTA' : 14, 'TTGTTG' : 16 } ,
        'LM' : { 'CTAATG' : 100, 'CTCATG' : 100, 'CTGATG' : 100, 'CTTATG' : 100, 'TTAATG' : 100, 'TTGATG' : 100 } ,
        'LN' : { 'CTAAAC' : 60, 'CTAAAT' : 40, 'CTCAAC' : 64, 'CTCAAT' : 36, 'CTGAAC' : 61, 'CTGAAT' : 39, 'CTTAAC' : 59, 'CTTAAT' : 41, 'TTAAAC' : 53, 'TTAAAT' : 47, 'TTGAAC' : 58, 'TTGAAT' : 42 } ,
        'LP' : { 'CTACCA' : 39, 'CTACCC' : 22, 'CTACCG' : 13, 'CTACCT' : 26, 'CTCCCA' : 27, 'CTCCCC' : 30, 'CTCCCG' : 15, 'CTCCCT' : 28, 'CTGCCA' : 27, 'CTGCCC' : 32, 'CTGCCG' : 18, 'CTGCCT' : 23, 'CTTCCA' : 39, 'CTTCCC' : 20, 'CTTCCG' : 13, 'CTTCCT' : 28, 'TTACCA' : 39, 'TTACCC' : 19, 'TTACCG' : 13, 'TTACCT' : 30, 'TTGCCA' : 38, 'TTGCCC' : 22, 'TTGCCG' : 14, 'TTGCCT' : 26 } ,
        'LQ' : { 'CTACAA' : 30, 'CTACAG' : 70, 'CTCCAA' : 22, 'CTCCAG' : 78, 'CTGCAA' : 26, 'CTGCAG' : 74, 'CTTCAA' : 32, 'CTTCAG' : 68, 'TTACAA' : 36, 'TTACAG' : 64, 'TTGCAA' : 32, 'TTGCAG' : 68 } ,
        'LR' : { 'CTAAGA' : 32, 'CTAAGG' : 23, 'CTACGA' : 12, 'CTACGC' : 7, 'CTACGG' : 16, 'CTACGT' : 10, 'CTCAGA' : 16, 'CTCAGG' : 23, 'CTCCGA' : 13, 'CTCCGC' : 18, 'CTCCGG' : 23, 'CTCCGT' : 8, 'CTGAGA' : 23, 'CTGAGG' : 21, 'CTGCGA' : 11, 'CTGCGC' : 12, 'CTGCGG' : 22, 'CTGCGT' : 11, 'CTTAGA' : 26, 'CTTAGG' : 17, 'CTTCGA' : 18, 'CTTCGC' : 8, 'CTTCGG' : 19, 'CTTCGT' : 13, 'TTAAGA' : 35, 'TTAAGG' : 23, 'TTACGA' : 14, 'TTACGC' : 4, 'TTACGG' : 14, 'TTACGT' : 10, 'TTGAGA' : 32, 'TTGAGG' : 24, 'TTGCGA' : 12, 'TTGCGC' : 6, 'TTGCGG' : 16, 'TTGCGT' : 10 } ,
        'LS' : { 'CTAAGC' : 28, 'CTAAGT' : 16, 'CTATCA' : 19, 'CTATCC' : 8, 'CTATCG' : 8, 'CTATCT' : 21, 'CTCAGC' : 30, 'CTCAGT' : 15, 'CTCTCA' : 12, 'CTCTCC' : 23, 'CTCTCG' : 6, 'CTCTCT' : 15, 'CTGAGC' : 28, 'CTGAGT' : 14, 'CTGTCA' : 16, 'CTGTCC' : 12, 'CTGTCG' : 12, 'CTGTCT' : 18, 'CTTAGC' : 21, 'CTTAGT' : 10, 'CTTTCA' : 21, 'CTTTCC' : 11, 'CTTTCG' : 9, 'CTTTCT' : 27, 'TTAAGC' : 21, 'TTAAGT' : 16, 'TTATCA' : 21, 'TTATCC' : 8, 'TTATCG' : 9, 'TTATCT' : 24, 'TTGAGC' : 24, 'TTGAGT' : 14, 'TTGTCA' : 20, 'TTGTCC' : 10, 'TTGTCG' : 10, 'TTGTCT' : 23 } ,
        'LT' : { 'CTAACA' : 31, 'CTAACC' : 13, 'CTAACG' : 15, 'CTAACT' : 42, 'CTCACA' : 21, 'CTCACC' : 43, 'CTCACG' : 13, 'CTCACT' : 23, 'CTGACA' : 27, 'CTGACC' : 20, 'CTGACG' : 20, 'CTGACT' : 33, 'CTTACA' : 28, 'CTTACC' : 14, 'CTTACG' : 14, 'CTTACT' : 45, 'TTAACA' : 34, 'TTAACC' : 11, 'TTAACG' : 14, 'TTAACT' : 41, 'TTGACA' : 30, 'TTGACC' : 14, 'TTGACG' : 16, 'TTGACT' : 40 } ,
        'LV' : { 'CTAGTA' : 20, 'CTAGTC' : 9, 'CTAGTG' : 53, 'CTAGTT' : 18, 'CTCGTA' : 6, 'CTCGTC' : 30, 'CTCGTG' : 52, 'CTCGTT' : 12, 'CTGGTA' : 17, 'CTGGTC' : 11, 'CTGGTG' : 56, 'CTGGTT' : 17, 'CTTGTA' : 19, 'CTTGTC' : 12, 'CTTGTG' : 49, 'CTTGTT' : 20, 'TTAGTA' : 23, 'TTAGTC' : 8, 'TTAGTG' : 46, 'TTAGTT' : 23, 'TTGGTA' : 19, 'TTGGTC' : 10, 'TTGGTG' : 51, 'TTGGTT' : 20 } ,
        'LW' : { 'CTATGG' : 100, 'CTCTGG' : 100, 'CTGTGG' : 100, 'CTTTGG' : 100, 'TTATGG' : 100, 'TTGTGG' : 100 } ,
        'LY' : { 'CTATAC' : 59, 'CTATAT' : 41, 'CTCTAC' : 63, 'CTCTAT' : 37, 'CTGTAC' : 60, 'CTGTAT' : 40, 'CTTTAC' : 56, 'CTTTAT' : 44, 'TTATAC' : 51, 'TTATAT' : 49, 'TTGTAC' : 57, 'TTGTAT' : 43 } ,
        'MA' : { 'ATGGCA' : 23, 'ATGGCC' : 36, 'ATGGCG' : 15, 'ATGGCT' : 26 } ,
        'MC' : { 'ATGTGC' : 51, 'ATGTGT' : 49 } ,
        'MU' : { 'ATGTGC' : 51, 'ATGTGT' : 49 } ,
        'MD' : { 'ATGGAC' : 51, 'ATGGAT' : 49 } ,
        'ME' : { 'ATGGAA' : 42, 'ATGGAG' : 58 } ,
        'MF' : { 'ATGTTC' : 49, 'ATGTTT' : 51 } ,
        'MG' : { 'ATGGGA' : 25, 'ATGGGC' : 34, 'ATGGGG' : 24, 'ATGGGT' : 17 } ,
        'MH' : { 'ATGCAC' : 58, 'ATGCAT' : 42 } ,
        'MI' : { 'ATGATA' : 19, 'ATGATC' : 43, 'ATGATT' : 38 } ,
        'MK' : { 'ATGAAA' : 41, 'ATGAAG' : 59 } ,
        'ML' : { 'ATGCTA' : 8, 'ATGCTC' : 15, 'ATGCTG' : 43, 'ATGCTT' : 13, 'ATGTTA' : 8, 'ATGTTG' : 13 } ,
        'MM' : { 'ATGATG' : 100 } ,
        'MN' : { 'ATGAAC' : 52, 'ATGAAT' : 48 } ,
        'MP' : { 'ATGCCA' : 27, 'ATGCCC' : 33, 'ATGCCG' : 11, 'ATGCCT' : 29 } ,
        'MQ' : { 'ATGCAA' : 27, 'ATGCAG' : 73 } ,
        'MR' : { 'ATGAGA' : 26, 'ATGAGG' : 27, 'ATGCGA' : 8, 'ATGCGC' : 14, 'ATGCGG' : 19, 'ATGCGT' : 7 } ,
        'MS' : { 'ATGAGC' : 22, 'ATGAGT' : 15, 'ATGTCA' : 15, 'ATGTCC' : 21, 'ATGTCG' : 7, 'ATGTCT' : 20 } ,
        'MT' : { 'ATGACA' : 28, 'ATGACC' : 35, 'ATGACG' : 10, 'ATGACT' : 26 } ,
        'MV' : { 'ATGGTA' : 12, 'ATGGTC' : 21, 'ATGGTG' : 49, 'ATGGTT' : 18 } ,
        'MW' : { 'ATGTGG' : 100 } ,
        'MY' : { 'ATGTAC' : 56, 'ATGTAT' : 44 } ,
        'NA' : { 'AACGCA' : 25, 'AACGCC' : 24, 'AACGCG' : 23, 'AACGCT' : 28, 'AATGCA' : 28, 'AATGCC' : 18, 'AATGCG' : 15, 'AATGCT' : 40 } ,
        'NC' : { 'AACTGC' : 60, 'AACTGT' : 40, 'AATTGC' : 71, 'AATTGT' : 29 } ,
        'NU' : { 'AACTGC' : 60, 'AACTGT' : 40, 'AATTGC' : 71, 'AATTGT' : 29 } ,
        'ND' : { 'AACGAC' : 83, 'AACGAT' : 17, 'AATGAC' : 75, 'AATGAT' : 26 } ,
        'NE' : { 'AACGAA' : 37, 'AACGAG' : 63, 'AATGAA' : 50, 'AATGAG' : 50 } ,
        'NF' : { 'AACTTC' : 80, 'AACTTT' : 20, 'AATTTC' : 72, 'AATTTT' : 28 } ,
        'NG' : { 'AACGGA' : 16, 'AACGGC' : 45, 'AACGGG' : 34, 'AACGGT' : 6, 'AATGGA' : 19, 'AATGGC' : 43, 'AATGGG' : 29, 'AATGGT' : 9 } ,
        'NH' : { 'AACCAC' : 82, 'AACCAT' : 19, 'AATCAC' : 73, 'AATCAT' : 27 } ,
        'NI' : { 'AACATA' : 22, 'AACATC' : 39, 'AACATT' : 39, 'AATATA' : 20, 'AATATC' : 17, 'AATATT' : 63 } ,
        'NK' : { 'AACAAA' : 46, 'AACAAG' : 54, 'AATAAA' : 54, 'AATAAG' : 46 } ,
        'NL' : { 'AACCTA' : 12, 'AACCTC' : 11, 'AACCTG' : 32, 'AACCTT' : 15, 'AACTTA' : 14, 'AACTTG' : 17, 'AATCTA' : 10, 'AATCTC' : 8, 'AATCTG' : 30, 'AATCTT' : 16, 'AATTTA' : 18, 'AATTTG' : 18 } ,
        'NM' : { 'AACATG' : 100, 'AATATG' : 100 } ,
        'NN' : { 'AACAAC' : 63, 'AACAAT' : 37, 'AATAAC' : 71, 'AATAAT' : 29 } ,
        'NP' : { 'AACCCA' : 31, 'AACCCC' : 28, 'AACCCG' : 17, 'AACCCT' : 24, 'AATCCA' : 56, 'AATCCC' : 12, 'AATCCG' : 9, 'AATCCT' : 23 } ,
        'NQ' : { 'AACCAA' : 31, 'AACCAG' : 69, 'AATCAA' : 32, 'AATCAG' : 68 } ,
        'NR' : { 'AACAGA' : 29, 'AACAGG' : 22, 'AACCGA' : 13, 'AACCGC' : 10, 'AACCGG' : 17, 'AACCGT' : 11, 'AATAGA' : 35, 'AATAGG' : 20, 'AATCGA' : 14, 'AATCGC' : 7, 'AATCGG' : 13, 'AATCGT' : 11 } ,
        'NS' : { 'AACAGC' : 29, 'AACAGT' : 15, 'AACTCA' : 17, 'AACTCC' : 9, 'AACTCG' : 11, 'AACTCT' : 18, 'AATAGC' : 28, 'AATAGT' : 7, 'AATTCA' : 23, 'AATTCC' : 9, 'AATTCG' : 6, 'AATTCT' : 26 } ,
        'NT' : { 'AACACA' : 30, 'AACACC' : 17, 'AACACG' : 19, 'AACACT' : 34, 'AATACA' : 30, 'AATACC' : 13, 'AATACG' : 12, 'AATACT' : 44 } ,
        'NV' : { 'AACGTA' : 16, 'AACGTC' : 14, 'AACGTG' : 55, 'AACGTT' : 16, 'AATGTA' : 18, 'AATGTC' : 12, 'AATGTG' : 53, 'AATGTT' : 18 } ,
        'NW' : { 'AACTGG' : 100, 'AATTGG' : 100 } ,
        'NY' : { 'AACTAC' : 61, 'AACTAT' : 39, 'AATTAC' : 72, 'AATTAT' : 28 } ,
        'PA' : { 'CCAGCA' : 29, 'CCAGCC' : 19, 'CCAGCG' : 19, 'CCAGCT' : 33, 'CCCGCA' : 21, 'CCCGCC' : 23, 'CCCGCG' : 26, 'CCCGCT' : 31, 'CCGGCA' : 21, 'CCGGCC' : 24, 'CCGGCG' : 23, 'CCGGCT' : 32, 'CCTGCA' : 27, 'CCTGCC' : 19, 'CCTGCG' : 17, 'CCTGCT' : 37 } ,
        'PC' : { 'CCATGC' : 57, 'CCATGT' : 43, 'CCCTGC' : 83, 'CCCTGT' : 17, 'CCGTGC' : 74, 'CCGTGT' : 26, 'CCTTGC' : 67, 'CCTTGT' : 33 } ,
        'PU' : { 'CCATGC' : 57, 'CCATGT' : 43, 'CCCTGC' : 83, 'CCCTGT' : 17, 'CCGTGC' : 74, 'CCGTGT' : 26, 'CCTTGC' : 67, 'CCTTGT' : 33 } ,
        'PD' : { 'CCAGAC' : 74, 'CCAGAT' : 26, 'CCCGAC' : 87, 'CCCGAT' : 14, 'CCGGAC' : 83, 'CCGGAT' : 17, 'CCTGAC' : 75, 'CCTGAT' : 25 } ,
        'PE' : { 'CCAGAA' : 46, 'CCAGAG' : 54, 'CCCGAA' : 32, 'CCCGAG' : 68, 'CCGGAA' : 35, 'CCGGAG' : 65, 'CCTGAA' : 43, 'CCTGAG' : 57 } ,
        'PF' : { 'CCATTC' : 70, 'CCATTT' : 31, 'CCCTTC' : 83, 'CCCTTT' : 17, 'CCGTTC' : 80, 'CCGTTT' : 20, 'CCTTTC' : 76, 'CCTTTT' : 24 } ,
        'PG' : { 'CCAGGA' : 21, 'CCAGGC' : 41, 'CCAGGG' : 30, 'CCAGGT' : 8, 'CCCGGA' : 11, 'CCCGGC' : 50, 'CCCGGG' : 35, 'CCCGGT' : 4, 'CCGGGA' : 13, 'CCGGGC' : 49, 'CCGGGG' : 33, 'CCGGGT' : 6, 'CCTGGA' : 17, 'CCTGGC' : 43, 'CCTGGG' : 32, 'CCTGGT' : 8 } ,
        'PH' : { 'CCACAC' : 76, 'CCACAT' : 24, 'CCCCAC' : 83, 'CCCCAT' : 17, 'CCGCAC' : 89, 'CCGCAT' : 12, 'CCTCAC' : 75, 'CCTCAT' : 25 } ,
        'PI' : { 'CCAATA' : 27, 'CCAATC' : 26, 'CCAATT' : 47, 'CCCATA' : 14, 'CCCATC' : 32, 'CCCATT' : 54, 'CCGATA' : 19, 'CCGATC' : 28, 'CCGATT' : 53, 'CCTATA' : 21, 'CCTATC' : 22, 'CCTATT' : 57 } ,
        'PK' : { 'CCAAAA' : 47, 'CCAAAG' : 53, 'CCCAAA' : 43, 'CCCAAG' : 57, 'CCGAAA' : 41, 'CCGAAG' : 59, 'CCTAAA' : 51, 'CCTAAG' : 49 } ,
        'PL' : { 'CCACTA' : 12, 'CCACTC' : 9, 'CCACTG' : 32, 'CCACTT' : 16, 'CCATTA' : 15, 'CCATTG' : 16, 'CCCCTA' : 8, 'CCCCTC' : 10, 'CCCCTG' : 41, 'CCCCTT' : 11, 'CCCTTA' : 13, 'CCCTTG' : 17, 'CCGCTA' : 8, 'CCGCTC' : 11, 'CCGCTG' : 46, 'CCGCTT' : 11, 'CCGTTA' : 12, 'CCGTTG' : 13, 'CCTCTA' : 10, 'CCTCTC' : 9, 'CCTCTG' : 34, 'CCTCTT' : 15, 'CCTTTA' : 15, 'CCTTTG' : 18 } ,
        'PM' : { 'CCAATG' : 100, 'CCCATG' : 100, 'CCGATG' : 100, 'CCTATG' : 100 } ,
        'PN' : { 'CCAAAC' : 55, 'CCAAAT' : 45, 'CCCAAC' : 82, 'CCCAAT' : 18, 'CCGAAC' : 71, 'CCGAAT' : 29, 'CCTAAC' : 66, 'CCTAAT' : 34 } ,
        'PP' : { 'CCACCA' : 34, 'CCACCC' : 24, 'CCACCG' : 16, 'CCACCT' : 26, 'CCCCCA' : 55, 'CCCCCC' : 10, 'CCCCCG' : 17, 'CCCCCT' : 18, 'CCGCCA' : 38, 'CCGCCC' : 22, 'CCGCCG' : 24, 'CCGCCT' : 16, 'CCTCCA' : 47, 'CCTCCC' : 17, 'CCTCCG' : 13, 'CCTCCT' : 24 } ,
        'PQ' : { 'CCACAA' : 31, 'CCACAG' : 69, 'CCCCAA' : 21, 'CCCCAG' : 79, 'CCGCAA' : 20, 'CCGCAG' : 80, 'CCTCAA' : 27, 'CCTCAG' : 73 } ,
        'PR' : { 'CCAAGA' : 31, 'CCAAGG' : 23, 'CCACGA' : 12, 'CCACGC' : 6, 'CCACGG' : 17, 'CCACGT' : 10, 'CCCAGA' : 29, 'CCCAGG' : 26, 'CCCCGA' : 11, 'CCCCGC' : 9, 'CCCCGG' : 18, 'CCCCGT' : 7, 'CCGAGA' : 24, 'CCGAGG' : 21, 'CCGCGA' : 10, 'CCGCGC' : 15, 'CCGCGG' : 22, 'CCGCGT' : 8, 'CCTAGA' : 27, 'CCTAGG' : 17, 'CCTCGA' : 16, 'CCTCGC' : 8, 'CCTCGG' : 21, 'CCTCGT' : 11 } ,
        'PS' : { 'CCAAGC' : 22, 'CCAAGT' : 15, 'CCATCA' : 20, 'CCATCC' : 11, 'CCATCG' : 11, 'CCATCT' : 22, 'CCCAGC' : 39, 'CCCAGT' : 8, 'CCCTCA' : 17, 'CCCTCC' : 9, 'CCCTCG' : 8, 'CCCTCT' : 19, 'CCGAGC' : 31, 'CCGAGT' : 8, 'CCGTCA' : 16, 'CCGTCC' : 14, 'CCGTCG' : 12, 'CCGTCT' : 19, 'CCTAGC' : 24, 'CCTAGT' : 7, 'CCTTCA' : 22, 'CCTTCC' : 13, 'CCTTCG' : 8, 'CCTTCT' : 27 } ,
        'PT' : { 'CCAACA' : 32, 'CCAACC' : 14, 'CCAACG' : 17, 'CCAACT' : 37, 'CCCACA' : 24, 'CCCACC' : 21, 'CCCACG' : 17, 'CCCACT' : 38, 'CCGACA' : 24, 'CCGACC' : 19, 'CCGACG' : 19, 'CCGACT' : 38, 'CCTACA' : 29, 'CCTACC' : 15, 'CCTACG' : 14, 'CCTACT' : 42 } ,
        'PV' : { 'CCAGTA' : 20, 'CCAGTC' : 10, 'CCAGTG' : 49, 'CCAGTT' : 21, 'CCCGTA' : 13, 'CCCGTC' : 13, 'CCCGTG' : 64, 'CCCGTT' : 10, 'CCGGTA' : 14, 'CCGGTC' : 12, 'CCGGTG' : 62, 'CCGGTT' : 12, 'CCTGTA' : 17, 'CCTGTC' : 12, 'CCTGTG' : 55, 'CCTGTT' : 17 } ,
        'PW' : { 'CCATGG' : 100, 'CCCTGG' : 100, 'CCGTGG' : 100, 'CCTTGG' : 100 } ,
        'PY' : { 'CCATAC' : 53, 'CCATAT' : 47, 'CCCTAC' : 81, 'CCCTAT' : 19, 'CCGTAC' : 70, 'CCGTAT' : 30, 'CCTTAC' : 66, 'CCTTAT' : 34 } ,
        'QA' : { 'CAAGCA' : 30, 'CAAGCC' : 17, 'CAAGCG' : 17, 'CAAGCT' : 36, 'CAGGCA' : 26, 'CAGGCC' : 21, 'CAGGCG' : 22, 'CAGGCT' : 31 } ,
        'QC' : { 'CAATGC' : 52, 'CAATGT' : 48, 'CAGTGC' : 55, 'CAGTGT' : 45 } ,
        'QU' : { 'CAATGC' : 52, 'CAATGT' : 48, 'CAGTGC' : 55, 'CAGTGT' : 45 } ,
        'QD' : { 'CAAGAC' : 72, 'CAAGAT' : 28, 'CAGGAC' : 80, 'CAGGAT' : 21 } ,
        'QE' : { 'CAAGAA' : 51, 'CAAGAG' : 49, 'CAGGAA' : 41, 'CAGGAG' : 59 } ,
        'QF' : { 'CAATTC' : 70, 'CAATTT' : 30, 'CAGTTC' : 76, 'CAGTTT' : 25 } ,
        'QG' : { 'CAAGGA' : 23, 'CAAGGC' : 40, 'CAAGGG' : 28, 'CAAGGT' : 9, 'CAGGGA' : 19, 'CAGGGC' : 45, 'CAGGGG' : 28, 'CAGGGT' : 9 } ,
        'QH' : { 'CAACAC' : 70, 'CAACAT' : 30, 'CAGCAC' : 81, 'CAGCAT' : 20 } ,
        'QI' : { 'CAAATA' : 27, 'CAAATC' : 26, 'CAAATT' : 47, 'CAGATA' : 24, 'CAGATC' : 38, 'CAGATT' : 38 } ,
        'QK' : { 'CAAAAA' : 49, 'CAAAAG' : 51, 'CAGAAA' : 45, 'CAGAAG' : 55 } ,
        'QL' : { 'CAACTA' : 14, 'CAACTC' : 9, 'CAACTG' : 29, 'CAACTT' : 19, 'CAATTA' : 15, 'CAATTG' : 14, 'CAGCTA' : 12, 'CAGCTC' : 11, 'CAGCTG' : 34, 'CAGCTT' : 15, 'CAGTTA' : 13, 'CAGTTG' : 14 } ,
        'QM' : { 'CAAATG' : 100, 'CAGATG' : 100 } ,
        'QN' : { 'CAAAAC' : 57, 'CAAAAT' : 43, 'CAGAAC' : 57, 'CAGAAT' : 43 } ,
        'QP' : { 'CAACCA' : 38, 'CAACCC' : 22, 'CAACCG' : 14, 'CAACCT' : 27, 'CAGCCA' : 28, 'CAGCCC' : 29, 'CAGCCG' : 18, 'CAGCCT' : 25 } ,
        'QQ' : { 'CAACAA' : 35, 'CAACAG' : 65, 'CAGCAA' : 31, 'CAGCAG' : 69 } ,
        'QR' : { 'CAAAGA' : 38, 'CAAAGG' : 25, 'CAACGA' : 11, 'CAACGC' : 4, 'CAACGG' : 12, 'CAACGT' : 10, 'CAGAGA' : 26, 'CAGAGG' : 23, 'CAGCGA' : 12, 'CAGCGC' : 9, 'CAGCGG' : 18, 'CAGCGT' : 11 } ,
        'QS' : { 'CAAAGC' : 29, 'CAAAGT' : 19, 'CAATCA' : 18, 'CAATCC' : 6, 'CAATCG' : 8, 'CAATCT' : 19, 'CAGAGC' : 27, 'CAGAGT' : 17, 'CAGTCA' : 17, 'CAGTCC' : 9, 'CAGTCG' : 11, 'CAGTCT' : 19 } ,
        'QT' : { 'CAAACA' : 32, 'CAAACC' : 13, 'CAAACG' : 15, 'CAAACT' : 40, 'CAGACA' : 29, 'CAGACC' : 18, 'CAGACG' : 21, 'CAGACT' : 32 } ,
        'QV' : { 'CAAGTA' : 20, 'CAAGTC' : 11, 'CAAGTG' : 49, 'CAAGTT' : 21, 'CAGGTA' : 18, 'CAGGTC' : 12, 'CAGGTG' : 51, 'CAGGTT' : 20 } ,
        'QW' : { 'CAATGG' : 100, 'CAGTGG' : 100 } ,
        'QY' : { 'CAATAC' : 54, 'CAATAT' : 46, 'CAGTAC' : 56, 'CAGTAT' : 44 } ,
        'RA' : { 'AGAGCA' : 31, 'AGAGCC' : 17, 'AGAGCG' : 19, 'AGAGCT' : 33, 'AGGGCA' : 28, 'AGGGCC' : 20, 'AGGGCG' : 20, 'AGGGCT' : 31, 'CGAGCA' : 28, 'CGAGCC' : 20, 'CGAGCG' : 19, 'CGAGCT' : 33, 'CGCGCA' : 9, 'CGCGCC' : 53, 'CGCGCG' : 27, 'CGCGCT' : 11, 'CGGGCA' : 26, 'CGGGCC' : 22, 'CGGGCG' : 23, 'CGGGCT' : 29, 'CGTGCA' : 25, 'CGTGCC' : 23, 'CGTGCG' : 19, 'CGTGCT' : 33 } ,
        'RC' : { 'AGATGC' : 48, 'AGATGT' : 52, 'AGGTGC' : 56, 'AGGTGT' : 44, 'CGATGC' : 55, 'CGATGT' : 45, 'CGCTGC' : 68, 'CGCTGT' : 32, 'CGGTGC' : 62, 'CGGTGT' : 38, 'CGTTGC' : 58, 'CGTTGT' : 42 } ,
        'RU' : { 'AGATGC' : 48, 'AGATGT' : 52, 'AGGTGC' : 56, 'AGGTGT' : 44, 'CGATGC' : 55, 'CGATGT' : 45, 'CGCTGC' : 68, 'CGCTGT' : 32, 'CGGTGC' : 62, 'CGGTGT' : 38, 'CGTTGC' : 58, 'CGTTGT' : 42 } ,
        'RD' : { 'AGAGAC' : 73, 'AGAGAT' : 27, 'AGGGAC' : 79, 'AGGGAT' : 21, 'CGAGAC' : 77, 'CGAGAT' : 24, 'CGCGAC' : 77, 'CGCGAT' : 23, 'CGGGAC' : 84, 'CGGGAT' : 17, 'CGTGAC' : 78, 'CGTGAT' : 22 } ,
        'RE' : { 'AGAGAA' : 51, 'AGAGAG' : 49, 'AGGGAA' : 43, 'AGGGAG' : 57, 'CGAGAA' : 45, 'CGAGAG' : 55, 'CGCGAA' : 12, 'CGCGAG' : 88, 'CGGGAA' : 37, 'CGGGAG' : 63, 'CGTGAA' : 43, 'CGTGAG' : 57 } ,
        'RF' : { 'AGATTC' : 69, 'AGATTT' : 31, 'AGGTTC' : 75, 'AGGTTT' : 26, 'CGATTC' : 73, 'CGATTT' : 28, 'CGCTTC' : 72, 'CGCTTT' : 28, 'CGGTTC' : 78, 'CGGTTT' : 22, 'CGTTTC' : 77, 'CGTTTT' : 24 } ,
        'RG' : { 'AGAGGA' : 28, 'AGAGGC' : 36, 'AGAGGG' : 27, 'AGAGGT' : 9, 'AGGGGA' : 21, 'AGGGGC' : 44, 'AGGGGG' : 26, 'AGGGGT' : 9, 'CGAGGA' : 19, 'CGAGGC' : 43, 'CGAGGG' : 30, 'CGAGGT' : 8, 'CGCGGA' : 10, 'CGCGGC' : 51, 'CGCGGG' : 31, 'CGCGGT' : 8, 'CGGGGA' : 15, 'CGGGGC' : 50, 'CGGGGG' : 27, 'CGGGGT' : 8, 'CGTGGA' : 16, 'CGTGGC' : 44, 'CGTGGG' : 30, 'CGTGGT' : 9 } ,
        'RH' : { 'AGACAC' : 75, 'AGACAT' : 26, 'AGGCAC' : 79, 'AGGCAT' : 21, 'CGACAC' : 75, 'CGACAT' : 25, 'CGCCAC' : 70, 'CGCCAT' : 30, 'CGGCAC' : 85, 'CGGCAT' : 15, 'CGTCAC' : 75, 'CGTCAT' : 25 } ,
        'RI' : { 'AGAATA' : 28, 'AGAATC' : 28, 'AGAATT' : 44, 'AGGATA' : 25, 'AGGATC' : 35, 'AGGATT' : 39, 'CGAATA' : 25, 'CGAATC' : 30, 'CGAATT' : 45, 'CGCATA' : 7, 'CGCATC' : 71, 'CGCATT' : 22, 'CGGATA' : 22, 'CGGATC' : 39, 'CGGATT' : 39, 'CGTATA' : 22, 'CGTATC' : 32, 'CGTATT' : 46 } ,
        'RK' : { 'AGAAAA' : 53, 'AGAAAG' : 47, 'AGGAAA' : 47, 'AGGAAG' : 53, 'CGAAAA' : 49, 'CGAAAG' : 51, 'CGCAAA' : 29, 'CGCAAG' : 71, 'CGGAAA' : 43, 'CGGAAG' : 57, 'CGTAAA' : 50, 'CGTAAG' : 50 } ,
        'RL' : { 'AGACTA' : 14, 'AGACTC' : 9, 'AGACTG' : 25, 'AGACTT' : 19, 'AGATTA' : 17, 'AGATTG' : 16, 'AGGCTA' : 13, 'AGGCTC' : 10, 'AGGCTG' : 33, 'AGGCTT' : 16, 'AGGTTA' : 14, 'AGGTTG' : 14, 'CGACTA' : 13, 'CGACTC' : 10, 'CGACTG' : 30, 'CGACTT' : 17, 'CGATTA' : 15, 'CGATTG' : 15, 'CGCCTA' : 5, 'CGCCTC' : 25, 'CGCCTG' : 49, 'CGCCTT' : 9, 'CGCTTA' : 3, 'CGCTTG' : 8, 'CGGCTA' : 11, 'CGGCTC' : 12, 'CGGCTG' : 40, 'CGGCTT' : 13, 'CGGTTA' : 12, 'CGGTTG' : 13, 'CGTCTA' : 12, 'CGTCTC' : 9, 'CGTCTG' : 31, 'CGTCTT' : 16, 'CGTTTA' : 15, 'CGTTTG' : 17 } ,
        'RM' : { 'AGAATG' : 100, 'AGGATG' : 100, 'CGAATG' : 100, 'CGCATG' : 100, 'CGGATG' : 100, 'CGTATG' : 100 } ,
        'RN' : { 'AGAAAC' : 50, 'AGAAAT' : 50, 'AGGAAC' : 58, 'AGGAAT' : 42, 'CGAAAC' : 57, 'CGAAAT' : 43, 'CGCAAC' : 70, 'CGCAAT' : 30, 'CGGAAC' : 65, 'CGGAAT' : 35, 'CGTAAC' : 60, 'CGTAAT' : 40 } ,
        'RP' : { 'AGACCA' : 30, 'AGACCC' : 26, 'AGACCG' : 15, 'AGACCT' : 29, 'AGGCCA' : 29, 'AGGCCC' : 28, 'AGGCCG' : 18, 'AGGCCT' : 25, 'CGACCA' : 34, 'CGACCC' : 25, 'CGACCG' : 15, 'CGACCT' : 26, 'CGCCCA' : 23, 'CGCCCC' : 37, 'CGCCCG' : 19, 'CGCCCT' : 21, 'CGGCCA' : 27, 'CGGCCC' : 32, 'CGGCCG' : 20, 'CGGCCT' : 21, 'CGTCCA' : 34, 'CGTCCC' : 25, 'CGTCCG' : 15, 'CGTCCT' : 26 } ,
        'RQ' : { 'AGACAA' : 42, 'AGACAG' : 58, 'AGGCAA' : 34, 'AGGCAG' : 66, 'CGACAA' : 33, 'CGACAG' : 67, 'CGCCAA' : 16, 'CGCCAG' : 84, 'CGGCAA' : 25, 'CGGCAG' : 75, 'CGTCAA' : 30, 'CGTCAG' : 70 } ,
        'RR' : { 'AGAAGA' : 34, 'AGAAGG' : 23, 'AGACGA' : 13, 'AGACGC' : 5, 'AGACGG' : 14, 'AGACGT' : 11, 'AGGAGA' : 29, 'AGGAGG' : 26, 'AGGCGA' : 12, 'AGGCGC' : 7, 'AGGCGG' : 17, 'AGGCGT' : 10, 'CGAAGA' : 32, 'CGAAGG' : 24, 'CGACGA' : 12, 'CGACGC' : 7, 'CGACGG' : 15, 'CGACGT' : 10, 'CGCAGA' : 11, 'CGCAGG' : 21, 'CGCCGA' : 8, 'CGCCGC' : 31, 'CGCCGG' : 23, 'CGCCGT' : 7, 'CGGAGA' : 23, 'CGGAGG' : 22, 'CGGCGA' : 11, 'CGGCGC' : 12, 'CGGCGG' : 22, 'CGGCGT' : 10, 'CGTAGA' : 25, 'CGTAGG' : 18, 'CGTCGA' : 14, 'CGTCGC' : 11, 'CGTCGG' : 20, 'CGTCGT' : 13 } ,
        'RS' : { 'AGAAGC' : 24, 'AGAAGT' : 21, 'AGATCA' : 19, 'AGATCC' : 6, 'AGATCG' : 11, 'AGATCT' : 19, 'AGGAGC' : 29, 'AGGAGT' : 16, 'AGGTCA' : 17, 'AGGTCC' : 9, 'AGGTCG' : 11, 'AGGTCT' : 18, 'CGAAGC' : 28, 'CGAAGT' : 18, 'CGATCA' : 17, 'CGATCC' : 8, 'CGATCG' : 9, 'CGATCT' : 19, 'CGCAGC' : 34, 'CGCAGT' : 10, 'CGCTCA' : 10, 'CGCTCC' : 26, 'CGCTCG' : 10, 'CGCTCT' : 11, 'CGGAGC' : 33, 'CGGAGT' : 15, 'CGGTCA' : 16, 'CGGTCC' : 9, 'CGGTCG' : 11, 'CGGTCT' : 17, 'CGTAGC' : 21, 'CGTAGT' : 11, 'CGTTCA' : 19, 'CGTTCC' : 14, 'CGTTCG' : 10, 'CGTTCT' : 25 } ,
        'RT' : { 'AGAACA' : 34, 'AGAACC' : 12, 'AGAACG' : 18, 'AGAACT' : 36, 'AGGACA' : 30, 'AGGACC' : 17, 'AGGACG' : 20, 'AGGACT' : 34, 'CGAACA' : 30, 'CGAACC' : 15, 'CGAACG' : 16, 'CGAACT' : 38, 'CGCACA' : 22, 'CGCACC' : 42, 'CGCACG' : 21, 'CGCACT' : 14, 'CGGACA' : 27, 'CGGACC' : 20, 'CGGACG' : 21, 'CGGACT' : 33, 'CGTACA' : 27, 'CGTACC' : 18, 'CGTACG' : 16, 'CGTACT' : 39 } ,
        'RV' : { 'AGAGTA' : 23, 'AGAGTC' : 12, 'AGAGTG' : 39, 'AGAGTT' : 26, 'AGGGTA' : 19, 'AGGGTC' : 16, 'AGGGTG' : 45, 'AGGGTT' : 20, 'CGAGTA' : 20, 'CGAGTC' : 12, 'CGAGTG' : 49, 'CGAGTT' : 19, 'CGCGTA' : 4, 'CGCGTC' : 27, 'CGCGTG' : 62, 'CGCGTT' : 7, 'CGGGTA' : 17, 'CGGGTC' : 14, 'CGGGTG' : 53, 'CGGGTT' : 16, 'CGTGTA' : 18, 'CGTGTC' : 12, 'CGTGTG' : 53, 'CGTGTT' : 18 } ,
        'RW' : { 'AGATGG' : 100, 'AGGTGG' : 100, 'CGATGG' : 100, 'CGCTGG' : 100, 'CGGTGG' : 100, 'CGTTGG' : 100 } ,
        'RY' : { 'AGATAC' : 54, 'AGATAT' : 46, 'AGGTAC' : 57, 'AGGTAT' : 43, 'CGATAC' : 58, 'CGATAT' : 42, 'CGCTAC' : 68, 'CGCTAT' : 32, 'CGGTAC' : 63, 'CGGTAT' : 37, 'CGTTAC' : 59, 'CGTTAT' : 41 } ,
        'SA' : { 'AGCGCA' : 22, 'AGCGCC' : 27, 'AGCGCG' : 24, 'AGCGCT' : 28, 'AGTGCA' : 26, 'AGTGCC' : 20, 'AGTGCG' : 15, 'AGTGCT' : 39, 'TCAGCA' : 30, 'TCAGCC' : 18, 'TCAGCG' : 18, 'TCAGCT' : 35, 'TCCGCA' : 17, 'TCCGCC' : 45, 'TCCGCG' : 21, 'TCCGCT' : 17, 'TCGGCA' : 21, 'TCGGCC' : 25, 'TCGGCG' : 22, 'TCGGCT' : 32, 'TCTGCA' : 30, 'TCTGCC' : 17, 'TCTGCG' : 18, 'TCTGCT' : 36 } ,
        'SC' : { 'AGCTGC' : 62, 'AGCTGT' : 38, 'AGTTGC' : 70, 'AGTTGT' : 30, 'TCATGC' : 59, 'TCATGT' : 41, 'TCCTGC' : 62, 'TCCTGT' : 38, 'TCGTGC' : 70, 'TCGTGT' : 30, 'TCTTGC' : 56, 'TCTTGT' : 44 } ,
        'SU' : { 'AGCTGC' : 62, 'AGCTGT' : 38, 'AGTTGC' : 70, 'AGTTGT' : 30, 'TCATGC' : 59, 'TCATGT' : 41, 'TCCTGC' : 62, 'TCCTGT' : 38, 'TCGTGC' : 70, 'TCGTGT' : 30, 'TCTTGC' : 56, 'TCTTGT' : 44 } ,
        'SD' : { 'AGCGAC' : 85, 'AGCGAT' : 16, 'AGTGAC' : 76, 'AGTGAT' : 24, 'TCAGAC' : 73, 'TCAGAT' : 28, 'TCCGAC' : 67, 'TCCGAT' : 33, 'TCGGAC' : 82, 'TCGGAT' : 19, 'TCTGAC' : 73, 'TCTGAT' : 27 } ,
        'SE' : { 'AGCGAA' : 34, 'AGCGAG' : 66, 'AGTGAA' : 48, 'AGTGAG' : 52, 'TCAGAA' : 48, 'TCAGAG' : 52, 'TCCGAA' : 26, 'TCCGAG' : 74, 'TCGGAA' : 37, 'TCGGAG' : 63, 'TCTGAA' : 47, 'TCTGAG' : 53 } ,
        'SF' : { 'AGCTTC' : 82, 'AGCTTT' : 19, 'AGTTTC' : 73, 'AGTTTT' : 28, 'TCATTC' : 68, 'TCATTT' : 32, 'TCCTTC' : 63, 'TCCTTT' : 37, 'TCGTTC' : 75, 'TCGTTT' : 25, 'TCTTTC' : 74, 'TCTTTT' : 26 } ,
        'SG' : { 'AGCGGA' : 13, 'AGCGGC' : 49, 'AGCGGG' : 32, 'AGCGGT' : 6, 'AGTGGA' : 17, 'AGTGGC' : 44, 'AGTGGG' : 30, 'AGTGGT' : 9, 'TCAGGA' : 21, 'TCAGGC' : 40, 'TCAGGG' : 30, 'TCAGGT' : 9, 'TCCGGA' : 18, 'TCCGGC' : 36, 'TCCGGG' : 35, 'TCCGGT' : 12, 'TCGGGA' : 13, 'TCGGGC' : 49, 'TCGGGG' : 33, 'TCGGGT' : 5, 'TCTGGA' : 21, 'TCTGGC' : 39, 'TCTGGG' : 31, 'TCTGGT' : 9 } ,
        'SH' : { 'AGCCAC' : 83, 'AGCCAT' : 17, 'AGTCAC' : 75, 'AGTCAT' : 25, 'TCACAC' : 76, 'TCACAT' : 25, 'TCCCAC' : 62, 'TCCCAT' : 38, 'TCGCAC' : 87, 'TCGCAT' : 13, 'TCTCAC' : 73, 'TCTCAT' : 27 } ,
        'SI' : { 'AGCATA' : 21, 'AGCATC' : 40, 'AGCATT' : 39, 'AGTATA' : 19, 'AGTATC' : 18, 'AGTATT' : 63, 'TCAATA' : 28, 'TCAATC' : 23, 'TCAATT' : 49, 'TCCATA' : 13, 'TCCATC' : 57, 'TCCATT' : 30, 'TCGATA' : 21, 'TCGATC' : 26, 'TCGATT' : 53, 'TCTATA' : 25, 'TCTATC' : 27, 'TCTATT' : 48 } ,
        'SK' : { 'AGCAAA' : 46, 'AGCAAG' : 54, 'AGTAAA' : 55, 'AGTAAG' : 45, 'TCAAAA' : 50, 'TCAAAG' : 50, 'TCCAAA' : 41, 'TCCAAG' : 59, 'TCGAAA' : 44, 'TCGAAG' : 56, 'TCTAAA' : 51, 'TCTAAG' : 49 } ,
        'SL' : { 'AGCCTA' : 11, 'AGCCTC' : 13, 'AGCCTG' : 34, 'AGCCTT' : 15, 'AGCTTA' : 13, 'AGCTTG' : 15, 'AGTCTA' : 9, 'AGTCTC' : 8, 'AGTCTG' : 31, 'AGTCTT' : 15, 'AGTTTA' : 18, 'AGTTTG' : 19, 'TCACTA' : 11, 'TCACTC' : 7, 'TCACTG' : 31, 'TCACTT' : 16, 'TCATTA' : 17, 'TCATTG' : 18, 'TCCCTA' : 6, 'TCCCTC' : 19, 'TCCCTG' : 41, 'TCCCTT' : 10, 'TCCTTA' : 7, 'TCCTTG' : 16, 'TCGCTA' : 9, 'TCGCTC' : 10, 'TCGCTG' : 45, 'TCGCTT' : 11, 'TCGTTA' : 12, 'TCGTTG' : 13, 'TCTCTA' : 12, 'TCTCTC' : 8, 'TCTCTG' : 29, 'TCTCTT' : 16, 'TCTTTA' : 16, 'TCTTTG' : 19 } ,
        'SM' : { 'AGCATG' : 100, 'AGTATG' : 100, 'TCAATG' : 100, 'TCCATG' : 100, 'TCGATG' : 100, 'TCTATG' : 100 } ,
        'SN' : { 'AGCAAC' : 64, 'AGCAAT' : 36, 'AGTAAC' : 73, 'AGTAAT' : 27, 'TCAAAC' : 57, 'TCAAAT' : 43, 'TCCAAC' : 62, 'TCCAAT' : 38, 'TCGAAC' : 68, 'TCGAAT' : 32, 'TCTAAC' : 58, 'TCTAAT' : 42 } ,
        'SP' : { 'AGCCCA' : 29, 'AGCCCC' : 29, 'AGCCCG' : 18, 'AGCCCT' : 24, 'AGTCCA' : 52, 'AGTCCC' : 15, 'AGTCCG' : 9, 'AGTCCT' : 24, 'TCACCA' : 39, 'TCACCC' : 22, 'TCACCG' : 14, 'TCACCT' : 25, 'TCCCCA' : 33, 'TCCCCC' : 23, 'TCCCCG' : 18, 'TCCCCT' : 25, 'TCGCCA' : 37, 'TCGCCC' : 26, 'TCGCCG' : 19, 'TCGCCT' : 18, 'TCTCCA' : 38, 'TCTCCC' : 21, 'TCTCCG' : 14, 'TCTCCT' : 27 } ,
        'SQ' : { 'AGCCAA' : 30, 'AGCCAG' : 70, 'AGTCAA' : 30, 'AGTCAG' : 70, 'TCACAA' : 34, 'TCACAG' : 66, 'TCCCAA' : 22, 'TCCCAG' : 78, 'TCGCAA' : 22, 'TCGCAG' : 78, 'TCTCAA' : 33, 'TCTCAG' : 67 } ,
        'SR' : { 'AGCAGA' : 26, 'AGCAGG' : 22, 'AGCCGA' : 12, 'AGCCGC' : 12, 'AGCCGG' : 18, 'AGCCGT' : 11, 'AGTAGA' : 32, 'AGTAGG' : 19, 'AGTCGA' : 14, 'AGTCGC' : 9, 'AGTCGG' : 16, 'AGTCGT' : 11, 'TCAAGA' : 36, 'TCAAGG' : 23, 'TCACGA' : 11, 'TCACGC' : 5, 'TCACGG' : 14, 'TCACGT' : 10, 'TCCAGA' : 21, 'TCCAGG' : 26, 'TCCCGA' : 10, 'TCCCGC' : 17, 'TCCCGG' : 21, 'TCCCGT' : 6, 'TCGAGA' : 26, 'TCGAGG' : 22, 'TCGCGA' : 9, 'TCGCGC' : 13, 'TCGCGG' : 21, 'TCGCGT' : 8, 'TCTAGA' : 27, 'TCTAGG' : 19, 'TCTCGA' : 16, 'TCTCGC' : 7, 'TCTCGG' : 19, 'TCTCGT' : 12 } ,
        'SS' : { 'AGCAGC' : 31, 'AGCAGT' : 15, 'AGCTCA' : 15, 'AGCTCC' : 10, 'AGCTCG' : 11, 'AGCTCT' : 17, 'AGTAGC' : 29, 'AGTAGT' : 7, 'AGTTCA' : 22, 'AGTTCC' : 10, 'AGTTCG' : 6, 'AGTTCT' : 26, 'TCAAGC' : 22, 'TCAAGT' : 13, 'TCATCA' : 21, 'TCATCC' : 10, 'TCATCG' : 9, 'TCATCT' : 25, 'TCCAGC' : 25, 'TCCAGT' : 14, 'TCCTCA' : 15, 'TCCTCC' : 23, 'TCCTCG' : 7, 'TCCTCT' : 15, 'TCGAGC' : 26, 'TCGAGT' : 8, 'TCGTCA' : 17, 'TCGTCC' : 16, 'TCGTCG' : 13, 'TCGTCT' : 20, 'TCTAGC' : 19, 'TCTAGT' : 10, 'TCTTCA' : 23, 'TCTTCC' : 12, 'TCTTCG' : 10, 'TCTTCT' : 25 } ,
        'ST' : { 'AGCACA' : 29, 'AGCACC' : 20, 'AGCACG' : 18, 'AGCACT' : 33, 'AGTACA' : 29, 'AGTACC' : 14, 'AGTACG' : 12, 'AGTACT' : 45, 'TCAACA' : 32, 'TCAACC' : 12, 'TCAACG' : 15, 'TCAACT' : 42, 'TCCACA' : 27, 'TCCACC' : 37, 'TCCACG' : 15, 'TCCACT' : 21, 'TCGACA' : 25, 'TCGACC' : 17, 'TCGACG' : 17, 'TCGACT' : 41, 'TCTACA' : 31, 'TCTACC' : 14, 'TCTACG' : 15, 'TCTACT' : 40 } ,
        'SV' : { 'AGCGTA' : 15, 'AGCGTC' : 16, 'AGCGTG' : 54, 'AGCGTT' : 15, 'AGTGTA' : 17, 'AGTGTC' : 13, 'AGTGTG' : 53, 'AGTGTT' : 17, 'TCAGTA' : 20, 'TCAGTC' : 10, 'TCAGTG' : 51, 'TCAGTT' : 19, 'TCCGTA' : 7, 'TCCGTC' : 25, 'TCCGTG' : 58, 'TCCGTT' : 10, 'TCGGTA' : 14, 'TCGGTC' : 11, 'TCGGTG' : 63, 'TCGGTT' : 12, 'TCTGTA' : 18, 'TCTGTC' : 11, 'TCTGTG' : 51, 'TCTGTT' : 19 } ,
        'SW' : { 'AGCTGG' : 100, 'AGTTGG' : 100, 'TCATGG' : 100, 'TCCTGG' : 100, 'TCGTGG' : 100, 'TCTTGG' : 100 } ,
        'SY' : { 'AGCTAC' : 62, 'AGCTAT' : 38, 'AGTTAC' : 76, 'AGTTAT' : 24, 'TCATAC' : 56, 'TCATAT' : 44, 'TCCTAC' : 61, 'TCCTAT' : 39, 'TCGTAC' : 70, 'TCGTAT' : 30, 'TCTTAC' : 58, 'TCTTAT' : 42 } ,
        'TA' : { 'ACAGCA' : 29, 'ACAGCC' : 18, 'ACAGCG' : 18, 'ACAGCT' : 34, 'ACCGCA' : 16, 'ACCGCC' : 50, 'ACCGCG' : 17, 'ACCGCT' : 17, 'ACGGCA' : 24, 'ACGGCC' : 24, 'ACGGCG' : 22, 'ACGGCT' : 30, 'ACTGCA' : 30, 'ACTGCC' : 17, 'ACTGCG' : 19, 'ACTGCT' : 34 } ,
        'TC' : { 'ACATGC' : 55, 'ACATGT' : 45, 'ACCTGC' : 62, 'ACCTGT' : 38, 'ACGTGC' : 63, 'ACGTGT' : 37, 'ACTTGC' : 51, 'ACTTGT' : 49 } ,
        'TU' : { 'ACATGC' : 55, 'ACATGT' : 45, 'ACCTGC' : 62, 'ACCTGT' : 38, 'ACGTGC' : 63, 'ACGTGT' : 37, 'ACTTGC' : 51, 'ACTTGT' : 49 } ,
        'TD' : { 'ACAGAC' : 74, 'ACAGAT' : 27, 'ACCGAC' : 70, 'ACCGAT' : 30, 'ACGGAC' : 81, 'ACGGAT' : 19, 'ACTGAC' : 74, 'ACTGAT' : 26 } ,
        'TE' : { 'ACAGAA' : 48, 'ACAGAG' : 52, 'ACCGAA' : 25, 'ACCGAG' : 75, 'ACGGAA' : 38, 'ACGGAG' : 62, 'ACTGAA' : 49, 'ACTGAG' : 51 } ,
        'TF' : { 'ACATTC' : 69, 'ACATTT' : 31, 'ACCTTC' : 66, 'ACCTTT' : 34, 'ACGTTC' : 75, 'ACGTTT' : 25, 'ACTTTC' : 74, 'ACTTTT' : 26 } ,
        'TG' : { 'ACAGGA' : 21, 'ACAGGC' : 40, 'ACAGGG' : 28, 'ACAGGT' : 10, 'ACCGGA' : 17, 'ACCGGC' : 42, 'ACCGGG' : 30, 'ACCGGT' : 11, 'ACGGGA' : 14, 'ACGGGC' : 48, 'ACGGGG' : 31, 'ACGGGT' : 7, 'ACTGGA' : 27, 'ACTGGC' : 37, 'ACTGGG' : 28, 'ACTGGT' : 8 } ,
        'TH' : { 'ACACAC' : 76, 'ACACAT' : 24, 'ACCCAC' : 65, 'ACCCAT' : 35, 'ACGCAC' : 87, 'ACGCAT' : 14, 'ACTCAC' : 77, 'ACTCAT' : 24 } ,
        'TI' : { 'ACAATA' : 27, 'ACAATC' : 25, 'ACAATT' : 47, 'ACCATA' : 12, 'ACCATC' : 60, 'ACCATT' : 28, 'ACGATA' : 23, 'ACGATC' : 31, 'ACGATT' : 46, 'ACTATA' : 26, 'ACTATC' : 29, 'ACTATT' : 45 } ,
        'TK' : { 'ACAAAA' : 50, 'ACAAAG' : 50, 'ACCAAA' : 37, 'ACCAAG' : 63, 'ACGAAA' : 44, 'ACGAAG' : 56, 'ACTAAA' : 54, 'ACTAAG' : 46 } ,
        'TL' : { 'ACACTA' : 12, 'ACACTC' : 8, 'ACACTG' : 31, 'ACACTT' : 16, 'ACATTA' : 17, 'ACATTG' : 17, 'ACCCTA' : 5, 'ACCCTC' : 21, 'ACCCTG' : 41, 'ACCCTT' : 11, 'ACCTTA' : 6, 'ACCTTG' : 15, 'ACGCTA' : 10, 'ACGCTC' : 9, 'ACGCTG' : 43, 'ACGCTT' : 11, 'ACGTTA' : 13, 'ACGTTG' : 14, 'ACTCTA' : 13, 'ACTCTC' : 7, 'ACTCTG' : 27, 'ACTCTT' : 16, 'ACTTTA' : 16, 'ACTTTG' : 19 } ,
        'TM' : { 'ACAATG' : 100, 'ACCATG' : 100, 'ACGATG' : 100, 'ACTATG' : 100 } ,
        'TN' : { 'ACAAAC' : 56, 'ACAAAT' : 44, 'ACCAAC' : 63, 'ACCAAT' : 37, 'ACGAAC' : 64, 'ACGAAT' : 36, 'ACTAAC' : 53, 'ACTAAT' : 47 } ,
        'TP' : { 'ACACCA' : 36, 'ACACCC' : 23, 'ACACCG' : 15, 'ACACCT' : 26, 'ACCCCA' : 33, 'ACCCCC' : 26, 'ACCCCG' : 14, 'ACCCCT' : 27, 'ACGCCA' : 32, 'ACGCCC' : 29, 'ACGCCG' : 20, 'ACGCCT' : 20, 'ACTCCA' : 34, 'ACTCCC' : 23, 'ACTCCG' : 15, 'ACTCCT' : 28 } ,
        'TQ' : { 'ACACAA' : 33, 'ACACAG' : 67, 'ACCCAA' : 22, 'ACCCAG' : 78, 'ACGCAA' : 23, 'ACGCAG' : 77, 'ACTCAA' : 36, 'ACTCAG' : 64 } ,
        'TR' : { 'ACAAGA' : 33, 'ACAAGG' : 23, 'ACACGA' : 12, 'ACACGC' : 5, 'ACACGG' : 15, 'ACACGT' : 10, 'ACCAGA' : 20, 'ACCAGG' : 26, 'ACCCGA' : 9, 'ACCCGC' : 18, 'ACCCGG' : 20, 'ACCCGT' : 7, 'ACGAGA' : 25, 'ACGAGG' : 21, 'ACGCGA' : 10, 'ACGCGC' : 13, 'ACGCGG' : 22, 'ACGCGT' : 9, 'ACTAGA' : 25, 'ACTAGG' : 17, 'ACTCGA' : 17, 'ACTCGC' : 8, 'ACTCGG' : 19, 'ACTCGT' : 14 } ,
        'TS' : { 'ACAAGC' : 22, 'ACAAGT' : 14, 'ACATCA' : 21, 'ACATCC' : 9, 'ACATCG' : 10, 'ACATCT' : 23, 'ACCAGC' : 30, 'ACCAGT' : 14, 'ACCTCA' : 14, 'ACCTCC' : 20, 'ACCTCG' : 7, 'ACCTCT' : 14, 'ACGAGC' : 24, 'ACGAGT' : 10, 'ACGTCA' : 19, 'ACGTCC' : 14, 'ACGTCG' : 13, 'ACGTCT' : 20, 'ACTAGC' : 17, 'ACTAGT' : 13, 'ACTTCA' : 23, 'ACTTCC' : 11, 'ACTTCG' : 11, 'ACTTCT' : 25 } ,
        'TT' : { 'ACAACA' : 32, 'ACAACC' : 14, 'ACAACG' : 15, 'ACAACT' : 39, 'ACCACA' : 25, 'ACCACC' : 42, 'ACCACG' : 13, 'ACCACT' : 20, 'ACGACA' : 27, 'ACGACC' : 18, 'ACGACG' : 19, 'ACGACT' : 36, 'ACTACA' : 33, 'ACTACC' : 14, 'ACTACG' : 18, 'ACTACT' : 36 } ,
        'TV' : { 'ACAGTA' : 19, 'ACAGTC' : 11, 'ACAGTG' : 51, 'ACAGTT' : 19, 'ACCGTA' : 6, 'ACCGTC' : 28, 'ACCGTG' : 56, 'ACCGTT' : 10, 'ACGGTA' : 15, 'ACGGTC' : 12, 'ACGGTG' : 60, 'ACGGTT' : 14, 'ACTGTA' : 20, 'ACTGTC' : 12, 'ACTGTG' : 47, 'ACTGTT' : 21 } ,
        'TW' : { 'ACATGG' : 100, 'ACCTGG' : 100, 'ACGTGG' : 100, 'ACTTGG' : 100 } ,
        'TY' : { 'ACATAC' : 57, 'ACATAT' : 43, 'ACCTAC' : 67, 'ACCTAT' : 33, 'ACGTAC' : 65, 'ACGTAT' : 35, 'ACTTAC' : 57, 'ACTTAT' : 43 } ,
        'VA' : { 'GTAGCA' : 31, 'GTAGCC' : 16, 'GTAGCG' : 17, 'GTAGCT' : 36, 'GTCGCA' : 16, 'GTCGCC' : 48, 'GTCGCG' : 15, 'GTCGCT' : 21, 'GTGGCA' : 26, 'GTGGCC' : 23, 'GTGGCG' : 21, 'GTGGCT' : 31, 'GTTGCA' : 29, 'GTTGCC' : 15, 'GTTGCG' : 16, 'GTTGCT' : 39 } ,
        'VC' : { 'GTATGC' : 50, 'GTATGT' : 50, 'GTCTGC' : 59, 'GTCTGT' : 41, 'GTGTGC' : 55, 'GTGTGT' : 45, 'GTTTGC' : 55, 'GTTTGT' : 45 } ,
        'VU' : { 'GTATGC' : 50, 'GTATGT' : 50, 'GTCTGC' : 59, 'GTCTGT' : 41, 'GTGTGC' : 55, 'GTGTGT' : 45, 'GTTTGC' : 55, 'GTTTGT' : 45 } ,
        'VD' : { 'GTAGAC' : 71, 'GTAGAT' : 30, 'GTCGAC' : 59, 'GTCGAT' : 41, 'GTGGAC' : 79, 'GTGGAT' : 21, 'GTTGAC' : 72, 'GTTGAT' : 28 } ,
        'VE' : { 'GTAGAA' : 52, 'GTAGAG' : 48, 'GTCGAA' : 33, 'GTCGAG' : 67, 'GTGGAA' : 41, 'GTGGAG' : 59, 'GTTGAA' : 53, 'GTTGAG' : 47 } ,
        'VF' : { 'GTATTC' : 66, 'GTATTT' : 34, 'GTCTTC' : 64, 'GTCTTT' : 36, 'GTGTTC' : 75, 'GTGTTT' : 25, 'GTTTTC' : 75, 'GTTTTT' : 26 } ,
        'VG' : { 'GTAGGA' : 25, 'GTAGGC' : 37, 'GTAGGG' : 28, 'GTAGGT' : 10, 'GTCGGA' : 17, 'GTCGGC' : 34, 'GTCGGG' : 35, 'GTCGGT' : 14, 'GTGGGA' : 20, 'GTGGGC' : 43, 'GTGGGG' : 29, 'GTGGGT' : 8, 'GTTGGA' : 22, 'GTTGGC' : 39, 'GTTGGG' : 28, 'GTTGGT' : 11 } ,
        'VH' : { 'GTACAC' : 72, 'GTACAT' : 28, 'GTCCAC' : 62, 'GTCCAT' : 38, 'GTGCAC' : 84, 'GTGCAT' : 17, 'GTTCAC' : 73, 'GTTCAT' : 27 } ,
        'VI' : { 'GTAATA' : 28, 'GTAATC' : 24, 'GTAATT' : 48, 'GTCATA' : 11, 'GTCATC' : 58, 'GTCATT' : 31, 'GTGATA' : 24, 'GTGATC' : 38, 'GTGATT' : 38, 'GTTATA' : 24, 'GTTATC' : 26, 'GTTATT' : 50 } ,
        'VK' : { 'GTAAAA' : 52, 'GTAAAG' : 48, 'GTCAAA' : 38, 'GTCAAG' : 62, 'GTGAAA' : 44, 'GTGAAG' : 56, 'GTTAAA' : 55, 'GTTAAG' : 45 } ,
        'VL' : { 'GTACTA' : 13, 'GTACTC' : 6, 'GTACTG' : 26, 'GTACTT' : 18, 'GTATTA' : 18, 'GTATTG' : 19, 'GTCCTA' : 7, 'GTCCTC' : 22, 'GTCCTG' : 40, 'GTCCTT' : 12, 'GTCTTA' : 6, 'GTCTTG' : 13, 'GTGCTA' : 12, 'GTGCTC' : 10, 'GTGCTG' : 36, 'GTGCTT' : 14, 'GTGTTA' : 13, 'GTGTTG' : 15, 'GTTCTA' : 12, 'GTTCTC' : 8, 'GTTCTG' : 27, 'GTTCTT' : 18, 'GTTTTA' : 17, 'GTTTTG' : 18 } ,
        'VM' : { 'GTAATG' : 100, 'GTCATG' : 100, 'GTGATG' : 100, 'GTTATG' : 100 } ,
        'VN' : { 'GTAAAC' : 51, 'GTAAAT' : 49, 'GTCAAC' : 61, 'GTCAAT' : 39, 'GTGAAC' : 55, 'GTGAAT' : 45, 'GTTAAC' : 56, 'GTTAAT' : 44 } ,
        'VP' : { 'GTACCA' : 34, 'GTACCC' : 23, 'GTACCG' : 14, 'GTACCT' : 28, 'GTCCCA' : 25, 'GTCCCC' : 35, 'GTCCCG' : 11, 'GTCCCT' : 29, 'GTGCCA' : 25, 'GTGCCC' : 33, 'GTGCCG' : 17, 'GTGCCT' : 25, 'GTTCCA' : 36, 'GTTCCC' : 22, 'GTTCCG' : 13, 'GTTCCT' : 29 } ,
        'VQ' : { 'GTACAA' : 37, 'GTACAG' : 63, 'GTCCAA' : 22, 'GTCCAG' : 78, 'GTGCAA' : 31, 'GTGCAG' : 69, 'GTTCAA' : 35, 'GTTCAG' : 65 } ,
        'VR' : { 'GTAAGA' : 31, 'GTAAGG' : 21, 'GTACGA' : 14, 'GTACGC' : 6, 'GTACGG' : 16, 'GTACGT' : 12, 'GTCAGA' : 20, 'GTCAGG' : 20, 'GTCCGA' : 12, 'GTCCGC' : 20, 'GTCCGG' : 20, 'GTCCGT' : 9, 'GTGAGA' : 21, 'GTGAGG' : 21, 'GTGCGA' : 13, 'GTGCGC' : 12, 'GTGCGG' : 21, 'GTGCGT' : 12, 'GTTAGA' : 26, 'GTTAGG' : 16, 'GTTCGA' : 19, 'GTTCGC' : 8, 'GTTCGG' : 17, 'GTTCGT' : 15 } ,
        'VS' : { 'GTAAGC' : 21, 'GTAAGT' : 13, 'GTATCA' : 21, 'GTATCC' : 10, 'GTATCG' : 10, 'GTATCT' : 26, 'GTCAGC' : 28, 'GTCAGT' : 15, 'GTCTCA' : 11, 'GTCTCC' : 25, 'GTCTCG' : 4, 'GTCTCT' : 16, 'GTGAGC' : 22, 'GTGAGT' : 14, 'GTGTCA' : 17, 'GTGTCC' : 14, 'GTGTCG' : 12, 'GTGTCT' : 21, 'GTTAGC' : 19, 'GTTAGT' : 10, 'GTTTCA' : 22, 'GTTTCC' : 12, 'GTTTCG' : 9, 'GTTTCT' : 28 } ,
        'VT' : { 'GTAACA' : 33, 'GTAACC' : 12, 'GTAACG' : 16, 'GTAACT' : 39, 'GTCACA' : 22, 'GTCACC' : 45, 'GTCACG' : 9, 'GTCACT' : 24, 'GTGACA' : 29, 'GTGACC' : 19, 'GTGACG' : 21, 'GTGACT' : 32, 'GTTACA' : 30, 'GTTACC' : 14, 'GTTACG' : 14, 'GTTACT' : 41 } ,
        'VV' : { 'GTAGTA' : 23, 'GTAGTC' : 8, 'GTAGTG' : 45, 'GTAGTT' : 23, 'GTCGTA' : 8, 'GTCGTC' : 33, 'GTCGTG' : 47, 'GTCGTT' : 12, 'GTGGTA' : 19, 'GTGGTC' : 12, 'GTGGTG' : 48, 'GTGGTT' : 21, 'GTTGTA' : 20, 'GTTGTC' : 12, 'GTTGTG' : 46, 'GTTGTT' : 22 } ,
        'VW' : { 'GTATGG' : 100, 'GTCTGG' : 100, 'GTGTGG' : 100, 'GTTTGG' : 100 } ,
        'VY' : { 'GTATAC' : 49, 'GTATAT' : 51, 'GTCTAC' : 62, 'GTCTAT' : 38, 'GTGTAC' : 56, 'GTGTAT' : 44, 'GTTTAC' : 54, 'GTTTAT' : 46 } ,
        'WA' : { 'TGGGCA' : 27, 'TGGGCC' : 20, 'TGGGCG' : 20, 'TGGGCT' : 32 } ,
        'WC' : { 'TGGTGC' : 59, 'TGGTGT' : 41 } ,
        'WU' : { 'TGGTGC' : 59, 'TGGTGT' : 41 } ,
        'WD' : { 'TGGGAC' : 78, 'TGGGAT' : 23 } ,
        'WE' : { 'TGGGAA' : 44, 'TGGGAG' : 56 } ,
        'WF' : { 'TGGTTC' : 76, 'TGGTTT' : 24 } ,
        'WG' : { 'TGGGGA' : 19, 'TGGGGC' : 44, 'TGGGGG' : 27, 'TGGGGT' : 9 } ,
        'WH' : { 'TGGCAC' : 79, 'TGGCAT' : 22 } ,
        'WI' : { 'TGGATA' : 23, 'TGGATC' : 34, 'TGGATT' : 43 } ,
        'WK' : { 'TGGAAA' : 46, 'TGGAAG' : 54 } ,
        'WL' : { 'TGGCTA' : 12, 'TGGCTC' : 10, 'TGGCTG' : 35, 'TGGCTT' : 16, 'TGGTTA' : 13, 'TGGTTG' : 14 } ,
        'WM' : { 'TGGATG' : 100 } ,
        'WN' : { 'TGGAAC' : 59, 'TGGAAT' : 41 } ,
        'WP' : { 'TGGCCA' : 32, 'TGGCCC' : 27, 'TGGCCG' : 17, 'TGGCCT' : 24 } ,
        'WQ' : { 'TGGCAA' : 31, 'TGGCAG' : 69 } ,
        'WR' : { 'TGGAGA' : 29, 'TGGAGG' : 25, 'TGGCGA' : 11, 'TGGCGC' : 8, 'TGGCGG' : 16, 'TGGCGT' : 10 } ,
        'WS' : { 'TGGAGC' : 28, 'TGGAGT' : 16, 'TGGTCA' : 17, 'TGGTCC' : 9, 'TGGTCG' : 11, 'TGGTCT' : 19 } ,
        'WT' : { 'TGGACA' : 30, 'TGGACC' : 16, 'TGGACG' : 19, 'TGGACT' : 35 } ,
        'WV' : { 'TGGGTA' : 18, 'TGGGTC' : 13, 'TGGGTG' : 51, 'TGGGTT' : 18 } ,
        'WW' : { 'TGGTGG' : 100 } ,
        'WY' : { 'TGGTAC' : 59, 'TGGTAT' : 41 } ,
        'YA' : { 'TACGCA' : 24, 'TACGCC' : 23, 'TACGCG' : 25, 'TACGCT' : 28, 'TATGCA' : 28, 'TATGCC' : 18, 'TATGCG' : 15, 'TATGCT' : 39 } ,
        'YC' : { 'TACTGC' : 63, 'TACTGT' : 37, 'TATTGC' : 74, 'TATTGT' : 26 } ,
        'YU' : { 'TACTGC' : 63, 'TACTGT' : 37, 'TATTGC' : 74, 'TATTGT' : 26 } ,
        'YD' : { 'TACGAC' : 83, 'TACGAT' : 17, 'TATGAC' : 75, 'TATGAT' : 25 } ,
        'YE' : { 'TACGAA' : 36, 'TACGAG' : 64, 'TATGAA' : 49, 'TATGAG' : 51 } ,
        'YF' : { 'TACTTC' : 81, 'TACTTT' : 19, 'TATTTC' : 73, 'TATTTT' : 28 } ,
        'YG' : { 'TACGGA' : 15, 'TACGGC' : 47, 'TACGGG' : 32, 'TACGGT' : 6, 'TATGGA' : 19, 'TATGGC' : 44, 'TATGGG' : 29, 'TATGGT' : 9 } ,
        'YH' : { 'TACCAC' : 82, 'TACCAT' : 19, 'TATCAC' : 75, 'TATCAT' : 25 } ,
        'YI' : { 'TACATA' : 20, 'TACATC' : 38, 'TACATT' : 42, 'TATATA' : 18, 'TATATC' : 19, 'TATATT' : 63 } ,
        'YK' : { 'TACAAA' : 46, 'TACAAG' : 54, 'TATAAA' : 53, 'TATAAG' : 47 } ,
        'YL' : { 'TACCTA' : 11, 'TACCTC' : 9, 'TACCTG' : 37, 'TACCTT' : 14, 'TACTTA' : 13, 'TACTTG' : 16, 'TATCTA' : 9, 'TATCTC' : 8, 'TATCTG' : 31, 'TATCTT' : 15, 'TATTTA' : 18, 'TATTTG' : 19 } ,
        'YM' : { 'TACATG' : 100, 'TATATG' : 100 } ,
        'YN' : { 'TACAAC' : 67, 'TACAAT' : 33, 'TATAAC' : 71, 'TATAAT' : 29 } ,
        'YP' : { 'TACCCA' : 36, 'TACCCC' : 24, 'TACCCG' : 18, 'TACCCT' : 22, 'TATCCA' : 56, 'TATCCC' : 12, 'TATCCG' : 9, 'TATCCT' : 23 } ,
        'YQ' : { 'TACCAA' : 28, 'TACCAG' : 72, 'TATCAA' : 29, 'TATCAG' : 71 } ,
        'YR' : { 'TACAGA' : 29, 'TACAGG' : 21, 'TACCGA' : 12, 'TACCGC' : 11, 'TACCGG' : 17, 'TACCGT' : 10, 'TATAGA' : 33, 'TATAGG' : 20, 'TATCGA' : 14, 'TATCGC' : 7, 'TATCGG' : 14, 'TATCGT' : 10 } ,
        'YS' : { 'TACAGC' : 32, 'TACAGT' : 14, 'TACTCA' : 17, 'TACTCC' : 10, 'TACTCG' : 11, 'TACTCT' : 17, 'TATAGC' : 29, 'TATAGT' : 7, 'TATTCA' : 22, 'TATTCC' : 11, 'TATTCG' : 6, 'TATTCT' : 26 } ,
        'YT' : { 'TACACA' : 28, 'TACACC' : 18, 'TACACG' : 20, 'TACACT' : 34, 'TATACA' : 28, 'TATACC' : 15, 'TATACG' : 12, 'TATACT' : 44 } ,
        'YV' : { 'TACGTA' : 15, 'TACGTC' : 12, 'TACGTG' : 59, 'TACGTT' : 14, 'TATGTA' : 17, 'TATGTC' : 12, 'TATGTG' : 55, 'TATGTT' : 16 } ,
        'YW' : { 'TACTGG' : 100, 'TATTGG' : 100 } ,
        'YY' : { 'TACTAC' : 64, 'TACTAT' : 36, 'TATTAC' : 75, 'TATTAT' : 25 }}

    #Bicodon Usage corrected for highly expressed genes and tRNA abundance in B-cells.
    #Filtered to include only the bicodons with a statistically significant difference in their frequencies
    #(highly expressed genes vs low-expressed and transcriptome).
    elif ex_sys == '3':
        CC_dict = {
        'AI' : { 'GCTATT' : 50, 'GCTATC' : 28, 'GCTATA' : 22} ,
        'AL' : { 'GCATTA' : 13, 'GCATTG' : 17, 'GCACTT' : 15, 'GCACTC' : 11, 'GCACTA' : 12, 'GCACTG' : 31, 'GCGTTA' : 13, 'GCGTTG' : 14, 'GCGCTT' : 12, 'GCGCTC' : 9, 'GCGCTA' : 11, 'GCGCTG' : 42} ,
        'AR' : { 'GCCCGT' : 11, 'GCCCGC' : 14, 'GCCCGA' : 9, 'GCCCGG' : 22, 'GCCAGA' : 21, 'GCCAGG' : 22, 'GCTCGT' : 12, 'GCTCGC' : 11, 'GCTCGA' : 19, 'GCTCGG' : 17, 'GCTAGA' : 24, 'GCTAGG' : 16} ,
        'AT' : { 'GCGACT' : 41, 'GCGACC' : 16, 'GCGACA' : 25, 'GCGACG' : 18} ,
        'CA' : { 'TGCGCT' : 33, 'TGCGCC' : 23, 'TGCGCA' : 22, 'TGCGCG' : 23} ,
        'CG' : { 'TGCGGG' : 32, 'TGCGGC' : 47, 'TGCGGT' : 10, 'TGCGGA' : 11} ,
        'CL' : { 'TGCTTA' : 12, 'TGCTTG' : 15, 'TGCCTT' : 14, 'TGCCTC' : 16, 'TGCCTA' : 8, 'TGCCTG' : 35, 'TGTTTA' : 16, 'TGTTTG' : 19, 'TGTCTT' : 15, 'TGTCTC' : 10, 'TGTCTA' : 6, 'TGTCTG' : 34} ,
        'CR' : { 'TGTCGT' : 17, 'TGTCGC' : 7, 'TGTCGA' : 12, 'TGTCGG' : 14, 'TGTAGA' : 32, 'TGTAGG' : 19} ,
        'UA' : { 'TGCGCT' : 33, 'TGCGCC' : 23, 'TGCGCA' : 22, 'TGCGCG' : 23} ,
        'UG' : { 'TGCGGG' : 32, 'TGCGGC' : 47, 'TGCGGT' : 10, 'TGCGGA' : 11} ,
        'UL' : { 'TGCTTA' : 12, 'TGCTTG' : 15, 'TGCCTT' : 14, 'TGCCTC' : 16, 'TGCCTA' : 8, 'TGCCTG' : 35, 'TGTTTA' : 16, 'TGTTTG' : 19, 'TGTCTT' : 15, 'TGTCTC' : 10, 'TGTCTA' : 6, 'TGTCTG' : 34} ,
        'UR' : { 'TGTCGT' : 17, 'TGTCGC' : 7, 'TGTCGA' : 12, 'TGTCGG' : 14, 'TGTAGA' : 32, 'TGTAGG' : 19} ,
        'DA' : { 'GACGCT' : 27, 'GACGCC' : 28, 'GACGCA' : 21, 'GACGCG' : 24} ,
        'DP' : { 'GATCCT' : 30, 'GATCCC' : 21, 'GATCCA' : 38, 'GATCCG' : 11} ,
        'DR' : { 'GATCGT' : 11, 'GATCGC' : 7, 'GATCGA' : 14, 'GATCGG' : 20, 'GATAGA' : 32, 'GATAGG' : 17} ,
        'ER' : { 'GAACGT' : 14, 'GAACGC' : 10, 'GAACGA' : 13, 'GAACGG' : 13, 'GAAAGA' : 29, 'GAAAGG' : 21} ,
        'ES' : { 'GAGTCT' : 16, 'GAGTCC' : 17, 'GAGTCA' : 14, 'GAGTCG' : 12, 'GAGAGT' : 16, 'GAGAGC' : 24} ,
        'FV' : { 'TTCGTT' : 17, 'TTCGTC' : 13, 'TTCGTA' : 17, 'TTCGTG' : 53} ,
        'GA' : { 'GGTGCT' : 37, 'GGTGCC' : 35, 'GGTGCA' : 17, 'GGTGCG' : 11} ,
        'GG' : { 'GGCGGG' : 26, 'GGCGGC' : 41, 'GGCGGT' : 17, 'GGCGGA' : 16, 'GGGGGG' : 20, 'GGGGGC' : 46, 'GGGGGT' : 16, 'GGGGGA' : 18} ,
        'GI' : { 'GGTATT' : 41, 'GGTATC' : 44, 'GGTATA' : 14} ,
        'GK' : { 'GGTAAA' : 49, 'GGTAAG' : 51} ,
        'GL' : { 'GGGTTA' : 12, 'GGGTTG' : 13, 'GGGCTT' : 13, 'GGGCTC' : 18, 'GGGCTA' : 12, 'GGGCTG' : 33, 'GGTTTA' : 12, 'GGTTTG' : 16, 'GGTCTT' : 14, 'GGTCTC' : 12, 'GGTCTA' : 11, 'GGTCTG' : 35} ,
        'GP' : { 'GGGCCT' : 25, 'GGGCCC' : 32, 'GGGCCA' : 26, 'GGGCCG' : 17} ,
        'GR' : { 'GGACGT' : 17, 'GGACGC' : 9, 'GGACGA' : 11, 'GGACGG' : 13, 'GGAAGA' : 36, 'GGAAGG' : 15, 'GGCCGT' : 12, 'GGCCGC' : 20, 'GGCCGA' : 15, 'GGCCGG' : 18, 'GGCAGA' : 18, 'GGCAGG' : 17, 'GGGCGT' : 12, 'GGGCGC' : 16, 'GGGCGA' : 13, 'GGGCGG' : 19, 'GGGAGA' : 20, 'GGGAGG' : 21, 'GGTCGT' : 22, 'GGTCGC' : 16, 'GGTCGA' : 11, 'GGTCGG' : 13, 'GGTAGA' : 24, 'GGTAGG' : 13} ,
        'GS' : { 'GGCTCT' : 16, 'GGCTCC' : 21, 'GGCTCA' : 12, 'GGCTCG' : 11, 'GGCAGT' : 15, 'GGCAGC' : 25, 'GGTTCT' : 33, 'GGTTCC' : 18, 'GGTTCA' : 18, 'GGTTCG' : 5, 'GGTAGT' : 8, 'GGTAGC' : 18} ,
        'GT' : { 'GGGACT' : 28, 'GGGACC' : 29, 'GGGACA' : 26, 'GGGACG' : 17, 'GGTACT' : 42, 'GGTACC' : 27, 'GGTACA' : 19, 'GGTACG' : 12} ,
        'GV' : { 'GGAGTT' : 22, 'GGAGTC' : 20, 'GGAGTA' : 14, 'GGAGTG' : 43} ,
        'HA' : { 'CATGCT' : 36, 'CATGCC' : 29, 'CATGCA' : 27, 'CATGCG' : 8} ,
        'HR' : { 'CATCGT' : 18, 'CATCGC' : 8, 'CATCGA' : 13, 'CATCGG' : 15, 'CATAGA' : 28, 'CATAGG' : 18} ,
        'HT' : { 'CATACT' : 47, 'CATACC' : 11, 'CATACA' : 24, 'CATACG' : 18} ,
        'IA' : { 'ATCGCT' : 32, 'ATCGCC' : 31, 'ATCGCA' : 18, 'ATCGCG' : 19} ,
        'II' : { 'ATAATT' : 63, 'ATAATC' : 19, 'ATAATA' : 17} ,
        'IL' : { 'ATATTA' : 16, 'ATATTG' : 20, 'ATACTT' : 12, 'ATACTC' : 8, 'ATACTA' : 13, 'ATACTG' : 32} ,
        'IS' : { 'ATATCT' : 26, 'ATATCC' : 11, 'ATATCA' : 17, 'ATATCG' : 9, 'ATAAGT' : 12, 'ATAAGC' : 25} ,
        'KL' : { 'AAGTTA' : 11, 'AAGTTG' : 14, 'AAGCTT' : 14, 'AAGCTC' : 16, 'AAGCTA' : 11, 'AAGCTG' : 33} ,
        'KR' : { 'AAACGT' : 16, 'AAACGC' : 9, 'AAACGA' : 12, 'AAACGG' : 13, 'AAAAGA' : 29, 'AAAAGG' : 20, 'AAGCGT' : 12, 'AAGCGC' : 13, 'AAGCGA' : 11, 'AAGCGG' : 16, 'AAGAGA' : 25, 'AAGAGG' : 23} ,
        'LA' : { 'CTCGCT' : 40, 'CTCGCC' : 14, 'CTCGCA' : 24, 'CTCGCG' : 22, 'CTTGCT' : 41, 'CTTGCC' : 19, 'CTTGCA' : 25, 'CTTGCG' : 14, 'TTAGCT' : 36, 'TTAGCC' : 18, 'TTAGCA' : 25, 'TTAGCG' : 20, 'TTGGCT' : 35, 'TTGGCC' : 23, 'TTGGCA' : 27, 'TTGGCG' : 15} ,
        'LE' : { 'CTAGAA' : 52, 'CTAGAG' : 48, 'CTGGAA' : 34, 'CTGGAG' : 66} ,
        'LG' : { 'TTAGGG' : 31, 'TTAGGC' : 46, 'TTAGGT' : 11, 'TTAGGA' : 12, 'TTGGGG' : 27, 'TTGGGC' : 44, 'TTGGGT' : 11, 'TTGGGA' : 18} ,
        'LI' : { 'CTTATT' : 57, 'CTTATC' : 26, 'CTTATA' : 17, 'TTAATT' : 54, 'TTAATC' : 24, 'TTAATA' : 22, 'TTGATT' : 43, 'TTGATC' : 36, 'TTGATA' : 22} ,
        'LL' : { 'CTCTTA' : 9, 'CTCTTG' : 13, 'CTCCTT' : 15, 'CTCCTC' : 15, 'CTCCTA' : 7, 'CTCCTG' : 41, 'CTTTTA' : 15, 'CTTTTG' : 18, 'CTTCTT' : 17, 'CTTCTC' : 14, 'CTTCTA' : 12, 'CTTCTG' : 25} ,
        'LP' : { 'CTTCCT' : 27, 'CTTCCC' : 25, 'CTTCCA' : 36, 'CTTCCG' : 11, 'TTACCT' : 21, 'TTACCC' : 19, 'TTACCA' : 42, 'TTACCG' : 18, 'TTGCCT' : 27, 'TTGCCC' : 22, 'TTGCCA' : 35, 'TTGCCG' : 16} ,
        'LQ' : { 'TTACAA' : 30, 'TTACAG' : 70} ,
        'LR' : { 'CTCCGT' : 17, 'CTCCGC' : 18, 'CTCCGA' : 11, 'CTCCGG' : 12, 'CTCAGA' : 22, 'CTCAGG' : 20, 'CTTCGT' : 18, 'CTTCGC' : 9, 'CTTCGA' : 17, 'CTTCGG' : 14, 'CTTAGA' : 26, 'CTTAGG' : 16, 'TTACGT' : 12, 'TTACGC' : 6, 'TTACGA' : 14, 'TTACGG' : 13, 'TTAAGA' : 31, 'TTAAGG' : 23} ,
        'LS' : { 'CTATCT' : 20, 'CTATCC' : 11, 'CTATCA' : 20, 'CTATCG' : 10, 'CTAAGT' : 13, 'CTAAGC' : 26, 'TTATCT' : 19, 'TTATCC' : 13, 'TTATCA' : 17, 'TTATCG' : 11, 'TTAAGT' : 14, 'TTAAGC' : 25} ,
        'LV' : { 'TTAGTT' : 20, 'TTAGTC' : 11, 'TTAGTA' : 18, 'TTAGTG' : 51, 'TTGGTT' : 18, 'TTGGTC' : 17, 'TTGGTA' : 18, 'TTGGTG' : 48} ,
        'NG' : { 'AACGGG' : 31, 'AACGGC' : 43, 'AACGGT' : 9, 'AACGGA' : 18} ,
        'NR' : { 'AATCGT' : 10, 'AATCGC' : 11, 'AATCGA' : 13, 'AATCGG' : 13, 'AATAGA' : 34, 'AATAGG' : 19} ,
        'PA' : { 'CCCGCT' : 32, 'CCCGCC' : 25, 'CCCGCA' : 22, 'CCCGCG' : 20, 'CCGGCT' : 33, 'CCGGCC' : 25, 'CCGGCA' : 22, 'CCGGCG' : 19} ,
        'PG' : { 'CCGGGG' : 29, 'CCGGGC' : 53, 'CCGGGT' : 6, 'CCGGGA' : 12} ,
        'PH' : { 'CCGCAT' : 22, 'CCGCAC' : 78} ,
        'PI' : { 'CCTATT' : 62, 'CCTATC' : 24, 'CCTATA' : 13} ,
        'PL' : { 'CCATTA' : 14, 'CCATTG' : 16, 'CCACTT' : 15, 'CCACTC' : 13, 'CCACTA' : 12, 'CCACTG' : 29, 'CCCTTA' : 10, 'CCCTTG' : 16, 'CCCCTT' : 15, 'CCCCTC' : 13, 'CCCCTA' : 7, 'CCCCTG' : 40} ,
        'PR' : { 'CCACGT' : 13, 'CCACGC' : 9, 'CCACGA' : 12, 'CCACGG' : 15, 'CCAAGA' : 29, 'CCAAGG' : 22, 'CCTCGT' : 18, 'CCTCGC' : 8, 'CCTCGA' : 14, 'CCTCGG' : 18, 'CCTAGA' : 26, 'CCTAGG' : 16} ,
        'PS' : { 'CCGTCT' : 21, 'CCGTCC' : 12, 'CCGTCA' : 18, 'CCGTCG' : 10, 'CCGAGT' : 7, 'CCGAGC' : 33} ,
        'PT' : { 'CCGACT' : 45, 'CCGACC' : 13, 'CCGACA' : 25, 'CCGACG' : 17} ,
        'QL' : { 'CAATTA' : 13, 'CAATTG' : 15, 'CAACTT' : 18, 'CAACTC' : 13, 'CAACTA' : 13, 'CAACTG' : 28} ,
        'RA' : { 'AGGGCT' : 26, 'AGGGCC' : 26, 'AGGGCA' : 31, 'AGGGCG' : 17, 'CGAGCT' : 40, 'CGAGCC' : 23, 'CGAGCA' : 22, 'CGAGCG' : 15, 'CGGGCT' : 33, 'CGGGCC' : 30, 'CGGGCA' : 20, 'CGGGCG' : 17, 'CGTGCT' : 49, 'CGTGCC' : 24, 'CGTGCA' : 15, 'CGTGCG' : 13} ,
        'RG' : { 'AGAGGG' : 23, 'AGAGGC' : 31, 'AGAGGT' : 22, 'AGAGGA' : 24, 'AGGGGG' : 22, 'AGGGGC' : 45, 'AGGGGT' : 12, 'AGGGGA' : 21, 'CGTGGG' : 24, 'CGTGGC' : 35, 'CGTGGT' : 20, 'CGTGGA' : 22} ,
        'RI' : { 'CGTATT' : 44, 'CGTATC' : 34, 'CGTATA' : 22} ,
        'RL' : { 'CGTTTA' : 13, 'CGTTTG' : 15, 'CGTCTT' : 14, 'CGTCTC' : 12, 'CGTCTA' : 10, 'CGTCTG' : 36} ,
        'RP' : { 'CGCCCT' : 16, 'CGCCCC' : 27, 'CGCCCA' : 41, 'CGCCCG' : 16, 'CGTCCT' : 21, 'CGTCCC' : 34, 'CGTCCA' : 30, 'CGTCCG' : 15} ,
        'RR' : { 'CGACGT' : 15, 'CGACGC' : 9, 'CGACGA' : 11, 'CGACGG' : 14, 'CGAAGA' : 30, 'CGAAGG' : 21, 'CGCCGT' : 14, 'CGCCGC' : 21, 'CGCCGA' : 8, 'CGCCGG' : 19, 'CGCAGA' : 18, 'CGCAGG' : 20} ,
        'RS' : { 'AGATCT' : 21, 'AGATCC' : 14, 'AGATCA' : 16, 'AGATCG' : 11, 'AGAAGT' : 18, 'AGAAGC' : 21, 'AGGTCT' : 22, 'AGGTCC' : 14, 'AGGTCA' : 15, 'AGGTCG' : 10, 'AGGAGT' : 14, 'AGGAGC' : 25} ,
        'RT' : { 'AGGACT' : 31, 'AGGACC' : 27, 'AGGACA' : 28, 'AGGACG' : 14, 'CGAACT' : 32, 'CGAACC' : 28, 'CGAACA' : 25, 'CGAACG' : 15, 'CGTACT' : 33, 'CGTACC' : 31, 'CGTACA' : 22, 'CGTACG' : 15} ,
        'RV' : { 'AGGGTT' : 16, 'AGGGTC' : 23, 'AGGGTA' : 18, 'AGGGTG' : 43, 'CGAGTT' : 18, 'CGAGTC' : 17, 'CGAGTA' : 18, 'CGAGTG' : 47} ,
        'SA' : { 'AGTGCT' : 36, 'AGTGCC' : 28, 'AGTGCA' : 25, 'AGTGCG' : 10} ,
        'SC' : { 'TCGTGT' : 26, 'TCGTGC' : 74} ,
        'SU' : { 'TCGTGT' : 26, 'TCGTGC' : 74} ,
        'SE' : { 'TCCGAA' : 35, 'TCCGAG' : 65} ,
        'SI' : { 'AGTATT' : 66, 'AGTATC' : 23, 'AGTATA' : 11, 'TCAATT' : 55, 'TCAATC' : 24, 'TCAATA' : 21, 'TCGATT' : 65, 'TCGATC' : 15, 'TCGATA' : 20, 'TCTATT' : 56, 'TCTATC' : 23, 'TCTATA' : 22} ,
        'SK' : { 'TCGAAA' : 43, 'TCGAAG' : 57} ,
        'SL' : { 'TCATTA' : 14, 'TCATTG' : 17, 'TCACTT' : 16, 'TCACTC' : 11, 'TCACTA' : 12, 'TCACTG' : 30, 'TCCTTA' : 10, 'TCCTTG' : 16, 'TCCCTT' : 11, 'TCCCTC' : 18, 'TCCCTA' : 7, 'TCCCTG' : 38, 'TCTTTA' : 12, 'TCTTTG' : 17, 'TCTCTT' : 15, 'TCTCTC' : 17, 'TCTCTA' : 11, 'TCTCTG' : 28} ,
        'SP' : { 'TCCCCT' : 24, 'TCCCCC' : 20, 'TCCCCA' : 49, 'TCCCCG' : 8} ,
        'SR' : { 'TCCCGT' : 12, 'TCCCGC' : 7, 'TCCCGA' : 10, 'TCCCGG' : 18, 'TCCAGA' : 28, 'TCCAGG' : 25, 'TCGCGT' : 10, 'TCGCGC' : 3, 'TCGCGA' : 11, 'TCGCGG' : 18, 'TCGAGA' : 34, 'TCGAGG' : 23, 'TCTCGT' : 13, 'TCTCGC' : 11, 'TCTCGA' : 19, 'TCTCGG' : 15, 'TCTAGA' : 28, 'TCTAGG' : 15} ,
        'SS' : { 'AGTTCT' : 28, 'AGTTCC' : 12, 'AGTTCA' : 16, 'AGTTCG' : 6, 'AGTAGT' : 8, 'AGTAGC' : 29, 'TCCTCT' : 23, 'TCCTCC' : 17, 'TCCTCA' : 16, 'TCCTCG' : 7, 'TCCAGT' : 10, 'TCCAGC' : 28, 'TCGTCT' : 21, 'TCGTCC' : 11, 'TCGTCA' : 18, 'TCGTCG' : 16, 'TCGAGT' : 6, 'TCGAGC' : 29, 'TCTTCT' : 22, 'TCTTCC' : 20, 'TCTTCA' : 21, 'TCTTCG' : 9, 'TCTAGT' : 9, 'TCTAGC' : 18} ,
        'ST' : { 'AGTACT' : 42, 'AGTACC' : 13, 'AGTACA' : 27, 'AGTACG' : 17, 'TCAACT' : 43, 'TCAACC' : 10, 'TCAACA' : 32, 'TCAACG' : 16} ,
        'SV' : { 'TCCGTT' : 12, 'TCCGTC' : 6, 'TCCGTA' : 15, 'TCCGTG' : 67} ,
        'TA' : { 'ACCGCT' : 30, 'ACCGCC' : 25, 'ACCGCA' : 20, 'ACCGCG' : 25} ,
        'TI' : { 'ACGATT' : 50, 'ACGATC' : 25, 'ACGATA' : 25, 'ACTATT' : 46, 'ACTATC' : 29, 'ACTATA' : 25} ,
        'TL' : { 'ACATTA' : 13, 'ACATTG' : 17, 'ACACTT' : 11, 'ACACTC' : 13, 'ACACTA' : 13, 'ACACTG' : 33, 'ACTTTA' : 13, 'ACTTTG' : 18, 'ACTCTT' : 16, 'ACTCTC' : 14, 'ACTCTA' : 13, 'ACTCTG' : 26} ,
        'TR' : { 'ACACGT' : 10, 'ACACGC' : 9, 'ACACGA' : 12, 'ACACGG' : 11, 'ACAAGA' : 34, 'ACAAGG' : 24, 'ACCCGT' : 12, 'ACCCGC' : 13, 'ACCCGA' : 9, 'ACCCGG' : 17, 'ACCAGA' : 25, 'ACCAGG' : 24} ,
        'VD' : { 'GTCGAT' : 17, 'GTCGAC' : 83} ,
        'VF' : { 'GTATTT' : 45, 'GTATTC' : 55} ,
        'VI' : { 'GTAATT' : 51, 'GTAATC' : 26, 'GTAATA' : 23, 'GTCATT' : 39, 'GTCATC' : 54, 'GTCATA' : 7} ,
        'VL' : { 'GTGTTA' : 11, 'GTGTTG' : 14, 'GTGCTT' : 14, 'GTGCTC' : 17, 'GTGCTA' : 12, 'GTGCTG' : 33, 'GTTTTA' : 12, 'GTTTTG' : 18, 'GTTCTT' : 18, 'GTTCTC' : 14, 'GTTCTA' : 11, 'GTTCTG' : 26} ,
        'VP' : { 'GTACCT' : 20, 'GTACCC' : 24, 'GTACCA' : 38, 'GTACCG' : 18} ,
        'VR' : { 'GTACGT' : 12, 'GTACGC' : 9, 'GTACGA' : 14, 'GTACGG' : 15, 'GTAAGA' : 31, 'GTAAGG' : 18, 'GTTCGT' : 21, 'GTTCGC' : 10, 'GTTCGA' : 17, 'GTTCGG' : 12, 'GTTAGA' : 24, 'GTTAGG' : 15} ,
        'VT' : { 'GTAACT' : 38, 'GTAACC' : 16, 'GTAACA' : 30, 'GTAACG' : 17} ,
        'VV' : { 'GTCGTT' : 14, 'GTCGTC' : 19, 'GTCGTA' : 13, 'GTCGTG' : 53} ,
        'YA' : { 'TATGCT' : 36, 'TATGCC' : 25, 'TATGCA' : 29, 'TATGCG' : 10}}
        #this dictionary serves as reference to check if CC for B-cells should be used.
        CC_evaluation_dict = {
        'AL' : ['GCG', 'GCA'],
        'AR' : ['GCT', 'GCC'],
        'AI' : ['GCT'],
        'AT' : ['GCG'],
        'RG' : ['AGA', 'CGT', 'AGG'],
        'RL' : ['CGT'],
        'RS' : ['AGA', 'AGG'],
        'RP' : ['CGT', 'CGC'],
        'RR' : ['CGA', 'CGC'],
        'RI' : ['CGT'],
        'RT' : ['CGA', 'CGT', 'AGG'],
        'RV' : ['CGA', 'AGG'],
        'RA' : ['CGA', 'AGG', 'CGT', 'CGG'],
        'NG' : ['AAC'],
        'NR' : ['AAT'],
        'DP' : ['GAT'],
        'DR' : ['GAT'],
        'DA' : ['GAC'],
        'CG' : ['TGC'],
        'CL' : ['TGT', 'TGC'],
        'CR' : ['TGT'],
        'CA' : ['TGC'],
        'UG' : ['TGC'],
        'UL' : ['TGT', 'TGC'],
        'UR' : ['TGT'],
        'UA' : ['TGC'],
        'ES' : ['GAG'],
        'ER' : ['GAA'],
        'QL' : ['CAA'],
        'GG' : ['GGG', 'GGC'],
        'GL' : ['GGG', 'GGT'],
        'GS' : ['GGT', 'GGC'],
        'GP' : ['GGG'],
        'GR' : ['GGA', 'GGG', 'GGT', 'GGC'],
        'GI' : ['GGT'],
        'GT' : ['GGG', 'GGT'],
        'GK' : ['GGT'],
        'GV' : ['GGA'],
        'GA' : ['GGT'],
        'HR' : ['CAT'],
        'HT' : ['CAT'],
        'HA' : ['CAT'],
        'IL' : ['ATA'],
        'IS' : ['ATA'],
        'II' : ['ATA'],
        'IA' : ['ATC'],
        'LG' : ['TTG', 'TTA'],
        'LL' : ['CTC', 'CTT'],
        'LS' : ['CTA', 'TTA'],
        'LP' : ['TTG', 'TTA', 'CTT'],
        'LQ' : ['TTA'],
        'LR' : ['CTC', 'TTA', 'CTT'],
        'LI' : ['TTG', 'TTA', 'CTT'],
        'LV' : ['TTG', 'TTA'],
        'LA' : ['TTG', 'CTC', 'TTA', 'CTT'],
        'LE' : ['CTG', 'CTA'],
        'KL' : ['AAG'],
        'KR' : ['AAA', 'AAG'],
        'FV' : ['TTC'],
        'PG' : ['CCG'],
        'PL' : ['CCC', 'CCA'],
        'PS' : ['CCG'],
        'PH' : ['CCG'],
        'PR' : ['CCA', 'CCT'],
        'PI' : ['CCT'],
        'PT' : ['CCG'],
        'PA' : ['CCC', 'CCG'],
        'SL' : ['TCA', 'TCC', 'TCT'],
        'SS' : ['AGT', 'TCC', 'TCT', 'TCG'],
        'SC' : ['TCG'],
        'SU' : ['TCG'],
        'SP' : ['TCC'],
        'SR' : ['TCC', 'TCT', 'TCG'],
        'SI' : ['AGT', 'TCA', 'TCT', 'TCG'],
        'ST' : ['AGT', 'TCA'],
        'SK' : ['TCG'],
        'SV' : ['TCC'],
        'SA' : ['AGT'],
        'SE' : ['TCC'],
        'TL' : ['ACA', 'ACT'],
        'TR' : ['ACC', 'ACA'],
        'TI' : ['ACG', 'ACT'],
        'TA' : ['ACC'],
        'YA' : ['TAT'],
        'VF' : ['GTA'],
        'VL' : ['GTT', 'GTG'],
        'VP' : ['GTA'],
        'VR' : ['GTT', 'GTA'],
        'VI' : ['GTC', 'GTA'],
        'VT' : ['GTA'],
        'VV' : ['GTC'],
        'VD' : ['GTC']}

    #Bicodon Usage corrected for highly expressed genes and tRNA abundance in HEK293-cells.
    #Filtered to include only the bicodons with a statistically significant difference in their frequencies
    #(highly expressed genes vs low-expressed and transcriptome).
    else:
        CC_dict = {
        'AA' : { 'GCCGCA' : 17, 'GCCGCC' : 40, 'GCCGCG' : 8, 'GCCGCT' : 36, 'GCTGCA' : 24, 'GCTGCC' : 50, 'GCTGCG' : 5, 'GCTGCT' : 21} ,
        'AD' : { 'GCAGAC' : 74, 'GCAGAT' : 26, 'GCTGAC' : 75, 'GCTGAT' : 25} ,
        'AE' : { 'GCAGAA' : 61, 'GCAGAG' : 39, 'GCTGAA' : 46, 'GCTGAG' : 54} ,
        'AF' : { 'GCATTC' : 41, 'GCATTT' : 59, 'GCGTTC' : 46, 'GCGTTT' : 54} ,
        'AG' : { 'GCAGGA' : 41, 'GCAGGC' : 14, 'GCAGGG' : 20, 'GCAGGT' : 25, 'GCTGGA' : 14, 'GCTGGC' : 22, 'GCTGGG' : 29, 'GCTGGT' : 35} ,
        'AI' : { 'GCAATA' : 9, 'GCAATC' : 44, 'GCAATT' : 47, 'GCCATA' : 13, 'GCCATC' : 49, 'GCCATT' : 38, 'GCGATA' : 66, 'GCGATC' : 17, 'GCGATT' : 17, 'GCTATA' : 1, 'GCTATC' : 37, 'GCTATT' : 61} ,
        'AK' : { 'GCAAAA' : 59, 'GCAAAG' : 41, 'GCCAAA' : 30, 'GCCAAG' : 70} ,
        'AL' : { 'GCACTA' : 26, 'GCACTC' : 8, 'GCACTG' : 44, 'GCACTT' : 12, 'GCATTA' : 3, 'GCATTG' : 6, 'GCCCTA' : 6, 'GCCCTC' : 17, 'GCCCTG' : 40, 'GCCCTT' : 16, 'GCCTTA' : 12, 'GCCTTG' : 10, 'GCGCTA' : 24, 'GCGCTC' : 16, 'GCGCTG' : 27, 'GCGCTT' : 1, 'GCGTTA' : 10, 'GCGTTG' : 22} ,
        'AM' : { 'GCAATG' : 13, 'GCCATG' : 87} ,
        'AN' : { 'GCAAAC' : 22, 'GCAAAT' : 78, 'GCCAAC' : 57, 'GCCAAT' : 43} ,
        'AP' : { 'GCACCA' : 39, 'GCACCC' : 26, 'GCACCG' : 16, 'GCACCT' : 19, 'GCTCCA' : 26, 'GCTCCC' : 32, 'GCTCCG' : 18, 'GCTCCT' : 24} ,
        'AQ' : { 'GCACAA' : 62, 'GCACAG' : 38, 'GCCCAA' : 10, 'GCCCAG' : 90} ,
        'AR' : { 'GCAAGA' : 36, 'GCAAGG' : 14, 'GCACGA' : 20, 'GCACGC' : 11, 'GCACGG' : 8, 'GCACGT' : 11, 'GCCAGA' : 17, 'GCCAGG' : 15, 'GCCCGA' : 15, 'GCCCGC' : 25, 'GCCCGG' : 23, 'GCCCGT' : 6, 'GCTAGA' : 4, 'GCTAGG' : 18, 'GCTCGA' : 19, 'GCTCGC' : 21, 'GCTCGG' : 19, 'GCTCGT' : 20} ,
        'AS' : { 'GCGAGC' : 25, 'GCGAGT' : 9, 'GCGTCA' : 20, 'GCGTCC' : 9, 'GCGTCG' : 15, 'GCGTCT' : 23} ,
        'AT' : { 'GCCACA' : 15, 'GCCACC' : 46, 'GCCACG' : 10, 'GCCACT' : 29, 'GCGACA' : 27, 'GCGACC' : 27, 'GCGACG' : 18, 'GCGACT' : 28, 'GCTACA' : 38, 'GCTACC' : 11, 'GCTACG' : 6, 'GCTACT' : 45} ,
        'AV' : { 'GCTGTA' : 6, 'GCTGTC' : 20, 'GCTGTG' : 56, 'GCTGTT' : 18} ,
        'CE' : { 'TGTGAA' : 43, 'TGTGAG' : 57} ,
        'CH' : { 'TGCCAC' : 25, 'TGCCAT' : 75} ,
        'CK' : { 'TGTAAA' : 44, 'TGTAAG' : 56} ,
        'CP' : { 'TGTCCA' : 39, 'TGTCCC' : 16, 'TGTCCG' : 29, 'TGTCCT' : 16} ,
        'CR' : { 'TGTAGA' : 15, 'TGTAGG' : 10, 'TGTCGA' : 3, 'TGTCGC' : 19, 'TGTCGG' : 25, 'TGTCGT' : 28} ,
        'CT' : { 'TGTACA' : 24, 'TGTACC' : 21, 'TGTACG' : 26, 'TGTACT' : 29} ,
        'CW' : { 'TGCTGG' : 100, 'TGTTGG' : 100} ,
        'UE' : { 'TGTGAA' : 43, 'TGTGAG' : 57} ,
        'UH' : { 'TGCCAC' : 25, 'TGCCAT' : 75} ,
        'UK' : { 'TGTAAA' : 44, 'TGTAAG' : 56} ,
        'UP' : { 'TGTCCA' : 39, 'TGTCCC' : 16, 'TGTCCG' : 29, 'TGTCCT' : 16} ,
        'UR' : { 'TGTAGA' : 15, 'TGTAGG' : 10, 'TGTCGA' : 3, 'TGTCGC' : 19, 'TGTCGG' : 25, 'TGTCGT' : 28} ,
        'UT' : { 'TGTACA' : 24, 'TGTACC' : 21, 'TGTACG' : 26, 'TGTACT' : 29} ,
        'UW' : { 'TGCTGG' : 100, 'TGTTGG' : 100} ,
        'DA' : { 'GACGCA' : 29, 'GACGCC' : 43, 'GACGCG' : 25, 'GACGCT' : 3, 'GATGCA' : 23, 'GATGCC' : 41, 'GATGCG' : 6, 'GATGCT' : 30} ,
        'DC' : { 'GATTGC' : 15, 'GATTGT' : 85} ,
        'DU' : { 'GATTGC' : 15, 'GATTGT' : 85} ,
        'DE' : { 'GATGAA' : 67, 'GATGAG' : 33} ,
        'DF' : { 'GACTTC' : 87, 'GACTTT' : 13} ,
        'DG' : { 'GATGGA' : 44, 'GATGGC' : 27, 'GATGGG' : 15, 'GATGGT' : 14} ,
        'DI' : { 'GACATA' : 29, 'GACATC' : 53, 'GACATT' : 18, 'GATATA' : 63, 'GATATC' : 18, 'GATATT' : 19} ,
        'DK' : { 'GATAAA' : 79, 'GATAAG' : 21} ,
        'DL' : { 'GATCTA' : 37, 'GATCTC' : 7, 'GATCTG' : 14, 'GATCTT' : 8, 'GATTTA' : 26, 'GATTTG' : 9} ,
        'DN' : { 'GACAAC' : 83, 'GACAAT' : 17, 'GATAAC' : 23, 'GATAAT' : 77} ,
        'DQ' : { 'GACCAA' : 31, 'GACCAG' : 69, 'GATCAA' : 26, 'GATCAG' : 74} ,
        'DR' : { 'GACAGA' : 10, 'GACAGG' : 23, 'GACCGA' : 22, 'GACCGC' : 13, 'GACCGG' : 26, 'GACCGT' : 6, 'GATAGA' : 34, 'GATAGG' : 6, 'GATCGA' : 8, 'GATCGC' : 40, 'GATCGG' : 8, 'GATCGT' : 4} ,
        'DT' : { 'GATACA' : 37, 'GATACC' : 49, 'GATACG' : 2, 'GATACT' : 12} ,
        'DV' : { 'GATGTA' : 17, 'GATGTC' : 11, 'GATGTG' : 62, 'GATGTT' : 9} ,
        'DY' : { 'GATTAC' : 22, 'GATTAT' : 78} ,
        'EA' : { 'GAAGCA' : 26, 'GAAGCC' : 17, 'GAAGCG' : 21, 'GAAGCT' : 36} ,
        'EC' : { 'GAATGC' : 50, 'GAATGT' : 50, 'GAGTGC' : 54, 'GAGTGT' : 46} ,
        'EU' : { 'GAATGC' : 50, 'GAATGT' : 50, 'GAGTGC' : 54, 'GAGTGT' : 46} ,
        'EE' : { 'GAAGAA' : 49, 'GAAGAG' : 51} ,
        'EF' : { 'GAGTTC' : 78, 'GAGTTT' : 22} ,
        'EG' : { 'GAAGGA' : 15, 'GAAGGC' : 21, 'GAAGGG' : 30, 'GAAGGT' : 34, 'GAGGGA' : 23, 'GAGGGC' : 34, 'GAGGGG' : 21, 'GAGGGT' : 21} ,
        'EH' : { 'GAACAC' : 23, 'GAACAT' : 77} ,
        'EI' : { 'GAAATA' : 26, 'GAAATC' : 31, 'GAAATT' : 43, 'GAGATA' : 19, 'GAGATC' : 48, 'GAGATT' : 33} ,
        'EK' : { 'GAAAAA' : 43, 'GAAAAG' : 57, 'GAGAAA' : 39, 'GAGAAG' : 61} ,
        'EL' : { 'GAGCTA' : 3, 'GAGCTC' : 19, 'GAGCTG' : 42, 'GAGCTT' : 13, 'GAGTTA' : 9, 'GAGTTG' : 14} ,
        'EM' : { 'GAAATG' : 31, 'GAGATG' : 69} ,
        'EN' : { 'GAAAAC' : 49, 'GAAAAT' : 51, 'GAGAAC' : 56, 'GAGAAT' : 44} ,
        'EP' : { 'GAACCA' : 29, 'GAACCC' : 27, 'GAACCG' : 22, 'GAACCT' : 23} ,
        'EQ' : { 'GAACAA' : 48, 'GAACAG' : 52, 'GAGCAA' : 27, 'GAGCAG' : 73} ,
        'ER' : { 'GAAAGA' : 23, 'GAAAGG' : 18, 'GAACGA' : 5, 'GAACGC' : 18, 'GAACGG' : 19, 'GAACGT' : 18, 'GAGAGA' : 15, 'GAGAGG' : 19, 'GAGCGA' : 15, 'GAGCGC' : 22, 'GAGCGG' : 23, 'GAGCGT' : 6} ,
        'ES' : { 'GAAAGC' : 21, 'GAAAGT' : 18, 'GAATCA' : 15, 'GAATCC' : 7, 'GAATCG' : 15, 'GAATCT' : 24, 'GAGAGC' : 25, 'GAGAGT' : 9, 'GAGTCA' : 14, 'GAGTCC' : 19, 'GAGTCG' : 16, 'GAGTCT' : 17} ,
        'ET' : { 'GAAACA' : 16, 'GAAACC' : 37, 'GAAACG' : 3, 'GAAACT' : 44, 'GAGACA' : 25, 'GAGACC' : 32, 'GAGACG' : 20, 'GAGACT' : 23} ,
        'EY' : { 'GAGTAC' : 60, 'GAGTAT' : 40} ,
        'FA' : { 'TTCGCA' : 31, 'TTCGCC' : 59, 'TTCGCG' : 5, 'TTCGCT' : 5, 'TTTGCA' : 45, 'TTTGCC' : 30, 'TTTGCG' : 4, 'TTTGCT' : 22} ,
        'FD' : { 'TTTGAC' : 23, 'TTTGAT' : 77} ,
        'FE' : { 'TTCGAA' : 45, 'TTCGAG' : 55} ,
        'FF' : { 'TTCTTC' : 85, 'TTCTTT' : 15, 'TTTTTC' : 61, 'TTTTTT' : 39} ,
        'FG' : { 'TTCGGA' : 41, 'TTCGGC' : 9, 'TTCGGG' : 8, 'TTCGGT' : 42, 'TTTGGA' : 14, 'TTTGGC' : 18, 'TTTGGG' : 16, 'TTTGGT' : 52} ,
        'FH' : { 'TTTCAC' : 8, 'TTTCAT' : 92} ,
        'FI' : { 'TTCATA' : 19, 'TTCATC' : 41, 'TTCATT' : 40} ,
        'FK' : { 'TTCAAA' : 18, 'TTCAAG' : 82, 'TTTAAA' : 75, 'TTTAAG' : 25} ,
        'FL' : { 'TTCCTA' : 12, 'TTCCTC' : 25, 'TTCCTG' : 48, 'TTCCTT' : 6, 'TTCTTA' : 2, 'TTCTTG' : 8, 'TTTCTA' : 27, 'TTTCTC' : 6, 'TTTCTG' : 15, 'TTTCTT' : 8, 'TTTTTA' : 36, 'TTTTTG' : 8} ,
        'FM' : { 'TTCATG' : 81, 'TTTATG' : 19} ,
        'FN' : { 'TTCAAC' : 62, 'TTCAAT' : 38, 'TTTAAC' : 48, 'TTTAAT' : 52} ,
        'FP' : { 'TTCCCA' : 11, 'TTCCCC' : 19, 'TTCCCG' : 29, 'TTCCCT' : 42, 'TTTCCA' : 30, 'TTTCCC' : 13, 'TTTCCG' : 40, 'TTTCCT' : 17} ,
        'FQ' : { 'TTCCAA' : 9, 'TTCCAG' : 91, 'TTTCAA' : 10, 'TTTCAG' : 90} ,
        'FR' : { 'TTCAGA' : 6, 'TTCAGG' : 20, 'TTCCGA' : 6, 'TTCCGC' : 25, 'TTCCGG' : 24, 'TTCCGT' : 19, 'TTTAGA' : 14, 'TTTAGG' : 61, 'TTTCGA' : 6, 'TTTCGC' : 9, 'TTTCGG' : 7, 'TTTCGT' : 3} ,
        'FS' : { 'TTCAGC' : 18, 'TTCAGT' : 19, 'TTCTCA' : 4, 'TTCTCC' : 21, 'TTCTCG' : 18, 'TTCTCT' : 20} ,
        'FT' : { 'TTCACA' : 36, 'TTCACC' : 45, 'TTCACG' : 5, 'TTCACT' : 14, 'TTTACA' : 67, 'TTTACC' : 11, 'TTTACG' : 3, 'TTTACT' : 19} ,
        'FV' : { 'TTCGTA' : 1, 'TTCGTC' : 19, 'TTCGTG' : 27, 'TTCGTT' : 53} ,
        'FY' : { 'TTCTAC' : 78, 'TTCTAT' : 22, 'TTTTAC' : 18, 'TTTTAT' : 82} ,
        'GA' : { 'GGAGCA' : 24, 'GGAGCC' : 50, 'GGAGCG' : 3, 'GGAGCT' : 22, 'GGGGCA' : 8, 'GGGGCC' : 38, 'GGGGCG' : 22, 'GGGGCT' : 32} ,
        'GD' : { 'GGGGAC' : 77, 'GGGGAT' : 23} ,
        'GE' : { 'GGAGAA' : 49, 'GGAGAG' : 51, 'GGCGAA' : 59, 'GGCGAG' : 41} ,
        'GF' : { 'GGATTC' : 20, 'GGATTT' : 80, 'GGCTTC' : 60, 'GGCTTT' : 40} ,
        'GG' : { 'GGAGGA' : 21, 'GGAGGC' : 28, 'GGAGGG' : 10, 'GGAGGT' : 41, 'GGCGGA' : 29, 'GGCGGC' : 39, 'GGCGGG' : 9, 'GGCGGT' : 23, 'GGGGGA' : 43, 'GGGGGC' : 23, 'GGGGGG' : 3, 'GGGGGT' : 30, 'GGTGGA' : 25, 'GGTGGC' : 14, 'GGTGGG' : 32, 'GGTGGT' : 29} ,
        'GH' : { 'GGCCAC' : 38, 'GGCCAT' : 62} ,
        'GI' : { 'GGAATA' : 22, 'GGAATC' : 48, 'GGAATT' : 30, 'GGCATA' : 18, 'GGCATC' : 51, 'GGCATT' : 32, 'GGGATA' : 3, 'GGGATC' : 82, 'GGGATT' : 15, 'GGTATA' : 23, 'GGTATC' : 36, 'GGTATT' : 42} ,
        'GK' : { 'GGAAAA' : 41, 'GGAAAG' : 59, 'GGCAAA' : 26, 'GGCAAG' : 74, 'GGGAAA' : 17, 'GGGAAG' : 83} ,
        'GL' : { 'GGCCTA' : 3, 'GGCCTC' : 21, 'GGCCTG' : 29, 'GGCCTT' : 16, 'GGCTTA' : 13, 'GGCTTG' : 18, 'GGGCTA' : 2, 'GGGCTC' : 20, 'GGGCTG' : 37, 'GGGCTT' : 23, 'GGGTTA' : 1, 'GGGTTG' : 18, 'GGTCTA' : 3, 'GGTCTC' : 7, 'GGTCTG' : 32, 'GGTCTT' : 5, 'GGTTTA' : 19, 'GGTTTG' : 34} ,
        'GM' : { 'GGCATG' : 73, 'GGGATG' : 27} ,
        'GN' : { 'GGAAAC' : 32, 'GGAAAT' : 68, 'GGCAAC' : 51, 'GGCAAT' : 49} ,
        'GP' : { 'GGGCCA' : 26, 'GGGCCC' : 31, 'GGGCCG' : 22, 'GGGCCT' : 22, 'GGTCCA' : 8, 'GGTCCC' : 52, 'GGTCCG' : 2, 'GGTCCT' : 38} ,
        'GQ' : { 'GGCCAA' : 32, 'GGCCAG' : 68} ,
        'GR' : { 'GGAAGA' : 24, 'GGAAGG' : 19, 'GGACGA' : 4, 'GGACGC' : 27, 'GGACGG' : 6, 'GGACGT' : 20, 'GGCAGA' : 7, 'GGCAGG' : 22, 'GGCCGA' : 9, 'GGCCGC' : 25, 'GGCCGG' : 24, 'GGCCGT' : 13, 'GGGAGA' : 5, 'GGGAGG' : 7, 'GGGCGA' : 23, 'GGGCGC' : 16, 'GGGCGG' : 27, 'GGGCGT' : 22, 'GGTAGA' : 23, 'GGTAGG' : 20, 'GGTCGA' : 7, 'GGTCGC' : 19, 'GGTCGG' : 21, 'GGTCGT' : 10} ,
        'GS' : { 'GGAAGC' : 29, 'GGAAGT' : 28, 'GGATCA' : 23, 'GGATCC' : 9, 'GGATCG' : 1, 'GGATCT' : 10, 'GGCAGC' : 22, 'GGCAGT' : 18, 'GGCTCA' : 6, 'GGCTCC' : 21, 'GGCTCG' : 14, 'GGCTCT' : 19, 'GGTAGC' : 20, 'GGTAGT' : 5, 'GGTTCA' : 9, 'GGTTCC' : 12, 'GGTTCG' : 1, 'GGTTCT' : 53} ,
        'GT' : { 'GGCACA' : 11, 'GGCACC' : 42, 'GGCACG' : 25, 'GGCACT' : 22, 'GGTACA' : 11, 'GGTACC' : 13, 'GGTACG' : 25, 'GGTACT' : 51} ,
        'GV' : { 'GGAGTA' : 6, 'GGAGTC' : 12, 'GGAGTG' : 28, 'GGAGTT' : 53, 'GGGGTA' : 22, 'GGGGTC' : 24, 'GGGGTG' : 26, 'GGGGTT' : 28, 'GGTGTA' : 6, 'GGTGTC' : 13, 'GGTGTG' : 33, 'GGTGTT' : 48} ,
        'GY' : { 'GGATAC' : 63, 'GGATAT' : 37, 'GGCTAC' : 53, 'GGCTAT' : 47} ,
        'HA' : { 'CATGCA' : 38, 'CATGCC' : 31, 'CATGCG' : 7, 'CATGCT' : 25} ,
        'HD' : { 'CATGAC' : 24, 'CATGAT' : 76} ,
        'HE' : { 'CACGAA' : 11, 'CACGAG' : 89, 'CATGAA' : 68, 'CATGAG' : 32} ,
        'HF' : { 'CATTTC' : 12, 'CATTTT' : 88} ,
        'HI' : { 'CACATA' : 19, 'CACATC' : 45, 'CACATT' : 36} ,
        'HK' : { 'CACAAA' : 22, 'CACAAG' : 78, 'CATAAA' : 10, 'CATAAG' : 90} ,
        'HL' : { 'CACCTA' : 12, 'CACCTC' : 20, 'CACCTG' : 39, 'CACCTT' : 12, 'CACTTA' : 2, 'CACTTG' : 16, 'CATCTA' : 3, 'CATCTC' : 20, 'CATCTG' : 25, 'CATCTT' : 21, 'CATTTA' : 7, 'CATTTG' : 24} ,
        'HN' : { 'CACAAC' : 70, 'CACAAT' : 30} ,
        'HQ' : { 'CACCAA' : 7, 'CACCAG' : 93, 'CATCAA' : 16, 'CATCAG' : 84} ,
        'HR' : { 'CATAGA' : 9, 'CATAGG' : 3, 'CATCGA' : 4, 'CATCGC' : 10, 'CATCGG' : 60, 'CATCGT' : 14} ,
        'HS' : { 'CACAGC' : 25, 'CACAGT' : 19, 'CACTCA' : 21, 'CACTCC' : 27, 'CACTCG' : 4, 'CACTCT' : 3, 'CATAGC' : 13, 'CATAGT' : 6, 'CATTCA' : 36, 'CATTCC' : 20, 'CATTCG' : 4, 'CATTCT' : 20} ,
        'HT' : { 'CACACA' : 19, 'CACACC' : 55, 'CACACG' : 11, 'CACACT' : 14, 'CATACA' : 3, 'CATACC' : 8, 'CATACG' : 41, 'CATACT' : 49} ,
        'IA' : { 'ATAGCA' : 23, 'ATAGCC' : 40, 'ATAGCG' : 4, 'ATAGCT' : 33, 'ATCGCA' : 24, 'ATCGCC' : 41, 'ATCGCG' : 7, 'ATCGCT' : 28, 'ATTGCA' : 15, 'ATTGCC' : 47, 'ATTGCG' : 3, 'ATTGCT' : 35} ,
        'IC' : { 'ATATGC' : 37, 'ATATGT' : 63} ,
        'IU' : { 'ATATGC' : 37, 'ATATGT' : 63} ,
        'ID' : { 'ATAGAC' : 84, 'ATAGAT' : 16, 'ATCGAC' : 60, 'ATCGAT' : 40, 'ATTGAC' : 52, 'ATTGAT' : 48} ,
        'IE' : { 'ATAGAA' : 55, 'ATAGAG' : 45, 'ATCGAA' : 41, 'ATCGAG' : 59, 'ATTGAA' : 52, 'ATTGAG' : 48} ,
        'IF' : { 'ATATTC' : 43, 'ATATTT' : 57, 'ATCTTC' : 60, 'ATCTTT' : 40} ,
        'IG' : { 'ATAGGA' : 43, 'ATAGGC' : 8, 'ATAGGG' : 12, 'ATAGGT' : 37, 'ATCGGA' : 23, 'ATCGGC' : 24, 'ATCGGG' : 25, 'ATCGGT' : 28, 'ATTGGA' : 16, 'ATTGGC' : 28, 'ATTGGG' : 26, 'ATTGGT' : 31} ,
        'IH' : { 'ATACAC' : 9, 'ATACAT' : 91, 'ATCCAC' : 77, 'ATCCAT' : 23, 'ATTCAC' : 30, 'ATTCAT' : 70} ,
        'II' : { 'ATAATA' : 26, 'ATAATC' : 55, 'ATAATT' : 19, 'ATCATA' : 19, 'ATCATC' : 57, 'ATCATT' : 24, 'ATTATA' : 2, 'ATTATC' : 38, 'ATTATT' : 60} ,
        'IK' : { 'ATAAAA' : 57, 'ATAAAG' : 43, 'ATCAAA' : 40, 'ATCAAG' : 60, 'ATTAAA' : 41, 'ATTAAG' : 59} ,
        'IL' : { 'ATACTA' : 8, 'ATACTC' : 5, 'ATACTG' : 50, 'ATACTT' : 16, 'ATATTA' : 17, 'ATATTG' : 4, 'ATCCTA' : 5, 'ATCCTC' : 21, 'ATCCTG' : 35, 'ATCCTT' : 10, 'ATCTTA' : 17, 'ATCTTG' : 12, 'ATTCTA' : 4, 'ATTCTC' : 17, 'ATTCTG' : 22, 'ATTCTT' : 21, 'ATTTTA' : 15, 'ATTTTG' : 20} ,
        'IM' : { 'ATCATG' : 85, 'ATTATG' : 15} ,
        'IN' : { 'ATAAAC' : 17, 'ATAAAT' : 83, 'ATCAAC' : 47, 'ATCAAT' : 53, 'ATTAAC' : 52, 'ATTAAT' : 48} ,
        'IP' : { 'ATACCA' : 33, 'ATACCC' : 28, 'ATACCG' : 17, 'ATACCT' : 22, 'ATCCCA' : 19, 'ATCCCC' : 46, 'ATCCCG' : 17, 'ATCCCT' : 18} ,
        'IQ' : { 'ATACAA' : 60, 'ATACAG' : 40, 'ATCCAA' : 37, 'ATCCAG' : 63, 'ATTCAA' : 36, 'ATTCAG' : 64} ,
        'IR' : { 'ATAAGA' : 21, 'ATAAGG' : 22, 'ATACGA' : 18, 'ATACGC' : 19, 'ATACGG' : 16, 'ATACGT' : 3, 'ATCAGA' : 12, 'ATCAGG' : 12, 'ATCCGA' : 17, 'ATCCGC' : 20, 'ATCCGG' : 25, 'ATCCGT' : 14, 'ATTAGA' : 18, 'ATTAGG' : 7, 'ATTCGA' : 25, 'ATTCGC' : 7, 'ATTCGG' : 23, 'ATTCGT' : 21} ,
        'IS' : { 'ATAAGC' : 23, 'ATAAGT' : 28, 'ATATCA' : 9, 'ATATCC' : 5, 'ATATCG' : 22, 'ATATCT' : 14, 'ATCAGC' : 22, 'ATCAGT' : 13, 'ATCTCA' : 12, 'ATCTCC' : 31, 'ATCTCG' : 2, 'ATCTCT' : 19, 'ATTAGC' : 16, 'ATTAGT' : 19, 'ATTTCA' : 6, 'ATTTCC' : 9, 'ATTTCG' : 20, 'ATTTCT' : 29} ,
        'IT' : { 'ATAACA' : 29, 'ATAACC' : 27, 'ATAACG' : 12, 'ATAACT' : 33, 'ATCACA' : 23, 'ATCACC' : 36, 'ATCACG' : 7, 'ATCACT' : 34} ,
        'IV' : { 'ATCGTA' : 3, 'ATCGTC' : 24, 'ATCGTG' : 64, 'ATCGTT' : 9} ,
        'IW' : { 'ATCTGG' : 100} ,
        'IY' : { 'ATATAC' : 42, 'ATATAT' : 58, 'ATCTAC' : 47, 'ATCTAT' : 53} ,
        'KA' : { 'AAAGCA' : 14, 'AAAGCC' : 38, 'AAAGCG' : 3, 'AAAGCT' : 46, 'AAGGCA' : 27, 'AAGGCC' : 30, 'AAGGCG' : 20, 'AAGGCT' : 24} ,
        'KC' : { 'AAATGC' : 54, 'AAATGT' : 46, 'AAGTGC' : 48, 'AAGTGT' : 52} ,
        'KU' : { 'AAATGC' : 54, 'AAATGT' : 46, 'AAGTGC' : 48, 'AAGTGT' : 52} ,
        'KD' : { 'AAAGAC' : 78, 'AAAGAT' : 22, 'AAGGAC' : 73, 'AAGGAT' : 27} ,
        'KE' : { 'AAAGAA' : 47, 'AAAGAG' : 53, 'AAGGAA' : 45, 'AAGGAG' : 55} ,
        'KF' : { 'AAATTC' : 56, 'AAATTT' : 44, 'AAGTTC' : 46, 'AAGTTT' : 54} ,
        'KG' : { 'AAAGGA' : 29, 'AAAGGC' : 26, 'AAAGGG' : 16, 'AAAGGT' : 29, 'AAGGGA' : 34, 'AAGGGC' : 26, 'AAGGGG' : 14, 'AAGGGT' : 26} ,
        'KH' : { 'AAACAC' : 20, 'AAACAT' : 80, 'AAGCAC' : 73, 'AAGCAT' : 27} ,
        'KI' : { 'AAAATA' : 25, 'AAAATC' : 32, 'AAAATT' : 43, 'AAGATA' : 17, 'AAGATC' : 44, 'AAGATT' : 38} ,
        'KK' : { 'AAAAAA' : 51, 'AAAAAG' : 49, 'AAGAAA' : 41, 'AAGAAG' : 59} ,
        'KL' : { 'AAACTA' : 3, 'AAACTC' : 22, 'AAACTG' : 16, 'AAACTT' : 26, 'AAATTA' : 11, 'AAATTG' : 22, 'AAGCTA' : 10, 'AAGCTC' : 16, 'AAGCTG' : 36, 'AAGCTT' : 18, 'AAGTTA' : 6, 'AAGTTG' : 13} ,
        'KM' : { 'AAAATG' : 32, 'AAGATG' : 68} ,
        'KN' : { 'AAAAAC' : 61, 'AAAAAT' : 39, 'AAGAAC' : 53, 'AAGAAT' : 47} ,
        'KP' : { 'AAACCA' : 28, 'AAACCC' : 20, 'AAACCG' : 19, 'AAACCT' : 33, 'AAGCCA' : 28, 'AAGCCC' : 29, 'AAGCCG' : 14, 'AAGCCT' : 30} ,
        'KQ' : { 'AAACAA' : 49, 'AAACAG' : 51, 'AAGCAA' : 24, 'AAGCAG' : 76} ,
        'KR' : { 'AAAAGA' : 18, 'AAAAGG' : 21, 'AAACGA' : 18, 'AAACGC' : 5, 'AAACGG' : 19, 'AAACGT' : 19, 'AAGAGA' : 18, 'AAGAGG' : 15, 'AAGCGA' : 15, 'AAGCGC' : 18, 'AAGCGG' : 20, 'AAGCGT' : 13} ,
        'KS' : { 'AAAAGC' : 6, 'AAAAGT' : 18, 'AAATCA' : 25, 'AAATCC' : 7, 'AAATCG' : 17, 'AAATCT' : 27, 'AAGAGC' : 23, 'AAGAGT' : 15, 'AAGTCA' : 16, 'AAGTCC' : 19, 'AAGTCG' : 8, 'AAGTCT' : 19} ,
        'KT' : { 'AAAACA' : 11, 'AAAACC' : 27, 'AAAACG' : 30, 'AAAACT' : 32, 'AAGACA' : 13, 'AAGACC' : 34, 'AAGACG' : 25, 'AAGACT' : 28} ,
        'KV' : { 'AAAGTA' : 30, 'AAAGTC' : 14, 'AAAGTG' : 29, 'AAAGTT' : 27, 'AAGGTA' : 26, 'AAGGTC' : 17, 'AAGGTG' : 42, 'AAGGTT' : 14} ,
        'KY' : { 'AAATAC' : 60, 'AAATAT' : 40, 'AAGTAC' : 51, 'AAGTAT' : 49} ,
        'LA' : { 'CTGGCA' : 13, 'CTGGCC' : 39, 'CTGGCG' : 20, 'CTGGCT' : 29, 'TTAGCA' : 10, 'TTAGCC' : 4, 'TTAGCG' : 36, 'TTAGCT' : 50} ,
        'LC' : { 'TTGTGC' : 67, 'TTGTGT' : 33} ,
        'LU' : { 'TTGTGC' : 67, 'TTGTGT' : 33} ,
        'LD' : { 'TTAGAC' : 7, 'TTAGAT' : 93} ,
        'LE' : { 'CTAGAA' : 79, 'CTAGAG' : 21, 'CTGGAA' : 31, 'CTGGAG' : 69, 'TTAGAA' : 47, 'TTAGAG' : 53} ,
        'LF' : { 'CTCTTC' : 59, 'CTCTTT' : 41, 'CTGTTC' : 78, 'CTGTTT' : 22, 'TTATTC' : 90, 'TTATTT' : 10} ,
        'LG' : { 'CTGGGA' : 25, 'CTGGGC' : 35, 'CTGGGG' : 18, 'CTGGGT' : 21, 'TTAGGA' : 5, 'TTAGGC' : 53, 'TTAGGG' : 36, 'TTAGGT' : 6} ,
        'LH' : { 'CTGCAC' : 87, 'CTGCAT' : 13} ,
        'LI' : { 'CTCATA' : 8, 'CTCATC' : 61, 'CTCATT' : 31, 'CTGATA' : 21, 'CTGATC' : 55, 'CTGATT' : 24, 'CTTATA' : 44, 'CTTATC' : 30, 'CTTATT' : 26, 'TTAATA' : 21, 'TTAATC' : 62, 'TTAATT' : 16} ,
        'LK' : { 'CTAAAA' : 58, 'CTAAAG' : 42, 'CTCAAA' : 44, 'CTCAAG' : 56, 'CTGAAA' : 32, 'CTGAAG' : 68, 'CTTAAA' : 52, 'CTTAAG' : 48, 'TTAAAA' : 72, 'TTAAAG' : 28, 'TTGAAA' : 21, 'TTGAAG' : 79} ,
        'LL' : { 'CTACTA' : 16, 'CTACTC' : 2, 'CTACTG' : 16, 'CTACTT' : 34, 'CTATTA' : 2, 'CTATTG' : 30, 'CTGCTA' : 11, 'CTGCTC' : 20, 'CTGCTG' : 42, 'CTGCTT' : 13, 'CTGTTA' : 3, 'CTGTTG' : 11, 'CTTCTA' : 4, 'CTTCTC' : 15, 'CTTCTG' : 23, 'CTTCTT' : 18, 'CTTTTA' : 21, 'CTTTTG' : 19} ,
        'LM' : { 'CTCATG' : 91, 'TTAATG' : 9} ,
        'LN' : { 'CTCAAC' : 56, 'CTCAAT' : 44, 'CTGAAC' : 59, 'CTGAAT' : 41, 'TTAAAC' : 24, 'TTAAAT' : 76, 'TTGAAC' : 62, 'TTGAAT' : 38} ,
        'LP' : { 'CTTCCA' : 29, 'CTTCCC' : 19, 'CTTCCG' : 30, 'CTTCCT' : 22} ,
        'LQ' : { 'CTCCAA' : 34, 'CTCCAG' : 66, 'CTGCAA' : 10, 'CTGCAG' : 90, 'CTTCAA' : 40, 'CTTCAG' : 60, 'TTACAA' : 13, 'TTACAG' : 87} ,
        'LR' : { 'CTCAGA' : 13, 'CTCAGG' : 10, 'CTCCGA' : 12, 'CTCCGC' : 29, 'CTCCGG' : 21, 'CTCCGT' : 15, 'CTGAGA' : 10, 'CTGAGG' : 16, 'CTGCGA' : 9, 'CTGCGC' : 19, 'CTGCGG' : 26, 'CTGCGT' : 19, 'CTTAGA' : 12, 'CTTAGG' : 4, 'CTTCGA' : 19, 'CTTCGC' : 18, 'CTTCGG' : 22, 'CTTCGT' : 25, 'TTAAGA' : 5, 'TTAAGG' : 22, 'TTACGA' : 28, 'TTACGC' : 2, 'TTACGG' : 22, 'TTACGT' : 20, 'TTGAGA' : 27, 'TTGAGG' : 22, 'TTGCGA' : 15, 'TTGCGC' : 17, 'TTGCGG' : 6, 'TTGCGT' : 13} ,
        'LS' : { 'CTGAGC' : 23, 'CTGAGT' : 7, 'CTGTCA' : 11, 'CTGTCC' : 31, 'CTGTCG' : 8, 'CTGTCT' : 20, 'CTTAGC' : 6, 'CTTAGT' : 14, 'CTTTCA' : 21, 'CTTTCC' : 20, 'CTTTCG' : 18, 'CTTTCT' : 22, 'TTAAGC' : 20, 'TTAAGT' : 4, 'TTATCA' : 3, 'TTATCC' : 21, 'TTATCG' : 25, 'TTATCT' : 28} ,
        'LT' : { 'CTAACA' : 78, 'CTAACC' : 7, 'CTAACG' : 1, 'CTAACT' : 14, 'CTCACA' : 21, 'CTCACC' : 51, 'CTCACG' : 14, 'CTCACT' : 14, 'CTTACA' : 26, 'CTTACC' : 26, 'CTTACG' : 4, 'CTTACT' : 44, 'TTAACA' : 37, 'TTAACC' : 4, 'TTAACG' : 14, 'TTAACT' : 45} ,
        'LV' : { 'TTAGTA' : 3, 'TTAGTC' : 16, 'TTAGTG' : 20, 'TTAGTT' : 61} ,
        'LY' : { 'CTTTAC' : 25, 'CTTTAT' : 75} ,
        'MA' : { 'ATGGCA' : 17, 'ATGGCC' : 38, 'ATGGCG' : 16, 'ATGGCT' : 28} ,
        'MF' : { 'ATGTTC' : 60, 'ATGTTT' : 40} ,
        'MI' : { 'ATGATA' : 9, 'ATGATC' : 52, 'ATGATT' : 39} ,
        'MK' : { 'ATGAAA' : 26, 'ATGAAG' : 74} ,
        'MR' : { 'ATGAGA' : 19, 'ATGAGG' : 28, 'ATGCGA' : 8, 'ATGCGC' : 20, 'ATGCGG' : 19, 'ATGCGT' : 7} ,
        'MS' : { 'ATGAGC' : 20, 'ATGAGT' : 14, 'ATGTCA' : 8, 'ATGTCC' : 29, 'ATGTCG' : 10, 'ATGTCT' : 19} ,
        'MT' : { 'ATGACA' : 19, 'ATGACC' : 47, 'ATGACG' : 10, 'ATGACT' : 24} ,
        'MV' : { 'ATGGTA' : 10, 'ATGGTC' : 18, 'ATGGTG' : 58, 'ATGGTT' : 15} ,
        'NA' : { 'AATGCA' : 28, 'AATGCC' : 26, 'AATGCG' : 13, 'AATGCT' : 34} ,
        'ND' : { 'AACGAC' : 18, 'AACGAT' : 82} ,
        'NE' : { 'AACGAA' : 34, 'AACGAG' : 66, 'AATGAA' : 51, 'AATGAG' : 49} ,
        'NF' : { 'AACTTC' : 31, 'AACTTT' : 69, 'AATTTC' : 17, 'AATTTT' : 83} ,
        'NG' : { 'AACGGA' : 24, 'AACGGC' : 30, 'AACGGG' : 16, 'AACGGT' : 30, 'AATGGA' : 20, 'AATGGC' : 34, 'AATGGG' : 13, 'AATGGT' : 33} ,
        'NI' : { 'AACATA' : 6, 'AACATC' : 58, 'AACATT' : 36, 'AATATA' : 48, 'AATATC' : 34, 'AATATT' : 18} ,
        'NK' : { 'AACAAA' : 42, 'AACAAG' : 58, 'AATAAA' : 42, 'AATAAG' : 58} ,
        'NL' : { 'AACCTA' : 10, 'AACCTC' : 22, 'AACCTG' : 34, 'AACCTT' : 17, 'AACTTA' : 3, 'AACTTG' : 15, 'AATCTA' : 21, 'AATCTC' : 14, 'AATCTG' : 20, 'AATCTT' : 22, 'AATTTA' : 1, 'AATTTG' : 21} ,
        'NM' : { 'AACATG' : 100, 'AATATG' : 100} ,
        'NN' : { 'AACAAC' : 67, 'AACAAT' : 33, 'AATAAC' : 26, 'AATAAT' : 74} ,
        'NP' : { 'AACCCA' : 28, 'AACCCC' : 37, 'AACCCG' : 18, 'AACCCT' : 17, 'AATCCA' : 51, 'AATCCC' : 14, 'AATCCG' : 20, 'AATCCT' : 15} ,
        'NQ' : { 'AACCAA' : 26, 'AACCAG' : 74, 'AATCAA' : 54, 'AATCAG' : 46} ,
        'NR' : { 'AACAGA' : 23, 'AACAGG' : 19, 'AACCGA' : 7, 'AACCGC' : 22, 'AACCGG' : 22, 'AACCGT' : 7, 'AATAGA' : 29, 'AATAGG' : 12, 'AATCGA' : 19, 'AATCGC' : 8, 'AATCGG' : 19, 'AATCGT' : 13} ,
        'NS' : { 'AACAGC' : 34, 'AACAGT' : 10, 'AACTCA' : 12, 'AACTCC' : 21, 'AACTCG' : 17, 'AACTCT' : 6, 'AATAGC' : 9, 'AATAGT' : 13, 'AATTCA' : 25, 'AATTCC' : 20, 'AATTCG' : 2, 'AATTCT' : 31} ,
        'NT' : { 'AATACA' : 15, 'AATACC' : 33, 'AATACG' : 16, 'AATACT' : 36} ,
        'NV' : { 'AACGTA' : 2, 'AACGTC' : 14, 'AACGTG' : 51, 'AACGTT' : 32, 'AATGTA' : 20, 'AATGTC' : 21, 'AATGTG' : 39, 'AATGTT' : 21} ,
        'NY' : { 'AACTAC' : 43, 'AACTAT' : 57, 'AATTAC' : 32, 'AATTAT' : 68} ,
        'PA' : { 'CCTGCA' : 12, 'CCTGCC' : 42, 'CCTGCG' : 16, 'CCTGCT' : 30} ,
        'PE' : { 'CCCGAA' : 38, 'CCCGAG' : 62} ,
        'PF' : { 'CCATTC' : 50, 'CCATTT' : 50, 'CCCTTC' : 75, 'CCCTTT' : 25} ,
        'PG' : { 'CCAGGA' : 9, 'CCAGGC' : 34, 'CCAGGG' : 20, 'CCAGGT' : 36, 'CCGGGA' : 25, 'CCGGGC' : 19, 'CCGGGG' : 23, 'CCGGGT' : 33, 'CCTGGA' : 17, 'CCTGGC' : 36, 'CCTGGG' : 21, 'CCTGGT' : 27} ,
        'PH' : { 'CCACAC' : 47, 'CCACAT' : 53, 'CCTCAC' : 42, 'CCTCAT' : 58} ,
        'PI' : { 'CCCATA' : 20, 'CCCATC' : 57, 'CCCATT' : 23, 'CCTATA' : 10, 'CCTATC' : 17, 'CCTATT' : 73} ,
        'PK' : { 'CCAAAA' : 46, 'CCAAAG' : 54, 'CCCAAA' : 39, 'CCCAAG' : 61} ,
        'PL' : { 'CCACTA' : 21, 'CCACTC' : 26, 'CCACTG' : 17, 'CCACTT' : 27, 'CCATTA' : 2, 'CCATTG' : 7, 'CCCCTA' : 5, 'CCCCTC' : 18, 'CCCCTG' : 31, 'CCCCTT' : 15, 'CCCTTA' : 15, 'CCCTTG' : 16, 'CCGCTA' : 26, 'CCGCTC' : 18, 'CCGCTG' : 24, 'CCGCTT' : 2, 'CCGTTA' : 17, 'CCGTTG' : 13, 'CCTCTA' : 5, 'CCTCTC' : 16, 'CCTCTG' : 24, 'CCTCTT' : 21, 'CCTTTA' : 13, 'CCTTTG' : 20} ,
        'PN' : { 'CCAAAC' : 43, 'CCAAAT' : 57, 'CCCAAC' : 67, 'CCCAAT' : 33, 'CCGAAC' : 48, 'CCGAAT' : 52} ,
        'PP' : { 'CCGCCA' : 47, 'CCGCCC' : 21, 'CCGCCG' : 19, 'CCGCCT' : 14} ,
        'PR' : { 'CCAAGA' : 18, 'CCAAGG' : 25, 'CCACGA' : 13, 'CCACGC' : 17, 'CCACGG' : 11, 'CCACGT' : 16, 'CCCAGA' : 12, 'CCCAGG' : 19, 'CCCCGA' : 19, 'CCCCGC' : 17, 'CCCCGG' : 24, 'CCCCGT' : 9, 'CCGAGA' : 18, 'CCGAGG' : 10, 'CCGCGA' : 19, 'CCGCGC' : 28, 'CCGCGG' : 23, 'CCGCGT' : 2, 'CCTAGA' : 14, 'CCTAGG' : 11, 'CCTCGA' : 10, 'CCTCGC' : 19, 'CCTCGG' : 23, 'CCTCGT' : 23} ,
        'PT' : { 'CCGACA' : 1, 'CCGACC' : 32, 'CCGACG' : 26, 'CCGACT' : 41} ,
        'PV' : { 'CCTGTA' : 12, 'CCTGTC' : 12, 'CCTGTG' : 63, 'CCTGTT' : 13} ,
        'PY' : { 'CCCTAC' : 52, 'CCCTAT' : 48} ,
        'QA' : { 'CAGGCA' : 20, 'CAGGCC' : 41, 'CAGGCG' : 7, 'CAGGCT' : 32} ,
        'QC' : { 'CAATGC' : 32, 'CAATGT' : 68, 'CAGTGC' : 59, 'CAGTGT' : 41} ,
        'QU' : { 'CAATGC' : 32, 'CAATGT' : 68, 'CAGTGC' : 59, 'CAGTGT' : 41} ,
        'QD' : { 'CAAGAC' : 37, 'CAAGAT' : 63, 'CAGGAC' : 47, 'CAGGAT' : 53} ,
        'QE' : { 'CAAGAA' : 46, 'CAAGAG' : 54, 'CAGGAA' : 35, 'CAGGAG' : 65} ,
        'QF' : { 'CAATTC' : 8, 'CAATTT' : 92} ,
        'QG' : { 'CAGGGA' : 23, 'CAGGGC' : 25, 'CAGGGG' : 26, 'CAGGGT' : 26} ,
        'QH' : { 'CAGCAC' : 81, 'CAGCAT' : 19} ,
        'QI' : { 'CAAATA' : 26, 'CAAATC' : 48, 'CAAATT' : 26, 'CAGATA' : 19, 'CAGATC' : 41, 'CAGATT' : 40} ,
        'QK' : { 'CAAAAA' : 52, 'CAAAAG' : 48, 'CAGAAA' : 36, 'CAGAAG' : 64} ,
        'QL' : { 'CAACTA' : 3, 'CAACTC' : 16, 'CAACTG' : 19, 'CAACTT' : 25, 'CAATTA' : 18, 'CAATTG' : 19, 'CAGCTA' : 12, 'CAGCTC' : 20, 'CAGCTG' : 39, 'CAGCTT' : 13, 'CAGTTA' : 2, 'CAGTTG' : 14} ,
        'QN' : { 'CAAAAC' : 38, 'CAAAAT' : 62, 'CAGAAC' : 61, 'CAGAAT' : 39} ,
        'QP' : { 'CAACCA' : 32, 'CAACCC' : 18, 'CAACCG' : 33, 'CAACCT' : 17} ,
        'QQ' : { 'CAACAA' : 45, 'CAACAG' : 55, 'CAGCAA' : 22, 'CAGCAG' : 78} ,
        'QR' : { 'CAAAGA' : 7, 'CAAAGG' : 30, 'CAACGA' : 22, 'CAACGC' : 4, 'CAACGG' : 19, 'CAACGT' : 18, 'CAGAGA' : 16, 'CAGAGG' : 20, 'CAGCGA' : 14, 'CAGCGC' : 20, 'CAGCGG' : 16, 'CAGCGT' : 14} ,
        'QS' : { 'CAAAGC' : 22, 'CAAAGT' : 24, 'CAATCA' : 17, 'CAATCC' : 6, 'CAATCG' : 11, 'CAATCT' : 20, 'CAGAGC' : 26, 'CAGAGT' : 11, 'CAGTCA' : 16, 'CAGTCC' : 19, 'CAGTCG' : 16, 'CAGTCT' : 12} ,
        'QT' : { 'CAAACA' : 33, 'CAAACC' : 28, 'CAAACG' : 4, 'CAAACT' : 35} ,
        'QV' : { 'CAAGTA' : 7, 'CAAGTC' : 51, 'CAAGTG' : 26, 'CAAGTT' : 15} ,
        'QY' : { 'CAATAC' : 29, 'CAATAT' : 71} ,
        'RA' : { 'AGAGCA' : 27, 'AGAGCC' : 37, 'AGAGCG' : 16, 'AGAGCT' : 20, 'AGGGCA' : 9, 'AGGGCC' : 31, 'AGGGCG' : 24, 'AGGGCT' : 35, 'CGGGCA' : 20, 'CGGGCC' : 27, 'CGGGCG' : 25, 'CGGGCT' : 27, 'CGTGCA' : 6, 'CGTGCC' : 34, 'CGTGCG' : 24, 'CGTGCT' : 37} ,
        'RD' : { 'AGAGAC' : 36, 'AGAGAT' : 64, 'AGGGAC' : 55, 'AGGGAT' : 45, 'CGAGAC' : 24, 'CGAGAT' : 76, 'CGCGAC' : 82, 'CGCGAT' : 18} ,
        'RE' : { 'AGAGAA' : 50, 'AGAGAG' : 50, 'AGGGAA' : 27, 'AGGGAG' : 73, 'CGAGAA' : 26, 'CGAGAG' : 74, 'CGCGAA' : 43, 'CGCGAG' : 57, 'CGGGAA' : 38, 'CGGGAG' : 62, 'CGTGAA' : 50, 'CGTGAG' : 50} ,
        'RF' : { 'AGATTC' : 52, 'AGATTT' : 48, 'CGCTTC' : 68, 'CGCTTT' : 32} ,
        'RG' : { 'AGAGGA' : 30, 'AGAGGC' : 20, 'AGAGGG' : 21, 'AGAGGT' : 29, 'AGGGGA' : 42, 'AGGGGC' : 23, 'AGGGGG' : 4, 'AGGGGT' : 31, 'CGAGGA' : 44, 'CGAGGC' : 27, 'CGAGGG' : 14, 'CGAGGT' : 15, 'CGCGGA' : 22, 'CGCGGC' : 17, 'CGCGGG' : 40, 'CGCGGT' : 21, 'CGTGGA' : 16, 'CGTGGC' : 28, 'CGTGGG' : 22, 'CGTGGT' : 33} ,
        'RH' : { 'AGACAC' : 77, 'AGACAT' : 23} ,
        'RI' : { 'AGAATA' : 24, 'AGAATC' : 48, 'AGAATT' : 28, 'AGGATA' : 19, 'AGGATC' : 43, 'AGGATT' : 38, 'CGAATA' : 4, 'CGAATC' : 73, 'CGAATT' : 23, 'CGTATA' : 40, 'CGTATC' : 34, 'CGTATT' : 26} ,
        'RK' : { 'AGAAAA' : 46, 'AGAAAG' : 54, 'AGGAAA' : 40, 'AGGAAG' : 60, 'CGCAAA' : 42, 'CGCAAG' : 58} ,
        'RL' : { 'AGACTA' : 7, 'AGACTC' : 4, 'AGACTG' : 26, 'AGACTT' : 26, 'AGATTA' : 13, 'AGATTG' : 24, 'AGGCTA' : 10, 'AGGCTC' : 16, 'AGGCTG' : 30, 'AGGCTT' : 16, 'AGGTTA' : 7, 'AGGTTG' : 22, 'CGCCTA' : 16, 'CGCCTC' : 27, 'CGCCTG' : 31, 'CGCCTT' : 12, 'CGCTTA' : 8, 'CGCTTG' : 6, 'CGTCTA' : 21, 'CGTCTC' : 16, 'CGTCTG' : 27, 'CGTCTT' : 7, 'CGTTTA' : 13, 'CGTTTG' : 16} ,
        'RM' : { 'CGTATG' : 100} ,
        'RN' : { 'AGAAAC' : 53, 'AGAAAT' : 47, 'CGGAAC' : 64, 'CGGAAT' : 36} ,
        'RP' : { 'AGACCA' : 28, 'AGACCC' : 32, 'AGACCG' : 11, 'AGACCT' : 28, 'AGGCCA' : 23, 'AGGCCC' : 41, 'AGGCCG' : 17, 'AGGCCT' : 19, 'CGTCCA' : 20, 'CGTCCC' : 26, 'CGTCCG' : 29, 'CGTCCT' : 24} ,
        'RQ' : { 'AGACAA' : 45, 'AGACAG' : 55, 'AGGCAA' : 28, 'AGGCAG' : 72, 'CGCCAA' : 40, 'CGCCAG' : 60} ,
        'RR' : { 'AGAAGA' : 24, 'AGAAGG' : 17, 'AGACGA' : 21, 'AGACGC' : 21, 'AGACGG' : 15, 'AGACGT' : 2, 'CGAAGA' : 40, 'CGAAGG' : 6, 'CGACGA' : 35, 'CGACGC' : 9, 'CGACGG' : 6, 'CGACGT' : 5, 'CGCAGA' : 8, 'CGCAGG' : 16, 'CGCCGA' : 6, 'CGCCGC' : 33, 'CGCCGG' : 18, 'CGCCGT' : 19, 'CGGAGA' : 15, 'CGGAGG' : 15, 'CGGCGA' : 11, 'CGGCGC' : 21, 'CGGCGG' : 24, 'CGGCGT' : 14, 'CGTAGA' : 20, 'CGTAGG' : 17, 'CGTCGA' : 5, 'CGTCGC' : 19, 'CGTCGG' : 21, 'CGTCGT' : 18} ,
        'RS' : { 'AGAAGC' : 16, 'AGAAGT' : 15, 'AGATCA' : 22, 'AGATCC' : 29, 'AGATCG' : 3, 'AGATCT' : 15, 'AGGAGC' : 24, 'AGGAGT' : 19, 'AGGTCA' : 10, 'AGGTCC' : 14, 'AGGTCG' : 14, 'AGGTCT' : 18, 'CGCAGC' : 32, 'CGCAGT' : 12, 'CGCTCA' : 9, 'CGCTCC' : 18, 'CGCTCG' : 12, 'CGCTCT' : 17, 'CGGAGC' : 30, 'CGGAGT' : 13, 'CGGTCA' : 6, 'CGGTCC' : 19, 'CGGTCG' : 13, 'CGGTCT' : 19} ,
        'RT' : { 'CGCACA' : 15, 'CGCACC' : 39, 'CGCACG' : 15, 'CGCACT' : 31, 'CGGACA' : 29, 'CGGACC' : 34, 'CGGACG' : 13, 'CGGACT' : 24, 'CGTACA' : 9, 'CGTACC' : 17, 'CGTACG' : 29, 'CGTACT' : 45} ,
        'RV' : { 'AGAGTA' : 19, 'AGAGTC' : 15, 'AGAGTG' : 36, 'AGAGTT' : 30, 'AGGGTA' : 3, 'AGGGTC' : 17, 'AGGGTG' : 42, 'AGGGTT' : 38, 'CGGGTA' : 14, 'CGGGTC' : 13, 'CGGGTG' : 71, 'CGGGTT' : 2, 'CGTGTA' : 39, 'CGTGTC' : 20, 'CGTGTG' : 32, 'CGTGTT' : 10} ,
        'RW' : { 'CGCTGG' : 100} ,
        'RY' : { 'AGATAC' : 58, 'AGATAT' : 42, 'CGATAC' : 54, 'CGATAT' : 46} ,
        'SA' : { 'AGTGCA' : 24, 'AGTGCC' : 21, 'AGTGCG' : 11, 'AGTGCT' : 44, 'TCAGCA' : 38, 'TCAGCC' : 20, 'TCAGCG' : 28, 'TCAGCT' : 14, 'TCCGCA' : 7, 'TCCGCC' : 64, 'TCCGCG' : 8, 'TCCGCT' : 21, 'TCTGCA' : 12, 'TCTGCC' : 34, 'TCTGCG' : 17, 'TCTGCT' : 38} ,
        'SC' : { 'TCATGC' : 15, 'TCATGT' : 85, 'TCCTGC' : 37, 'TCCTGT' : 63} ,
        'SU' : { 'TCATGC' : 15, 'TCATGT' : 85, 'TCCTGC' : 37, 'TCCTGT' : 63} ,
        'SD' : { 'TCCGAC' : 55, 'TCCGAT' : 45, 'TCTGAC' : 50, 'TCTGAT' : 50} ,
        'SE' : { 'TCAGAA' : 50, 'TCAGAG' : 50, 'TCGGAA' : 39, 'TCGGAG' : 61} ,
        'SF' : { 'AGTTTC' : 85, 'AGTTTT' : 15, 'TCATTC' : 60, 'TCATTT' : 40, 'TCCTTC' : 79, 'TCCTTT' : 21} ,
        'SG' : { 'TCAGGA' : 15, 'TCAGGC' : 45, 'TCAGGG' : 34, 'TCAGGT' : 6, 'TCGGGA' : 28, 'TCGGGC' : 34, 'TCGGGG' : 24, 'TCGGGT' : 14} ,
        'SI' : { 'TCAATA' : 69, 'TCAATC' : 17, 'TCAATT' : 15, 'TCCATA' : 13, 'TCCATC' : 39, 'TCCATT' : 48} ,
        'SK' : { 'AGCAAA' : 41, 'AGCAAG' : 59, 'AGTAAA' : 83, 'AGTAAG' : 17, 'TCAAAA' : 49, 'TCAAAG' : 51, 'TCCAAA' : 39, 'TCCAAG' : 61, 'TCGAAA' : 10, 'TCGAAG' : 90, 'TCTAAA' : 40, 'TCTAAG' : 60} ,
        'SL' : { 'TCACTA' : 19, 'TCACTC' : 6, 'TCACTG' : 28, 'TCACTT' : 9, 'TCATTA' : 33, 'TCATTG' : 5, 'TCCCTA' : 3, 'TCCCTC' : 13, 'TCCCTG' : 28, 'TCCCTT' : 24, 'TCCTTA' : 10, 'TCCTTG' : 23, 'TCTCTA' : 7, 'TCTCTC' : 19, 'TCTCTG' : 19, 'TCTCTT' : 23, 'TCTTTA' : 13, 'TCTTTG' : 18} ,
        'SN' : { 'AGTAAC' : 47, 'AGTAAT' : 53, 'TCCAAC' : 57, 'TCCAAT' : 43, 'TCGAAC' : 61, 'TCGAAT' : 39, 'TCTAAC' : 78, 'TCTAAT' : 22} ,
        'SP' : { 'TCCCCA' : 16, 'TCCCCC' : 56, 'TCCCCG' : 9, 'TCCCCT' : 18, 'TCGCCA' : 9, 'TCGCCC' : 50, 'TCGCCG' : 14, 'TCGCCT' : 26} ,
        'SQ' : { 'AGCCAA' : 39, 'AGCCAG' : 61, 'AGTCAA' : 35, 'AGTCAG' : 65, 'TCACAA' : 54, 'TCACAG' : 46, 'TCCCAA' : 16, 'TCCCAG' : 84, 'TCTCAA' : 34, 'TCTCAG' : 66} ,
        'SR' : { 'AGTAGA' : 18, 'AGTAGG' : 17, 'AGTCGA' : 12, 'AGTCGC' : 21, 'AGTCGG' : 4, 'AGTCGT' : 27, 'TCAAGA' : 38, 'TCAAGG' : 34, 'TCACGA' : 4, 'TCACGC' : 5, 'TCACGG' : 11, 'TCACGT' : 7, 'TCCAGA' : 10, 'TCCAGG' : 15, 'TCCCGA' : 11, 'TCCCGC' : 27, 'TCCCGG' : 20, 'TCCCGT' : 17, 'TCGAGA' : 17, 'TCGAGG' : 28, 'TCGCGA' : 1, 'TCGCGC' : 28, 'TCGCGG' : 23, 'TCGCGT' : 3, 'TCTAGA' : 18, 'TCTAGG' : 17, 'TCTCGA' : 5, 'TCTCGC' : 17, 'TCTCGG' : 25, 'TCTCGT' : 18} ,
        'SS' : { 'AGTAGC' : 10, 'AGTAGT' : 11, 'AGTTCA' : 30, 'AGTTCC' : 18, 'AGTTCG' : 2, 'AGTTCT' : 28, 'TCAAGC' : 15, 'TCAAGT' : 15, 'TCATCA' : 23, 'TCATCC' : 7, 'TCATCG' : 12, 'TCATCT' : 29} ,
        'ST' : { 'AGCACA' : 34, 'AGCACC' : 46, 'AGCACG' : 9, 'AGCACT' : 11, 'AGTACA' : 22, 'AGTACC' : 22, 'AGTACG' : 22, 'AGTACT' : 34, 'TCCACA' : 18, 'TCCACC' : 58, 'TCCACG' : 7, 'TCCACT' : 17, 'TCGACA' : 40, 'TCGACC' : 13, 'TCGACG' : 6, 'TCGACT' : 41} ,
        'SV' : { 'AGTGTA' : 9, 'AGTGTC' : 37, 'AGTGTG' : 46, 'AGTGTT' : 9, 'TCCGTA' : 4, 'TCCGTC' : 43, 'TCCGTG' : 51, 'TCCGTT' : 2, 'TCGGTA' : 16, 'TCGGTC' : 49, 'TCGGTG' : 24, 'TCGGTT' : 11} ,
        'SW' : { 'TCTTGG' : 100} ,
        'SY' : { 'AGTTAC' : 42, 'AGTTAT' : 58, 'TCGTAC' : 79, 'TCGTAT' : 21} ,
        'TC' : { 'ACGTGC' : 35, 'ACGTGT' : 65} ,
        'TU' : { 'ACGTGC' : 35, 'ACGTGT' : 65} ,
        'TD' : { 'ACAGAC' : 31, 'ACAGAT' : 69, 'ACTGAC' : 62, 'ACTGAT' : 38} ,
        'TE' : { 'ACCGAA' : 44, 'ACCGAG' : 56, 'ACTGAA' : 47, 'ACTGAG' : 53} ,
        'TF' : { 'ACATTC' : 18, 'ACATTT' : 82, 'ACCTTC' : 74, 'ACCTTT' : 26} ,
        'TG' : { 'ACAGGA' : 38, 'ACAGGC' : 32, 'ACAGGG' : 17, 'ACAGGT' : 13, 'ACCGGA' : 7, 'ACCGGC' : 42, 'ACCGGG' : 9, 'ACCGGT' : 43, 'ACTGGA' : 18, 'ACTGGC' : 26, 'ACTGGG' : 24, 'ACTGGT' : 32} ,
        'TH' : { 'ACACAC' : 29, 'ACACAT' : 71} ,
        'TI' : { 'ACAATA' : 38, 'ACAATC' : 34, 'ACAATT' : 28, 'ACCATA' : 7, 'ACCATC' : 52, 'ACCATT' : 42, 'ACTATA' : 18, 'ACTATC' : 28, 'ACTATT' : 53} ,
        'TK' : { 'ACAAAA' : 64, 'ACAAAG' : 36, 'ACCAAA' : 32, 'ACCAAG' : 68, 'ACTAAA' : 39, 'ACTAAG' : 61} ,
        'TL' : { 'ACACTA' : 6, 'ACACTC' : 13, 'ACACTG' : 23, 'ACACTT' : 19, 'ACATTA' : 21, 'ACATTG' : 18} ,
        'TN' : { 'ACAAAC' : 21, 'ACAAAT' : 79, 'ACCAAC' : 71, 'ACCAAT' : 29} ,
        'TQ' : { 'ACCCAA' : 34, 'ACCCAG' : 66, 'ACTCAA' : 30, 'ACTCAG' : 70} ,
        'TR' : { 'ACAAGA' : 25, 'ACAAGG' : 26, 'ACACGA' : 6, 'ACACGC' : 7, 'ACACGG' : 30, 'ACACGT' : 6, 'ACCAGA' : 12, 'ACCAGG' : 16, 'ACCCGA' : 17, 'ACCCGC' : 23, 'ACCCGG' : 15, 'ACCCGT' : 16, 'ACTAGA' : 20, 'ACTAGG' : 16, 'ACTCGA' : 8, 'ACTCGC' : 24, 'ACTCGG' : 25, 'ACTCGT' : 8} ,
        'TS' : { 'ACAAGC' : 8, 'ACAAGT' : 23, 'ACATCA' : 15, 'ACATCC' : 30, 'ACATCG' : 11, 'ACATCT' : 12, 'ACCAGC' : 39, 'ACCAGT' : 11, 'ACCTCA' : 9, 'ACCTCC' : 19, 'ACCTCG' : 8, 'ACCTCT' : 15} ,
        'TT' : { 'ACAACA' : 38, 'ACAACC' : 14, 'ACAACG' : 3, 'ACAACT' : 44, 'ACCACA' : 16, 'ACCACC' : 40, 'ACCACG' : 11, 'ACCACT' : 34} ,
        'TY' : { 'ACATAC' : 44, 'ACATAT' : 56} ,
        'VA' : { 'GTAGCA' : 9, 'GTAGCC' : 24, 'GTAGCG' : 57, 'GTAGCT' : 10, 'GTTGCA' : 9, 'GTTGCC' : 55, 'GTTGCG' : 10, 'GTTGCT' : 26} ,
        'VD' : { 'GTAGAC' : 54, 'GTAGAT' : 46, 'GTGGAC' : 38, 'GTGGAT' : 62} ,
        'VE' : { 'GTAGAA' : 82, 'GTAGAG' : 18, 'GTCGAA' : 3, 'GTCGAG' : 97} ,
        'VF' : { 'GTGTTC' : 31, 'GTGTTT' : 69, 'GTTTTC' : 59, 'GTTTTT' : 41} ,
        'VG' : { 'GTGGGA' : 26, 'GTGGGC' : 36, 'GTGGGG' : 12, 'GTGGGT' : 26, 'GTTGGA' : 18, 'GTTGGC' : 38, 'GTTGGG' : 21, 'GTTGGT' : 24} ,
        'VH' : { 'GTCCAC' : 90, 'GTCCAT' : 10} ,
        'VI' : { 'GTAATA' : 68, 'GTAATC' : 10, 'GTAATT' : 22, 'GTCATA' : 7, 'GTCATC' : 71, 'GTCATT' : 22, 'GTGATA' : 48, 'GTGATC' : 27, 'GTGATT' : 25, 'GTTATA' : 21, 'GTTATC' : 39, 'GTTATT' : 40} ,
        'VK' : { 'GTAAAA' : 86, 'GTAAAG' : 14, 'GTCAAA' : 36, 'GTCAAG' : 64, 'GTGAAA' : 19, 'GTGAAG' : 81, 'GTTAAA' : 50, 'GTTAAG' : 50} ,
        'VL' : { 'GTACTA' : 1, 'GTACTC' : 10, 'GTACTG' : 59, 'GTACTT' : 14, 'GTATTA' : 2, 'GTATTG' : 15, 'GTGCTA' : 15, 'GTGCTC' : 19, 'GTGCTG' : 34, 'GTGCTT' : 8, 'GTGTTA' : 3, 'GTGTTG' : 21, 'GTTCTA' : 4, 'GTTCTC' : 21, 'GTTCTG' : 26, 'GTTCTT' : 23, 'GTTTTA' : 3, 'GTTTTG' : 23} ,
        'VM' : { 'GTAATG' : 100} ,
        'VN' : { 'GTAAAC' : 18, 'GTAAAT' : 82, 'GTGAAC' : 71, 'GTGAAT' : 29} ,
        'VP' : { 'GTACCA' : 13, 'GTACCC' : 19, 'GTACCG' : 59, 'GTACCT' : 9, 'GTGCCA' : 33, 'GTGCCC' : 21, 'GTGCCG' : 27, 'GTGCCT' : 19, 'GTTCCA' : 48, 'GTTCCC' : 15, 'GTTCCG' : 12, 'GTTCCT' : 25} ,
        'VQ' : { 'GTACAA' : 68, 'GTACAG' : 32} ,
        'VR' : { 'GTAAGA' : 46, 'GTAAGG' : 15, 'GTACGA' : 3, 'GTACGC' : 14, 'GTACGG' : 18, 'GTACGT' : 4, 'GTGAGA' : 17, 'GTGAGG' : 17, 'GTGCGA' : 6, 'GTGCGC' : 29, 'GTGCGG' : 25, 'GTGCGT' : 6, 'GTTAGA' : 15, 'GTTAGG' : 12, 'GTTCGA' : 27, 'GTTCGC' : 9, 'GTTCGG' : 15, 'GTTCGT' : 21} ,
        'VS' : { 'GTGAGC' : 25, 'GTGAGT' : 3, 'GTGTCA' : 10, 'GTGTCC' : 31, 'GTGTCG' : 19, 'GTGTCT' : 12, 'GTTAGC' : 6, 'GTTAGT' : 8, 'GTTTCA' : 6, 'GTTTCC' : 15, 'GTTTCG' : 29, 'GTTTCT' : 37} ,
        'WA' : { 'TGGGCA' : 20, 'TGGGCC' : 40, 'TGGGCG' : 5, 'TGGGCT' : 35} ,
        'WE' : { 'TGGGAA' : 41, 'TGGGAG' : 59} ,
        'WK' : { 'TGGAAA' : 41, 'TGGAAG' : 59} ,
        'WN' : { 'TGGAAC' : 67, 'TGGAAT' : 33} ,
        'WR' : { 'TGGAGA' : 16, 'TGGAGG' : 25, 'TGGCGA' : 1, 'TGGCGC' : 20, 'TGGCGG' : 22, 'TGGCGT' : 16} ,
        'WY' : { 'TGGTAC' : 69, 'TGGTAT' : 31} ,
        'YA' : { 'TACGCA' : 20, 'TACGCC' : 48, 'TACGCG' : 13, 'TACGCT' : 18, 'TATGCA' : 9, 'TATGCC' : 52, 'TATGCG' : 3, 'TATGCT' : 37} ,
        'YE' : { 'TATGAA' : 49, 'TATGAG' : 51} ,
        'YG' : { 'TACGGA' : 6, 'TACGGC' : 34, 'TACGGG' : 26, 'TACGGT' : 35} ,
        'YH' : { 'TATCAC' : 15, 'TATCAT' : 85} ,
        'YI' : { 'TACATA' : 29, 'TACATC' : 49, 'TACATT' : 22, 'TATATA' : 6, 'TATATC' : 55, 'TATATT' : 39} ,
        'YK' : { 'TACAAA' : 45, 'TACAAG' : 55, 'TATAAA' : 36, 'TATAAG' : 64} ,
        'YL' : { 'TATCTA' : 2, 'TATCTC' : 10, 'TATCTG' : 20, 'TATCTT' : 16, 'TATTTA' : 31, 'TATTTG' : 21} ,
        'YN' : { 'TATAAC' : 21, 'TATAAT' : 79} ,
        'YP' : { 'TACCCA' : 21, 'TACCCC' : 29, 'TACCCG' : 28, 'TACCCT' : 23} ,
        'YQ' : { 'TACCAA' : 27, 'TACCAG' : 73, 'TATCAA' : 45, 'TATCAG' : 55} ,
        'YR' : { 'TACAGA' : 19, 'TACAGG' : 16, 'TACCGA' : 8, 'TACCGC' : 24, 'TACCGG' : 20, 'TACCGT' : 13} ,
        'YS' : { 'TATAGC' : 7, 'TATAGT' : 31, 'TATTCA' : 21, 'TATTCC' : 13, 'TATTCG' : 6, 'TATTCT' : 22} ,
        'YT' : { 'TATACA' : 34, 'TATACC' : 14, 'TATACG' : 17, 'TATACT' : 35} ,
        'YY' : { 'TACTAC' : 68, 'TACTAT' : 32, 'TATTAC' : 64, 'TATTAT' : 36}}
        #this dictionary serves as reference to check if CC for HEK293-cells should be used.
        CC_evaluation_dict = {
        'AA': {'GCC', 'GCT'},
        'AD': {'GCA', 'GCT'},
        'AE': {'GCA', 'GCT'},
        'AF': {'GCA', 'GCG'},
        'AG': {'GCA', 'GCT'},
        'AI': {'GCA', 'GCC', 'GCG', 'GCT'},
        'AK': {'GCA', 'GCC'},
        'AL': {'GCA', 'GCC', 'GCG'},
        'AM': {'GCA', 'GCC'},
        'AN': {'GCA', 'GCC'},
        'AP': {'GCA', 'GCT'},
        'AQ': {'GCA', 'GCC'},
        'AR': {'GCA', 'GCC', 'GCT'},
        'AS': {'GCG'},
        'AT': {'GCC', 'GCG', 'GCT'},
        'AV': {'GCT'},
        'CE': {'TGT'},
        'CH': {'TGC'},
        'CK': {'TGT'},
        'CP': {'TGT'},
        'CR': {'TGT'},
        'CT': {'TGT'},
        'CW': {'TGC', 'TGT'},
        'UE': {'TGT'},
        'UH': {'TGC'},
        'UK': {'TGT'},
        'UP': {'TGT'},
        'UR': {'TGT'},
        'UT': {'TGT'},
        'UW': {'TGC', 'TGT'},
        'DA': {'GAC', 'GAT'},
        'DC': {'GAT'},
        'DU': {'GAT'},
        'DE': {'GAT'},
        'DF': {'GAC'},
        'DG': {'GAT'},
        'DI': {'GAC', 'GAT'},
        'DK': {'GAT'},
        'DL': {'GAT'},
        'DN': {'GAC', 'GAT'},
        'DQ': {'GAC', 'GAT'},
        'DR': {'GAC', 'GAT'},
        'DT': {'GAT'},
        'DV': {'GAT'},
        'DY': {'GAT'},
        'EA': {'GAA'},
        'EC': {'GAA', 'GAG'},
        'EU': {'GAA', 'GAG'},
        'EE': {'GAA'},
        'EF': {'GAG'},
        'EG': {'GAA', 'GAG'},
        'EH': {'GAA'},
        'EI': {'GAA', 'GAG'},
        'EK': {'GAA', 'GAG'},
        'EL': {'GAG'},
        'EM': {'GAA', 'GAG'},
        'EN': {'GAA', 'GAG'},
        'EP': {'GAA'},
        'EQ': {'GAA', 'GAG'},
        'ER': {'GAA', 'GAG'},
        'ES': {'GAA', 'GAG'},
        'ET': {'GAA', 'GAG'},
        'EY': {'GAG'},
        'FA': {'TTC', 'TTT'},
        'FD': {'TTT'},
        'FE': {'TTC'},
        'FF': {'TTC', 'TTT'},
        'FG': {'TTC', 'TTT'},
        'FH': {'TTT'},
        'FI': {'TTC'},
        'FK': {'TTC', 'TTT'},
        'FL': {'TTC', 'TTT'},
        'FM': {'TTC', 'TTT'},
        'FN': {'TTC', 'TTT'},
        'FP': {'TTC', 'TTT'},
        'FQ': {'TTC', 'TTT'},
        'FR': {'TTC', 'TTT'},
        'FS': {'TTC'},
        'FT': {'TTC', 'TTT'},
        'FV': {'TTC'},
        'FY': {'TTC', 'TTT'},
        'GA': {'GGA', 'GGG'},
        'GD': {'GGG'},
        'GE': {'GGA', 'GGC'},
        'GF': {'GGA', 'GGC'},
        'GG': {'GGA', 'GGC', 'GGG', 'GGT'},
        'GH': {'GGC'},
        'GI': {'GGA', 'GGC', 'GGG', 'GGT'},
        'GK': {'GGA', 'GGC', 'GGG'},
        'GL': {'GGC', 'GGG', 'GGT'},
        'GM': {'GGC', 'GGG'},
        'GN': {'GGA', 'GGC'},
        'GP': {'GGG', 'GGT'},
        'GQ': {'GGC'},
        'GR': {'GGA', 'GGC', 'GGG', 'GGT'},
        'GS': {'GGA', 'GGC', 'GGT'},
        'GT': {'GGC', 'GGT'},
        'GV': {'GGA', 'GGG', 'GGT'},
        'GY': {'GGA', 'GGC'},
        'HA': {'CAT'},
        'HD': {'CAT'},
        'HE': {'CAC', 'CAT'},
        'HF': {'CAT'},
        'HI': {'CAC'},
        'HK': {'CAC', 'CAT'},
        'HL': {'CAC', 'CAT'},
        'HN': {'CAC'},
        'HQ': {'CAC', 'CAT'},
        'HR': {'CAT'},
        'HS': {'CAC', 'CAT'},
        'HT': {'CAC', 'CAT'},
        'IA': {'ATA', 'ATC', 'ATT'},
        'IC': {'ATA'},
        'IU': {'ATA'},
        'ID': {'ATA', 'ATC', 'ATT'},
        'IE': {'ATA', 'ATC', 'ATT'},
        'IF': {'ATA', 'ATC'},
        'IG': {'ATA', 'ATC', 'ATT'},
        'IH': {'ATA', 'ATC', 'ATT'},
        'II': {'ATA', 'ATC', 'ATT'},
        'IK': {'ATA', 'ATC', 'ATT'},
        'IL': {'ATA', 'ATC', 'ATT'},
        'IM': {'ATC', 'ATT'},
        'IN': {'ATA', 'ATC', 'ATT'},
        'IP': {'ATA', 'ATC'},
        'IQ': {'ATA', 'ATC', 'ATT'},
        'IR': {'ATA', 'ATC', 'ATT'},
        'IS': {'ATA', 'ATC', 'ATT'},
        'IT': {'ATA', 'ATC'},
        'IV': {'ATC'},
        'IW': {'ATC'},
        'IY': {'ATA', 'ATC'},
        'KA': {'AAA', 'AAG'},
        'KC': {'AAA', 'AAG'},
        'KU': {'AAA', 'AAG'},
        'KD': {'AAA', 'AAG'},
        'KE': {'AAA', 'AAG'},
        'KF': {'AAA', 'AAG'},
        'KG': {'AAA', 'AAG'},
        'KH': {'AAA', 'AAG'},
        'KI': {'AAA', 'AAG'},
        'KK': {'AAA', 'AAG'},
        'KL': {'AAA', 'AAG'},
        'KM': {'AAA', 'AAG'},
        'KN': {'AAA', 'AAG'},
        'KP': {'AAA', 'AAG'},
        'KQ': {'AAA', 'AAG'},
        'KR': {'AAA', 'AAG'},
        'KS': {'AAA', 'AAG'},
        'KT': {'AAA', 'AAG'},
        'KV': {'AAA', 'AAG'},
        'KY': {'AAA', 'AAG'},
        'LA': {'CTG', 'TTA'},
        'LC': {'TTG'},
        'LU': {'TTG'},
        'LD': {'TTA'},
        'LE': {'CTA', 'CTG', 'TTA'},
        'LF': {'CTC', 'CTG', 'TTA'},
        'LG': {'CTG', 'TTA'},
        'LH': {'CTG'},
        'LI': {'CTC', 'CTG', 'CTT', 'TTA'},
        'LK': {'CTA', 'CTC', 'CTG', 'CTT', 'TTA', 'TTG'},
        'LL': {'CTA', 'CTG', 'CTT'},
        'LM': {'CTC', 'TTA'},
        'LN': {'CTC', 'CTG', 'TTA', 'TTG'},
        'LP': {'CTT'},
        'LQ': {'CTC', 'CTG', 'CTT', 'TTA'},
        'LR': {'CTC', 'CTG', 'CTT', 'TTA', 'TTG'},
        'LS': {'CTG', 'CTT', 'TTA'},
        'LT': {'CTA', 'CTC', 'CTT', 'TTA'},
        'LV': {'TTA'},
        'LY': {'CTT'},
        'MA': {'ATG'},
        'MF': {'ATG'},
        'MI': {'ATG'},
        'MK': {'ATG'},
        'MR': {'ATG'},
        'MS': {'ATG'},
        'MT': {'ATG'},
        'MV': {'ATG'},
        'NA': {'AAT'},
        'ND': {'AAC'},
        'NE': {'AAC', 'AAT'},
        'NF': {'AAC', 'AAT'},
        'NG': {'AAC', 'AAT'},
        'NI': {'AAC', 'AAT'},
        'NK': {'AAC', 'AAT'},
        'NL': {'AAC', 'AAT'},
        'NM': {'AAC', 'AAT'},
        'NN': {'AAC', 'AAT'},
        'NP': {'AAC', 'AAT'},
        'NQ': {'AAC', 'AAT'},
        'NR': {'AAC', 'AAT'},
        'NS': {'AAC', 'AAT'},
        'NT': {'AAT'},
        'NV': {'AAC', 'AAT'},
        'NY': {'AAC', 'AAT'},
        'PA': {'CCT'},
        'PE': {'CCC'},
        'PF': {'CCA', 'CCC'},
        'PG': {'CCA', 'CCG', 'CCT'},
        'PH': {'CCA', 'CCT'},
        'PI': {'CCC', 'CCT'},
        'PK': {'CCA', 'CCC'},
        'PL': {'CCA', 'CCC', 'CCG', 'CCT'},
        'PN': {'CCA', 'CCC', 'CCG'},
        'PP': {'CCG'},
        'PR': {'CCA', 'CCC', 'CCG', 'CCT'},
        'PT': {'CCG'},
        'PV': {'CCT'},
        'PY': {'CCC'},
        'QA': {'CAG'},
        'QC': {'CAA', 'CAG'},
        'QU': {'CAA', 'CAG'},
        'QD': {'CAA', 'CAG'},
        'QE': {'CAA', 'CAG'},
        'QF': {'CAA'},
        'QG': {'CAG'},
        'QH': {'CAG'},
        'QI': {'CAA', 'CAG'},
        'QK': {'CAA', 'CAG'},
        'QL': {'CAA', 'CAG'},
        'QN': {'CAA', 'CAG'},
        'QP': {'CAA'},
        'QQ': {'CAA', 'CAG'},
        'QR': {'CAA', 'CAG'},
        'QS': {'CAA', 'CAG'},
        'QT': {'CAA'},
        'QV': {'CAA'},
        'QY': {'CAA'},
        'RA': {'AGA', 'AGG', 'CGG', 'CGT'},
        'RD': {'AGA', 'AGG', 'CGA', 'CGC'},
        'RE': {'AGA', 'AGG', 'CGA', 'CGC', 'CGG', 'CGT'},
        'RF': {'AGA', 'CGC'},
        'RG': {'AGA', 'AGG', 'CGA', 'CGC', 'CGT'},
        'RH': {'AGA'},
        'RI': {'AGA', 'AGG', 'CGA', 'CGT'},
        'RK': {'AGA', 'AGG', 'CGC'},
        'RL': {'AGA', 'AGG', 'CGC', 'CGT'},
        'RM': {'CGT'},
        'RN': {'AGA', 'CGG'},
        'RP': {'AGA', 'AGG', 'CGT'},
        'RQ': {'AGA', 'AGG', 'CGC'},
        'RR': {'AGA', 'CGA', 'CGC', 'CGG', 'CGT'},
        'RS': {'AGA', 'AGG', 'CGC', 'CGG'},
        'RT': {'CGC', 'CGG', 'CGT'},
        'RV': {'AGA', 'AGG', 'CGG', 'CGT'},
        'RW': {'CGC'},
        'RY': {'AGA', 'CGA'},
        'SA': {'AGT', 'TCA', 'TCC', 'TCT'},
        'SC': {'TCA', 'TCC'},
        'SU': {'TCA', 'TCC'},
        'SD': {'TCC', 'TCT'},
        'SE': {'TCA', 'TCG'},
        'SF': {'AGT', 'TCA', 'TCC'},
        'SG': {'TCA', 'TCG'},
        'SI': {'TCA', 'TCC'},
        'SK': {'AGC', 'AGT', 'TCA', 'TCC', 'TCG', 'TCT'},
        'SL': {'TCA', 'TCC', 'TCT'},
        'SN': {'AGT', 'TCC', 'TCG', 'TCT'},
        'SP': {'TCC', 'TCG'},
        'SQ': {'AGC', 'AGT', 'TCA', 'TCC', 'TCT'},
        'SR': {'AGT', 'TCA', 'TCC', 'TCG', 'TCT'},
        'SS': {'AGT', 'TCA'},
        'ST': {'AGC', 'AGT', 'TCC', 'TCG'},
        'SV': {'AGT', 'TCC', 'TCG'},
        'SW': {'TCT'},
        'SY': {'AGT', 'TCG'},
        'TC': {'ACG'},
        'TU': {'ACG'},
        'TD': {'ACA', 'ACT'},
        'TE': {'ACC', 'ACT'},
        'TF': {'ACA', 'ACC'},
        'TG': {'ACA', 'ACC', 'ACT'},
        'TH': {'ACA'},
        'TI': {'ACA', 'ACC', 'ACT'},
        'TK': {'ACA', 'ACC', 'ACT'},
        'TL': {'ACA'},
        'TN': {'ACA', 'ACC'},
        'TQ': {'ACC', 'ACT'},
        'TR': {'ACA', 'ACC', 'ACT'},
        'TS': {'ACA', 'ACC'},
        'TT': {'ACA', 'ACC'},
        'TY': {'ACA'},
        'VA': {'GTA', 'GTT'},
        'VD': {'GTA', 'GTG'},
        'VE': {'GTA', 'GTC'},
        'VF': {'GTG', 'GTT'},
        'VG': {'GTG', 'GTT'},
        'VH': {'GTC'},
        'VI': {'GTA', 'GTC', 'GTG', 'GTT'},
        'VK': {'GTA', 'GTC', 'GTG', 'GTT'},
        'VL': {'GTA', 'GTG', 'GTT'},
        'VM': {'GTA'},
        'VN': {'GTA', 'GTG'},
        'VP': {'GTA', 'GTG', 'GTT'},
        'VQ': {'GTA'},
        'VR': {'GTA', 'GTG', 'GTT'},
        'VS': {'GTG', 'GTT'},
        'WA': {'TGG'},
        'WE': {'TGG'},
        'WK': {'TGG'},
        'WN': {'TGG'},
        'WR': {'TGG'},
        'WY': {'TGG'},
        'YA': {'TAC', 'TAT'},
        'YE': {'TAT'},
        'YG': {'TAC'},
        'YH': {'TAT'},
        'YI': {'TAC', 'TAT'},
        'YK': {'TAC', 'TAT'},
        'YL': {'TAT'},
        'YN': {'TAT'},
        'YP': {'TAC'},
        'YQ': {'TAC', 'TAT'},
        'YR': {'TAC'},
        'YS': {'TAT'},
        'YT': {'TAT'},
        'YY': {'TAC', 'TAT'}}

#------Autocorrelation Bias------------
#A dictionary for Autocorrelation Bias. 'W' and 'M' were ommitted. Only codons correlated with
#more than one codon are included. Only correlations > 3 Standard Deviations from expected were
#taken into account. In that order, codons for the following AAs are included: L, S, R, G, A, P, V and T.
#(Data from Cannarozzi et al., 2010)
CoBias_dict = {'CTC': ['CTC', 'CTG'], 'CTT': ['CTT', 'CTC'], 'CTG': ['CTG', 'CTC'],
'TCC': ['TCC', 'TCT'], 'TCT': ['TCT', 'TCC'], 'AGC': ['AGC', 'AGT'], 'AGT': ['AGT', 'AGC'],
'CGA': ['CGA', 'CGC', 'CGG', 'AGG'], 'CGC': ['CGC', 'CGA', 'CGT', 'CGG', 'AGG'], 'CGG': ['CGA', 'CGC', 'CGG'],
'GGC': ['GGA', 'GGC', 'GGG'], 'GGA': ['GGA', 'GGC', 'GGG'], 'GGG': ['GGG', 'GGA', 'GGC'],
'GCA': ['GCA', 'GCG'], 'GCG': ['GCG', 'GCA'],
'CCC': ['CCC', 'CCG'],
'GTA': ['GTA', 'GTG'], 'GTG': ['GTG', 'GTA'],
'ACC': ['ACC', 'ACT'], 'ACT': ['ACT', 'ACC'], 'ACA': ['ACA', 'ACG'], 'ACG': ['ACG', 'ACA']}

#
##---------------Defining the functions that will be needed--------------------------------
#

#GC_correction uses the 4-parameter model. Thus, we first need to optimize the parameters A, B, C, D. The correction
#ratio will be callibrated according to the the desired GC%, chosen by the user. This means that the correction ratio
# y == 0 when GC% == desired GC (i.e any deviation from desired GC% will be countered by the correction ratio. The higher
#the deviation, the stronger the correction).
def logistic4(X, A, B, C, D):
    """4PL logistic equation."""
    return ((A-D)/(1.0+((X/C)**B))) + D
def residuals(p, Y, X):
    """Deviations of data from fitted 4PL curve"""
    A,B,C,D = p
    err = Y-logistic4(X, A, B, C, D)
    return err


#A function that returns a list of GC_content-corrected codon weights. It corrects according to the 4-parameter logistic model.
def Correct4_GCcontent(AA, CodonChoices, CodonWeights, GCcontent, lenCdnSeq, lst_Parameters, GC_aim):
    """This function receives as input a list with codons, a list with their weights, the GCcontent,
    the length of the growing codon seq and a list with the fitted parameters A, B, C, D. As a function of
    GC_content, it returns a corrected list of Weights to be used in the random codon selection.
    If CdnSeq is < 10, it will return the unmodified CodonWeights"""
    newWeights = []; cdns_n_vals = {}; aa=AA; des_GC = GC_aim
    #correct weights if the growing seq has incorporated at least 10 nucleotides
    if lenCdnSeq >= 10:
        #if the aa is coded by only one codon, return its weight: [100]
        if aa in ['M', 'W']:
            newWeights = [100]
        #if it has several codons, correct the weights
        else:
            #with optimized parameters, calculate correction ration 'y' according to GC_content:
            x=GCcontent; a, b, c, d = lst_Parameters[0]
            y = abs(logistic4(x, a, b, c, d))
            #define other necessary variables to calculate the new GC_corrected weights
            #(will be filled below in loop)
            prevATtotal = 0; prevGCtotal = 0; newATtotal = 0; newGCtotal = 0
            #create a new dictionary with the codons and values(weights) to be modified
            #also fill in values needed.
            for Num, Cdn in zip(CodonWeights, CodonChoices):
                cdns_n_vals[Cdn]=Num
                if Cdn[2] in ['A', 'T']:
                    prevATtotal+=Num; newATtotal+=Num; newGCtotal+=(y*Num)
                else:
                    prevGCtotal+=Num; newATtotal+=(y*Num); newGCtotal+=Num
            #if GC is high, correct in favor of AT wobbles
            if GCcontent > des_GC:
                for pair in cdns_n_vals.items():
                    if pair[0][2] in ['C', 'G']:
                        newWeights.append(round((1-y)*pair[1],2))
                    else:
                        newWeights.append(round((newATtotal*pair[1])/prevATtotal,2))
            #if GC too low, correct in favor of GC wobbles
            elif GCcontent < des_GC:
                for pair in cdns_n_vals.items():
                    if pair[0][2] in ['A', 'T']:
                        newWeights.append(round((1-y)*pair[1], 2))
                    else:
                        newWeights.append(round((newGCtotal*pair[1])/prevGCtotal, 2))
            #if GC == desired GC, don't modify weights
            else:
                newWeights = CodonWeights
    #if the growing codon seq hasn't incorporated at least 10 codons, don't modify weights
    else:
        newWeights = CodonWeights
    return newWeights


#A function for Autocorrelation Bias. Based on data for Homo sapiens from Cannarozzi et al., 2010.
#Returns a list of Autocorrelation-corrected codon weights to be used in the random codon selection.
def Correct4_Autocorr_Bias(AA, Index, AAseq, NAseq, LstWeights, dict_codons):
    codons_dict = dict_codons
    #continue inside if the current AA is not one of the exceptions.
    if AA not in ['M', 'W']:
        #Inner function to set the scanning window of maximum 25AA, after which scanning makes no sense
        #since autocorrelation bias effect is considered 0.
        def relevantLen(input_AAlist):
            element = 0; output_AAlist = []
            while len(output_AAlist) < 25:
                if len(output_AAlist) == len(input_AAlist):
                    break
                output_AAlist.append(input_AAlist[element])
                element += 1
            return output_AAlist
        #Prepare our AAseq for previous-instance scanning:
        #First shorten the AAseq from its first element to the element we're in, then, reverse it to scan it backwards,
        #and last, make sure that it has max 25 AAs.
        AAseqShort = AAseq[:Index]; AAseqShort_rv = AAseqShort[::-1]; AAseqShorter_rv = relevantLen(AAseqShort_rv)
        #if our AA has previously occurred within the last 25 AAs.
        if AA in AAseqShorter_rv:
            #Scan for previous instance of AA (which position?). Store the index.
            prev = AAseqShorter_rv.index(AA)
            #based on that distance, calculate the effect of codon bias
            Wght = round(((-0.1601*prev)+11.247), 2)
            #Using the distance (index stored in 'prev'), find which synonymous codon was used.
            cdn_start= len(NAseq)-((prev*3)+3); cdn_end = cdn_start+3
            Cdn_used = NAseq[cdn_start:cdn_end]
            #Create a list of codons that will be favored (= toBias). For this, it needs to check the CoBias_dict.
            #if the codon appears in CoBias_dict, it is correlated with more codons than itself.
            if Cdn_used in CoBias_dict.keys():
                toBias = CoBias_dict[Cdn_used] #also correlated with other codons
            else:
                toBias = [Cdn_used] #only correlated with itself
            #define some other variables necessary for calculating the final weights.
            toDistribute = Wght*len(toBias); prevTotal = 0; Cdns = codons_dict[AA][1]; Vals = LstWeights
            #fill in prevTotal (= sum(weights of disfavored codons before being disfavored))
            for codon, weight in zip(Cdns, Vals):
                if codon not in toBias:
                    prevTotal += weight
            newTotal = prevTotal - toDistribute #(= sum(weights of disfavored codons after being disfavored))
            #create the list of autocorrelation bias-corrected weights and return it as output.
            IndexCdns_inDict = [Cdns.index(x) for x in Cdns if x in toBias]
            newWeights = [x+Wght if Vals.index(x) in IndexCdns_inDict else (newTotal*x)/prevTotal for x in Vals]
            return newWeights
        #if first instance of our current AA, no codon bias effect; return a list with unmodified Weights.
        else:
            return LstWeights
    #if it's one of the exceptions, no codon bias effect; just return the unmodified Weight: [100]
    else:
        return LstWeights


#A function that calculates the GC content of the nascent NA seq.
def GCcont(yourSeq):
    """Calculates the GC content in a string and returns the value."""
    if len(yourSeq) == 0:
        return 0
    else:
        GCcontent = ((yourSeq.count('G') + yourSeq.count('C'))/len(yourSeq))*100
        return round(GCcontent, 1)


#A function that inspects a string for restriction sites or other motifs.
def Motifs(yourSeq, RS=False, CpG=False, HP=False, ATs=False, Pyr=False):
    """Inspects a string for restriction sites or other motifs. Input = DNA sequence and the desired option."""
    #check for restriction sites (RS == True). If found returns matching object (m.ob) where the first one locates (left to right). If not returns None.
    if RS:
        RSite = re.search("(GGATCC)|(ACTAGT)", yourSeq) # Only BamHI and SpeI will be controlled for. (TGTACA)|(GAATTC)|(AAGCTT)|(GGTACC)|(TGGCCA)|(CCATGG)|(GCTAGC)|(TCTAGA)|(CTCGAG)
        return RSite
    #or count the number of CGs (CpG == True) in the seq.
    elif CpG:
        CpGs = len(re.findall('CG', yourSeq))
        return CpGs
    #or homopolymer stretches >= 6 (HP == True) and return m.ob/None.
    elif HP:
        homopolymer = re.search('A{6,}|T{6,}|C{6,}|G{6,}', yourSeq)
        return homopolymer
    #or A/T/AT stretches >= 8 (ATs == True) and return m.ob/None:
    elif ATs:
        AT = re.search('[AT]{8,}', yourSeq)
        return AT
    #or Pyrimidine stretches >= 10 (Pyr == True) and return m.ob/None:
    elif Pyr:
        Pyrimidines = re.search('[CT]{10,}', yourSeq)
        return Pyrimidines


#Function to convert a string of NA seq into a list of its codons.
def toCodonList(NA_seq):
    CodonList = []; i=0
    while i != len(NA_seq):
        CodonList.append(NA_seq[i:i+3])
        i += 3
    return CodonList


#Function to avoid overflow when calculating the geometric mean (last step of the Codon Adaptation Index calculation).
def geomean(xs):
    return math.exp(math.fsum(math.log(x) for x in xs) / len(xs))


#Since a bottleneck in translation lies at the initiation step, the first codons (20) have to be as unstructured as
#possible. For this, the following function generates 10 candidates (first 20 codons), and returns the string with the
#highest minimum free energy (MFE) (i.e. a "...less stable structure contributes to the increase of mRNA expression levels." in
# Jia, M, and Li, Y. 2005. https://doi.org/10.1016/j.febslet.2005.08.059).
#The MFE is calculated with seqfold package, developed by JJTimmons (https://pypi.org/project/seqfold/).
def highest_MFE_start(AminoAcid_Seq, tuple_inherited):
    ex_sys, des_GC, codons_dict, CC_dict, CC_evaluation_dict, CoBias_dict, lst_parameters, seq_fold = tuple_inherited
    aaSeq = AminoAcid_Seq[:20] ; lenAASeq = len(aaSeq) ; candidates = {}
    for Round in range(10):
        ATruns_Off = 0 ; PyrRuns_Off = 0; rSite_counter = 0 #avoids looping infinitely
        newSeq = '' ; lenNewSeq = len(newSeq)/3 ; i=0
        while lenNewSeq < lenAASeq:
            counter = 0
            #build the seq before checking for rSites and other motifs
            while counter < 20:
                #In case there are < 10 aa left, this avoids index errors adjusting the number of iterations.
                if lenAASeq-i < 20:
                    counter = i
                aa = aaSeq[i]
                Choices = codons_dict[aa][1]
                #No codon context influence for first codon (i.e. usually ATG), weights are taken from the single codons.
                if i == 0:
                    Wghts = codons_dict[aa][0]
                #Assess if codon weights will be taken from Codon Context dictionary or single codon dictionary.
                else:
                    aa_ = aaSeq[i-1]
                    #Weights from CC if:
                    #the 'expression system' option 1 was chosen (always codon context for this option).
                    if ex_sys == '1' or ex_sys == '2':
                        if aa == "*": #(no need for codon context for Stop codon)
                            Wghts = codons_dict[aa][0]
                        else:
                            Wghts = []
                            for _codon in Choices:
                                bicodon_weight = CC_dict[aa_+aa][newSeq[-3:]+_codon]
                                Wghts.append(bicodon_weight)
                    #(for B-cells or HEK, there is codon context influence only when needed)
                    #or if the AA combination (previous (aa_) + current (aa)) AND previously used codon
                    #are in CC_evaluation_dictionary)
                    elif aa_+aa in CC_evaluation_dict and newSeq[-3:] in CC_evaluation_dict[aa_+aa]:
                        Wghts = []
                        for _codon in Choices:
                            bicodon_weight = CC_dict[aa_+aa][newSeq[-3:]+_codon]
                            Wghts.append(bicodon_weight)
                    #(the weights from the single codon usage table are used in B-cells or HEK when CC isn't needed)
                    #If not, weights are take from single codon usage dictionary.
                    else:
                        Wghts = codons_dict[aa][0]
                Wghts_GC = Correct4_GCcontent(aa, Choices, Wghts, GCcont(newSeq), len(newSeq), lst_parameters, des_GC) #correction of weights according to GC%
                Wghts_CoBias = Correct4_Autocorr_Bias(aa, i, aaSeq, newSeq, Wghts_GC, codons_dict) #correction of weights according to Autocorrelation Bias
                codon = random.choices(Choices, weights=Wghts_CoBias, k=1)
                newSeq += codon[0]
                i += 1
                counter +=1
            #Once the seq is finished, assess various motifs:
            #---Restriction sites:
            rSite = Motifs(newSeq, RS=True) ; rSiteBool = not rSite == None
            if rSiteBool:#if restriction site found:
                if rSite_counter >= 150:
                    newSeq = ''; i = 0; rSite_counter = 0; ATruns_Off = 0; PyrRuns_Off = 0
                else:
                    newSeq = newSeq[:rSite.start()-rSite.start()%3] #slice the seq at the beginning of the codon containing the start of the rSite
                    i = int(len(newSeq)/3)#update the aa position to continue backtranslating in the correct site
                    rSite_counter += 1
            #---Homopolymers >= 6:
            HPoly = Motifs(newSeq, HP=True) ; HPBool = not HPoly == None
            if HPBool:
                newSeq = newSeq[:HPoly.start()-HPoly.start()%3]
                i = int(len(newSeq)/3)
            #---A/T/AT stretches >= 8:
            #a limit of 100 corrections of A/T/AT stretches per sequence is set.
            if not ATruns_Off > 100:
                ATruns = Motifs(newSeq, ATs=True) ; ATrunsBool = not ATruns == None
                if ATrunsBool:
                    newSeq = newSeq[:ATruns.start()-ATruns.start()%3]
                    i = int(len(newSeq)/3)
                    ATruns_Off +=1
            #---Pyrimidine stretches >= 10:
            #a limit of 100 corrections of Pyrimidine stretches per sequence is set.
            if not PyrRuns_Off > 100:
                PyrRuns = Motifs(newSeq, Pyr=True) ; PyrRunsBool = not PyrRuns == None
                if PyrRunsBool:
                    newSeq = newSeq[:PyrRuns.start()-PyrRuns.start()%3]
                    i = int(len(newSeq)/3)
                    PyrRuns_Off +=1

            lenNewSeq = len(newSeq)/3
        #Once candidate finished, calculate MFE and save in dictionary
        MFE = dg(newSeq)
        candidates[MFE]= newSeq
    #Once all candidates finished, return the one with the highest MFE
    return candidates[max(candidates.keys())]


#In order to apply multi-processing, the main while loop for backtranslation had to be converted into a function.
#Arguments needed: the name of the gene, aminoacid sequence to backtranslate, the maximum threshold (set at the beginning of the script).
#and a tuple with the variables that need to be inherited to the parallel child processes.
def back_translate(geneName, AminoAcid_Seq, Max_threshold, tuple_inherited):
    #unpack values from tuple
    ex_sys, des_GC, codons_dict, CC_dict, CC_evaluation_dict, CoBias_dict, lst_parameters, seq_fold = tuple_inherited
    #Defining all the parameters that are needed for the backtranslation
    candidates_dict = {} #to store the 10 candidates.
    Gene_Name = geneName ; aaSeq = AminoAcid_Seq ; lenAASeq = len(aaSeq)
    #Generate the seq start with the highes MFE
    if seq_fold:
        Seq_start = highest_MFE_start(AminoAcid_Seq, tuple_inherited)
    else:
        Seq_start = ''
    #run the backtranslation 10 times to create 10 candidates
    for Round in range(10):
        newSeq = Seq_start ; lenNewSeq = len(newSeq)/3
        MinThreshold = 48 ; MaxThreshold = Max_threshold
        relaxMax = 0 ; relaxMin = 0 ; ATruns_Off = 0 ; PyrRuns_Off = 0 #The counters to relax thresholds and to turn off some motif checkups (avoids getting infinitely stuck).
        rSite_counter = 0 # Restart the seq after 200 to avoid getting stuck infinitely growing and cutting fragments with restriction sites.
        i = int(len(newSeq)/3) #specifying the index position to retrieve the aa to backtranslate from Seq
        #
        #this loop will continue until the seq is completely backtranslated (assessed by size)
        while lenNewSeq < lenAASeq:
        #
            #For speed, build 10 codons before checking rSites (controlled with 'counter').
            counter = 0
            while counter < 10:
                #In case there are < 10 aa left, this avoids index errors adjusting the number of iterations.
                if lenAASeq-i < 10:
                    counter = 10 - (lenAASeq-i)
                #For aminoacid 'aa', randomly select a 'codon' from 'Choices' according to 'Wghts_*' and add it to the 'newSeq' of codons.
                aa = aaSeq[i]
                Choices = codons_dict[aa][1]
                #No codon context influence for first codon (i.e. usually ATG), weights are taken from the single codons.
                if i == 0:
                    Wghts = codons_dict[aa][0]
                else:
                    aa_ = aaSeq[i-1]
                    #Assess if codon weights will be taken from Codon Context dictionary or single codon dictionary.
                    #Weights from CC if:
                    #
                    #the 'expression system' option 1 or 2 was chosen (always codon context for these options).
                    if ex_sys == '1' or ex_sys == '2':
                        if aa == "*": #(no need for codon context for Stop codon)
                            Wghts = codons_dict[aa][0]
                        else:
                            Wghts = []
                            for _codon in Choices:
                                bicodon_weight = CC_dict[aa_+aa][newSeq[-3:]+_codon]
                                Wghts.append(bicodon_weight)
                    #
                    #(for B-cells or HEK, there is codon context influence only when needed)
                    #or if the AA combination (previous (aa_) + current (aa)) AND previously used codon
                    #are in CC_evaluation_dictionary)
                    elif aa_+aa in CC_evaluation_dict and newSeq[-3:] in CC_evaluation_dict[aa_+aa]:
                        Wghts = []
                        for _codon in Choices:
                            bicodon_weight = CC_dict[aa_+aa][newSeq[-3:]+_codon]
                            Wghts.append(bicodon_weight)

                    #Otherwise, weights are take from single codon usage dictionary.
                    #This applies to B-cells or HEK when CC isn't needed.
                    else:
                        Wghts = codons_dict[aa][0]
                    #
                Wghts_GC = Correct4_GCcontent(aa, Choices, Wghts, GCcont(newSeq), len(newSeq), lst_parameters, des_GC) #correction of weights according to GC%
                Wghts_CoBias = Correct4_Autocorr_Bias(aa, i, aaSeq, newSeq, Wghts_GC, codons_dict) #correction of weights according to Autocorrelation Bias
                codon = random.choices(Choices, weights=Wghts_CoBias, k=1)
                newSeq += codon[0]
                i += 1
                counter += 1
                #
            #Once 10 codons have been added to the newSeq, assess various motifs:
            #---Restriction sites:
            rSite = Motifs(newSeq, RS=True); rSiteBool = not rSite == None
            if rSiteBool:#if restriction site found:
                if rSite_counter >= 200:
                    newSeq = Seq_start; i = int(len(newSeq)/3); rSite_counter = 0; ATruns_Off = 0; PyrRuns_Off = 0
                else:
                #Whether the restriction site starts at the beginning of a codon or in the middle, this will slice
                #the seq in the correct site (always at the beginning of the codon containing the start of the rSite).
                    newSeq = newSeq[:rSite.start()-rSite.start()%3]
                    i = int(len(newSeq)/3)#update the aa position to continue backtranslating in the correct site
                    rSite_counter += 1
            #---Homopolymers >= 6:
            HPoly = Motifs(newSeq, HP=True); HPBool = not HPoly == None
            if HPBool:
                newSeq = newSeq[:HPoly.start()-HPoly.start()%3]
                i = int(len(newSeq)/3)
            #---A/T/AT stretches >= 8:
            #a limit of 100 corrections of A/T/AT stretches per sequence is set.
            if not ATruns_Off > 100:
                ATruns = Motifs(newSeq, ATs=True); ATrunsBool = not ATruns == None
                if ATrunsBool:
                    newSeq = newSeq[:ATruns.start()-ATruns.start()%3]
                    i = int(len(newSeq)/3)
                    ATruns_Off += 1
            #---Pyrimidine stretches >= 10:
            #a limit of 100 corrections of Pyrimidine stretches per sequence is set.
            if not PyrRuns_Off > 100:
                PyrRuns = Motifs(newSeq, Pyr=True); PyrRunsBool = not PyrRuns == None
                if PyrRunsBool:
                    newSeq = newSeq[:PyrRuns.start()-PyrRuns.start()%3]
                    i = int(len(newSeq)/3)
                    PyrRuns_Off += 1

            #Before it reaches the while loop above again, if the seq is finished but the GC content is off,
            #start over. Independently, each 10 times it starts over due to GC > MaxThreshold or GC < MinThreshold,
            #the respective threshold is relaxed by 0.5%. Repeat this until sequence passess this test and can exit the main while loop.
            #Also, restart the ATruns_Off and PyrRuns_Off counters if the backtranslation has to start over.
            lenNewSeq = len(newSeq)/3 #length in codons
            if lenNewSeq == lenAASeq:
                GC_content = GCcont(newSeq)
                if GC_content > MaxThreshold:
                    newSeq = Seq_start; relaxMax +=1; i = int(len(newSeq)/3)
                    ATruns_Off = 0; PyrRuns_Off = 0; rSite_counter = 0
                    if relaxMax % 10 == 0:
                        MaxThreshold += 0.5
                elif GC_content < MinThreshold:
                    newSeq = Seq_start; relaxMin +=1; i = int(len(newSeq)/3)
                    ATruns_Off = 0; PyrRuns_Off = 0; rSite_counter = 0
                    if relaxMin % 10 == 0:
                        MinThreshold -= 0.5
                lenNewSeq = len(newSeq)/3
            #
        #AT this point the sequence is finished. Calculate the values for Tournament Selection.
        #Specific weight for GC%:
        Weight_GC = 2
        #---GC_content---
        GC_score = -(abs(des_GC-GC_content)**Weight_GC) #calculated previously
        #---Codon Adaptation Index---
        CAI = 1; RA_list = []
        for AA, codon in zip(aaSeq, toCodonList(newSeq)):
            Index = codons_dict[AA][1].index(codon)
            weight_codon = codons_dict[AA][0][Index]
            weight_maxCodon = max(codons_dict[AA][0])
            codon_RA = weight_codon/weight_maxCodon #relative adaptiveness of the codon
            RA_list.append(codon_RA)
        CAI = geomean(RA_list)*100 #Codon Adaptation Index of our candidate sequence, expressed in %
        #---CpG motifs---
        CpG_score = -((Motifs(newSeq, CpG=True)/lenAASeq)*100) #Number of CGs / length of Seq in codons, expressed in %
        #---FINAL SCORE----
        SeqScore = sum([CAI, GC_score, CpG_score])
        #save the candidate in the candidates_dict
        candidates_dict[SeqScore] = newSeq
        #
    #select the candidate with the highest score from the candidates_dict
    winner_seq = candidates_dict[max(candidates_dict.keys())] #Find the highest key (i.e score) and get the stored seq.
    #
    #--END OF THE FUNCTION--
    #Return the GeneName with the winner NAseq.
    #
    print(f"\n{a_line*30}\n{Gene_Name} SUCCESSFULLY backtranslated!\nLength = {len(winner_seq)}\nGC% = {GCcont(winner_seq)}\n{a_line*30}\n")
    return (Gene_Name, winner_seq)


#this generator controls for empty lines
def nonblank_lines(f):
    for l in f:
        line = l.rstrip()
        if line:
            yield line

#
#This chunk of code runs ONLY in the MAIN script (not in child parallel processes).
#
if __name__ == '__main__':
    #
    #--------------Callibration of the 4 parameters (A, B, C, D) according to --------------
    #                                user-defined desired GC.
    #
    #Data for calculating the 4 parameters (A, B, C, D)
    x_vals = [0.000000001,40,des_GC,70,100] #represents GC_content (to avoid errors, first val isn't 0).
    y_vals = [-1, -0.9, 0, 0.9, 1] #correction ratio as a function of GC_content, when x == desired GC, y == 0
    p0 = [0, 1, 1, 1] #initial guess for parameters A, B, C, D (arbitrary)
    #Optimize A, B, C and D using least squares method.
    lst_parameters = leastsq(residuals, p0, args=(y_vals, x_vals))

    #
    #------------------------Getting the entries from the input file -------------------------
    #
    # The input should be a text file with FASTA format or NAMEtabSEQUENCE for each
    # amino acid sequence to be backtranslated. Empty lines are skipped.
    with open(InFilename, 'r') as f: #check if the file is FASTA or NAMEtabSEQUENCE format.
        first_line = f.readline()
        is_fasta = re.match('>.*[^\t]\n', first_line)
        is_NAMEtabSEQ = re.match('.*\t.*', first_line)

    entries_dict = {}
    with open(InFilename, 'r') as f_in:
        if is_fasta:
            contents = f_in.read()
            entries = re.finditer('(>.*)\n([A-Z\n\*]+)', contents)
            for entry in entries:
                GeneName = entry[1]
                aaSeq = entry[2].replace('\n', '')
                entries_dict[GeneName]=aaSeq
        elif is_NAMEtabSEQ:
            for line in nonblank_lines(f_in):
                GeneName, aaSeq = line.split('\t')
                entries_dict[GeneName] = aaSeq
        else:
            print(f"\nSorry, I cannot process your file: unsupported format.\nMake sure your file is either fasta or each line is 'NAMEtabSEQUENCE'\n")

    #
    #-----------------------Backtranslating and saving the output-----------------------
    #
    # General flow: Create parallel processes to loop over the entries dictionary 'entries_dict'.
    # For each parallel process (e.g. each entry), backtranslate the aaSeq 10 times. Choose the candidate with the
    # best score and save it in the output dictionary 'out_dict'.
    t1 = time.perf_counter() #start time
    out_dict = {} #for the output

    #Create a tuple with the variables that need to be inherited to the child processes
    inherited_tuple = (ex_sys, des_GC, codons_dict, CC_dict, CC_evaluation_dict, CoBias_dict, lst_parameters, seq_fold)

    #Backtranslation in parallel
    with concurrent.futures.ProcessPoolExecutor() as executor:
        results = [executor.submit(back_translate, gene_name, aaSeq, MaxThreshold, inherited_tuple) for gene_name, aaSeq in entries_dict.items()]
        #the output of the function "back_translate" is a tuple = (GeneName, winner_seq).
        #The tuple contains the name of the gene backtranslated and the seq that obtained the highest score (score according to GC%, Codon Adaptation Index (CAI) and CG dinucleotide counts).
        #Below, the sequences are saved in the output dictionary as they are being completed. Stored as GeneName:winner_seq (key:value).
        for process in concurrent.futures.as_completed(results):
            result = process.result()
            out_dict[result[0]] = result[1]
    
    #Save all the sequences in the output file and print status.
    with open (OutFilename, 'a') as f_out:
        for GeneName, NAseq in out_dict.items():
            f_out.write(f"{GeneName}\t{NAseq}\n")

    t2 = time.perf_counter() #stop time
    print(f"\nFinished in {round((t2-t1)/60, 2)} minutes (in secs: {round(t2-t1,2)})\n")

    # All the sequences are backtranslated and saved in the desired output file
    print(f"{a_space*30}ALL THE SEQUENCES HAVE BEEN SUCCESSFULLY BACKTRANSLATED AND SAVED!!\n")

    #
    #------------------------Wanna see the results printed on the screen?---------------------
    #
    answer = input("Would you like to print your new NA sequences on the screen? (y/n)\n>>> "); options = ['y', 'n']
    while answer not in options:
        answer = input("Please choose a valid option:\n(y/n)\n>>> ")
    if answer == 'y':
        with open(OutFilename, 'r') as f:
            for line in f:
                print(f"{line}\n")
            print(f"END\n{a_space}")
