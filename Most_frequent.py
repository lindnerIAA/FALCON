#!/anaconda3/bin/python

#This is a script that only uses the MOST common codons for backtranslation.

import re, pathlib, time, random

a_line = '-' ; a_space = ' ' #for output aesthetics.

#-----------------------Dialogue with the user to set desired options------------
#
print(f"\n{a_space*30}Hi, I am not FALCON but I wish I was! :c\n\n")

#Input file and check if exists in working directory
InFilename = input('Anyways, please write the name of your INPUT file (e.g. "sequences.txt", "AAseqs.fasta")\n>>> ')
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

#Codon usage dictionary for B-cells cells (tRNA-corrected).
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

#A function that inspects a string for restriction sites.
def Motifs(yourSeq):
    """Inspects a string for restriction sites. Input = DNA sequence."""
    #If restriction site found, returns matching object (m.ob) where the first one locates (left to right). If not returns None.
    RSite = re.search("(GGATCC)|(ACTAGT)", yourSeq) # Only BamHI and SpeI will be controlled for. (TGTACA)|(GAATTC)|(AAGCTT)|(GGTACC)|(TGGCCA)|(CCATGG)|(GCTAGC)|(TCTAGA)|(CTCGAG)
    return RSite

#A function that calculates the GC content of the nascent NA seq.
def GCcont(yourSeq):
    """Calculates the GC content in a string and returns the value."""
    if len(yourSeq) == 0:
        return 0
    else:
        GCcontent = ((yourSeq.count('G') + yourSeq.count('C'))/len(yourSeq))*100
        return round(GCcontent, 1)

#this generator controls for empty lines
def nonblank_lines(f):
    for l in f:
        line = l.rstrip()
        if line:
            yield line

#returns the second most common codon
def Second(Weights, choices, index):
    Wghts_2nd = Weights.copy() ; Choices_2nd = choices.copy()
    del Wghts_2nd[index] ; del Choices_2nd[index]
    new_Idx = Wghts_2nd.index(max(Wghts_2nd))

    return Choices_2nd[new_Idx]


#function for backtranslation
def back_translate(geneName, AminoAcid_Seq):
    #Defining all the parameters that are needed for the backtranslation
    Gene_Name = geneName ; aaSeq = AminoAcid_Seq ; lenAASeq = len(aaSeq)
    newSeq = '' ; lenNewSeq = len(newSeq)/3
    second_most_common = False; use_random = 1 #
    i = 0 #specifying the index position to retrieve the aa to backtranslate from Seq
    #
    #this loop will continue until the seq is completely backtranslated (assessed by size)
    while lenNewSeq < lenAASeq:
        #For aminoacid 'aa', select most common codon from 'Choices' and add it to the 'newSeq' of codons.
        aa = aaSeq[i] ; Choices = codons_dict[aa][1] ; Wghts = codons_dict[aa][0] ; Index = Wghts.index(max(Wghts))
        #
        if aaSeq[i] in ['M', 'W']:
            codon = Choices[Index] ; newSeq += codon
        #avoids getting stuck in restriction sites by introducing a randomly chosen codon if needed
        #(when the second most common codon also creates a restriction site).
        elif use_random % 11 == 0:
            codon = random.choices(Choices, weights=Wghts, k=1) ; newSeq += codon[0] ; use_random = 1
        #if a rSite found, the second most common codon has to be used.
        elif second_most_common:
            codon = Second(Wghts, Choices, Index) ; newSeq += codon ; second_most_common = False
        else:
            codon = Choices[Index] ; newSeq += codon
        i += 1
        #
        #Assess restriction sites:
        rSite = Motifs(newSeq); rSiteBool = not rSite == None
        #if restriction site found:
        if rSiteBool:
            second_most_common = True ; use_random += 1
            #slice the seq at the beginning of the codon containing the start of the rSite.
            newSeq = newSeq[:rSite.start()-rSite.start()%3]; i = int(len(newSeq)/3)#update the aa position
            #
        #Refresh the variable storing len(newSeq) before going back up to the next round.
        lenNewSeq = len(newSeq)/3 #length in codons
        #
    #Return the finished sequence.
    print(f"\n{a_line*30}\n{Gene_Name} SUCCESSFULLY backtranslated!\nLength = {lenNewSeq*3}\nGC% = {GCcont(newSeq)}\n{a_line*30}\n")
    return (Gene_Name, newSeq)

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
t1 = time.perf_counter() #start time
out_dict = {} #for the output
#Backtranslate
for Gene_name, AAseq in entries_dict.items():
    NAseq = back_translate(Gene_name, AAseq)[1]
    out_dict[Gene_name] = NAseq
#Save all the sequences in the output file and print status.
with open (OutFilename, 'a') as f_out:
    for GeneName, NAseq in out_dict.items():
        f_out.write(f"{GeneName}\t{NAseq}\n")

t2 = time.perf_counter() #stop time
print(f"\nFinished in {round((t2-t1)/60, 2)} minutes (in secs: {round(t2-t1,2)})\n")

# All the sequences are backtranslated and saved in the desired output file
print(f"{a_space*30}ALL THE SEQUENCES HAVE BEEN SUCCESSFULLY BACKTRANSLATED AND SAVED!!\n")
