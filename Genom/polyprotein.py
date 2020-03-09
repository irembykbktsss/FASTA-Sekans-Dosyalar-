from Bio import pairwise2
from Bio.SubsMat import MatrixInfo as matlist
import os

codon_table = {"TTT" : "F", "CTT" : "L", "ATT" : "I", "GTT" : "V",
           "TTC" : "F", "CTC" : "L", "ATC" : "I", "GTC" : "V",
           "TTA" : "L", "CTA" : "L", "ATA" : "I", "GTA" : "V",
           "TTG" : "L", "CTG" : "L", "ATG" : "M", "GTG" : "V",
           "TCT" : "S", "CCT" : "P", "ACT" : "T", "GCT" : "A",
           "TCC" : "S", "CCC" : "P", "ACC" : "T", "GCC" : "A",
           "TCA" : "S", "CCA" : "P", "ACA" : "T", "GCA" : "A",
           "TCG" : "S", "CCG" : "P", "ACG" : "T", "GCG" : "A",
           "TAT" : "Y", "CAT" : "H", "AAT" : "N", "GAT" : "D",
           "TAC" : "Y", "CAC" : "H", "AAC" : "N", "GAC" : "D",
           "TAA" : "STOP", "CAA" : "Q", "AAA" : "K", "GAA" : "E",
           "TAG" : "STOP", "CAG" : "Q", "AAG" : "K", "GAG" : "E",
           "TGT" : "C", "CGT" : "R", "AGT" : "S", "GGT" : "G",
           "TGC" : "C", "CGC" : "R", "AGC" : "S", "GGC" : "G",
           "TGA" : "STOP", "CGA" : "R", "AGA" : "R", "GGA" : "G",
           "TGG" : "W", "CGG" : "R", "AGG" : "R", "GGG" : "G" 
           }

file1 = open('C:/Users/90531/Desktop/Genom/sars.fasta', 'r')
file1.seek(69)
dna1=file1.read()
s1= dna1.replace('\n',"") 
sequence1 = ""                        
for i in range(0, len(s1)-(3+len(s1)%3), 3):                      
    
    sequence1 += codon_table[s1[i:i+3]]


file2 = open('C:/Users/90531/Desktop/Genom/nCov.fasta', 'r')
file2.seek(100)
dna2=file2.read()
s2 = dna2.replace('\n',"") 
sequence2 = "" 
for j in range(0, len(s2)-(3+len(s2)%3), 3):                      
    
    sequence2 += codon_table[s2[i:i+3]]


def format_alignment(align1, align2, score, begin, end, full_sequences=False):
    align_begin = begin 
    align_end = end 
    start1 = start2 = "" 
    start_m = begin  # Begin of match line (how many spaces to include) 
    # For local alignments: 
    if not full_sequences and (begin != 0 or end != len(align1)): 
    # Calculate the actual start positions in the un-aligned sequences 
        # This will only work if the gap symbol is '-' or ['-']! 
        start1 = str(len(align1[:begin]) - align1[:begin].count("-") + 1) + " " 
        start2 = str(len(align2[:begin]) - align2[:begin].count("-") + 1) + " " 
        start_m = max(len(start1), len(start2)) 
    elif full_sequences: 
        start_m = 0 
        begin = 0 
        end = len(align1) 
   
    if isinstance(align1, list): 
    # List elements will be separated by spaces, since they can be 
    # of different lengths 
        align1 = [a + " " for a in align1] 
        align2 = [a + " " for a in align2] 
   
    s1_line = ["{:>{width}}".format(start1, width=start_m)]  # seq1 line 
    m_line = [" " * start_m]  # match line 
    s2_line = ["{:>{width}}".format(start2, width=start_m)]  # seq2 line 
   
    for n, (a, b) in enumerate(zip(align1[begin:end],  align2[begin:end])): 
        # Since list elements can be of different length, we center them, 
        # using the maximum length of the two compared elements as width 
        m_len = max(len(a), len(b)) 
        s1_line.append("{:^{width}}".format(a, width=m_len)) 
        s2_line.append("{:^{width}}".format(b, width=m_len)) 
        if full_sequences and (n < align_begin or n >= align_end): 
            m_line.append("{:^{width}}".format(" ", width=m_len))  # space 
            continue 
        if a == b: 
            m_line.append("{:^{width}}".format("|", width=m_len))  # match 
        elif a.strip() == "-" or b.strip() == "-": 
            m_line.append("{:^{width}}".format(" ", width=m_len))  # gap 
        else: 
            m_line.append("{:^{width}}".format(".", width=m_len))  # mismatch 
   
    s2_line.append("\n  Score=%g\n" % score) 
    return "\n".join(["".join(s1_line), "".join(m_line), "".join(s2_line)]) 
    
sars = s1[225:21445]
nCov = s2[265:21554]
matrix = matlist.blosum62
for a in pairwise2.align.globaldx(sars, nCov, matrix):
    print (format_alignment(*a))