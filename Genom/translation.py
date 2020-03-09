import re
from Bio import SeqIO
from Bio import pairwise2
from Bio.SubsMat import MatrixInfo as matlist
import os

file = open('C:/Users/90531/Desktop/Genom/nCov.fasta', 'r')
#file.seek(69)   #sars
file.seek(100)  #nCov
dna=file.read()        

print ("DNA Sequence: ", dna)
a = dna.replace('\n',"")                   #boşlukları kaldır

# DNA codon table
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

print("-----------------------------------------------------------------------------------------------------------------")
aminoacid_sequence = ""                                   #aminoasit sentezleme 
 
for i in range(0, len(a)-(3+len(a)%3), 3):                      
    
    aminoacid_sequence += codon_table[a[i:i+3]] 
       
print ("Aminoacid Sequence: ", aminoacid_sequence )


print("----------------------------------------------------------------------------------------------------------------")

print("ORF")
records = SeqIO.parse('C:/Users/90531/Desktop/sars.fasta', 'fasta')
for record in records:
    for strand, seq in (1, record.seq), (-1, record.seq.reverse_complement()):
        for frame in range(3):
            length = 3 * ((len(seq)-frame) // 3)
            for pro in seq[frame:frame+length].translate(table = 1).split("*")[:-1]:
                if 'M' in pro:
                    orf = pro[pro.find('M'):]
                    pos = seq[frame:frame+length].translate(table=1).find(orf)*3 + frame +1
                    if len(orf)*3 +3 > 7:
                        print("{}...{} - length {}, strand {}, frame {}, pos {}, name {}".format\
                               (orf[:3], orf[-3:], len(orf)*3+3, strand, frame, pos, record.id))


