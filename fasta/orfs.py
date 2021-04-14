import collections

# read file
read_file = "dogma.fasta.txt"

# DNA formatted to tuple 
dna = collections.namedtuple("DNA", ['description', 'sequence'])

# DNA_list in list format 
DNA_list = list()

# Initialize
description = ""
dna_sequence = ""

# Given in problem 1 (Assuming that we use it)
permitted_code = ["A", "B", "C", "D", "E", "F", "G", "H", "I", "K", "L", "M", "N","P",
                  "Q", "R", "S", "T", "U", "V", "W", "X", "Y", "Z", "*", "-"]

def check_line(line):
    '''
    Check line to see whether not permitted code exists
    '''
    for l in line:
        if (l not in permitted_code):
            print("Warning")
            print(line+ "has not permitted_code, the code is " + l)

with open(read_file, "r") as ifp:
    for line in ifp:
        if (line[0:].find(">") != -1):
            if (dna_sequence != ""):
                DNA_list.append(dna(description, dna_sequence))
            temp = line[1:]
            while (temp.find(">") != -1):
                DNA_list.append(dna(temp[temp.find(" "):temp.find(">")].strip(), ""))
                temp = temp[temp.find(">")+1:] if temp.find(">") != -1 else temp
            description = temp[temp.find(" ")+1:-1] if temp.find(" ") != -1 else ""
            dna_sequence = ""
        else:
            check_line(line.strip())
            dna_sequence += line.strip()
    DNA_list.append(dna(description, dna_sequence))

class Biology():
    '''
    conduct transcription & translation
    input: sequence, description
    output: * frame_position [-3,-2,-1,1,2,3] | starting_base_number | ending_base_number | 
            protein_sequence_length | protein sequence
            (ex: * 1 | 1  | 45 | 7 | ABCDEFG
                 * 2 | 34 | 5 | 3 | NCD)
    '''
    
    def __init__(self, dna):
        '''
        Initializing 
        '''
        self.dna = dna
        self.transcribe(dna)
        # print(self.rna)
        
        self.protein = []
        self.translate(self.rna, 1)
        self.translate(self.c_rna, -1)
        
    def print_protein(self):
        '''
        Output: print in given format (protein)
        '''
        for idx, p in enumerate(self.protein):
            print("*" + " | ".join([str(p[0]), str(p[1]), str(p[2]), str(p[3]), p[4]]))
    
    def transcribe(self, dna):
        '''
        DNA -> RNA
        '''
        transcribe_map = {'A': 'U', 'T': 'A', 'C': 'G', 'G': 'C'} # Transcription
        comp_map = {'U': 'A', 'A': 'U', 'G': 'C', 'C': 'G'} # Complementary
        
        rna = ""
        for d in dna:
            rna += transcribe_map[d]
        self.rna = rna
        
        c_rna = ""
        for d in self.rna:
            c_rna += comp_map[d]
        self.c_rna = c_rna[::-1]
    
    
    def translate(self, rna, direction):
        '''
        RNA -> Protein
        '''
        # Codon table for reference in dictionary w/ val and key 
        translate_map = {"UUU":"F", "UUC":"F", "UUA":"L", "UUG":"L",
            "UCU":"S", "UCC":"s", "UCA":"S", "UCG":"S",
            "UAU":"Y", "UAC":"Y", "UAA":"STOP", "UAG":"STOP",
            "UGU":"C", "UGC":"C", "UGA":"STOP", "UGG":"W",
            "CUU":"L", "CUC":"L", "CUA":"L", "CUG":"L",
            "CCU":"P", "CCC":"P", "CCA":"P", "CCG":"P",
            "CAU":"H", "CAC":"H", "CAA":"Q", "CAG":"Q",
            "CGU":"R", "CGC":"R", "CGA":"R", "CGG":"R",
            "AUU":"I", "AUC":"I", "AUA":"I", "AUG":"M",
            "ACU":"T", "ACC":"T", "ACA":"T", "ACG":"T",
            "AAU":"N", "AAC":"N", "AAA":"K", "AAG":"K",
            "AGU":"S", "AGC":"S", "AGA":"R", "AGG":"R",
            "GUU":"V", "GUC":"V", "GUA":"V", "GUG":"V",
            "GCU":"A", "GCC":"A", "GCA":"A", "GCG":"A",
            "GAU":"D", "GAC":"D", "GAA":"E", "GAG":"E",
            "GGU":"G", "GGC":"G", "GGA":"G", "GGG":"G",}
        
        now = rna
        protein_list = []
        
        starting_base_number = 0
        
        for starting_base_number in range(0, 3):
            frame_position = starting_base_number + 1
            
            while (rna[starting_base_number+3:].find("AUG") != -1): # START CODON
                starting_base_number += rna[starting_base_number+3:].find("AUG")
                ending_base_number = starting_base_number + 3
                protein_sequence = ""
                checking = starting_base_number

                is_finished = 0
                
                while(ending_base_number + 2 < len(rna)):
                    codon = rna[ending_base_number:ending_base_number + 3]
                    
                    if codon == "UAG" and "UAA" and "UGA": # END CODON 
                        is_finished = 1
                        break
                        
                    ending_base_number += 3
                    protein_sequence += translate_map[codon]
                
                if (is_finished == 1):
                    protein_list.append((direction * frame_position, starting_base_number, ending_base_number, len(protein_sequence), protein_sequence))
                starting_base_number += 3
            
        self.protein.extend(protein_list)

for d in DNA_list:
    b = Biology(d.sequence)
    b.print_protein()