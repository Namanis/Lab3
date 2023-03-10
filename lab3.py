import re

SP33 = r'[!-I]'
SS64 = r'[;-h]'
IP64 = r'[@-h]'
IP2_64 = r'[C-h]'
IP33 = r'[!-I]'
seq = "GCCGGGCCCGGCGCACGAGCCCCACCCCGAGGACGAGCCCCCGGCGGAGGAGGTCCCAAGGGGGCCTTCCCCGCCAGAGGTCCACCAGGGCGAGGGCGGCGAGGAGGGGCAAGGCGAGGCCGGGCCAGCCGAGCCAGGCCTCCAGCCCGGC"
sequences = []
# Open the Fastq file
with open("reads_for_analysis.fastq", "r") as f:
  # Read the first line of the file
    line = f.readline()
  # Initialize a variable to store the quality encoding
    encoding = "unknown"

    # determin how many symbols to check
    check = 0 

    matchesSP33 = 0
    matchesSS64 = 0
    matchesIP64 = 0
    matchesIP2_64 = 0
    matchesIP33 = 0
    readCounter = 0
    bestSoFar = 0

  
  # Read the rest of the file, one line at a time
    while line:

        if line[0] == "@":
      # Read the next line, which is the sequence data
      #jump few lines
            line=f.readline().strip()
            countC = line.count("C")
            countG = line.count("G")
            countLine = len(line)
            percentage = ((countC+countG)/countLine)*100

            if(percentage>=bestSoFar):
                bestLine = line
                bestSoFar = percentage

            if(percentage>=85):
                print(line)
                print(percentage)
                sequences.append(line)
            line = f.readline().strip()
        elif line[0] == "+":
      # Read the next line, which is the quality scores
            quality_line = f.readline().strip()

        
            matchesSP33 = matchesSP33 + len(re.findall(SP33, quality_line))
            matchesSS64 = matchesSS64 + len(re.findall(SS64, quality_line))
            matchesIP64 = matchesIP64 + len(re.findall(IP64, quality_line))
            matchesIP2_64 = matchesIP2_64 + len(re.findall(IP2_64, quality_line))
            matchesIP33 = matchesIP33 + len(re.findall(IP33, quality_line))
            check = check + len(quality_line)
            

            line=f.readline()
        readCounter = readCounter + 1

      
    # Read the next line
    line = f.readline()
# print(matchesSP33)
# print(matchesSS64)
# print(matchesIP64)
# print(matchesIP2_64)
# print(matchesIP33)


if matchesSP33 >= max(matchesSS64, matchesIP64, matchesIP2_64, matchesIP33):
    encoding = "Sanger Phred+33"
elif matchesSS64 > max(matchesSP33, matchesIP64, matchesIP2_64, matchesIP33):
    encoding = "Solexa Solexa+64"
elif matchesIP64 >= max(matchesSP33, matchesSS64, matchesIP2_64, matchesIP33):
    encoding = "Illumina 1.3+ Phred+64"
elif matchesIP2_64 >= max(matchesSP33, matchesIP64, matchesIP64, matchesIP33):
    encoding = "Illumina 1.5+ Phred+64"
elif matchesIP33 >= max(matchesSP33, matchesIP64, matchesIP64, matchesIP2_64):
    encoding = "Illumina 1.8+ Phred+33"
    
# Print the encoding
print(encoding)
print(readCounter)
print(bestLine)


print(seq)
countC = seq.count("C")
countG = seq.count("G")
countLine = len(seq)
percentage = ((countC+countG)/countLine)*100

print(len(bestLine))
print(len(seq))
print(percentage)



from Bio.Blast import NCBIWWW
from Bio.Blast import NCBIXML

# Paimkite po 5 kiekvieno piko vir????n??s sekas ir sukurkite s??ra???? su sekomis
#sequences = ['ATGTCAGTTAATGCTCAGTTAA', 'ATGTCAGTTAATGCTCAGTTAA', 'ATGTCAGTTAATGCTCAGTTAA']

# Sukurkite tu????i?? s??ra????, kuris bus naudojamas lentel??s duomenims saugoti
results = []

# Eikite per kiekvien?? sek?? ir atlikite BLAST paie??k??
for sequence in sequences:
    result_handle = NCBIWWW.qblast("blastn", "nt", sequence, hitlist_size=1)
    blast_results = NCBIXML.read(result_handle)
    
    # I??traukite read'o id ir rasto mikroorganizmo r?????? i?? BLAST rezultat??
    for alignment in blast_results.alignments:
        for hsp in alignment.hsps:
            read_id = alignment.title.split('|')[3]
            organism = alignment.hit_def.split(' ')[0]
            
            # Prid??kite duomenis ?? lentel??
            results.append([read_id, organism])

print(results)# Atspausdinkite lentel?? su read'o id ir r