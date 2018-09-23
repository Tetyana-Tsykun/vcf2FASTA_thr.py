import gzip
import csv
import argparse
import sys

parser = argparse.ArgumentParser(description="script to convert an all sites vcf to FASTA format. FASTA description will be the sample name in the VCF header.Only does one chromosome/region at a time.")
parser.add_argument("-v", "--vcf", action="store", required=True, help="Input VCF file. Should be a multisample vcf, though it should theoretically work with a single sample.")
parser.add_argument("-o", "--out", action="store", required=True, help="Output filename")
parser.add_argument("-c", "--chromosome", action="store", required=True, help="Chromosome to output. Should be something in the first column of the vcf.")
parser.add_argument("-g", "--gzip", action="store_true", required=False, help="Set if the VCF is gzipped.")
parser.add_argument("-thr", "--threshold", action="store", required=False, help="Specify threshold in percent")
parser.add_argument("--noiupac", action="store_true", required=False, help="Change to IUPAC (default - yes)")


args = parser.parse_args()

vcf_in = args.vcf
out_name = args.out
out_chr = args.chromosome

sample_names = []
sample_seqs = []

iupac = {
	'AA' : 'A',
	'CC' : 'C',
	'TT' : 'T',
	'GG' : 'G',
	'AG' : 'R',
	'CT' : 'Y',
	'CG' : 'S',
	'AT' : 'W',
	'GT' : 'K',
	'AC' : 'M',
	'ACG' : 'V',
	'ACT' : 'H',
	'AGT' : 'D',
	'CGT' : 'B',
	'ACGT' : 'N'
}

iupac_reverse = {
	'R' : 'AG',
	'Y' : 'CT',
	'S' : 'CG',
	'W' : 'AT',
	'K' : 'GT',
	'M' : 'AC',
	'V' : 'ACG',
	'H' : 'ACT',
	'D' : 'AGT',
	'B' : 'CGT',
	'N' :'ACGT'
}

#print iupac_reverse.keys()
#default value for threshold
threshold = 0.33
#iupac_switch="yes"
counter_none=0
counter_ok = 0
counter_fail = 0
if args.gzip:
    opener = gzip.open
else:
    opener = open

if args.threshold:
    threshold=float(args.threshold)/100.
	
print "Threshold=", threshold
with opener(vcf_in, 'r') as tsvin:
    tsvin = csv.reader(tsvin, delimiter='\t')

    for row in tsvin:
        if any('##' in strings for strings in row):
            continue
        if any('#CHROM' in strings for strings in row):
            sample_names = row[9:]
            for sample in sample_names:
                sample_seqs.append([""])
            continue

        chrom,pos,id,ref,alt,qual,filter,info,format=row[0:9]
        haplotypes = row[9:]        
        fullset = "%s,%s" %(ref,alt)
#        print fullset
        if chrom == out_chr:
            alt_list = alt.split(",")
            for index,haplotype in enumerate(haplotypes):
				diplo=haplotype.split(":")	
				#Get frequencies 
				total_freq = float(diplo[1])
				ref_freq = float(diplo[3])
				alt1_freq = float(diplo[5].split(",")[0])
				alt2_freq = 0
				#Get alt1_freq only if it exists
				if (len(diplo[5].split(",")[0])!=len(diplo[5])):
					alt2_freq=float(diplo[5].split(",")[1])
				total_sum = ref_freq+alt1_freq+alt2_freq 
#				print "Haplotypes= ", diplo, " Total= ", total_sum, " Ref= ", ref_freq, "Alt1= ", alt1_freq, "Alt2= ", alt2_freq 	
#				print " diplo", diplo[0]
#				Next check basically not needed, just to check how many reads were dropped out by filter
				if (total_freq!=ref_freq+alt1_freq+alt2_freq):
					counter_fail+=1.
#					print "Haplotypes= ", diplo
#					print "Normalisation failed!", "Total= ", total_freq, " Ref= ", ref_freq, "Alt1= ", alt1_freq, "Alt2= ", alt2_freq 
#					print "Difference = ", total_freq-ref_freq-alt1_freq-alt2_freq
				else:
					counter_ok+=1.
#					print "Normalisation OK!"
				result=""
				if diplo[0] == "."or total_sum==0:
					result = "N"*len(ref)
#					print result
				else:
#					base1 = fullset.split(",")[int(diplo[0].split("/")[0])]
#					base2 = fullset.split(",")[int(diplo[0].split("/")[1])]				
#					print diplo[0], " ", fullset, " ", base1, " ", base2 
					
					for i,c in enumerate(ref):
						temp=""
#						print i, " ", c, " ", base2[i]
						if (ref_freq/total_sum>threshold):
							if (iupac_reverse.has_key(c)):
								temp=iupac_reverse[c]
							else:
								temp+=c
						if (alt1_freq/total_sum>threshold):
							#print "alt1", alt.split(",")[0], "temp ", temp
							if (iupac_reverse.has_key(alt.split(",")[0][i])):
								temp+=iupac_reverse[alt.split(",")[0][i]]
							else:
								#print alt1_freq/total_sum, " ", alt, "i", i, "temp", temp
								temp+=alt.split(",")[0][i]
						
						if (alt2_freq/total_sum>threshold):
							if (iupac_reverse.has_key(alt.split(",")[1][i])):
								temp+=iupac_reverse[alt.split(",")[1][i]]
							else:
								temp+=alt.split(",")[1][i]
						#set() to eliminate repeated characters
						temp= ''.join(set(temp)) 
						# condition to avoid looking for single characters produced by set()
						if (len(temp)>1):
							if(args.noiupac):
								result+="N"
							else:
								result+=iupac[''.join(sorted(temp))]
						else:
							result+=temp				
# if none of the variants crossed the threshold result should be filled with N*len(ref) 
				if (len(result)==0):
					result="N"*len(ref)
					#print result
					counter_none+=1		
				sample_seqs[index][0] = sample_seqs[index][0]+result
				#print result
				
        else:
            continue
print "Sum mismatch ratio (fail/OK) = ", float(counter_fail/(counter_fail+counter_ok))
print "Number of cases with no threshold crossed = ", counter_none 
fasta_out = open(out_name, 'w')
for i in range(len(sample_seqs)):
    fasta_out.write(">"+str(sample_names[i])+"\n"+sample_seqs[i][0]+"\n")
fasta_out.close()
