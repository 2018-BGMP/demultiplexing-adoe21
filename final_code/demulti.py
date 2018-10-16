import argparse
import gzip

def get_arguments():
	parser = argparse.ArgumentParser(description="Program to demultiplex sequencing reads")
	parser.add_argument("-f1", "--filename1", help="sequence file", required=True, type=str)
	parser.add_argument("-f2", "--filename2", help="index file", required=True, type=str)
	parser.add_argument("-f3", "--filename3", help="Reverse index file", required=True, type=str)
	parser.add_argument("-f4", "--filename4", help="Reverse sequence file", required=True, type=str)
	parser.add_argument("-i", "--indexfile", help="File of correct barcodes", required=True, type=str)
	return parser.parse_args()
	
args = get_arguments()

def reverseseq(string):
	'''This function inputs a string of DNA and outputs the reverse compliment of that DNA string'''
	'''Example: input of ACCCGTT will return AACGGGT'''
	revstring = string[::-1]
	listrev = list(revstring)
	for i in range(len(listrev)):
		if listrev[i] == "A":
			listrev[i] = "T"
			continue
		if listrev[i] == "T":
			listrev[i] = "A"
			continue
		if listrev[i] == "C":
			listrev[i] = "G"
			continue
		if listrev[i] == "G":
			listrev[i] = "C"
			continue 
	revcomp = "".join(listrev)
	return revcomp

def convert_phred(letter):
    "Converts a single character into a phred score"
    score = ord(letter)
    score = score - 33
    return score
indexdict = {}
fhctr = 1
with open(args.indexfile,"r") as indexes:
	for lines in indexes:
		lines=lines.strip("\r\n")                                             #Creates list of good indexes
		indexdict[lines] = open(lines + "forwardout.fastq", "a")
		indexdict[reverseseq(lines)] = open(lines + "reverseout.fastq","a")   #Opens output files for writing good outputs to correct files
with gzip.open(args.filename1, "rt") as file1:
	with gzip.open(args.filename2, "rt") as file2:
		with gzip.open(args.filename3, "rt") as file3:
			with gzip.open(args.filename4, "rt") as file4:                    #Opens 4 input files
				with open("badqscoreforward.fastq", "a") as badqc1:
					with open("badqscorereverse.fastq", "a") as badqc2:
						with open("forwardbad.fastq", "a") as bad1:
							with open("reversebad.fastq", "a") as bad2:        #Opens 4 ouput files to write the reads that either have a bad qscore or
								ctr = 0                                        #index hopping to
								hopctr = 0
								goodctrdict = {}
								var = "start"
								while var != "":                              #Checks if we are at the end of the input files
									header1=file1.readline()
									var = header1
									if var == "":                             #Ends loop if we parsed through all of the input files
										continue
									sequence1=file1.readline()
									comment1=file1.readline()
									qscore1=file1.readline()
									header2=file2.readline()
									sequence2=file2.readline()
									comment2=file2.readline()
									qscore2=file2.readline()
									header3=file3.readline()
									sequence3=file3.readline()
									comment3=file3.readline()
									qscore3=file3.readline()
									header4=file4.readline()
									sequence4=file4.readline()                #Saves first read from each file and saves them as variables
									comment4=file4.readline()                 #This rewrites the variables for each read so we don't save too 
									qscore4=file4.readline()                  #much to memory
									if len(sequence2) < 2:
										continue
									if sequence2[:-1] == reverseseq(sequence3[:-1]) and sequence2[:-1] in indexdict:    #Checks that the index didn't hop and is
										qscore_list1 = []                                                               #one of the good indexes we are looking for
										qscore_list2 = []
										for elem in qscore1[:-1]:
											qscore_list1.append(convert_phred(elem))
										for item in qscore4[:-1]:
											qscore_list2.append(convert_phred(item))
										if sum(qscore_list1)/len(qscore_list1) < 25 or sum(qscore_list2)/len(qscore_list2) < 25:    #Checks that read has a good qscore
											badqc1.write(header1[:-1] + ":" + sequence2)
											badqc1.write(sequence1)
											badqc1.write(comment1)
											badqc1.write(qscore1)
											badqc2.write(header4[:-1] + ":" + sequence3)
											badqc2.write(sequence4)
											badqc2.write(comment4)
											badqc2.write(qscore4)                           #Writes reads to bad qscore output file if qscore of read is bad
										else:
											if sequence2[:-1] not in goodctrdict:
												goodctrdict[sequence2[:-1]] = 1
											if sequence2[:-1] in goodctrdict:
												goodctrdict[sequence2[:-1]] += 1 
											indexdict[sequence2[:-1]].write(header1[:-1] + ":" + sequence2)
											indexdict[sequence2[:-1]].write(sequence1)
											indexdict[sequence2[:-1]].write(comment1)
											indexdict[sequence2[:-1]].write(qscore1)
											indexdict[sequence3[:-1]].write(header4[:-1] + ":" + sequence3)
											indexdict[sequence3[:-1]].write(sequence4)
											indexdict[sequence3[:-1]].write(comment4)
											indexdict[sequence3[:-1]].write(qscore4)      #Writes good reads to proper output files in proper format
												
									else:
										hopctr += 1 
										bad1.write(header1[:-1] + ":" + sequence2)
										bad1.write(sequence1)
										bad1.write(comment1)
										bad1.write(qscore1)
										bad2.write(header4[:-1] + ":" + sequence3)
										bad2.write(sequence4)
										bad2.write(comment4)
										bad2.write(qscore4)                              #Writes reads that index hopped to the proper output files
for elem in indexdict:
	indexdict[elem].close()
	indexdict[reverseseq(elem)].close()							#Closes ouput files I wrote good outputs to 		
indexhoppercent = hopctr/(1452986940/4)*100		                    #Gives a percent of index hopping 
for key in goodctrdict:
	goodctrdict[key] = goodctrdict[key]/(1452986940/4)*100          #Gives a percent of index hopping per index
print("The indexes hopped " + str(hopctr) + " times") 
print("The index hopped " + str(indexhoppercent) + "% of the time")
print("Index" + "\t" + "Percent of Total Reads from Index")
for elem in goodctrdict:
	print(str(elem) + "\t" + str(goodctrdict[elem]))

	