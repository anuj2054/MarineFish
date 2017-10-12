#!/usr/bin/env python

Usage = """This is a program for to do an entire reciprocal BLAST search.

The nomenclature is as follows:

Sequences A vs Sequences B is the BLAST search

Step one = Seqs A vs BLAST database of Seqs B

Step two (this step) = Seqs B vs BLAST database of Seqs A

The columns of the results table MUST be as follows:

Column 1: query id
Column 2: query length
Column 3: subject id
Column 4: subject length
Column 5: bitscore
Column 6: e-value
Column 7 onwards can be however you like. 

"""

#import any modules needed
import subprocess

#Give information required
print "\nA program to run a full reciprocal BLAST program. \n\
If you are having any problems, it is worth reading the usage for the assumptions and lay out\n\n"

import sys
print sys.argv
#The inputs required are as follows:\n\n\
#\
#1. Project Name\n\n\
#\
#2. Name of the fasta files containing the first [A] set of sequences (must be in the working folder)\n\n\
#\
#3. Name of the fasta files containing the second [B] set of sequences (must be in the working folder)\n\n\
#\
#4. Name of the second database to search against [A] (this is must be in the BLAST db folder or the working folder)\n\n\
#\
#5. Name of the first database to search against [B] (this is must be in the BLAST db folder or the working folder)\n\n\
#\
#6. What kind of BLAST search are you doing?\n\n\
#\
#7. What e-value cut-off you would like?\n\n\"
#


#Get re required inputs
#ProjectName = raw_input("1. What is the project name? ")
#ProjectName = "RBHAbylopsisTetragona3"
ProjectName = sys.argv[1]

#SeqsFastaA = raw_input("2. What is the name of the fasta file where the first sequences [A] are kept? ")
#SeqsFastaA = "/media/dnarules/BioSurvey/dnarules/Anuj/perseus-data/genesUniprotdownloaded/Uniprot187UniformHeaders.fasta"
#SeqsFastaA = "/media/dnarules/BioSurvey/dnarules/Anuj/perseus-data/genesUniprotdownloaded/Uniprot187OriginalHeaders.fasta"
SeqsFastaA = sys.argv[2]

#SeqsFastaB = raw_input("3. What is the name of the fasta file where the second sequences [B] are kept? ")
#SeqsFastaB = "/media/dnarules/BioSurvey/dnarules/Anuj/perseus-data/genomesandgenes/AbylopsisTetragona.fasta"
#SeqsFastaB = "/media/dnarules/BioSurvey/dnarules/Anuj/perseus-data/genomesUniformHeaderWithLCLWithCnidaria/ThirdAbylopsisTetragona.fasta"
#SeqsFastaB = "/media/dnarules/BioSurvey/dnarules/Anuj/perseus-data/genomesAllFromCasey/AbylopsisTetragona.fasta"
SeqsFastaB = sys.argv[3]

#BlastDbA = raw_input("4. What is the name of the first BLAST database [A] for the search? ")
#BlastDbA = "/media/dnarules/BioSurvey/dnarules/Anuj/Uniprot187/Uniprot187"
#BlastDbA = "/media/dnarules/BioSurvey/dnarules/Anuj/Uniprot187_2/Uniprot187_2"
BlastDbA = sys.argv[4]

#BlastDbB = raw_input("5. What is the name of the second BLAST database [B] for the search? ")
#BlastDbB = "/media/dnarules/BioSurvey/dnarules/Anuj/RBHtest/ThirdAbylopsisTetragona"
#BlastDbB = "/media/dnarules/BioSurvey/dnarules/Anuj/transcriptome2/ThirdAbylopsisTetragona.fasta"
BlastDbB = sys.argv[5]

#BlastType = raw_input("6. What type of blast search are you using? ")
#BlastType = "blastn"

#Evalue = 1e-05raw_input("7. What e-value would you like to use? e.g. 1e-05 ")
Evalue = 1e-05



#####################################################################################3
#BLAST SEARCH ONE!

#Set up the results table
ResultsTableStep1 = ProjectName + "_resultstableA.txt"

#Make the Blastcommand
Blastcommand1 = "tblastn -query " + SeqsFastaA + " -db " + str(BlastDbB) +\
" -out " + ResultsTableStep1 + " -outfmt \"6 qseqid qlen sseqid slen qframe qstart qend " +\
"sframe sstart send evalue bitscore pident nident length\" -num_alignments 1 -evalue " + str(Evalue)


print "\nRunning the first BLAST search\n\
The following BLAST command is running.... \n" + Blastcommand1 + "\n"

ShellOutput = subprocess.call(Blastcommand1, shell=True)

#wait for the program to finish

#Check it has worked, the output should be 0

if ShellOutput == 0:
	print "\nBlast search 1 has worked. Now pulling out the reciprocal hits...\n"
	
else:
	print "\nThere has been an error with the BLAST search, you might need to run it manually.\n\
	then run the next part of this program by itself, using reciprocalblaststep2.py\n"


###################################################################################33
#BLAST SEARCH TWO!

#Get names of the BLAST hits and put into a list
#Set up the list of sequence names and the line counter
SeqsToBlast = []
LineCounter = 0

#Open the Blast results table
InFileA = open(ResultsTableStep1, 'r')
ResultsTableA = "".join(InFileA)

#Replace \r with \n if needed
ResultsTableA = ResultsTableA.replace("\r","\n")
BlastResultsA = ResultsTableA.split('\n')

WarningMessages = 0

#Filter for line e-value and get out unique seqs
for LineA in BlastResultsA:
	LineSplitA = LineA.split('\t')
	if len(LineSplitA) > 1: 
		print LineSplitA[10]
		#this was changed from LineSplitA[5] to LineSplitA[10]
		####################################################3
		if float(LineSplitA[10]) < Evalue and (LineSplitA[2]) not in SeqsToBlast:
			SeqsToBlast.append(LineSplitA[2])
			LineCounter += 1

#Now retrieve the BLAST hits from the BLAST search
print "Retrieving the BLAST hits from the FASTA file...\n"
print "Number of lines present in BLAST result table = " + str(len(BlastResultsA))
print "Number of unique sequences to find = " + str(len(SeqsToBlast))


#Open the fasta file:
FastaFile = open(SeqsFastaB, 'r')
FastaFile = "".join(FastaFile)

#Replace \r with \n if needed, the split the file
FastaFile = FastaFile.replace("\r","\n")
FastaSeqs = FastaFile.split('>')

#Go through fasta file, checking if sequences is in SeqsToBlast
#If it is, then the sequence is added to FastaOutput file
#Also count the number of sequences retrieved, so we know its working... 
FastaOutput = ""
SequenceCounter = 0
#print SeqsToBlast
for Seq in FastaSeqs:
	FastaSplit = Seq.split("\n")
	if len(FastaSplit) > 1:
		#print FastaSplit[0]
		Name = FastaSplit[0].split(" ")[0]
		if Name in SeqsToBlast:
			Sequence = "".join(FastaSplit[1:]) 
			FastaOutput = FastaOutput + ">" + str(FastaSplit[0]) + "\n" + Sequence + "\n"
			SequenceCounter += 1

#Name for fasta file for second blast search:

SeqsForBlast2 = ProjectName + "_seqsforsearch2.fasta"

#Write the output to file
WriteOutFile = True
if WriteOutFile == True:
	OutFileName = SeqsForBlast2
	OutFile = open(OutFileName, 'w')
	OutFile.write(FastaOutput)

print "Number of seqences added = " + str(SequenceCounter)

if SequenceCounter != len(SeqsToBlast):
	print "\nWARNING! WARNING! Not all sequences were found! WARNING! WARNING!\n"
	WarningMessages += 1	
	
OutFile.close()
InFileA.close()

#Make the name for the second results table

ResultsTableStep2 = ProjectName + "_resultstableB.txt"

#Now perform the reciprocal blast search. 
#Make the blastcommand
#At the moment its pretty simple, it needs to be made more complicated in future
Blastcommand2 = "blastx -query " + SeqsForBlast2 + " -db " + str(BlastDbA) +\
" -out " + ResultsTableStep2 + " -outfmt \"6 qseqid qlen sseqid slen qframe qstart qend " +\
"sframe sstart send evalue bitscore pident nident length\" -num_alignments 1 -evalue " + str(Evalue)

print "\nRunning the second BLAST search\n\
The following BLAST command is running.... \n" + Blastcommand2 + "\n"

ShellOutput = subprocess.call(Blastcommand2, shell=True)

#wait for the program to finish

#Check it has worked, the output should be 0

if ShellOutput == 0:
	print "Blast search 2 has worked. Now pulling out the reciprocal hits...\n"
	
else:
	print "There has been an error with the BLAST search, you might need to run it manually.\n\
	then run the third part of this program by itself, reciprocalblaststep3.py"

##########################################################################################
#Now pull out the reciprocal hits
#Open Blast results table and make a list out of the results.
#Open the Blast results table

InFileB = open(ResultsTableStep2, 'r')
ResultsTableB = "".join(InFileB)

#Replace \r with \n if needed, the split the file
ResultsTableB = ResultsTableB.replace("\r","\n")
BlastResultsB = ResultsTableB.split('\n')

#Go through results table and make a dictionary out of it,
#Only if evalue is appropriate
#First make a dictionary 
ResultsB = {}

for LineB in BlastResultsB:
	if len(LineB) > 0:
		LineSplitB = LineB.split("\t")
		#this was changed from LineSplitA[5] to LineSplitA[10]
		###############################################
		if float(LineSplitB[10]) < Evalue:
			ResultsB[LineSplitB[0]] = "\t".join(LineSplitB)
	
#print ResultsB


#Go through lines of first blast search, retrieve the reciprocal BLAST
#Then try and compare if the subject is the same as the original query

#Start a counter
ReciprocalHitCounter = 0

#Setup Final Output

FinalOutput = ""

for LineA in BlastResultsA:
	LineSplitA = LineA.split('\t')
	if len(LineA) > 0 and LineSplitA[2] in ResultsB.keys():
		ComparisonLineSplit = (ResultsB[(LineSplitA[2])]).split('\t')	
		if LineSplitA[0] == ComparisonLineSplit[2]:
			#print "Reciprocal hit found!\n"
			ReciprocalHitCounter += 1
			#Add hits to reciprocal output
			FinalOutput = FinalOutput + LineA + "\t" + str("\t".join(ComparisonLineSplit[4:])) + "\n"
			
		
			
#print "Program is almost finished.\n\n\
			
FinalOutFileName = ProjectName + "_ReciprocalHits.txt"

if WriteOutFile == True:
	FinalOutFile = open(FinalOutFileName, 'w')
	FinalOutFile.write(FinalOutput)

InFileB.close()
FinalOutFile.close()

print "\nProgram is finished\n"
print "The results are in: " + str(FinalOutFileName) + "\n"
print "The number of reciprocal hits was: " + str(ReciprocalHitCounter) + "\n"

#Move any intermediate files created into the project folder
#MakeDirectory first
MakeDirectoy = "mkdir " + ProjectName + "_intermediatefiles"

subprocess.call(MakeDirectoy, shell = True)

#Move generated files into the newly created folder

Move1 = "mv " + ResultsTableStep1 + " " + ProjectName + "_intermediatefiles"
Move2 = "mv " + ResultsTableStep2 + " " + ProjectName + "_intermediatefiles"
Move3 = "mv " + SeqsForBlast2 + " " + ProjectName + "_intermediatefiles"

subprocess.call(Move1, shell = True)
subprocess.call(Move2, shell = True)
subprocess.call(Move3, shell = True)

print "Intermediate files can be found in this folder: " + ProjectName + "_intermediatefiles\n"

if WarningMessages > 0:
	print "Check for warning messages above.\n" +\
	"There was " + str(WarningMessages) + " warning messages\n"
