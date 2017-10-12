###### fastqc
#### trimmmomatic
### trinity
## transdecoder

############################# ORTHOLOG DETECTION:  reciprocal best hits####################

#creating database for the OXPHOS genes
makeblastdb -in /home/dnarules/Documents/Sam_MarineFishProject/Zebrafish_OXPHOSgenes.fasta -dbtype 'prot' -out /home/dnarules/Documents/Sam_MarineFishProject/Zebrafish_OXPHOSgenes -parse_seqids
# creating the databases for the transcriptomes
echo " creating databases"
cd /home/dnarules/Documents/Sam_MarineFishProject/TransDecoder_outputs/
for transcriptome in *.fasta
do
makeblastdb -in $transcriptome -dbtype 'nucl' -out /home/dnarules/Documents/Sam_MarineFishProject/TransDecoder_outputs/$transcriptome -parse_seqids
echo "database created"
echo "moving onto next transcriptome"
done

#tabulating the reciprocal best hits for all the transcriptomes
mkdir /home/dnarules/Documents/Sam_MarineFishProject/reciprocalBLAST
cd /home/dnarules/Documents/Sam_MarineFishProject/reciprocalBLAST
for transcriptomePath in /home/dnarules/Documents/Sam_MarineFishProject/TransDecoder_outputs/*.fasta
do
transcriptome=`basename $transcriptomePath`
ProjectName="RBH_$transcriptome"
SeqsFastaA="/home/dnarules/Documents/Sam_MarineFishProject/Zebrafish_OXPHOSgenes.fasta"
SeqsFastaB="/home/dnarules/Documents/Sam_MarineFishProject/TransDecoder_outputs/$transcriptome"
BlastDbA="/home/dnarules/Documents/Sam_MarineFishProject/Zebrafish_OXPHOSgenes"
BlastDbB="/home/dnarules/Documents/Sam_MarineFishProject/TransDecoder_outputs/$transcriptome"
python /home/dnarules/Documents/Sam_MarineFishProject/Myreciprocalblast_allsteps.py $ProjectName $SeqsFastaA $SeqsFastaB $BlastDbA $BlastDbB 
echo "completing" $transcriptome 
done

#extract the names of the sequences and only the best hit portion that had the RBH


for transcriptomePath in /home/dnarules/Documents/Sam_MarineFishProject/TransDecoder_outputs/*.fasta
do
Transcriptome=`basename $transcriptomePath`
echo $Transcriptome
infile="/home/dnarules/Documents/Sam_MarineFishProject/reciprocalBLAST/RBH_$Transcriptome""_ReciprocalHits.txt"
exec 4<$infile
while read -u4 row 
do
ProteinSeqID=$(echo $row | awk '{ print $1}')  
FishSeqID=$(echo $row | awk '{ print $3}')
SStart=$(echo $row | awk '{ print $17}')  
SEnd=$(echo $row | awk '{ print $18}')
outfile="/home/dnarules/Documents/Sam_MarineFishProject/reciprocalBLAST/Blasted_$ProteinSeqID$Transcriptome"
echo $ProteinSeqID
echo $Transcriptome
if [ $SStart -gt $SEnd ] 
then 
strand="minus"
echo "extracting the reversed hit"
blastdbcmd -db /home/dnarules/Documents/Sam_MarineFishProject/TransDecoder_outputs/$Transcriptome -dbtype nucl -entry $FishSeqID -range $SEnd-$SStart -strand $strand -outfmt %f -out $outfile
fi
if [ $SEnd -gt $SStart ] 
then 
strand="plus"
echo "extracting the reversed hit"
blastdbcmd -db /home/dnarules/Documents/Sam_MarineFishProject/TransDecoder_outputs/$Transcriptome -dbtype nucl -entry $FishSeqID -range $SStart-$SEnd -strand $strand -outfmt %f -out $outfile
fi
#there is no .fasta at end because it already has a .fasta in name
#blastdbcmd -db /media/dnarules/BioSurvey/dnarules/Anuj/DatabaseForGenomesAllFromCaseyUniformHeaderWithLCLWithFish/$Transcriptome -dbtype nucl -entry $FishSeqID -range $SStart-$SEnd -strand $strand -outfmt %f -out $outfile
echo "rename the header of the fasta file to something that has only animal"
Rscript -e 'library(seqinr);arg=commandArgs(TRUE);sequence=read.fasta(arg[1],seqonly=TRUE); write.fasta(sequence,c(arg[2]),file.out=arg[1])' $outfile $Transcriptome
done
done
cd /home/dnarules/Documents/Sam_MarineFishProject/reciprocalBLAST/
find . -name "*" -size 0k |xargs rm

###############aligning the blast results#############################################

# go into geneious and split the big fasta file into multiple fasta files use the batch export function
mkdir /home/dnarules/Documents/Sam_MarineFishProject/separateOXPHOSproteins
perl /home/dnarules/Documents/Sam_MarineFishProject/split_multifasta.pl --input_file=/home/dnarules/Documents/Sam_MarineFishProject/Zebrafish_OXPHOSgenes.fasta --output_dir=/home/dnarules/Documents/Sam_MarineFishProject/separateOXPHOSproteins

mkdir /home/dnarules/Documents/Sam_MarineFishProject/separateOXPHOSgenes
perl /home/dnarules/Documents/Sam_MarineFishProject/split_multifasta.pl --input_file=/home/dnarules/Documents/Sam_MarineFishProject/Danio_rerio_OXPHOSgenes.fasta --output_dir=/home/dnarules/Documents/Sam_MarineFishProject/separateOXPHOSgenes


mkdir /home/dnarules/Documents/Sam_MarineFishProject/aligned_genes
for proteinPath in /home/dnarules/Documents/Sam_MarineFishProject/separateOXPHOSgenes/*.fasta 
do
for transcriptomePath in /home/dnarules/Documents/Sam_MarineFishProject/TransDecoder_outputs/*.fasta
do
Protein=`basename $proteinPath .fasta`
Transcriptome=`basename $transcriptomePath`
echo " concantenating and aligning the DNA sequences of seperate species into one file "
echo $Protein
echo $Transcriptome
cat /home/dnarules/Documents/Sam_MarineFishProject/reciprocalBLAST/Blasted_$Protein$Transcriptome >>  /home/dnarules/Documents/Sam_MarineFishProject/aligned_genes/Unaligned_$Protein.fasta
done
#clustalo -i /media/dnarules/BioSurvey/dnarules/Anuj/Aligned27/Unaligned$Protein.fasta -o /media/dnarules/BioSurvey/dnarules/Anuj/Aligned27/Aligned$Protein.fasta --force
cat /home/dnarules/Documents/Sam_MarineFishProject/separateOXPHOSgenes/$Protein.fasta >>  /home/dnarules/Documents/Sam_MarineFishProject/aligned_genes/Unaligned_$Protein.fasta
mafft --auto /home/dnarules/Documents/Sam_MarineFishProject/aligned_genes/Unaligned_$Protein.fasta > /home/dnarules/Documents/Sam_MarineFishProject/aligned_genes/Aligned_$Protein.fasta
done
#removing empty alignments
cd /home/dnarules/Documents/Sam_MarineFishProject/aligned_genes
find . -name "*" -size 0k |xargs rm

######################################################333



#####import the alignments into geneious or any other gene visualization software to check to make sure that the alignments are good
### if the alignments are not good, trim some of them manually
###RUN CODEML








