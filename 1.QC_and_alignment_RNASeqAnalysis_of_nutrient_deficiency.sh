#Download trim_galore version 0.4.2 software https://www.bioinformatics.babraham.ac.uk/projects/trim_galore/
#Download STAR version 2.5.2b https://github.com/alexdobin/STAR 
#Download Samtools http://samtools.sourceforge.net/
#Download RSEQC http://rseqc.sourceforge.net/
#Download sratoolkit.2.8.0-centos_linux64 https://www.ncbi.nlm.nih.gov/sra/docs/toolkitsoft/
#install softwares to a directory (for example /home/soft/)
#copy database to a directory (for example /home/database/)
#make project directory (for example /home/project/)


#Analyze RNASeq DataSet of  (QC trim and STAR alignment)

#All the 15 RNA-Seq data of Micronutrient deficiency  are available under the accession number SRP117202. https://www.ncbi.nlm.nih.gov/sra/?term=SRP117202


Sampleids=("ControlR1" "ControlR2" "ControlR3" "FeR1" "FeR2" "FeR3" "ZnR1" "ZnR2" "ZnR3" "CuR1" "CuR2" "CuR3" "MnR1" "MnR2" "MnR3")
Adaptor1s=("AGATCGGAAGAGCACACGTCTGAACTCCAGTCACGTCCGCACATCTCGTATGCCGTCTTCTGCTTG" "AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC GTGAAA ACATCTCGTATGCCGTCTTCTGCTTG" "AGATCGGAAGAGCACACGTCTGAACTCCAGTCACGTGGCCACATCTCGTATGCCGTCTTCTGCTTG" "AGATCGGAAGAGCACACGTCTGAACTCCAGTCACGTTTCGACATCTCGTATGCCGTCTTCTGCTTG" "AGATCGGAAGAGCACACGTCTGAACTCCAGTCACCGTACGACATCTCGTATGCCGTCTTCTGCTTG" "AGATCGGAAGAGCACACGTCTGAACTCCAGTCACGAGTGGACATCTCGTATGCCGTCTTCTGCTTG" "AGATCGGAAGAGCACACGTCTGAACTCCAGTCACATTCCTACATCTCGTATGCCGTCTTCTGCTTG" "AGATCGGAAGAGCACACGTCTGAACTCCAGTCACCAAAAGACATCTCGTATGCCGTCTTCTGCTTG" "AGATCGGAAGAGCACACGTCTGAACTCCAGTCACCGGATTACATCTCGTATGCCGTCTTCTGCTTG" "AGATCGGAAGAGCACACGTCTGAACTCCAGTCACCAACTAAGACATCTCGTATGCCGTCTTCTGCTTG" "AGATCGGAAGAGCACACGTCTGAACTCCAGTCACCACGATAGACATCTCGTATGCCGTCTTCTGCTTG" "AGATCGGAAGAGCACACGTCTGAACTCCAGTCACCACTCAATACATCTCGTATGCCGTCTTCTGCTTG" "AGATCGGAAGAGCACACGTCTGAACTCCAGTCACGGTAGCACATCTCGTATGCCGTCTTCTGCTTG" "AGATCGGAAGAGCACACGTCTGAACTCCAGTCACACTGATACATCTCGTATGCCGTCTTCTGCTTG" "AGATCGGAAGAGCACACGTCTGAACTCCAGTCACATTAGCACATCTCGTATGCCGTCTTCTGCTTG")
Adaptor2="AATGATACGGCGACCACCGAGATCTACACTCTTTCCCTACACGACGCTCTTCCGATCT"
for((i=0; i<=14; i++));
do
cd ${Sampleids[i]}
perl /home/soft/trim_galore -q 25 --stringency 5 --dont_gzip --fastqc --retain_unpaired -r1 31 -r2 31 --length 30 -o ./ --paired --phred33 -a $Adaptor1 -a2 $Adaptor2 ${Sampleids[i]}"_R1.fastq.gz" ${Sampleids[i]}"_R2.fastq.gz"
read1=*_val_1.fq
read2=*_val_2.fq
/home/soft/STAR --runThreadN 24 --alignMatesGapMax 100000 --alignIntronMin 35 --alignIntronMax 2000 --outFilterMismatchNmax 8 --outFilterMultimapNmax 20 --alignSJoverhangMin 8 --alignSJDBoverhangMin 8 --outSAMstrandField intronMotif --outSAMtype SAM --genomeDir /home/database/STARgenomedir --sjdbGTFfile /home/database/RAPDBMSU.gtf --readFilesIn $read1 $read2
/home/soft/samtools view -@ 16 -bS -o Aligned.out.bam Aligned.out.sam
/home/soft/samtools sort -@ 16 -o Aligned.out.sort.bam Aligned.out.bam
mv Aligned.out.sam accepted_hits.sam
rm -f Aligned.out.bam
rm -rf _STARtmp
rm -rf _STARgenome
rm -f Log*
mv Aligned.out.sort.bam accepted_hits.bam
python /home/soft/RSEQC/bam_stat.py -i accepted_hits.bam | cat > "RSEQC_alignment.stat"
read1=*_val_1.fq
read2=*_val_2.fq
trimlen=75
perl /home/pcd/trim_arg.pl -readlen $trimlen -fq1 $read1 -fq2 $read2
read1=$read1"."$trimlen
read2=$read2"."$trimlen
/home/soft/STAR --chimSegmentMin 2 --outFilterMismatchNmax 3 --alignEndsType EndToEnd --runThreadN 24 --outSAMstrandField intronMotif --outSAMtype SAM --alignSJDBoverhangMin 8 --alignIntronMax 300000 --genomeDir /home/database/STARgenomedir --sjdbGTFfile /home/database/RAPDBMSU.gtf --outFileNamePrefix rMATS --readFilesIn $read1 $read2
/home/soft/samtools view -@ 16 -bS -o rMATSAligned.out.bam rMATSAligned.out.sam
/home/soft/samtools sort -@ 16 -o rMATSAligned.out.sort.bam rMATSAligned.out.bam
mv rMATSAligned.out.sort.bam "accepted_hits"$trimlen".bam"
/home/soft/samtools index accepted_hits75.bam
cd ..
done




#Analyze RNASeq DataSet of Pi deficiency (QC trim and STAR alignment)

#All the 126 RNA-Seq data of the time course of Pi-starved root and shoot samples are available under the accession number SRA097415 (Secco et al., 2013). https://www.ncbi.nlm.nih.gov/sra/?term=SRA097415

RootSamples=("H1piPlusR1_SRR1005257" "H1piPlusR2_SRR1005258" "H1piPlusR3_SRR1005259" "H1piMinusR1_SRR1005260" "H1piMinusR2_SRR1005261" "H1piMinusR3_SRR1005262" "H6piPlusR1_SRR1005308" "H6piPlusR2_SRR1005309" "H6piPlusR3_SRR1005310" "H6piMinusR1_SRR1005311" "H6piMinusR2_SRR1005312" "H6piMinusR3_SRR1005313" "H24piPlusR1_SRR1005296" "H24piPlusR2_SRR1005297" "H24piPlusR3_SRR1005298" "H24piMinusR1_SRR1005299" "H24piMinusR2_SRR1005300" "H24piMinusR3_SRR1005301" "D3piPlusR1_SRR1005302" "D3piPlusR2_SRR1005303" "D3piPlusR3_SRR1005304" "D3piMinusR1_SRR1005305" "D3piMinusR2_SRR1005306" "D3piMinusR3_SRR1005307" "D7piPlusR1_SRR1005314" "D7piPlusR2_SRR1005315" "D7piPlusR3_SRR1005316" "D7piMinusR1_SRR1005317" "D7piMinusR2_SRR1005318" "D7piMinusR3_SRR1005319" "D21piPlusR1_SRR1005290" "D21piPlusR2_SRR1005291" "D21piPlusR3_SRR1005292" "D21piMinusR1_SRR1005293" "D21piMinusR2_SRR1005294" "D21piMinusR3_SRR1005295" "D21H1piPlusR1_SRR1005263" "D21H1piPlusR2_SRR1005264" "D21H1piPlusR3_SRR1005265" "D21H1piMinusR1_SRR1005266" "D21H1piMinusR2_SRR1005267" "D21H1piMinusR3_SRR1005268" "D21H1piPlusRecR1_SRR1005269" "D21H1piPlusRecR2_SRR1005270" "D21H1piPlusRecR3_SRR1005271" "D21H6piPlusR1_SRR1005281" "D21H6piPlusR2_SRR1005282" "D21H6piPlusR3_SRR1005283" "D21H6piMinusR1_SRR1005284" "D21H6piMinusR2_SRR1005285" "D21H6piMinusR3_SRR1005286" "D21H6piPlusRecR1_SRR1005287" "D21H6piPlusRecR2_SRR1005288" "D21H6piPlusRecR3_SRR1005289" "D21H24piPlusR1_SRR1005272" "D21H24piPlusR2_SRR1005273" "D21H24piPlusR3_SRR1005274" "D21H24piMinusR1_SRR1005275" "D21H24piMinusR2_SRR1005276" "D21H24piMinusR3_SRR1005277" "D21H24piPlusRecR1_SRR1005278" "D21H24piPlusRecR2_SRR1005279" "D21H24piPlusRecR3_SRR1005280")
ShootSamples=("H1piPlusR1_SRR1005320" "H1piPlusR2_SRR1005321" "H1piPlusR3_SRR1005322" "H1piMinusR1_SRR1005323" "H1piMinusR2_SRR1005324" "H1piMinusR3_SRR1005325" "H6piPlusR1_SRR1005371" "H6piPlusR2_SRR1005372" "H6piPlusR3_SRR1005373" "H6piMinusR1_SRR1005374" "H6piMinusR2_SRR1005375" "H6piMinusR3_SRR1005376" "H24piPlusR1_SRR1005359" "H24piPlusR2_SRR1005360" "H24piPlusR3_SRR1005361" "H24piMinusR1_SRR1005362" "H24piMinusR2_SRR1005363" "H24piMinusR3_SRR1005364" "D3piPlusR1_SRR1005365" "D3piPlusR2_SRR1005366" "D3piPlusR3_SRR1005367" "D3piMinusR1_SRR1005368" "D3piMinusR2_SRR1005369" "D3piMinusR3_SRR1005370" "D7piPlusR1_SRR1005377" "D7piPlusR2_SRR1005378" "D7piPlusR3_SRR1005379" "D7piMinusR1_SRR1005380" "D7piMinusR2_SRR1005381" "D7piMinusR3_SRR1005382" "D21piPlusR1_SRR1005353" "D21piPlusR2_SRR1005354" "D21piPlusR3_SRR1005355" "D21piMinusR1_SRR1005356" "D21piMinusR2_SRR1005357" "D21piMinusR3_SRR1005358" "D21H1piPlusR1_SRR1005326" "D21H1piPlusR2_SRR1005327" "D21H1piPlusR3_SRR1005328" "D21H1piMinusR1_SRR1005329" "D21H1piMinusR2_SRR1005330" "D21H1piMinusR3_SRR1005331" "D21H1piPlusRecR1_SRR1005332" "D21H1piPlusRecR2_SRR1005333" "D21H1piPlusRecR3_SRR1005334" "D21H6piPlusR1_SRR1005344" "D21H6piPlusR2_SRR1005345" "D21H6piPlusR3_SRR1005346" "D21H6piMinusR1_SRR1005347" "D21H6piMinusR2_SRR1005348" "D21H6piMinusR3_SRR1005349" "D21H6piPlusRecR1_SRR1005350" "D21H6piPlusRecR2_SRR1005351" "D21H6piPlusRecR3_SRR1005352" "D21H24piPlusR1_SRR1005335" "D21H24piPlusR2_SRR1005336" "D21H24piPlusR3_SRR1005337" "D21H24piMinusR1_SRR1005338" "D21H24piMinusR2_SRR1005339" "D21H24piMinusR3_SRR1005340" "D21H24piPlusRecR1_SRR1005341" "D21H24piPlusRecR2_SRR1005342" "D21H24piPlusRecR3_SRR1005343")
RootShootSamples=(${RootSamples[@]} ${ShootSamples[@]})
for((i=0; i<=125; i++));
do
#echo ${RootShootSamples[$i]}
cd ${RootShootSamples[$i]}
sname=`ls *.sra`
sname=`echo ${sname} | cut -d '.' -f 1`
/home/soft/sratoolkit.2.8.0-centos_linux64/fastq-dump --split-files $sname".sra"
perl /home/soft/trim_galore -q 25 --stringency 5 -o ./ --paired $sname"_1.fastq" $sname"_2.fastq"
read1=$sname"_1_val_1.fq"
read2=$sname"_2_val_2.fq"
/home/soft/STAR --runThreadN 24 --alignMatesGapMax 100000 --alignIntronMin 35 --alignIntronMax 2000 --outFilterMismatchNmax 8 --outFilterMultimapNmax 20 --alignSJoverhangMin 8 --alignSJDBoverhangMin 8 --outSAMstrandField intronMotif --outSAMtype SAM --genomeDir /home/database/STARgenomedir --sjdbGTFfile /home/database/RAPDBMSU.gtf --readFilesIn $read1 $read2
/home/soft/samtools view -@ 16 -bS -o Aligned.out.bam Aligned.out.sam
/home/soft/samtools sort -@ 16 -o Aligned.out.sort.bam Aligned.out.bam
mv Aligned.out.sam accepted_hits.sam
rm -f Aligned.out.bam
rm -rf _STARtmp
rm -rf _STARgenome
rm -f Log*
mv Aligned.out.sort.bam accepted_hits.bam
python /home/soft/RSEQC/bam_stat.py -i accepted_hits.bam | cat > "RSEQC_alignment.stat"
read1=*_val_1.fq
read2=*_val_2.fq
trimlen=75
perl /home/pcd/trim_arg.pl -readlen $trimlen -fq1 $read1 -fq2 $read2
read1=$read1"."$trimlen
read2=$read2"."$trimlen
/home/soft/STAR --chimSegmentMin 2 --outFilterMismatchNmax 3 --alignEndsType EndToEnd --runThreadN 24 --outSAMstrandField intronMotif --outSAMtype SAM --alignSJDBoverhangMin 8 --alignIntronMax 300000 --genomeDir /home/database/STARgenomedir --sjdbGTFfile /home/database/RAPDBMSU.gtf --outFileNamePrefix rMATS --readFilesIn $read1 $read2
/home/soft/samtools view -@ 16 -bS -o rMATSAligned.out.bam rMATSAligned.out.sam
/home/soft/samtools sort -@ 16 -o rMATSAligned.out.sort.bam rMATSAligned.out.bam
mv rMATSAligned.out.sort.bam "accepted_hits"$trimlen".bam"
/home/soft/samtools index accepted_hits75.bam
cd ..
done
