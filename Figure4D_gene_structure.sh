sampledir=/home/project/
ASEexprs=/home/project/DASG/ASE/PiFeZnCuMn_in_ex.xls
PSIexprs=/home/project/DASG/ASE/PiFeZnCuMn_psi.xls
cdsgtf=/home/project/merged_filter_multiexon_junction40_TPM_CDS.gtf
exongtf=/home/project/merged_filter_multiexon_junction40_TPM.gtf
fpkmfile=/home/project/DEG/FeZnCuMnPiRootShoot/PiFeZnCuMn_transcript_sample_TPM_filter.xls
SRdomain=/home/database/LOC_Os05g48390.1_domain.gtf
fafile=/home/project/merged_filter_multiexon_junction40_TPM.fa
genebed=/home/project/merged_filter_multiexon_junction40_TPM.bed
Gene="G39496;D21;D21piPlusR1_SRR1005353,D21piPlusR2_SRR1005354,D21piPlusR3_SRR1005355,D21piMinusR1_SRR1005356,D21piMinusR2_SRR1005357,D21piMinusR3_SRR1005358"
ASEs="RI&G39496&chr05&-&27740004&27741420&27740004&27740703&27741108&27741420;RI&G39496&chr05&-&27741108&27741657&27741108&27741420&27741517&27741657;RI&G39496&chr05&-&27742673&27743033&27742673&27742817&27742901&27743033;RI&G39496&chr05&-&27742901&27744618&27742901&27743033&27743203&27744618;RI&G39496&chr05&-&27743203&27744868&27743203&27744618&27744737&27744868"
arr2=(${Gene//;/ }) 
time=${arr2[1]}
Gene=${arr2[0]}
samplebams=${arr2[2]}
rm -rf $Gene"_"$time
mkdir $Gene"_"$time
cd $Gene"_"$time
perl /home/pcd/filter_gene_for_as_visual.pl -gene $Gene -bed $genebed > one.bed
arr2=(${samplebams//,/ })
unset arr3
k=0
for i in ${arr2[@]}  
do  
   echo ${arr2[$k]}
   label=`echo $i | awk '{split($0,a,"_" ); print a[1]}'`
   /home/soft/samtools view -L one.bed -o $label".bam" $sampledir"/"${arr2[$k]}"/accepted_hits.bam"
   arr3[$k]=$label   
   k=$k+1;    	
done
arr2=(${samplebams//,/ })
unset arr3
k=0
for i in ${arr2[@]}  
do  
   echo ${arr2[$k]}
   label=`echo $i | awk '{split($0,a,"_" ); print a[1]}'`
   arr3[$k]=$label   
   k=$k+1;    	
done
/home/soft/samtools merge -@ 16 -O BAM "Control.bam" ${arr3[0]}".bam" ${arr3[1]}".bam" ${arr3[2]}".bam"
rm -f ${arr3[0]}".bam" ${arr3[1]}".bam" ${arr3[2]}".bam"
/home/soft/samtools merge -@ 16 -O BAM "Pi.bam" ${arr3[3]}".bam" ${arr3[4]}".bam" ${arr3[5]}".bam"
rm -f ${arr3[3]}".bam" ${arr3[4]}".bam" ${arr3[5]}".bam"
/home/soft/samtools sort -o "Control_sort.bam" "Control.bam"
/home/soft/samtools sort -o "Pi_sort.bam" "Pi.bam"
mv Control_sort.bam Control.bam
mv Pi_sort.bam Pi.bam
/home/soft/samtools depth "Control.bam" > "Control.txt"
/home/soft/samtools depth "Pi.bam" > "Pi.txt"
perl /home/pcd/filter_ASE_for_as_visual_tt.pl -ASE $ASEs -inex $ASEexprs -PSI $PSIexprs -samplelist "D21piPlusR1_S,D21piPlusR2_S,D21piPlusR3_S" > "Control_junction.txt"
perl /home/pcd/filter_ASE_for_as_visual_tt.pl -ASE $ASEs -inex $ASEexprs -PSI $PSIexprs -samplelist "D21piMinusR1_S,D21piMinusR2_S,D21piMinusR3_S" > "Pi_junction.txt"
tid1="TU180856"
perl /home/pcd/filter_fa_by_tid.pl -tid $tid1 -trans $fafile > $tid1"_mRNA.fa"
perl /home/pcd/filter_fa_by_tid_len.pl -tid $tid1 -trans $fafile
perl /home/pcd/filter_gtf_by_tid.pl.pl -targettid $tid1 -allgtf $cdsgtf > $tid1"_CDS.gtf"
perl /home/pcd/filter_gtf_by_tid.pl.pl -targettid $tid1 -allgtf $exongtf > $tid1"_exon.gtf"
tid2="TU180857"
perl /home/pcd/filter_fa_by_tid.pl -tid $tid2 -trans $fafile > $tid2"_mRNA.fa"
perl /home/pcd/filter_fa_by_tid_len.pl -tid $tid2 -trans $fafile
perl /home/pcd/filter_gtf_by_tid.pl.pl -targettid $tid2 -allgtf $cdsgtf > $tid2"_CDS.gtf"
perl /home/pcd/filter_gtf_by_tid.pl.pl -targettid $tid2 -allgtf $exongtf > $tid2"_exon.gtf"
tid3="LOC_Os05g48390.1"
perl /home/pcd/filter_fa_by_tid.pl -tid $tid3 -trans $fafile > $tid3"_mRNA.fa"
perl /home/pcd/filter_fa_by_tid_len.pl -tid $tid3 -trans $fafile
perl /home/pcd/filter_gtf_by_tid.pl.pl -targettid $tid3 -allgtf $cdsgtf > $tid3"_CDS.gtf"
perl /home/pcd/filter_gtf_by_tid.pl.pl -targettid $tid3 -allgtf $exongtf > $tid3"_exon.gtf"



start_codon=`perl /home/pcd/get_CDS_len_gtf.pl -inGTF $tid1"_CDS.gtf" -exGTF $tid3"_CDS.gtf"`
trans_start=`perl /home/pcd/genomepos_mRNApos.pl -targetpos $start_codon -GTF $tid1"_exon.gtf" -mRNA $tid1"_mRNA.fa"`
perl /home/pcd/translate_mRNA.pl -fa $tid1"_mRNA.fa" -pos $trans_start -gtf $tid1"_exon.gtf" > tid1.gtf

start_codon=`perl /home/pcd/get_CDS_len_gtf.pl -inGTF $tid2"_CDS.gtf" -exGTF $tid3"_CDS.gtf"`
trans_start=`perl /home/pcd/genomepos_mRNApos.pl -targetpos $start_codon -GTF $tid2"_exon.gtf" -mRNA $tid2"_mRNA.fa"`
perl /home/pcd/translate_mRNA.pl -fa $tid2"_mRNA.fa" -pos $trans_start -gtf $tid2"_exon.gtf" > tid2.gtf

trans_start=`perl /home/pcd/genomepos_mRNApos.pl -targetpos $start_codon -GTF $tid3"_exon.gtf" -mRNA $tid3"_mRNA.fa"`
perl /home/pcd/translate_mRNA.pl -fa $tid3"_mRNA.fa" -pos $trans_start -gtf $tid3"_exon.gtf" > tid3.gtf



start_codon=`perl /home/pcd/get_CDS_len_gtf.pl -inGTF $incltid"_CDS.gtf" -exGTF $excltid"_CDS.gtf"`
trans_start=`perl /home/pcd/genomepos_mRNApos.pl -targetpos $start_codon -GTF $incltid"_exon.gtf" -mRNA $incltid"_mRNA.fa"`
perl /home/pcd/translate_mRNA.pl -fa $incltid"_mRNA.fa" -pos $trans_start -gtf $incltid"_exon.gtf" > include.gtf
trans_start=`perl /home/pcd/genomepos_mRNApos.pl -targetpos $start_codon -GTF $excltid"_exon.gtf" -mRNA $excltid"_mRNA.fa"`
perl /home/pcd/translate_mRNA.pl -fa $excltid"_mRNA.fa" -pos $trans_start -gtf $excltid"_exon.gtf" > exclude.gtf
perl /home/pcd/filter_SR_domain_gtf_by_gene.pl -gene $Gene -SRdomain $SRdomain > Domain.gtf
head -10 Domain.gtf
Rscript /home/pcd/plot_gene_exon_structure_whole_gene_PHO2.r "one.bed" /home/pcd/ "Control,Pi"
