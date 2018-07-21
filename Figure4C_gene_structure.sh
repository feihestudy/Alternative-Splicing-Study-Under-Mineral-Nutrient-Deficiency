sampledir=/home/project/
ASEexprs=/home/project/DASG/ASE/PiFeZnCuMn_in_ex.xls
PSIexprs=/home/project/DASG/ASE/PiFeZnCuMn_psi.xls
cdsgtf=/home/project/merged_filter_multiexon_junction40_TPM_CDS.gtf
exongtf=/home/project/merged_filter_multiexon_junction40_TPM.gtf
fpkmfile=/home/project/DEG/FeZnCuMnPiRootShoot/PiFeZnCuMn_transcript_sample_TPM_filter.xls
SRdomain=/home/database/LOC_Os02g55910.1_domain.gtf
fafile=/home/project/merged_filter_multiexon_junction40_TPM.fa
genebed=/home/project/merged_filter_multiexon_junction40_TPM.bed
Gene="G17200;D21;D21piPlusR1_SRR1005290,D21piPlusR2_SRR1005291,D21piPlusR3_SRR1005292,D21piMinusR1_SRR1005293,D21piMinusR2_SRR1005294,D21piMinusR3_SRR1005295"
ASEs="RI&G17200&chr02&+&34223856&34225167&34223856&34224029&34225005&34225167"
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
perl /home/pcd/filter_ASE_for_as_visual_tt.pl -ASE $ASEs -inex $ASEexprs -PSI $PSIexprs -samplelist "D21piPlusR1_R,D21piPlusR2_R,D21piPlusR3_R" > "Control_junction.txt"
perl /home/pcd/filter_ASE_for_as_visual_tt.pl -ASE $ASEs -inex $ASEexprs -PSI $PSIexprs -samplelist "D21piMinusR1_R,D21piMinusR2_R,D21piMinusR3_R" > "Pi_junction.txt"
incltid="TU82879"
perl /home/pcd/filter_fa_by_tid.pl -tid $incltid -trans $fafile > $incltid"_mRNA.fa"
perl /home/pcd/filter_fa_by_tid_len.pl -tid $incltid -trans $fafile
perl /home/pcd/filter_gtf_by_tid.pl.pl -targettid $incltid -allgtf $cdsgtf > $incltid"_CDS.gtf"
perl /home/pcd/filter_gtf_by_tid.pl.pl -targettid $incltid -allgtf $exongtf > $incltid"_exon.gtf"
excltid="LOC_Os02g55910.1"
perl /home/pcd/filter_fa_by_tid.pl -tid $excltid -trans $fafile > $excltid"_mRNA.fa"
perl /home/pcd/filter_gtf_by_tid.pl.pl -targettid $excltid -allgtf $cdsgtf > $excltid"_CDS.gtf"
perl /home/pcd/filter_gtf_by_tid.pl.pl -targettid $excltid -allgtf $exongtf > $excltid"_exon.gtf"
start_codon=`perl /home/pcd/get_CDS_len_gtf.pl -inGTF $incltid"_CDS.gtf" -exGTF $excltid"_CDS.gtf"`
trans_start=`perl /home/pcd/genomepos_mRNApos.pl -targetpos $start_codon -GTF $incltid"_exon.gtf" -mRNA $incltid"_mRNA.fa"`
perl /home/pcd/translate_mRNA.pl -fa $incltid"_mRNA.fa" -pos $trans_start -gtf $incltid"_exon.gtf" > include.gtf
trans_start=`perl /home/pcd/genomepos_mRNApos.pl -targetpos $start_codon -GTF $excltid"_exon.gtf" -mRNA $excltid"_mRNA.fa"`
perl /home/pcd/translate_mRNA.pl -fa $excltid"_mRNA.fa" -pos $trans_start -gtf $excltid"_exon.gtf" > exclude.gtf
perl /home/pcd/filter_SR_domain_gtf_by_gene.pl -gene $Gene -SRdomain $SRdomain > Domain.gtf
head -10 Domain.gtf
Rscript /home/pcd/plot_gene_exon_structure_AS_RI_run_pi_2sample_wholeG_MGD.r "one.bed" /home/pcd/ "Control,Pi"
