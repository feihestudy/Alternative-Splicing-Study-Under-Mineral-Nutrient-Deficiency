sampledir=/home/project/
ASEexprs=/home/project/DASG/ASE/PiFeZnCuMn_in_ex.xls
PSIexprs=/home/project/DASG/ASE/PiFeZnCuMn_psi.xls
sampledir=/home/project/
cdsgtf=/home/project/merged_filter_multiexon_junction40_TPM_CDS.gtf
exongtf=/home/project/merged_filter_multiexon_junction40_TPM.gtf
fpkmfile=/home/project/DEG/FeZnCuMnPiRootShoot/PiFeZnCuMn_transcript_sample_TPM_filter.xls
fafile=/home/project/merged_filter_multiexon_junction40_TPM.fa
ASEdir2=/home/project/DASG/ASE/
ASE="A3SS&G9853&chr02&-&1194628&1195251&1194628&1194741&1195962&1196094;D21;D21piPlusR1_SRR1005290,D21piPlusR2_SRR1005291,D21piPlusR3_SRR1005292,D21piMinusR1_SRR1005293,D21piMinusR2_SRR1005294,D21piMinusR3_SRR1005295"
arr2=(${ASE//;/ }) 
time=${arr2[1]}
ASE=${arr2[0]}
samplebams=${arr2[2]}
rm -rf $ASE"_"$time
mkdir $ASE"_"$time
cd $ASE"_"$time
perl /home/pcd/filter_ASE_for_as_visual.pl -ASE $ASE
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
arr2=(${samplebams//,/ }) 
unset gnames
gnames=`echo ${arr2[0]} | awk '{split($0,a,"_" ); print a[1]}'`
gnames=$gnames"_R"
k=0
for i in ${arr2[@]:1}  
do  
   label=`echo $i | awk '{split($0,a,"_" ); print a[1]}'`
   gnames=$gnames","$label"_R"
done
echo $gnames
perl /home/pcd/filter_ASE_for_as_visual_nx.pl -ASE $ASE -inex $ASEexprs -PSI $PSIexprs -samplelist $gnames
perl /home/pcd/merge_junction.pl -input1 ${arr3[0]}"_junction.txt" -input2 ${arr3[1]}"_junction.txt" -input3 ${arr3[2]}"_junction.txt" > "Control_junction.txt"
perl /home/pcd/merge_junction.pl -input1 ${arr3[3]}"_junction.txt" -input2 ${arr3[4]}"_junction.txt" -input3 ${arr3[5]}"_junction.txt" > "Pi_junction.txt"
rm -f ${arr3[0]}"_junction.txt" ${arr3[1]}"_junction.txt" ${arr3[2]}"_junction.txt"
rm -f ${arr3[3]}"_junction.txt" ${arr3[4]}"_junction.txt" ${arr3[5]}"_junction.txt"
ASEtype=`echo $ASE | awk '{split($0,a,"&" ); print a[1]}'`
perl $pcd/filter_major_tid_by_ASE.pl -ASE $ASE -dir $ASEdir2 > ASE2transcript_major.txt
incltid=`cut -f 2 ASE2transcript_major.txt`
perl /home/pcd/filter_fa_by_tid.pl -tid $incltid -trans $fafile > $incltid"_mRNA.fa"
perl /home/pcd/filter_fa_by_tid_len.pl -tid $incltid -trans $fafile
perl /home/pcd/filter_gtf_by_tid.pl -targettid $incltid -allgtf $cdsgtf > $incltid"_CDS.gtf"
perl /home/pcd/filter_gtf_by_tid.pl -targettid $incltid -allgtf $exongtf > $incltid"_exon.gtf"
excltid=`cut -f 3 ASE2transcript_major.txt`
perl /home/pcd/filter_fa_by_tid.pl -tid $excltid -trans $fafile > $excltid"_mRNA.fa"
perl /home/pcd/filter_gtf_by_tid.pl -targettid $excltid -allgtf $cdsgtf > $excltid"_CDS.gtf"
perl /home/pcd/filter_gtf_by_tid.pl -targettid $excltid -allgtf $exongtf > $excltid"_exon.gtf"
start_codon=`perl /home/pcd/get_CDS_len_gtf.pl -inGTF $incltid"_CDS.gtf" -exGTF $excltid"_CDS.gtf"`
trans_start=`perl /home/pcd/genomepos_mRNApos.pl -targetpos $start_codon -GTF $incltid"_exon.gtf" -mRNA $incltid"_mRNA.fa"`
perl /home/pcd/translate_mRNA.pl -fa $incltid"_mRNA.fa" -pos $trans_start -gtf $incltid"_exon.gtf" > include.gtf
trans_start=`perl /home/pcd/genomepos_mRNApos.pl -targetpos $start_codon -GTF $excltid"_exon.gtf" -mRNA $excltid"_mRNA.fa"`
perl /home/pcd/translate_mRNA.pl -fa $excltid"_mRNA.fa" -pos $trans_start -gtf $excltid"_exon.gtf" > exclude.gtf
perl /home/pcd/filter_SR_domain_gtf.pl -ASE $ASE -SRdomain $SRdomain > Domain.gtf
Rscript /home/pcd/plot_gene_exon_structure_DASE_SR02.r $ASE /hone/pcd/ "+Pi,-Pi"


