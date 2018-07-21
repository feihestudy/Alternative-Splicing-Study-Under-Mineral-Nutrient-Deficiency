sampledir=/home/project/
ASEexprs=/home/project/DASG/ASE/PiFeZnCuMn_in_ex.xls
PSIexprs=/home/project/DASG/ASE/PiFeZnCuMn_psi.xls
cdsgtf=/home/project/merged_filter_multiexon_junction40_TPM_CDS.gtf
exongtf=/home/project/merged_filter_multiexon_junction40_TPM.gtf
fpkmfile=/home/project/DEG/FeZnCuMnPiRootShoot/PiFeZnCuMn_transcript_sample_TPM_filter.xls
SRdomain=/home/database/LOC_Os03g06520.1_domain.gtf
fafile=/home/project/merged_filter_multiexon_junction40_TPM.fa
genebed=/home/project/merged_filter_multiexon_junction40_TPM.bed
Gene="G18448;Fe;ControlR1,ControlR2,ControlR3,FeR1,FeR2,FeR3"
ASEs="RI&G18448&chr03&-&3273065&3274928&3273065&3273179&3274862&3274928"
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
/home/soft/samtools merge -@ 16 -O BAM $time".bam" ${arr3[3]}".bam" ${arr3[4]}".bam" ${arr3[5]}".bam"
rm -f ${arr3[3]}".bam" ${arr3[4]}".bam" ${arr3[5]}".bam"
/home/soft/samtools sort -o "Control_sort.bam" "Control.bam"
/home/soft/samtools sort -o $time"_sort.bam" $time".bam"
mv Control_sort.bam Control.bam
mv $time"_sort.bam" $time".bam"
/home/soft/samtools depth "Control.bam" > "Control.txt"
/home/soft/samtools depth $time".bam" > $time".txt"
perl /home/pcd/filter_ASE_for_as_visual_tt.pl -ASE $ASEs -inex $ASEexprs -PSI $PSIexprs -samplelist "ControlR1,ControlR2,ControlR3" > "Control_junction.txt"
perl /home/pcd/filter_ASE_for_as_visual_tt.pl -ASE $ASEs -inex $ASEexprs -PSI $PSIexprs -samplelist "FeR1,FeR2,FeR3" > $time"_junction.txt"
incltid="TU90424"
perl /home/pcd/filter_fa_by_tid.pl -tid $incltid -trans $fafile > $incltid"_mRNA.fa"
perl /home/pcd/filter_fa_by_tid_len.pl -tid $incltid -trans $fafile
perl /home/pcd/filter_gtf_by_tid.pl.pl -targettid $incltid -allgtf $cdsgtf > $incltid"_CDS.gtf"
perl /home/pcd/filter_gtf_by_tid.pl.pl -targettid $incltid -allgtf $exongtf > $incltid"_exon.gtf"
excltid="LOC_Os03g06520.1"
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
Rscript /home/pcd/plot_gene_exon_structure_whole_gene_S_transport.r "one.bed" /home/pcd/ "Control,Fe"
