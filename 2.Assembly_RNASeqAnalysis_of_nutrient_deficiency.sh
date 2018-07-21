#Download Stringtie version 0.4.2 software https://ccb.jhu.edu/software/stringtie/
#Download TACO v0.7.3 https://tacorna.github.io/
#Download gffcompare v0.9.9c https://github.com/gpertea/gffcompare
#Download salmon 0.8.2 https://combine-lab.github.io/salmon/
#Download rMATS 3.2.5 http://rnaseq-mats.sourceforge.net/
#Download CPC2 1.0 http://cpc2.cbi.pku.edu.cn/download.php
#Download TransDecoder 3.0.1 https://github.com/TransDecoder/TransDecoder/releases
#install softwares to a directory (for example /home/soft/)
#copy database to a directory (for example /home/database/)
#make project directory (for example /home/project/)



#Splicing junction profile of multiple samples from STAR alignment result (SJ.out.tab)
samplenames=${AllSamples[0]}
for i in ${AllSamples[@]:1}  
do
   samplenames=$samplenames","${i};
done
perl /home/pcd/junction_profile.pl -dir "./" -samplelist $samplenames > SJ_matrix.txt



#Single-Sample StringTie Assembly)

FeZnCuMnSampleids=("ControlR1" "ControlR2" "ControlR3" "FeR1" "FeR2" "FeR3" "ZnR1" "ZnR2" "ZnR3" "CuR1" "CuR2" "CuR3" "MnR1" "MnR2" "MnR3")
PiRootSamples=("H1piPlusR1_SRR1005257" "H1piPlusR2_SRR1005258" "H1piPlusR3_SRR1005259" "H1piMinusR1_SRR1005260" "H1piMinusR2_SRR1005261" "H1piMinusR3_SRR1005262" "H6piPlusR1_SRR1005308" "H6piPlusR2_SRR1005309" "H6piPlusR3_SRR1005310" "H6piMinusR1_SRR1005311" "H6piMinusR2_SRR1005312" "H6piMinusR3_SRR1005313" "H24piPlusR1_SRR1005296" "H24piPlusR2_SRR1005297" "H24piPlusR3_SRR1005298" "H24piMinusR1_SRR1005299" "H24piMinusR2_SRR1005300" "H24piMinusR3_SRR1005301" "D3piPlusR1_SRR1005302" "D3piPlusR2_SRR1005303" "D3piPlusR3_SRR1005304" "D3piMinusR1_SRR1005305" "D3piMinusR2_SRR1005306" "D3piMinusR3_SRR1005307" "D7piPlusR1_SRR1005314" "D7piPlusR2_SRR1005315" "D7piPlusR3_SRR1005316" "D7piMinusR1_SRR1005317" "D7piMinusR2_SRR1005318" "D7piMinusR3_SRR1005319" "D21piPlusR1_SRR1005290" "D21piPlusR2_SRR1005291" "D21piPlusR3_SRR1005292" "D21piMinusR1_SRR1005293" "D21piMinusR2_SRR1005294" "D21piMinusR3_SRR1005295" "D21H1piPlusR1_SRR1005263" "D21H1piPlusR2_SRR1005264" "D21H1piPlusR3_SRR1005265" "D21H1piMinusR1_SRR1005266" "D21H1piMinusR2_SRR1005267" "D21H1piMinusR3_SRR1005268" "D21H1piPlusRecR1_SRR1005269" "D21H1piPlusRecR2_SRR1005270" "D21H1piPlusRecR3_SRR1005271" "D21H6piPlusR1_SRR1005281" "D21H6piPlusR2_SRR1005282" "D21H6piPlusR3_SRR1005283" "D21H6piMinusR1_SRR1005284" "D21H6piMinusR2_SRR1005285" "D21H6piMinusR3_SRR1005286" "D21H6piPlusRecR1_SRR1005287" "D21H6piPlusRecR2_SRR1005288" "D21H6piPlusRecR3_SRR1005289" "D21H24piPlusR1_SRR1005272" "D21H24piPlusR2_SRR1005273" "D21H24piPlusR3_SRR1005274" "D21H24piMinusR1_SRR1005275" "D21H24piMinusR2_SRR1005276" "D21H24piMinusR3_SRR1005277" "D21H24piPlusRecR1_SRR1005278" "D21H24piPlusRecR2_SRR1005279" "D21H24piPlusRecR3_SRR1005280")
PiShootSamples=("H1piPlusR1_SRR1005320" "H1piPlusR2_SRR1005321" "H1piPlusR3_SRR1005322" "H1piMinusR1_SRR1005323" "H1piMinusR2_SRR1005324" "H1piMinusR3_SRR1005325" "H6piPlusR1_SRR1005371" "H6piPlusR2_SRR1005372" "H6piPlusR3_SRR1005373" "H6piMinusR1_SRR1005374" "H6piMinusR2_SRR1005375" "H6piMinusR3_SRR1005376" "H24piPlusR1_SRR1005359" "H24piPlusR2_SRR1005360" "H24piPlusR3_SRR1005361" "H24piMinusR1_SRR1005362" "H24piMinusR2_SRR1005363" "H24piMinusR3_SRR1005364" "D3piPlusR1_SRR1005365" "D3piPlusR2_SRR1005366" "D3piPlusR3_SRR1005367" "D3piMinusR1_SRR1005368" "D3piMinusR2_SRR1005369" "D3piMinusR3_SRR1005370" "D7piPlusR1_SRR1005377" "D7piPlusR2_SRR1005378" "D7piPlusR3_SRR1005379" "D7piMinusR1_SRR1005380" "D7piMinusR2_SRR1005381" "D7piMinusR3_SRR1005382" "D21piPlusR1_SRR1005353" "D21piPlusR2_SRR1005354" "D21piPlusR3_SRR1005355" "D21piMinusR1_SRR1005356" "D21piMinusR2_SRR1005357" "D21piMinusR3_SRR1005358" "D21H1piPlusR1_SRR1005326" "D21H1piPlusR2_SRR1005327" "D21H1piPlusR3_SRR1005328" "D21H1piMinusR1_SRR1005329" "D21H1piMinusR2_SRR1005330" "D21H1piMinusR3_SRR1005331" "D21H1piPlusRecR1_SRR1005332" "D21H1piPlusRecR2_SRR1005333" "D21H1piPlusRecR3_SRR1005334" "D21H6piPlusR1_SRR1005344" "D21H6piPlusR2_SRR1005345" "D21H6piPlusR3_SRR1005346" "D21H6piMinusR1_SRR1005347" "D21H6piMinusR2_SRR1005348" "D21H6piMinusR3_SRR1005349" "D21H6piPlusRecR1_SRR1005350" "D21H6piPlusRecR2_SRR1005351" "D21H6piPlusRecR3_SRR1005352" "D21H24piPlusR1_SRR1005335" "D21H24piPlusR2_SRR1005336" "D21H24piPlusR3_SRR1005337" "D21H24piMinusR1_SRR1005338" "D21H24piMinusR2_SRR1005339" "D21H24piMinusR3_SRR1005340" "D21H24piPlusRecR1_SRR1005341" "D21H24piPlusRecR2_SRR1005342" "D21H24piPlusRecR3_SRR1005343")
AllSamples=(${FeZnCuMnSampleids[@]} ${PiRootSamples[@]} ${PiShootSamples[@]})
for((i=0; i<=140; i++));
do
echo ${AllSamples[$i]}
cd ${AllSamples[$i]}
/home/soft/stringtie -p 16 -G /home/database/RAPDBMSU.gtf -m 200 -o  "stringtie.gtf" "accepted_hits.bam"
mv "stringtie.gtf" ${arr[$i]}"_stringtie.gtf"
cd ..
done



#Multi-Sample StringTie Assembly)


rm -f gtf.list
for((i=0; i<=140; i++));
do
	echo "/home/project/"${AllSamples[$i]}"/"${AllSamples[$i]}"_stringtie.gtf" >> gtf.list
done

rm -rf TACO
/home/soft/taco-v0.7.3.Linux_x86_64/taco_run -p 8 -o ./TACO gtf.list

cp ./TACO/assembly.gtf ./



#Filter transcript based one length (>=200bp)

perl /home/pcd/filter_gtf_by_transcript_len.pl -GTF assembly.gtf -LEN 200 > assembly_len.gtf
mv assembly_len.gtf merged.gtf


#Compare merged.gtf with RAPDB and MSU GTF file and build GTF file for AS event

/home/soft/gffcompare -r /home/database/RAPDBMSU.gtf merged.gtf
perl /home/pcd/parse_cuffcompare.pl -tmap gffcmp.merged.gtf.tmap
perl /home/pcd/get_transcript_cuffcompare_gtf_fa_osa_equal_taco.pl -tmap gffcmp.merged.gtf.tmap -gtf /home/database/RAPDBMSU.gtf > equal.gtf
perl /home/pcd/get_transcript_cuffcompare_gtf_fa_osa_v2.pl -tmap gffcmp.merged.gtf.tmap -gtf merged.gtf > jou.gtf
perl /home/pcd/filter_singleexon_gtf.pl -gtf jou.gtf > jou_multiexon.gtf # Single-Exon Filter for novel transcript
perl /home/pcd/splice_junction_filter.pl -GTF jou_multiexon.gtf -SJ SJ_matrix.txt -sumcutoff 40 > jou_multiexon_junction40.gtf
perl /home/pcd/get_transcript_cuffcompare_gtf_fa_osa_overlap_taco.pl -tmap gffcmp.merged.gtf.tmap -gtf /home/database/RAPDBMSU.gtf > overlap.gtf
cat jou_multiexon_junction40.gtf equal.gtf overlap.gtf > merged_filter_multiexon_junction40.gtf
perl /home/pcd/get_gid2tid_from_gtf_osa.pl -gtf merged_filter_multiexon_junction40.gtf > gid2tid.txt #transcript_id and gene_id relationship
perl /home/pcd/get_seq_from_gtf.pl -database /home/database/genome.fa -GTF merged_filter_multiexon_junction40.gtf > merged_filter_multiexon_junction40.fa #Extract transcript sequence from GTF using "merged_filter_multiexon_junction40" GTF


#Transcript and gene quantification using Salmon 

/home/soft/Salmon-0.8.2_linux_x86_64/bin/salmon index -t merged_filter_multiexon_junction40.fa -i merged_filter_multiexon_junction40_salmon_index 

mkdir salmon
cd salmon
for((i=0; i<=140; i++));
do
cd "/home/project/"${AllSamples[$i]}"/"
read1=`ls *_val_1.fq`
read2=`ls *_val_2.fq`
cd /home/project/salmon
read1="/home/project/"${AllSamples[$i]}"/"$read1
read2="/home/project/"${AllSamples[$i]}"/"$read2
/home/soft/Salmon-0.8.2_linux_x86_64/bin/salmon quant -p 16 -l "IU" -i ../merged_filter_multiexon_junction40_salmon_index -o ${AllSamples[$i]} -1 $read1 -2 $read2
done;
cd ..

perl /home/pcd/get_sample_transcript_TPM_salmon.pl -dir /home/project/salmon/ -sample $samplenames > PiFeZnCuMn_transcript_sample_TPM.xls
perl /home/pcd/get_sample_gene_TPM_salmon.pl -dir /home/project/salmon/ -sample $samplenames > PiFeZnCuMn_gene_sample_TPM.xls


#Filter transcript using TPM >= 1 at least one sample
perl /home/pcd/filter_gtf_by_TPM.pl -GTF merged_filter_multiexon_junction40.gtf -TPM 1 -N 1 -TPMfile PiFeZnCuMn_transcript_sample_TPM.xls > merged_filter_multiexon_junction40_TPM.gtf


#AS event (RI SE A3SS A5SS MXE) detection from GTF
/home/soft/samtools view -b -H -o blank.bam /home/project/ControlR1/accepted_hits.bam
/home/soft/samtools index blank.bam
python /home/soft/rMATS.3.2.5/bin/processGTF.BAMs.py merged_filter_multiexon_junction40_TPM.gtf ./merged_filter_multiexon_junction40_TPM_AS blank.bam fr-unstranded PE 1
mv merged_filter_multiexon_junction40_TPM_AS.SE.txt AS.SE.txt
mv merged_filter_multiexon_junction40_TPM_AS.RI.txt AS.RI.txt
mv merged_filter_multiexon_junction40_TPM_AS.MXE.txt AS.MXE.txt
mv merged_filter_multiexon_junction40_TPM_AS.A5SS.txt AS.A5SS.txt
mv merged_filter_multiexon_junction40_TPM_AS.A3SS.txt AS.A3SS.txt


#### Gene locus annotation table from GTF

perl /home/pcd/gene_locus_annotation.pl -gtf merged_filter_multiexon_junction40_TPM.gtf -tmap gffcmp.merged.gtf.tmap > merged_filter_multiexon_junction40_TPM_annotation.xls

#### Intron list from GTF

perl /home/pcd/GTF_intron.pl -GTF merged_filter_multiexon_junction40_TPM.gtf > intron_list.txt

#### Gene feature from GTF
perl /home/pcd/select_tid_from_GTF.pl -GTF merged_filter_multiexon_junction40_TPM.gtf -Genome /home/database/genome.fa > gene_feature.xls

#### Intron length Distribution
perl /home/pcd/select_tid_from_GTF_intron_dis.pl -GTF merged_filter_multiexon_junction40_TPM.gtf


perl /home/pcd/filter_fa_by_GTF_tid.pl -GTF merged_filter_multiexon_junction40_TPM.gtf -input merged_filter_multiexon_junction40.fa > merged_filter_multiexon_junction40_TPM.fa


#### Coding potential prediction of novel assembled transcript sequences
perl /home/pcd/filter_gtf_by_novel.pl -GTF merged_filter_multiexon_junction40_TPM.gtf > novel.gtf
perl /home/pcd/filter_seq_by_novel.pl -fa merged_filter_multiexon_junction40_TPM.fa > novel.fa
python /home/soft/CPC2-beta/bin/CPC2.py -i novel.fa -o novel.fa.cpc2
/home/soft/TransDecoder-3.0.1/TransDecoder.LongOrfs -S -m 30 -t novel.fa
/home/soft/TransDecoder-3.0.1/TransDecoder.Predict --cpu 8 --single_best_orf --retain_long_orfs 100 -t novel.fa
perl /home/soft/select_longest_ORF.pl -fa novel.fa.transdecoder.pep > novel.fa.transdecoder.pep.aa
perl /home/soft/select_protein.pl
perl /home/soft/mRNA_pos_genome.pl > novel_AS_CDS.gtf
perl /home/pcd/filter_cds_gtf.pl > merged_filter_multiexon_junction40_TPM_CDS.gtf



#Known AS event (RI SE A3SS A5SS MXE) detection from GTF
cd /home/database/
/home/soft/samtools view -b -H -o blank.bam /home/project/ControlR1/accepted_hits.bam
/home/soft/samtools index blank.bam
python /public/analysis/hefei/software/rMATS.3.0.9/bin/processGTF.SAMs.py RAPDBMSU_mRNA.gtf ./AS blank.bam fr-unstranded PE 1
