#Download bedtools-2.17.0 software https://github.com/arq5x/bedtools2/
#Download STAR version 2.5.2b https://github.com/alexdobin/STAR
#install softwares to a directory (for example /home/soft/)
#copy database to a directory (for example /home/database/)


#Merge RAPDB and MSU databases
wget http://rapdb.dna.affrc.go.jp/download/archive/irgsp1/IRGSP-1.0_genome.fasta.gz
wget http://rapdb.dna.affrc.go.jp/download/archive/irgsp1/IRGSP-1.0_representative_2016-08-05.tar.gz
wget http://rapdb.dna.affrc.go.jp/download/archive/irgsp1/IRGSP-1.0_protein_2016-08-05.fasta.gz
wget http://rapdb.dna.affrc.go.jp/download/archive/irgsp1/IRGSP-1.0_predicted_2016-08-05.tar.gz
wget http://rapdb.dna.affrc.go.jp/download/archive/irgsp1/IRGSP-1.0_predicted-protein_2016-08-05.fasta.gz
wget ftp://ftp.plantbiology.msu.edu/pub/data/Eukaryotic_Projects/o_sativa/annotation_dbs/pseudomolecules/version_7.0/all.dir/all.locus_brief_info.7.0
wget ftp://ftp.plantbiology.msu.edu/pub/data/Eukaryotic_Projects/o_sativa/annotation_dbs/pseudomolecules/version_7.0/all.dir/all.pep
wget ftp://ftp.plantbiology.msu.edu/pub/data/Eukaryotic_Projects/o_sativa/annotation_dbs/pseudomolecules/version_7.0/all.dir/all.gff3
tar -xvf IRGSP-1.0_representative_2016-08-05.tar.gz
mv ./IRGSP-1.0_representative_2016-08-05/transcripts_exon.gff repre_transcripts_exon.gff
gunzip -c IRGSP-1.0_protein_2016-08-05.fasta.gz > repre_protein.fasta
tar -xvf IRGSP-1.0_predicted_2016-08-05.tar.gz
gunzip -c IRGSP-1.0_predicted-protein_2016-08-05.fasta.gz > predicted_protein.fasta
mv ./IRGSP-1.0_predicted_2016-08-05/transcripts.gff predicted_transcripts.gff

perl /home/pcd/change_gtf_chr_MSU.pl > all_MSU.gff3
perl /home/pcd/gff2gtf_MSU.pl > all_MSU.gtf
perl /home/pcd/gff2gtf_rapdb_predict.pl > predicted_transcripts.gtf
perl /home/pcd/gff2gtf_rapdb_repre.pl > repre_transcripts_exon.gtf
cat predicted_transcripts.gtf repre_transcripts_exon.gtf > all_RAPDB.gtf

perl /home/pcd/gene_bed.pl > RAPDB_MSU_mRNA.bed
/home/soft/bedtools-2.17.0/bin/sortBed -i RAPDB_MSU_mRNA.bed > RAPDB_MSU_mRNA_sort.bed
/home/soft/bedtools-2.17.0/bin/mergeBed -s -nms  -i RAPDB_MSU_mRNA_sort.bed > RAPDB_MSU_mRNA_merge.bed
rm -f RAPDB_MSU_mRNA_sort.bed
cat all_MSU.gtf repre_transcripts_exon.gtf predicted_transcripts.gtf > RAPDBMSU_mRNA.gtf

perl /home/pcd/replace_tid2gid_ref.pl > RAPDBMSU_mRNA_nx.gtf
perl /home/pcd/remove_uniq_transcript.pl -gtf RAPDBMSU_mRNA_nx.gtf > uniq.txt
perl /home/pcd/filter_gtf_by_tid_RAPDMSU.pl > RAPDBMSU_mRNA_nx2.gtf
mv RAPDBMSU_mRNA_nx2.gtf RAPDBMSU_mRNA.gtf


#STAR index genome for alignment
/home/soft/STAR --runMode genomeGenerate --runThreadN 24 --sjdbGTFfile RAPDBMSU_mRNA.gtf --genomeDir ./ --genomeFastaFiles genome.fa



