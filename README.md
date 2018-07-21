# Alternative-Splicing-Study-Under-Mineral-Nutrient-Deficiency

A critical role for alternative splicing in maintaining mineral nutrient homeostasis in rice (Oryza sativa). Chunlan Dong, Fei He, et al. (2018). Summary of the study was that genome-wide RNA-seq and mutant analyses reveal the importance of alternative splicing in the mineral deficiency response and in mineral uptake and remobilization in rice.

The repository has 2 folders, including pcd and database. The shell/R code in the main folder is used to replicate the alternative splicing analysis pipeline, including main figures. The shell script provide the step-by-step guide to analize mineral nutrient deficiency RNASeq dataset analysis. The perl/R in the pcd folder will be also called in each shell script. The genome sequence and several associated gene annotation files, which will be the intermediate result of the pipeline is deposited in the database fold. 
 
#1. 0.create_genome_database.sh
Download rice chromosome sequences and gene annotation files from RAPDB and MSU database. The non-redundant reference gene GTF file from two database will be merged. The merged gtf will be used for RNASeq read STAR alignment. These files will be deposited in the database fold.

2. 1.QC_and_alignment_RNASeqAnalysis_of_nutrient_deficiency.sh
Download all the 15 RNA-Seq data of Micronutrient deficiency, which are available under the accession number SRP117202(https://www.ncbi.nlm.nih.gov/sra/?term=SRP117202). Download all the 126 RNA-Seq data of the time course of Pi-starved root and shoot samples, which are available under the accession number SRA097415 (Secco et al., 2013)(https://www.ncbi.nlm.nih.gov/sra/?term=SRA097415). Create project folder, and each sample name folder under the project folder. FastQ trimming and alignment will be finished in each sample name folder.

3. 2.Assembly_RNASeqAnalysis_of_nutrient_deficiency.sh
Single sample assembly using StringTie, multiple sample assembly using TACO, a series of filters of assembled novel transcripts, and alternative splicing event detection from final GTF file.

4. 3.ASE_RNASeqAnalysis_of_nutrient_deficiency.sh
Alternative splicing event quantification for all samples using rMATS, splice Site distribution of known or novel ASEs, the correlation analysis between ASE frequency, GC percentage, expression level (TPM), exon number, intron number, exon length, intron length, and protein domain enrichment analysis of high AS frequency genes.

5. 4.DEG_RNASeqAnalysis_of_nutrient_deficiency.sh
Identify differential expressed genes under micronutrient deficiency or Pi deficiency, and Gene ontology, KEGG pathway enrichment analysis.

6. 5.DASG_RNASeqAnalysis_of_nutrient_deficiency.sh
Identify differential alternative splicing events under micronutrient deficiency or or Pi deficiency, and Gene ontology, KEGG pathway enrichment analysis.

7. Figure1A_ASE_structure_and_ratio.r
Overview of the five different types of AS and their frequency in the RNA-Seq datasets. The input file 'PiFeZnCuMn_ASE_type.xls' was the output file from '3.ASE_RNASeqAnalysis_of_nutrient_deficiency.sh'.

8. Figure1B-F_DASE.sh, SupplementalFigure4A-F_DASE.sh, Figure5C-D_DASE.sh
ASE visualization of RNASeq data. The input file 'PiFeZnCuMn_in_ex.xls' and 'PiFeZnCuMn_psi.xls' was the output file from '3.ASE_RNASeqAnalysis_of_nutrient_deficiency.sh'. The input file 'merged_filter_multiexon_junction40_TPM_CDS.gtf', 'merged_filter_multiexon_junction40_TPM.gtf' and 'merged_filter_multiexon_junction40_TPM.fa' was the output file from '2.Assembly_RNASeqAnalysis_of_nutrient_deficiency.sh'. The input file 'PiFeZnCuMn_transcript_sample_TPM_filter.xls' was the output file from '4.DEG_RNASeqAnalysis_of_nutrient_deficiency.sh'.

9. Figure3A_GO_enrichment_barplot.sh
Functional categorization (Gene Ontology) of DASGs and DEGs under Fe-deficiency and P-deficiency conditions in rice roots. The input file 'GO_enrichment.xls' was the output file from '4.DEG_RNASeqAnalysis_of_nutrient_deficiency.sh'.

10. Figure4A-D_gene_structure.sh
Alternative splicing of four genes under Fe-deficiency or P-deficiency conditions in rice roots. The input file 'PiFeZnCuMn_in_ex.xls' and 'PiFeZnCuMn_psi.xls' were the output file from '3.ASE_RNASeqAnalysis_of_nutrient_deficiency.sh'. The input file 'merged_filter_multiexon_junction40_TPM_CDS.gtf', 'merged_filter_multiexon_junction40_TPM.gtf' and 'merged_filter_multiexon_junction40_TPM.fa' were the output file from '2.Assembly_RNASeqAnalysis_of_nutrient_deficiency.sh'. The input file 'PiFeZnCuMn_transcript_sample_TPM_filter.xls' was the output file from '4.DEG_RNASeqAnalysis_of_nutrient_deficiency.sh'. 'LOC_Os03g06520.1_domain.gtf', 'LOC_Os07g30130.1_domain.gtf', 'LOC_Os02g55910.1_domain.gtf' and 'LOC_Os05g48390.1_domain.gtf' were deposited in the database folder.

11. Figure5A_high_ASG_domain_barplot_high.r
Protein domain barplot of high-frequency alternative splicing genes. The input file 'high_AS_genes_domain_enrichment_significant.xls' was the output file from '3.ASE_RNASeqAnalysis_of_nutrient_deficiency.sh'.

12. Figure5B_SR_DASG_heatmap.r
Alternative splicing of 22 rice SR genes under mineral nutrient deficiency. The input file 'SRs_FeZnCuMn_DASG_overlap.xls', 'SRs_PiRoot_DASG_overlap.xls', and 'SRs_PiShoot_DASG_overlap.xls' were the output file from '5.DASG_RNASeqAnalysis_of_nutrient_deficiency.sh'.

13. SupplementalFigure2A_and_B_intron_length_distribution.r
Intron length distribution plot. The input file 'gene_intron_distribution1.xls' and 'gene_intron_distribution2.xls' were the output file from '2.Assembly_RNASeqAnalysis_of_nutrient_deficiency.sh'.

14. SupplementalFigure2C_SpliceSite_proportion_known_novel.r
Proportion of canonical or non-canonical splice site usage in known and novel AS. The input file 'splice_site_stat.xls' was the output file from '3.ASE_RNASeqAnalysis_of_nutrient_deficiency.sh'.

15. SupplementalFigure2D_ASE_feature_boxplot.r
AS frequency of intron-containing genes associated with exon or intron number, exon or intron length, transcriptional abundance and GC content. All input files were the output file from '3.ASE_RNASeqAnalysis_of_nutrient_deficiency.sh'.

16. SupplementalFigure3_ASG_stat.r The input file 'ASG_condition.xls' was the output file from '3.ASE_RNASeqAnalysis_of_nutrient_deficiency.sh'.
 
