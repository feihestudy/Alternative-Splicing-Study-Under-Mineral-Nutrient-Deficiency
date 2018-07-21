# Identify differential expressed genes under micronutrient deficiency

mkdir DEG
rm -rf FeZnCuMn
mkdir FeZnCuMn
cd FeZnCuMn
rm -f FeZnCuMn_samplelist_condition
echo -e "ControlR1\0011ControlR1" > FeZnCuMn_samplelist_condition
echo -e "ControlR2\0011ControlR2" >> FeZnCuMn_samplelist_condition
echo -e "ControlR3\0011\ControlR3" >> FeZnCuMn_samplelist_condition
echo -e "FeR1\0011FeR1" >> FeZnCuMn_samplelist_condition
echo -e "FeR2\0011FeR2" >> FeZnCuMn_samplelist_condition
echo -e "FeR3\0011FeR3" >> FeZnCuMn_samplelist_condition
echo -e "ZnR1\0011ZnR1" >> FeZnCuMn_samplelist_condition
echo -e "ZnR2\0011ZnR2" >> FeZnCuMn_samplelist_condition
echo -e "ZnR3\0011ZnR3" >> FeZnCuMn_samplelist_condition
echo -e "CuR1\0011CuR1" >> FeZnCuMn_samplelist_condition
echo -e "CuR2\0011CuR2" >> FeZnCuMn_samplelist_condition
echo -e "CuR3\0011CuR3" >> FeZnCuMn_samplelist_condition
echo -e "MnR1\0011MnR1" >> FeZnCuMn_samplelist_condition
echo -e "MnR2\0011MnR2" >> FeZnCuMn_samplelist_condition
echo -e "MnR3\0011MnR3" >> FeZnCuMn_samplelist_condition
echo -e "MnR1\0011MnR1" >> FeZnCuMn_samplelist_condition
echo -e "MnR2\0011MnR2" >> FeZnCuMn_samplelist_condition
unset arr
k=1
for i in `cat FeZnCuMn_samplelist_condition | cut -f 1`
do
   arr[$k]=$i;
   k=$k+1;
done
samplenames=${arr[1]}
for i in ${arr[@]:2}  
do
   samplenames=$samplenames","${i};
done
cp /home/project/gid2tid.txt ./
perl /home/pcd/get_sample_gene_TPM_salmon.pl -dir /home/project/salmon/ -sample $samplenames > FeZnCuMn_gene_sample_TPM.xls
perl /home/pcd/filter_ex_by_GTF.pl -GTF /home/project/merged_filter_multiexon_junction40_TPM.gtf -input FeZnCuMn_gene_sample_TPM.xls > FeZnCuMn_gene_sample_TPM_filter.xls
Rscript /home/pcd/arrange_sample_arg.r FeZnCuMn_gene_sample_TPM_filter.xls FeZnCuMn_samplelist_condition
mv result.xls FeZnCuMn_gene_sample_TPM_filter.xls
Rscript /home/pcd/log2_arg.r FeZnCuMn_gene_sample_TPM_filter.xls
mv result.xls FeZnCuMn_gene_sample_log2TPM_filter.xls
perl /home/pcd/get_sample_gene_count_salmon.pl -dir /home/project/salmon/ -sample $samplenames > FeZnCuMn_gene_sample_count.xls
perl /home/pcd/filter_ex_by_GTF.pl -GTF /home/project/merged_filter_multiexon_junction40_TPM.gtf -input FeZnCuMn_gene_sample_count.xls > FeZnCuMn_gene_sample_count_filter.xls
Rscript /home/pcd/arrange_sample_arg.r FeZnCuMn_gene_sample_count_filter.xls FeZnCuMn_samplelist_condition
mv result.xls FeZnCuMn_gene_sample_count_filter.xls
Rscript /home/pcd/deseq2_FeZnCuMn.r
perl /home/pcd/add_anno_rice_DEG.pl -degene Fe_vs_Control_gene_exp_significant.xls -annotation /home/project/merged_filter_multiexon_junction40_TPM_annotation.xls > Fe_vs_Control_gene_exp_significant_with_annotation.xls
perl /home/pcd/add_anno_rice_DEG.pl -degene Zn_vs_Control_gene_exp_significant.xls -annotation /home/project/merged_filter_multiexon_junction40_TPM_annotation.xls > Zn_vs_Control_gene_exp_significant_with_annotation.xls
perl /home/pcd/add_anno_rice_DEG.pl -degene Cu_vs_Control_gene_exp_significant.xls -annotation /home/project/merged_filter_multiexon_junction40_TPM_annotation.xls > Cu_vs_Control_gene_exp_significant_with_annotation.xls
perl /home/pcd/add_anno_rice_DEG.pl -degene Mn_vs_Control_gene_exp_significant.xls -annotation /home/project/merged_filter_multiexon_junction40_TPM_annotation.xls > Mn_vs_Control_gene_exp_significant_with_annotation.xls
perl /home/pcd/add_anno_rice_DEG.pl -degene Fe_vs_Control_gene_exp.xls -annotation /home/project/merged_filter_multiexon_junction40_TPM_annotation.xls > Fe_vs_Control_gene_exp_with_annotation.xls
perl /home/pcd/add_anno_rice_DEG.pl -degene Zn_vs_Control_gene_exp.xls -annotation /home/project/merged_filter_multiexon_junction40_TPM_annotation.xls > Zn_vs_Control_gene_exp_with_annotation.xls
perl /home/pcd/add_anno_rice_DEG.pl -degene Cu_vs_Control_gene_exp.xls -annotation /home/project/merged_filter_multiexon_junction40_TPM_annotation.xls > Cu_vs_Control_gene_exp_with_annotation.xls
perl /home/pcd/add_anno_rice_DEG.pl -degene Mn_vs_Control_gene_exp.xls -annotation /home/project/merged_filter_multiexon_junction40_TPM_annotation.xls > Mn_vs_Control_gene_exp_with_annotation.xls
# GO and KEGG enrichment analysis of DEGs under micronutrient deficiency
Funs=(Fe Zn Cu Mn)
for((i=0; i<=3; i++))
do
rm -rf ${Funs[$i]}
mkdir ${Funs[$i]}
cd ${Funs[$i]}
ln "/home/project/DEG/FeZnCuMn/"${Funs[$i]}"_vs_Control_gene_exp_significant_with_annotation.xls" ./target.xls
perl /home/pcd/degene_for_DEG.pl -degene target.xls > diffgene.txt
perl /home/pcd/KEGG_enrich.pl -ispvalue "Pvalue" -bl /home/database/allgene.txt -gl diffgene.txt -symbol2KEGG /home/database/symbol2KEGG.txt
perl /home/pcd/GO_enrich.pl -ispvalue "Pvalue" -bl /home/database/allgene.txt -gl diffgene.txt -symbol2GO /home/database/symbol2GO.txt
cd ..
done


# Identify differential expressed genes under Pi deficiency in roots
rm -rf PiRoot
mkdir PiRoot
cd PiRoot

rm -f PiRoot_samplelist_condition
echo -e "H1piPlusR1_SRR1005257\0011H1piPlusR1_R" >> PiRoot_samplelist_condition
echo -e "H1piPlusR2_SRR1005258\0011H1piPlusR2_R" >> PiRoot_samplelist_condition
echo -e "H1piPlusR3_SRR1005259\0011H1piPlusR3_R" >> PiRoot_samplelist_condition
echo -e "H1piMinusR1_SRR1005260\0011H1piMinusR1_R" >> PiRoot_samplelist_condition
echo -e "H1piMinusR2_SRR1005261\0011H1piMinusR2_R" >> PiRoot_samplelist_condition
echo -e "H1piMinusR3_SRR1005262\0011H1piMinusR3_R" >> PiRoot_samplelist_condition
echo -e "H6piPlusR1_SRR1005308\0011H6piPlusR1_R" >> PiRoot_samplelist_condition
echo -e "H6piPlusR2_SRR1005309\0011H6piPlusR2_R" >> PiRoot_samplelist_condition
echo -e "H6piPlusR3_SRR1005310\0011H6piPlusR3_R" >> PiRoot_samplelist_condition
echo -e "H6piMinusR1_SRR1005311\0011H6piMinusR1_R" >> PiRoot_samplelist_condition
echo -e "H6piMinusR2_SRR1005312\0011H6piMinusR2_R" >> PiRoot_samplelist_condition
echo -e "H6piMinusR3_SRR1005313\0011H6piMinusR3_R" >> PiRoot_samplelist_condition
echo -e "H24piPlusR1_SRR1005296\0011H24piPlusR1_R" >> PiRoot_samplelist_condition
echo -e "H24piPlusR2_SRR1005297\0011H24piPlusR2_R" >> PiRoot_samplelist_condition
echo -e "H24piPlusR3_SRR1005298\0011H24piPlusR3_R" >> PiRoot_samplelist_condition
echo -e "H24piMinusR1_SRR1005299\0011H24piMinusR1_R" >> PiRoot_samplelist_condition
echo -e "H24piMinusR2_SRR1005300\0011H24piMinusR2_R" >> PiRoot_samplelist_condition
echo -e "H24piMinusR3_SRR1005301\0011H24piMinusR3_R" >> PiRoot_samplelist_condition
echo -e "D3piPlusR1_SRR1005302\0011D3piPlusR1_R" >> PiRoot_samplelist_condition
echo -e "D3piPlusR2_SRR1005303\0011D3piPlusR2_R" >> PiRoot_samplelist_condition
echo -e "D3piPlusR3_SRR1005304\0011D3piPlusR3_R" >> PiRoot_samplelist_condition
echo -e "D3piMinusR1_SRR1005305\0011D3piMinusR1_R" >> PiRoot_samplelist_condition
echo -e "D3piMinusR2_SRR1005306\0011D3piMinusR2_R" >> PiRoot_samplelist_condition
echo -e "D3piMinusR3_SRR1005307\0011D3piMinusR3_R" >> PiRoot_samplelist_condition
echo -e "D7piPlusR1_SRR1005314\0011D7piPlusR1_R" >> PiRoot_samplelist_condition
echo -e "D7piPlusR2_SRR1005315\0011D7piPlusR2_R" >> PiRoot_samplelist_condition
echo -e "D7piPlusR3_SRR1005316\0011D7piPlusR3_R" >> PiRoot_samplelist_condition
echo -e "D7piMinusR1_SRR1005317\0011D7piMinusR1_R" >> PiRoot_samplelist_condition
echo -e "D7piMinusR2_SRR1005318\0011D7piMinusR2_R" >> PiRoot_samplelist_condition
echo -e "D7piMinusR3_SRR1005319\0011D7piMinusR3_R" >> PiRoot_samplelist_condition
echo -e "D21piPlusR1_SRR1005290\0011D21piPlusR1_R" >> PiRoot_samplelist_condition
echo -e "D21piPlusR2_SRR1005291\0011D21piPlusR2_R" >> PiRoot_samplelist_condition
echo -e "D21piPlusR3_SRR1005292\0011D21piPlusR3_R" >> PiRoot_samplelist_condition
echo -e "D21piMinusR1_SRR1005293\0011D21piMinusR1_R" >> PiRoot_samplelist_condition
echo -e "D21piMinusR2_SRR1005294\0011D21piMinusR2_R" >> PiRoot_samplelist_condition
echo -e "D21piMinusR3_SRR1005295\0011D21piMinusR3_R" >> PiRoot_samplelist_condition
echo -e "D21H1piPlusR1_SRR1005263\0011D21H1piPlusR1_R" >> PiRoot_samplelist_condition
echo -e "D21H1piPlusR2_SRR1005264\0011D21H1piPlusR2_R" >> PiRoot_samplelist_condition
echo -e "D21H1piPlusR3_SRR1005265\0011D21H1piPlusR3_R" >> PiRoot_samplelist_condition
echo -e "D21H1piMinusR1_SRR1005266\0011D21H1piMinusR1_R" >> PiRoot_samplelist_condition
echo -e "D21H1piMinusR2_SRR1005267\0011D21H1piMinusR2_R" >> PiRoot_samplelist_condition
echo -e "D21H1piMinusR3_SRR1005268\0011D21H1piMinusR3_R" >> PiRoot_samplelist_condition
echo -e "D21H1piPlusRecR1_SRR1005269\0011D21H1piPlusRecR1_R" >> PiRoot_samplelist_condition
echo -e "D21H1piPlusRecR2_SRR1005270\0011D21H1piPlusRecR2_R" >> PiRoot_samplelist_condition
echo -e "D21H1piPlusRecR3_SRR1005271\0011D21H1piPlusRecR3_R" >> PiRoot_samplelist_condition
echo -e "D21H6piPlusR1_SRR1005281\0011D21H6piPlusR1_R" >> PiRoot_samplelist_condition
echo -e "D21H6piPlusR2_SRR1005282\0011D21H6piPlusR2_R" >> PiRoot_samplelist_condition
echo -e "D21H6piPlusR3_SRR1005283\0011D21H6piPlusR3_R" >> PiRoot_samplelist_condition
echo -e "D21H6piMinusR1_SRR1005284\0011D21H6piMinusR1_R" >> PiRoot_samplelist_condition
echo -e "D21H6piMinusR2_SRR1005285\0011D21H6piMinusR2_R" >> PiRoot_samplelist_condition
echo -e "D21H6piMinusR3_SRR1005286\0011D21H6piMinusR3_R" >> PiRoot_samplelist_condition
echo -e "D21H6piPlusRecR1_SRR1005287\0011D21H6piPlusRecR1_R" >> PiRoot_samplelist_condition
echo -e "D21H6piPlusRecR2_SRR1005288\0011D21H6piPlusRecR2_R" >> PiRoot_samplelist_condition
echo -e "D21H6piPlusRecR3_SRR1005289\0011D21H6piPlusRecR3_R" >> PiRoot_samplelist_condition
echo -e "D21H24piPlusR1_SRR1005272\0011D21H24piPlusR1_R" >> PiRoot_samplelist_condition
echo -e "D21H24piPlusR2_SRR1005273\0011D21H24piPlusR2_R" >> PiRoot_samplelist_condition
echo -e "D21H24piPlusR3_SRR1005274\0011D21H24piPlusR3_R" >> PiRoot_samplelist_condition
echo -e "D21H24piMinusR1_SRR1005275\0011D21H24piMinusR1_R" >> PiRoot_samplelist_condition
echo -e "D21H24piMinusR2_SRR1005276\0011D21H24piMinusR2_R" >> PiRoot_samplelist_condition
echo -e "D21H24piMinusR3_SRR1005277\0011D21H24piMinusR3_R" >> PiRoot_samplelist_condition
echo -e "D21H24piPlusRecR1_SRR1005278\0011D21H24piPlusRecR1_R" >> PiRoot_samplelist_condition
echo -e "D21H24piPlusRecR2_SRR1005279\0011D21H24piPlusRecR2_R" >> PiRoot_samplelist_condition
echo -e "D21H24piPlusRecR3_SRR1005280\0011D21H24piPlusRecR3_R" >> PiRoot_samplelist_condition

unset arr
k=1
for i in `cat PiRoot_samplelist_condition | cut -f 1`
do
   arr[$k]=$i;
   k=$k+1;
done;
samplenames=${arr[1]}
for i in ${arr[@]:2}  
do
   samplenames=$samplenames","${i};
done
cp /home/project/gid2tid.txt ./
perl /home/pcd/get_sample_gene_TPM_salmon.pl  -dir /home/project/salmon/ -sample $samplenames > PiRoot_gene_sample_TPM.xls
perl /home/pcd/filter_ex_by_GTF.pl -GTF /home/project/merged_filter_multiexon_junction40_TPM.gtf -input PiRoot_gene_sample_TPM.xls > PiRoot_gene_sample_TPM_filter.xls
Rscript /home/pcd/arrange_sample_arg.r PiRoot_gene_sample_TPM_filter.xls PiRoot_samplelist_condition
mv result.xls PiRoot_gene_sample_TPM_filter.xls
Rscript /home/pcd/log2_arg.r PiRoot_gene_sample_TPM_filter.xls
mv result.xls PiRoot_gene_sample_log2TPM_filter.xls
perl /home/pcd/get_sample_gene_count_salmon.pl -dir /home/project/salmon/ -sample $samplenames > PiRoot_gene_sample_count.xls
perl /home/pcd/filter_ex_by_GTF.pl -GTF /home/project/merged_filter_multiexon_junction40_TPM.gtf -input PiRoot_gene_sample_count.xls > PiRoot_gene_sample_count_filter.xls
Rscript /home/pcd/arrange_sample_arg.r PiRoot_gene_sample_count_filter.xls PiRoot_samplelist_condition
mv result.xls PiRoot_gene_sample_count_filter.xls

Rscript /home/pcd/deseq2_PiRoot.r

# union all DEGs form DEGs under -Pi at each time
perl /home/pcd/table2table_diffgene.pl -samplepair "H1piMinus_vs_H1piPlus_R,H6piMinus_vs_H6piPlus_R,H24piMinus_vs_H24piPlus_R,D3piMinus_vs_D3piPlus_R,D7piMinus_vs_D7piPlus_R,D21piMinus_vs_D21piPlus_R,D21H1piMinus_vs_D21H1piPlus_R,D21H6piMinus_vs_D21H6piPlus_R,D21H24piMinus_vs_D21H24piPlus_R" > PiRoot_DEG_overlap.xls
perl /home/pcd/table2table_diffgene_nx.pl -rawcount "PiRoot_gene_sample_count_filter.xls" -normcount "PiRoot_gene_sample_normalized_count.xls" -TPM "PiRoot_gene_sample_TPM_filter.xls" -target PiRoot_DEG_overlap.xls -samplepair "H1piMinus_vs_H1piPlus_R,H6piMinus_vs_H6piPlus_R,H24piMinus_vs_H24piPlus_R,D3piMinus_vs_D3piPlus_R,D7piMinus_vs_D7piPlus_R,D21piMinus_vs_D21piPlus_R,D21H1piMinus_vs_D21H1piPlus_R,D21H6piMinus_vs_D21H6piPlus_R,D21H24piMinus_vs_D21H24piPlus_R" > PiRoot_DEG_overlap_info.xls

perl /home/pcd/add_anno_DEG.pl -degene D21piMinus_vs_D21piPlus_R_gene_exp_significant.xls -annotation /home/project/merged_filter_multiexon_junction40_TPM_annotation.xls > D21piMinus_vs_D21piPlus_R_gene_exp_significant_with_annotation.xls

perl /home/pcd/add_anno_DEG.pl -degene PiRoot_DEG_overlap_info.xls -annotation /home/project/merged_filter_multiexon_junction40_TPM_annotation.xls > PiRoot_DEG_overlap_info_with_annotation.xls


# GO and KEGG enrichment analysis of DEGs under Pi deficiency in roots
perl /home/pcd/degene_for_DEG.pl -degene PiRoot_DEG_overlap_info_with_annotation.xls > diffgene.txt
perl /home/pcd/KEGG_enrich.pl -ispvalue "Pvalue" -bl /home/database/allgene.txt -gl diffgene.txt -symbol2KEGG /home/database/symbol2KEGG.txt
perl /home/pcd/GO_enrich.pl -ispvalue "Pvalue" -bl /home/database/allgene.txt -gl diffgene.txt -symbol2GO /home/database/symbol2GO.txt
perl /home/pcd/degene_for_DEG.pl -degene D21piMinus_vs_D21piPlus_R_gene_exp_significant_with_annotation.xls > diffgene.txt
perl /home/pcd/KEGG_enrich.pl -ispvalue "Pvalue" -bl /home/database/allgene.txt -gl diffgene.txt -symbol2KEGG /home/database/symbol2KEGG.txt
perl /home/pcd/GO_enrich.pl -ispvalue "Pvalue" -bl /home/database/allgene.txt -gl diffgene.txt -symbol2GO /home/database/symbol2GO.txt
mv GO_enrichment.xls D21_R_GO_enrichment.xls
mv GO_enrichment_significant.xls D21_R_GO_enrichment_significant.xls
mv KEGG_enrichment.xls D21_R_KEGG_enrichment.xls
mv KEGG_enrichment_significant.xls D21_R_KEGG_enrichment_significant.xls



# Identify differential expressed genes under Pi deficiency in shoots

rm -f PiShoot_samplelist_condition
echo -e "H1piPlusR1_SRR1005320\0011H1piPlusR1_S" >> PiShoot_samplelist_condition
echo -e "H1piPlusR2_SRR1005321\0011H1piPlusR2_S" >> PiShoot_samplelist_condition
echo -e "H1piPlusR3_SRR1005322\0011H1piPlusR3_S" >> PiShoot_samplelist_condition
echo -e "H1piMinusR1_SRR1005323\0011H1piMinusR1_S" >> PiShoot_samplelist_condition
echo -e "H1piMinusR2_SRR1005324\0011H1piMinusR2_S" >> PiShoot_samplelist_condition
echo -e "H1piMinusR3_SRR1005325\0011H1piMinusR3_S" >> PiShoot_samplelist_condition
echo -e "H6piPlusR1_SRR1005371\0011H6piPlusR1_S" >> PiShoot_samplelist_condition
echo -e "H6piPlusR2_SRR1005372\0011H6piPlusR2_S" >> PiShoot_samplelist_condition
echo -e "H6piPlusR3_SRR1005373\0011H6piPlusR3_S" >> PiShoot_samplelist_condition
echo -e "H6piMinusR1_SRR1005374\0011H6piMinusR1_S" >> PiShoot_samplelist_condition
echo -e "H6piMinusR2_SRR1005375\0011H6piMinusR2_S" >> PiShoot_samplelist_condition
echo -e "H6piMinusR3_SRR1005376\0011H6piMinusR3_S" >> PiShoot_samplelist_condition
echo -e "H24piPlusR1_SRR1005359\0011H24piPlusR1_S" >> PiShoot_samplelist_condition
echo -e "H24piPlusR2_SRR1005360\0011H24piPlusR2_S" >> PiShoot_samplelist_condition
echo -e "H24piPlusR3_SRR1005361\0011H24piPlusR3_S" >> PiShoot_samplelist_condition
echo -e "H24piMinusR1_SRR1005362\0011H24piMinusR1_S" >> PiShoot_samplelist_condition
echo -e "H24piMinusR2_SRR1005363\0011H24piMinusR2_S" >> PiShoot_samplelist_condition
echo -e "H24piMinusR3_SRR1005364\0011H24piMinusR3_S" >> PiShoot_samplelist_condition
echo -e "D3piPlusR1_SRR1005365\0011D3piPlusR1_S" >> PiShoot_samplelist_condition
echo -e "D3piPlusR2_SRR1005366\0011D3piPlusR2_S" >> PiShoot_samplelist_condition
echo -e "D3piPlusR3_SRR1005367\0011D3piPlusR3_S" >> PiShoot_samplelist_condition
echo -e "D3piMinusR1_SRR1005368\0011D3piMinusR1_S" >> PiShoot_samplelist_condition
echo -e "D3piMinusR2_SRR1005369\0011D3piMinusR2_S" >> PiShoot_samplelist_condition
echo -e "D3piMinusR3_SRR1005370\0011D3piMinusR3_S" >> PiShoot_samplelist_condition
echo -e "D7piPlusR1_SRR1005377\0011D7piPlusR1_S" >> PiShoot_samplelist_condition
echo -e "D7piPlusR2_SRR1005378\0011D7piPlusR2_S" >> PiShoot_samplelist_condition
echo -e "D7piPlusR3_SRR1005379\0011D7piPlusR3_S" >> PiShoot_samplelist_condition
echo -e "D7piMinusR1_SRR1005380\0011D7piMinusR1_S" >> PiShoot_samplelist_condition
echo -e "D7piMinusR2_SRR1005381\0011D7piMinusR2_S" >> PiShoot_samplelist_condition
echo -e "D7piMinusR3_SRR1005382\0011D7piMinusR3_S" >> PiShoot_samplelist_condition
echo -e "D21piPlusR1_SRR1005353\0011D21piPlusR1_S" >> PiShoot_samplelist_condition
echo -e "D21piPlusR2_SRR1005354\0011D21piPlusR2_S" >> PiShoot_samplelist_condition
echo -e "D21piPlusR3_SRR1005355\0011D21piPlusR3_S" >> PiShoot_samplelist_condition
echo -e "D21piMinusR1_SRR1005356\0011D21piMinusR1_S" >> PiShoot_samplelist_condition
echo -e "D21piMinusR2_SRR1005357\0011D21piMinusR2_S" >> PiShoot_samplelist_condition
echo -e "D21piMinusR3_SRR1005358\0011D21piMinusR3_S" >> PiShoot_samplelist_condition
echo -e "D21H1piPlusR1_SRR1005326\0011D21H1piPlusR1_S" >> PiShoot_samplelist_condition
echo -e "D21H1piPlusR2_SRR1005327\0011D21H1piPlusR2_S" >> PiShoot_samplelist_condition
echo -e "D21H1piPlusR3_SRR1005328\0011D21H1piPlusR3_S" >> PiShoot_samplelist_condition
echo -e "D21H1piMinusR1_SRR1005329\0011D21H1piMinusR1_S" >> PiShoot_samplelist_condition
echo -e "D21H1piMinusR2_SRR1005330\0011D21H1piMinusR2_S" >> PiShoot_samplelist_condition
echo -e "D21H1piMinusR3_SRR1005331\0011D21H1piMinusR3_S" >> PiShoot_samplelist_condition
echo -e "D21H1piPlusRecR1_SRR1005332\0011D21H1piPlusRecR1_S" >> PiShoot_samplelist_condition
echo -e "D21H1piPlusRecR2_SRR1005333\0011D21H1piPlusRecR2_S" >> PiShoot_samplelist_condition
echo -e "D21H1piPlusRecR3_SRR1005334\0011D21H1piPlusRecR3_S" >> PiShoot_samplelist_condition
echo -e "D21H6piPlusR1_SRR1005344\0011D21H6piPlusR1_S" >> PiShoot_samplelist_condition
echo -e "D21H6piPlusR2_SRR1005345\0011D21H6piPlusR2_S" >> PiShoot_samplelist_condition
echo -e "D21H6piPlusR3_SRR1005346\0011D21H6piPlusR3_S" >> PiShoot_samplelist_condition
echo -e "D21H6piMinusR1_SRR1005347\0011D21H6piMinusR1_S" >> PiShoot_samplelist_condition
echo -e "D21H6piMinusR2_SRR1005348\0011D21H6piMinusR2_S" >> PiShoot_samplelist_condition
echo -e "D21H6piMinusR3_SRR1005349\0011D21H6piMinusR3_S" >> PiShoot_samplelist_condition
echo -e "D21H6piPlusRecR1_SRR1005350\0011D21H6piPlusRecR1_S" >> PiShoot_samplelist_condition
echo -e "D21H6piPlusRecR2_SRR1005351\0011D21H6piPlusRecR2_S" >> PiShoot_samplelist_condition
echo -e "D21H6piPlusRecR3_SRR1005352\0011D21H6piPlusRecR3_S" >> PiShoot_samplelist_condition
echo -e "D21H24piPlusR1_SRR1005335\0011D21H24piPlusR1_S" >> PiShoot_samplelist_condition
echo -e "D21H24piPlusR2_SRR1005336\0011D21H24piPlusR2_S" >> PiShoot_samplelist_condition
echo -e "D21H24piPlusR3_SRR1005337\0011D21H24piPlusR3_S" >> PiShoot_samplelist_condition
echo -e "D21H24piMinusR1_SRR1005338\0011D21H24piMinusR1_S" >> PiShoot_samplelist_condition
echo -e "D21H24piMinusR2_SRR1005339\0011D21H24piMinusR2_S" >> PiShoot_samplelist_condition
echo -e "D21H24piMinusR3_SRR1005340\0011D21H24piMinusR3_S" >> PiShoot_samplelist_condition
echo -e "D21H24piPlusRecR1_SRR1005341\0011D21H24piPlusRecR1_S" >> PiShoot_samplelist_condition
echo -e "D21H24piPlusRecR2_SRR1005342\0011D21H24piPlusRecR2_S" >> PiShoot_samplelist_condition
echo -e "D21H24piPlusRecR3_SRR1005343\0011D21H24piPlusRecR3_S" >> PiShoot_samplelist_condition

rm -rf PiShoot
mkdir PiShoot
cd PiShoot
unset arr
k=1
for i in `cat PiShoot_samplelist_condition | cut -f 1`
do
   arr[$k]=$i;
   k=$k+1;
done;
samplenames=${arr[1]}
for i in ${arr[@]:2}  
do
   samplenames=$samplenames","${i};
done
cp /home/project/gid2tid.txt ./
perl /home/pcd/get_sample_gene_TPM_salmon.pl -dir /home/project/salmon/ -sample $samplenames > PiShoot_gene_sample_TPM.xls
perl /home/pcd/filter_ex_by_GTF.pl -GTF /home/project/merged_filter_multiexon_junction40_TPM.gtf -input PiShoot_gene_sample_TPM.xls > PiShoot_gene_sample_TPM_filter.xls
Rscript /home/pcd/arrange_sample_arg.r PiShoot_gene_sample_TPM_filter.xls PiShoot_samplelist_condition
mv result.xls PiShoot_gene_sample_TPM_filter.xls
Rscript /home/pcd/log2_arg.r PiShoot_gene_sample_TPM_filter.xls
mv result.xls PiShoot_gene_sample_log2TPM_filter.xls
perl /home/pcd/get_sample_gene_count_salmon.pl -dir /home/project/salmon/ -sample $samplenames > PiShoot_gene_sample_count.xls
perl /home/pcd/filter_ex_by_GTF.pl -GTF /home/project/merged_filter_multiexon_junction40_TPM.gtf -input PiShoot_gene_sample_count.xls > PiShoot_gene_sample_count_filter.xls
Rscript /home/pcd/arrange_sample_arg.r PiShoot_gene_sample_count_filter.xls PiShoot_samplelist_condition
mv result.xls PiShoot_gene_sample_count_filter.xls
Rscript /home/pcd/deseq2_PiShoot.r

perl /home/pcd/table2table_diffgene.pl -samplepair "H1piMinus_vs_H1piPlus_S,H6piMinus_vs_H6piPlus_S,H24piMinus_vs_H24piPlus_S,D3piMinus_vs_D3piPlus_S,D7piMinus_vs_D7piPlus_S,D21piMinus_vs_D21piPlus_S,D21H1piMinus_vs_D21H1piPlus_S,D21H6piMinus_vs_D21H6piPlus_S,D21H24piMinus_vs_D21H24piPlus_S" > PiShoot_DEG_overlap.xls
perl /home/pcd/table2table_diffgene_nx.pl -rawcount "PiShoot_gene_sample_count.xls" -normcount "PiShoot_gene_sample_normalized_count.xls" -TPM "PiShoot_gene_sample_TPM.xls" -target PiShoot_DEG_overlap.xls -samplepair "H1piMinus_vs_H1piPlus_S,H6piMinus_vs_H6piPlus_S,H24piMinus_vs_H24piPlus_S,D3piMinus_vs_D3piPlus_S,D7piMinus_vs_D7piPlus_S,D21piMinus_vs_D21piPlus_S,D21H1piMinus_vs_D21H1piPlus_S,D21H6piMinus_vs_D21H6piPlus_S,D21H24piMinus_vs_D21H24piPlus_S" > PiShoot_DEG_overlap_info.xls


perl /home/pcd/add_anno_DEG.pl -degene D21piMinus_vs_D21piPlus_S_gene_exp_significant.xls -annotation /home/project/merged_filter_multiexon_junction40_TPM.gtf > D21piMinus_vs_D21piPlus_S_gene_exp_significant_with_annotation.xls



# GO and KEGG enrichment analysis of DEGs under Pi deficiency in shoots
perl /home/pcd/degene_for_DEG_v2.pl -degene PiShoot_DEG_overlap_info_with_annotation.xls > diffgene.txt
perl /home/pcd/KEGG_enrich.pl -ispvalue "Pvalue" -bl /home/database/allgene.txt -gl diffgene.txt -symbol2KEGG /home/database/symbol2KEGG.txt
perl /home/pcd/GO_enrich.pl -ispvalue "Pvalue" -bl /home/database/allgene.txt -gl diffgene.txt -symbol2GO /home/database/symbol2GO.txt

perl /home/pcd/degene_for_DEG.pl -degene D21piMinus_vs_D21piPlus_S_gene_exp_significant_with_annotation.xls > diffgene.txt
perl /home/pcd/KEGG_enrich.pl -ispvalue "Pvalue" -bl /home/database/allgene.txt -gl diffgene.txt -symbol2KEGG /home/database/symbol2KEGG.txt
perl /home/pcd/GO_enrich.pl -ispvalue "Pvalue" -bl /home/database/allgene.txt -gl diffgene.txt -symbol2GO /home/database/symbol2GO.txt
mv GO_enrichment.xls D21_S_GO_enrichment.xls
mv GO_enrichment_significant.xls D21_S_GO_enrichment_significant.xls
mv KEGG_enrichment.xls D21_S_KEGG_enrichment.xls
mv KEGG_enrichment_significant.xls D21_S_KEGG_enrichment_significant.xls





# Gene and transcript quantification of all 141 samples

cat FeZnCuMn_samplelist_condition PiRoot_samplelist_condition PiShoot_samplelist_condition > PiFeZnCuMn_samplelist_condition

rm -rf FeZnCuMnPiRootShoot
mkdir FeZnCuMnPiRootShoot
cd FeZnCuMnPiRootShoot
unset arr　
k=1
for i in `cat PiFeZnCuMn_samplelist_condition | cut -f 1`
do
   #echo $i
   arr[$k]=$i;
   k=$k+1;
done;
samplenames=${arr[1]}
for i in ${arr[@]:2}  
do
   #echo ${i}
   samplenames=$samplenames","${i};
done
cp /home/project/gid2tid.txt ./
perl /home/pcd/get_sample_gene_TPM_salmon.pl  -dir /home/project/salmon/ -sample $samplenames > PiFeZnCuMn_gene_sample_TPM.xls
perl /home/pcd/filter_ex_by_GTF.pl -GTF /home/project/merged_filter_multiexon_junction40_TPM.gtf -input PiFeZnCuMn_gene_sample_TPM.xls > PiFeZnCuMn_gene_sample_TPM_filter.xls
Rscript /home/pcd/arrange_sample_arg.r PiFeZnCuMn_gene_sample_TPM_filter.xls PiFeZnCuMn_samplelist_condition
mv result.xls PiFeZnCuMn_gene_sample_TPM_filter.xls
Rscript /home/pcd/log2_arg.r PiFeZnCuMn_gene_sample_TPM_filter.xls
mv result.xls PiFeZnCuMn_gene_sample_log2TPM_filter.xls
perl /home/pcd/get_sample_gene_count_salmon.pl -dir /home/project/salmon/ -sample $samplenames > PiFeZnCuMn_gene_sample_count.xls
perl /home/pcd/filter_ex_by_GTF.pl -GTF /home/project/merged_filter_multiexon_junction40_TPM.gtf -input PiFeZnCuMn_gene_sample_count.xls > PiFeZnCuMn_gene_sample_count_filter.xls
Rscript /home/pcd/arrange_sample_arg.r PiFeZnCuMn_gene_sample_count_filter.xls PiFeZnCuMn_samplelist_condition
mv result.xls PiFeZnCuMn_gene_sample_count_filter.xls

perl /home/pcd/get_sample_transcript_TPM_salmon.pl -dir /home/project/salmon/ -sample $samplenames > PiFeZnCuMn_transcript_sample_TPM.xls
perl /home/pcd/filter_ex_by_GTF_tid.pl -GTF /home/project/merged_filter_multiexon_junction40_TPM.gtf -input PiFeZnCuMn_transcript_sample_TPM.xls > PiFeZnCuMn_transcript_sample_TPM_filter.xls
Rscript /home/pcd/arrange_sample_arg.r PiFeZnCuMn_transcript_sample_TPM_filter.xls PiFeZnCuMn_samplelist_condition
mv result.xls PiFeZnCuMn_transcript_sample_TPM_filter.xls
Rscript /home/pcd/log2_arg.r PiFeZnCuMn_transcript_sample_TPM_filter.xls
mv result.xls PiFeZnCuMn_transcript_sample_TPM_filter.xls
