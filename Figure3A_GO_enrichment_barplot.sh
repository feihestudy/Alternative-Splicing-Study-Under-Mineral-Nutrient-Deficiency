#1. GO enrichment analysis of DEGs and DASGs under Fe deficiency
Rscript /home/pcd/GO_Select_BP.r /home/project/DEG/FeZnCuMn/Fe/GO_enrichment.xls
perl /home/pcd/GO_remove_redundant.pl -GO GO_enrichment_BP.xls > Fe_DEG_all_GO_BP_simple.txt

Rscript /home/pcd/GO_Select_BP.r /home/project/DASG/Fe_vs_Control/GO_enrichment.xls
perl /home/pcd/GO_remove_redundant.pl -GO GO_enrichment_BP.xls > Fe_DASG_all_GO_BP_simple.txt

perl /home/pcd/GOBarPlot_top.pl -tissue "Fe"
mv GO_PVs_tmp.xls GO_PVs_Fe_top.xls
Rscript /home/pcd/GOBarPlot.r GO_PVs_Fe_top.xls
mv Plot.pdf Figure3A.pdf





#2. GO enrichment analysis of DEGs and DASGs under Pi deficiency at 21 days
Rscript /home/pcd/GO_Select_BP.r D21_R_GO_enrichment.xls
perl /home/pcd/GO_remove_redundant.pl -GO D21_R_GO_enrichment_BP.xls > Pi21DRoot_DEG_all_GO_BP_simple.txt

Rscript /home/pcd/GO_Select_BP.r /home/project/DASG/Root_D21piMinus_vs_D21piPlus/GO_enrichment.xls
perl /home/pcd/GO_remove_redundant.pl -GO GO_enrichment_BP.xls > Pi21DRoot_DASG_all_GO_BP_simple.txt

perl /home/pcd/GOBarPlot_top.pl -tissue "Pi21DRoot"
mv GO_PVs_tmp.xls GO_PVs_Pi21DRoot_top.xls
Rscript /home/pcd/GOBarPlot.r GO_PVs_Pi21DRoot_top.xls
mv Plot.pdf Figure3B.pdf






#3. GO enrichment analysis of DEGs and DASGs under Zn deficiency
Rscript /home/pcd/GO_Select_BP.r /home/project/DEG/FeZnCuMn/Zn/GO_enrichment.xls
perl /home/pcd/GO_remove_redundant.pl -GO GO_enrichment_BP.xls > Zn_DEG_all_GO_BP_simple.txt

Rscript /home/pcd/GO_Select_BP.r /home/project/DASG/Zn_vs_Control/GO_enrichment.xls
perl /home/pcd/GO_remove_redundant.pl -GO GO_enrichment_BP.xls > Zn_DASG_all_GO_BP_simple.txt

perl /home/pcd/GOBarPlot_top.pl -tissue "Zn"
mv GO_PVs_tmp.xls GO_PVs_Zn_top.xls
Rscript /home/pcd/GOBarPlot.r GO_PVs_Zn_top.xls
mv Plot.pdf SupplementalFigure5.pdf





#4. GO enrichment analysis of DEGs and DASGs under Cu deficiency
Rscript /home/pcd/GO_Select_BP.r /home/project/DEG/FeZnCuMn/Cu/GO_enrichment.xls
perl /home/pcd/GO_remove_redundant.pl -GO GO_enrichment_BP.xls > Cu_DEG_all_GO_BP_simple.txt

Rscript /home/pcd/GO_Select_BP.r /home/project/DASG/Cu_vs_Control/GO_enrichment.xls
perl /home/pcd/GO_remove_redundant.pl -GO GO_enrichment_BP.xls > Cu_DASG_all_GO_BP_simple.txt

perl /home/pcd/GOBarPlot_top.pl -tissue "Cu"
mv GO_PVs_tmp.xls GO_PVs_Cu_top.xls
Rscript /home/pcd/GOBarPlot.r GO_PVs_Cu_top.xls
mv Plot.pdf SupplementalFigure6.pdf





#5. GO enrichment analysis of DEGs and DASGs under Mn deficiency
Rscript /home/pcd/GO_Select_BP.r /home/project/DEG/FeZnCuMn/Mn/GO_enrichment.xls
perl /home/pcd/GO_remove_redundant.pl -GO GO_enrichment_BP.xls > Mn_DEG_all_GO_BP_simple.txt

Rscript /home/pcd/GO_Select_BP.r /home/project/DASG/Mn_vs_Control/GO_enrichment.xls
perl /home/pcd/GO_remove_redundant.pl -GO GO_enrichment_BP.xls > Mn_DASG_all_GO_BP_simple.txt

perl /home/pcd/GOBarPlot_top.pl -tissue "Mn"
mv GO_PVs_tmp.xls GO_PVs_Mn_top.xls
Rscript /home/pcd/GOBarPlot.r GO_PVs_Mn_top.xls
mv Plot.pdf SupplementalFigure7.pdf

