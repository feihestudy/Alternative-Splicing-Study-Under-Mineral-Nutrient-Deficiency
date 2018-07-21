# Identify differential alternative splicing events under micronutrient deficiency
mkdir DASG
cd DASG
rm -rf Fe_vs_Control
mkdir Fe_vs_Control
cd Fe_vs_Control
Rscript /home/pcd/merge_ASE_for_diff.r "ControlR1,ControlR2,ControlR3" "FeR1,FeR2,FeR3"
awk '{ split($2,ic_1,","); sum_ic_1=0; for (x in ic_1) sum_ic_1 += ic_1[x]; split($4,ic_2,","); sum_ic_2=0; for (x in ic_2) sum_ic_2 += ic_2[x]; split($3,sc_1,","); sum_sc_1=0; for (x in sc_1) sum_sc_1 += sc_1[x]; split($5,sc_2,","); sum_sc_2=0; for (x in sc_2) sum_sc_2 += sc_2[x]; if ( NR==1 || ( (sum_ic_2 + sum_sc_2 > 0) && (sum_ic_1 + sum_sc_1 > 0) && (sum_ic_1 != 0 || sum_ic_2 != 0) && (sum_sc_1 != 0 || sum_sc_2 != 0) && $6!=0 && $7!=0 ) ) {print $0}}' JCEC.RNASeq.RI.MATS.input.txt > RI.JCEC.input.txt
/home/soft/rMATS.3.2.5/MATS/rMATS.sh -d RI.JCEC.input.txt -o RIout_JCEC/ -c 0.0001 -p 8 -t U
awk '{ split($2,ic_1,","); sum_ic_1=0; for (x in ic_1) sum_ic_1 += ic_1[x]; split($4,ic_2,","); sum_ic_2=0; for (x in ic_2) sum_ic_2 += ic_2[x]; split($3,sc_1,","); sum_sc_1=0; for (x in sc_1) sum_sc_1 += sc_1[x]; split($5,sc_2,","); sum_sc_2=0; for (x in sc_2) sum_sc_2 += sc_2[x]; if ( NR==1 || ( (sum_ic_2 + sum_sc_2 > 0) && (sum_ic_1 + sum_sc_1 > 0) && (sum_ic_1 != 0 || sum_ic_2 != 0) && (sum_sc_1 != 0 || sum_sc_2 != 0) && $6!=0 && $7!=0 ) ) {print $0}}' JCEC.RNASeq.SE.MATS.input.txt > SE.JCEC.input.txt
/home/soft/rMATS.3.2.5/MATS/rMATS.sh -d SE.JCEC.input.txt -o SEout_JCEC/ -c 0.0001 -p 8 -t U
awk '{ split($2,ic_1,","); sum_ic_1=0; for (x in ic_1) sum_ic_1 += ic_1[x]; split($4,ic_2,","); sum_ic_2=0; for (x in ic_2) sum_ic_2 += ic_2[x]; split($3,sc_1,","); sum_sc_1=0; for (x in sc_1) sum_sc_1 += sc_1[x]; split($5,sc_2,","); sum_sc_2=0; for (x in sc_2) sum_sc_2 += sc_2[x]; if ( NR==1 || ( (sum_ic_2 + sum_sc_2 > 0) && (sum_ic_1 + sum_sc_1 > 0) && (sum_ic_1 != 0 || sum_ic_2 != 0) && (sum_sc_1 != 0 || sum_sc_2 != 0) && $6!=0 && $7!=0 ) ) {print $0}}' JCEC.RNASeq.MXE.MATS.input.txt > MXE.JCEC.input.txt
/home/soft/rMATS.3.2.5/MATS/rMATS.sh -d MXE.JCEC.input.txt -o MXEout_JCEC/ -c 0.0001 -p 8 -t U
awk '{ split($2,ic_1,","); sum_ic_1=0; for (x in ic_1) sum_ic_1 += ic_1[x]; split($4,ic_2,","); sum_ic_2=0; for (x in ic_2) sum_ic_2 += ic_2[x]; split($3,sc_1,","); sum_sc_1=0; for (x in sc_1) sum_sc_1 += sc_1[x]; split($5,sc_2,","); sum_sc_2=0; for (x in sc_2) sum_sc_2 += sc_2[x]; if ( NR==1 || ( (sum_ic_2 + sum_sc_2 > 0) && (sum_ic_1 + sum_sc_1 > 0) && (sum_ic_1 != 0 || sum_ic_2 != 0) && (sum_sc_1 != 0 || sum_sc_2 != 0) && $6!=0 && $7!=0 ) ) {print $0}}' JCEC.RNASeq.A5SS.MATS.input.txt > A5SS.JCEC.input.txt
/home/soft/rMATS.3.2.5/MATS/rMATS.sh -d A5SS.JCEC.input.txt -o A5SSout_JCEC/ -c 0.0001 -p 8 -t U
awk '{ split($2,ic_1,","); sum_ic_1=0; for (x in ic_1) sum_ic_1 += ic_1[x]; split($4,ic_2,","); sum_ic_2=0; for (x in ic_2) sum_ic_2 += ic_2[x]; split($3,sc_1,","); sum_sc_1=0; for (x in sc_1) sum_sc_1 += sc_1[x]; split($5,sc_2,","); sum_sc_2=0; for (x in sc_2) sum_sc_2 += sc_2[x]; if ( NR==1 || ( (sum_ic_2 + sum_sc_2 > 0) && (sum_ic_1 + sum_sc_1 > 0) && (sum_ic_1 != 0 || sum_ic_2 != 0) && (sum_sc_1 != 0 || sum_sc_2 != 0) && $6!=0 && $7!=0 ) ) {print $0}}' JCEC.RNASeq.A3SS.MATS.input.txt > A3SS.JCEC.input.txt
/home/soft/rMATS.3.2.5/MATS/rMATS.sh -d A3SS.JCEC.input.txt -o A3SSout_JCEC/ -c 0.0001 -p 8 -t U
cd ..
rm -rf Zn_vs_Control
mkdir Zn_vs_Control
cd Zn_vs_Control
Rscript /home/pcd/merge_ASE_for_diff.r "ControlR1,ControlR2,ControlR3" "ZnR1,ZnR2,ZnR3"
awk '{ split($2,ic_1,","); sum_ic_1=0; for (x in ic_1) sum_ic_1 += ic_1[x]; split($4,ic_2,","); sum_ic_2=0; for (x in ic_2) sum_ic_2 += ic_2[x]; split($3,sc_1,","); sum_sc_1=0; for (x in sc_1) sum_sc_1 += sc_1[x]; split($5,sc_2,","); sum_sc_2=0; for (x in sc_2) sum_sc_2 += sc_2[x]; if ( NR==1 || ( (sum_ic_2 + sum_sc_2 > 0) && (sum_ic_1 + sum_sc_1 > 0) && (sum_ic_1 != 0 || sum_ic_2 != 0) && (sum_sc_1 != 0 || sum_sc_2 != 0) && $6!=0 && $7!=0 ) ) {print $0}}' JCEC.RNASeq.RI.MATS.input.txt > RI.JCEC.input.txt
/home/soft/rMATS.3.2.5/MATS/rMATS.sh -d RI.JCEC.input.txt -o RIout_JCEC/ -c 0.0001 -p 8 -t U
awk '{ split($2,ic_1,","); sum_ic_1=0; for (x in ic_1) sum_ic_1 += ic_1[x]; split($4,ic_2,","); sum_ic_2=0; for (x in ic_2) sum_ic_2 += ic_2[x]; split($3,sc_1,","); sum_sc_1=0; for (x in sc_1) sum_sc_1 += sc_1[x]; split($5,sc_2,","); sum_sc_2=0; for (x in sc_2) sum_sc_2 += sc_2[x]; if ( NR==1 || ( (sum_ic_2 + sum_sc_2 > 0) && (sum_ic_1 + sum_sc_1 > 0) && (sum_ic_1 != 0 || sum_ic_2 != 0) && (sum_sc_1 != 0 || sum_sc_2 != 0) && $6!=0 && $7!=0 ) ) {print $0}}' JCEC.RNASeq.SE.MATS.input.txt > SE.JCEC.input.txt
/home/soft/rMATS.3.2.5/MATS/rMATS.sh -d SE.JCEC.input.txt -o SEout_JCEC/ -c 0.0001 -p 8 -t U
awk '{ split($2,ic_1,","); sum_ic_1=0; for (x in ic_1) sum_ic_1 += ic_1[x]; split($4,ic_2,","); sum_ic_2=0; for (x in ic_2) sum_ic_2 += ic_2[x]; split($3,sc_1,","); sum_sc_1=0; for (x in sc_1) sum_sc_1 += sc_1[x]; split($5,sc_2,","); sum_sc_2=0; for (x in sc_2) sum_sc_2 += sc_2[x]; if ( NR==1 || ( (sum_ic_2 + sum_sc_2 > 0) && (sum_ic_1 + sum_sc_1 > 0) && (sum_ic_1 != 0 || sum_ic_2 != 0) && (sum_sc_1 != 0 || sum_sc_2 != 0) && $6!=0 && $7!=0 ) ) {print $0}}' JCEC.RNASeq.MXE.MATS.input.txt > MXE.JCEC.input.txt
/home/soft/rMATS.3.2.5/MATS/rMATS.sh -d MXE.JCEC.input.txt -o MXEout_JCEC/ -c 0.0001 -p 8 -t U
awk '{ split($2,ic_1,","); sum_ic_1=0; for (x in ic_1) sum_ic_1 += ic_1[x]; split($4,ic_2,","); sum_ic_2=0; for (x in ic_2) sum_ic_2 += ic_2[x]; split($3,sc_1,","); sum_sc_1=0; for (x in sc_1) sum_sc_1 += sc_1[x]; split($5,sc_2,","); sum_sc_2=0; for (x in sc_2) sum_sc_2 += sc_2[x]; if ( NR==1 || ( (sum_ic_2 + sum_sc_2 > 0) && (sum_ic_1 + sum_sc_1 > 0) && (sum_ic_1 != 0 || sum_ic_2 != 0) && (sum_sc_1 != 0 || sum_sc_2 != 0) && $6!=0 && $7!=0 ) ) {print $0}}' JCEC.RNASeq.A5SS.MATS.input.txt > A5SS.JCEC.input.txt
/home/soft/rMATS.3.2.5/MATS/rMATS.sh -d A5SS.JCEC.input.txt -o A5SSout_JCEC/ -c 0.0001 -p 8 -t U
awk '{ split($2,ic_1,","); sum_ic_1=0; for (x in ic_1) sum_ic_1 += ic_1[x]; split($4,ic_2,","); sum_ic_2=0; for (x in ic_2) sum_ic_2 += ic_2[x]; split($3,sc_1,","); sum_sc_1=0; for (x in sc_1) sum_sc_1 += sc_1[x]; split($5,sc_2,","); sum_sc_2=0; for (x in sc_2) sum_sc_2 += sc_2[x]; if ( NR==1 || ( (sum_ic_2 + sum_sc_2 > 0) && (sum_ic_1 + sum_sc_1 > 0) && (sum_ic_1 != 0 || sum_ic_2 != 0) && (sum_sc_1 != 0 || sum_sc_2 != 0) && $6!=0 && $7!=0 ) ) {print $0}}' JCEC.RNASeq.A3SS.MATS.input.txt > A3SS.JCEC.input.txt
/home/soft/rMATS.3.2.5/MATS/rMATS.sh -d A3SS.JCEC.input.txt -o A3SSout_JCEC/ -c 0.0001 -p 8 -t U
cd ..
rm -rf Cu_vs_Control
mkdir Cu_vs_Control
cd Cu_vs_Control
Rscript /home/pcd/merge_ASE_for_diff.r "ControlR1,ControlR2,ControlR3" "CuR1,CuR2,CuR3"
awk '{ split($2,ic_1,","); sum_ic_1=0; for (x in ic_1) sum_ic_1 += ic_1[x]; split($4,ic_2,","); sum_ic_2=0; for (x in ic_2) sum_ic_2 += ic_2[x]; split($3,sc_1,","); sum_sc_1=0; for (x in sc_1) sum_sc_1 += sc_1[x]; split($5,sc_2,","); sum_sc_2=0; for (x in sc_2) sum_sc_2 += sc_2[x]; if ( NR==1 || ( (sum_ic_2 + sum_sc_2 > 0) && (sum_ic_1 + sum_sc_1 > 0) && (sum_ic_1 != 0 || sum_ic_2 != 0) && (sum_sc_1 != 0 || sum_sc_2 != 0) && $6!=0 && $7!=0 ) ) {print $0}}' JCEC.RNASeq.RI.MATS.input.txt > RI.JCEC.input.txt
/home/soft/rMATS.3.2.5/MATS/rMATS.sh -d RI.JCEC.input.txt -o RIout_JCEC/ -c 0.0001 -p 8 -t U
awk '{ split($2,ic_1,","); sum_ic_1=0; for (x in ic_1) sum_ic_1 += ic_1[x]; split($4,ic_2,","); sum_ic_2=0; for (x in ic_2) sum_ic_2 += ic_2[x]; split($3,sc_1,","); sum_sc_1=0; for (x in sc_1) sum_sc_1 += sc_1[x]; split($5,sc_2,","); sum_sc_2=0; for (x in sc_2) sum_sc_2 += sc_2[x]; if ( NR==1 || ( (sum_ic_2 + sum_sc_2 > 0) && (sum_ic_1 + sum_sc_1 > 0) && (sum_ic_1 != 0 || sum_ic_2 != 0) && (sum_sc_1 != 0 || sum_sc_2 != 0) && $6!=0 && $7!=0 ) ) {print $0}}' JCEC.RNASeq.SE.MATS.input.txt > SE.JCEC.input.txt
/home/soft/rMATS.3.2.5/MATS/rMATS.sh -d SE.JCEC.input.txt -o SEout_JCEC/ -c 0.0001 -p 8 -t U
awk '{ split($2,ic_1,","); sum_ic_1=0; for (x in ic_1) sum_ic_1 += ic_1[x]; split($4,ic_2,","); sum_ic_2=0; for (x in ic_2) sum_ic_2 += ic_2[x]; split($3,sc_1,","); sum_sc_1=0; for (x in sc_1) sum_sc_1 += sc_1[x]; split($5,sc_2,","); sum_sc_2=0; for (x in sc_2) sum_sc_2 += sc_2[x]; if ( NR==1 || ( (sum_ic_2 + sum_sc_2 > 0) && (sum_ic_1 + sum_sc_1 > 0) && (sum_ic_1 != 0 || sum_ic_2 != 0) && (sum_sc_1 != 0 || sum_sc_2 != 0) && $6!=0 && $7!=0 ) ) {print $0}}' JCEC.RNASeq.MXE.MATS.input.txt > MXE.JCEC.input.txt
/home/soft/rMATS.3.2.5/MATS/rMATS.sh -d MXE.JCEC.input.txt -o MXEout_JCEC/ -c 0.0001 -p 8 -t U
awk '{ split($2,ic_1,","); sum_ic_1=0; for (x in ic_1) sum_ic_1 += ic_1[x]; split($4,ic_2,","); sum_ic_2=0; for (x in ic_2) sum_ic_2 += ic_2[x]; split($3,sc_1,","); sum_sc_1=0; for (x in sc_1) sum_sc_1 += sc_1[x]; split($5,sc_2,","); sum_sc_2=0; for (x in sc_2) sum_sc_2 += sc_2[x]; if ( NR==1 || ( (sum_ic_2 + sum_sc_2 > 0) && (sum_ic_1 + sum_sc_1 > 0) && (sum_ic_1 != 0 || sum_ic_2 != 0) && (sum_sc_1 != 0 || sum_sc_2 != 0) && $6!=0 && $7!=0 ) ) {print $0}}' JCEC.RNASeq.A5SS.MATS.input.txt > A5SS.JCEC.input.txt
/home/soft/rMATS.3.2.5/MATS/rMATS.sh -d A5SS.JCEC.input.txt -o A5SSout_JCEC/ -c 0.0001 -p 8 -t U
awk '{ split($2,ic_1,","); sum_ic_1=0; for (x in ic_1) sum_ic_1 += ic_1[x]; split($4,ic_2,","); sum_ic_2=0; for (x in ic_2) sum_ic_2 += ic_2[x]; split($3,sc_1,","); sum_sc_1=0; for (x in sc_1) sum_sc_1 += sc_1[x]; split($5,sc_2,","); sum_sc_2=0; for (x in sc_2) sum_sc_2 += sc_2[x]; if ( NR==1 || ( (sum_ic_2 + sum_sc_2 > 0) && (sum_ic_1 + sum_sc_1 > 0) && (sum_ic_1 != 0 || sum_ic_2 != 0) && (sum_sc_1 != 0 || sum_sc_2 != 0) && $6!=0 && $7!=0 ) ) {print $0}}' JCEC.RNASeq.A3SS.MATS.input.txt > A3SS.JCEC.input.txt
/home/soft/rMATS.3.2.5/MATS/rMATS.sh -d A3SS.JCEC.input.txt -o A3SSout_JCEC/ -c 0.0001 -p 8 -t U
cd ..
rm -rf Mn_vs_Control
mkdir Mn_vs_Control
cd Mn_vs_Control
Rscript /home/pcd/merge_ASE_for_diff.r "ControlR1,ControlR2,ControlR3" "MnR1,MnR2,MnR3"
awk '{ split($2,ic_1,","); sum_ic_1=0; for (x in ic_1) sum_ic_1 += ic_1[x]; split($4,ic_2,","); sum_ic_2=0; for (x in ic_2) sum_ic_2 += ic_2[x]; split($3,sc_1,","); sum_sc_1=0; for (x in sc_1) sum_sc_1 += sc_1[x]; split($5,sc_2,","); sum_sc_2=0; for (x in sc_2) sum_sc_2 += sc_2[x]; if ( NR==1 || ( (sum_ic_2 + sum_sc_2 > 0) && (sum_ic_1 + sum_sc_1 > 0) && (sum_ic_1 != 0 || sum_ic_2 != 0) && (sum_sc_1 != 0 || sum_sc_2 != 0) && $6!=0 && $7!=0 ) ) {print $0}}' JCEC.RNASeq.RI.MATS.input.txt > RI.JCEC.input.txt
/home/soft/rMATS.3.2.5/MATS/rMATS.sh -d RI.JCEC.input.txt -o RIout_JCEC/ -c 0.0001 -p 8 -t U
awk '{ split($2,ic_1,","); sum_ic_1=0; for (x in ic_1) sum_ic_1 += ic_1[x]; split($4,ic_2,","); sum_ic_2=0; for (x in ic_2) sum_ic_2 += ic_2[x]; split($3,sc_1,","); sum_sc_1=0; for (x in sc_1) sum_sc_1 += sc_1[x]; split($5,sc_2,","); sum_sc_2=0; for (x in sc_2) sum_sc_2 += sc_2[x]; if ( NR==1 || ( (sum_ic_2 + sum_sc_2 > 0) && (sum_ic_1 + sum_sc_1 > 0) && (sum_ic_1 != 0 || sum_ic_2 != 0) && (sum_sc_1 != 0 || sum_sc_2 != 0) && $6!=0 && $7!=0 ) ) {print $0}}' JCEC.RNASeq.SE.MATS.input.txt > SE.JCEC.input.txt
/home/soft/rMATS.3.2.5/MATS/rMATS.sh -d SE.JCEC.input.txt -o SEout_JCEC/ -c 0.0001 -p 8 -t U
awk '{ split($2,ic_1,","); sum_ic_1=0; for (x in ic_1) sum_ic_1 += ic_1[x]; split($4,ic_2,","); sum_ic_2=0; for (x in ic_2) sum_ic_2 += ic_2[x]; split($3,sc_1,","); sum_sc_1=0; for (x in sc_1) sum_sc_1 += sc_1[x]; split($5,sc_2,","); sum_sc_2=0; for (x in sc_2) sum_sc_2 += sc_2[x]; if ( NR==1 || ( (sum_ic_2 + sum_sc_2 > 0) && (sum_ic_1 + sum_sc_1 > 0) && (sum_ic_1 != 0 || sum_ic_2 != 0) && (sum_sc_1 != 0 || sum_sc_2 != 0) && $6!=0 && $7!=0 ) ) {print $0}}' JCEC.RNASeq.MXE.MATS.input.txt > MXE.JCEC.input.txt
/home/soft/rMATS.3.2.5/MATS/rMATS.sh -d MXE.JCEC.input.txt -o MXEout_JCEC/ -c 0.0001 -p 8 -t U
awk '{ split($2,ic_1,","); sum_ic_1=0; for (x in ic_1) sum_ic_1 += ic_1[x]; split($4,ic_2,","); sum_ic_2=0; for (x in ic_2) sum_ic_2 += ic_2[x]; split($3,sc_1,","); sum_sc_1=0; for (x in sc_1) sum_sc_1 += sc_1[x]; split($5,sc_2,","); sum_sc_2=0; for (x in sc_2) sum_sc_2 += sc_2[x]; if ( NR==1 || ( (sum_ic_2 + sum_sc_2 > 0) && (sum_ic_1 + sum_sc_1 > 0) && (sum_ic_1 != 0 || sum_ic_2 != 0) && (sum_sc_1 != 0 || sum_sc_2 != 0) && $6!=0 && $7!=0 ) ) {print $0}}' JCEC.RNASeq.A5SS.MATS.input.txt > A5SS.JCEC.input.txt
/home/soft/rMATS.3.2.5/MATS/rMATS.sh -d A5SS.JCEC.input.txt -o A5SSout_JCEC/ -c 0.0001 -p 8 -t U
awk '{ split($2,ic_1,","); sum_ic_1=0; for (x in ic_1) sum_ic_1 += ic_1[x]; split($4,ic_2,","); sum_ic_2=0; for (x in ic_2) sum_ic_2 += ic_2[x]; split($3,sc_1,","); sum_sc_1=0; for (x in sc_1) sum_sc_1 += sc_1[x]; split($5,sc_2,","); sum_sc_2=0; for (x in sc_2) sum_sc_2 += sc_2[x]; if ( NR==1 || ( (sum_ic_2 + sum_sc_2 > 0) && (sum_ic_1 + sum_sc_1 > 0) && (sum_ic_1 != 0 || sum_ic_2 != 0) && (sum_sc_1 != 0 || sum_sc_2 != 0) && $6!=0 && $7!=0 ) ) {print $0}}' JCEC.RNASeq.A3SS.MATS.input.txt > A3SS.JCEC.input.txt
/home/soft/rMATS.3.2.5/MATS/rMATS.sh -d A3SS.JCEC.input.txt -o A3SSout_JCEC/ -c 0.0001 -p 8 -t U
cd ..



# Identify differential alternative splicing events under Pi deficiency in roots
rm -rf Root_H1piMinus_vs_H1piPlus
mkdir Root_H1piMinus_vs_H1piPlus
cd Root_H1piMinus_vs_H1piPlus
Rscript /home/pcd/merge_ASE_for_diff.r "H1piPlusR1_SRR1005257,H1piPlusR2_SRR1005258,H1piPlusR3_SRR1005259" "H1piMinusR1_SRR1005260,H1piMinusR2_SRR1005261,H1piMinusR3_SRR1005262"
awk '{ split($2,ic_1,","); sum_ic_1=0; for (x in ic_1) sum_ic_1 += ic_1[x]; split($4,ic_2,","); sum_ic_2=0; for (x in ic_2) sum_ic_2 += ic_2[x]; split($3,sc_1,","); sum_sc_1=0; for (x in sc_1) sum_sc_1 += sc_1[x]; split($5,sc_2,","); sum_sc_2=0; for (x in sc_2) sum_sc_2 += sc_2[x]; if ( NR==1 || ( (sum_ic_2 + sum_sc_2 > 0) && (sum_ic_1 + sum_sc_1 > 0) && (sum_ic_1 != 0 || sum_ic_2 != 0) && (sum_sc_1 != 0 || sum_sc_2 != 0) && $6!=0 && $7!=0 ) ) {print $0}}' JCEC.RNASeq.RI.MATS.input.txt > RI.JCEC.input.txt
/home/soft/rMATS.3.2.5/MATS/rMATS.sh -d RI.JCEC.input.txt -o RIout_JCEC/ -c 0.0001 -p 8 -t U
awk '{ split($2,ic_1,","); sum_ic_1=0; for (x in ic_1) sum_ic_1 += ic_1[x]; split($4,ic_2,","); sum_ic_2=0; for (x in ic_2) sum_ic_2 += ic_2[x]; split($3,sc_1,","); sum_sc_1=0; for (x in sc_1) sum_sc_1 += sc_1[x]; split($5,sc_2,","); sum_sc_2=0; for (x in sc_2) sum_sc_2 += sc_2[x]; if ( NR==1 || ( (sum_ic_2 + sum_sc_2 > 0) && (sum_ic_1 + sum_sc_1 > 0) && (sum_ic_1 != 0 || sum_ic_2 != 0) && (sum_sc_1 != 0 || sum_sc_2 != 0) && $6!=0 && $7!=0 ) ) {print $0}}' JCEC.RNASeq.SE.MATS.input.txt > SE.JCEC.input.txt
/home/soft/rMATS.3.2.5/MATS/rMATS.sh -d SE.JCEC.input.txt -o SEout_JCEC/ -c 0.0001 -p 8 -t U
awk '{ split($2,ic_1,","); sum_ic_1=0; for (x in ic_1) sum_ic_1 += ic_1[x]; split($4,ic_2,","); sum_ic_2=0; for (x in ic_2) sum_ic_2 += ic_2[x]; split($3,sc_1,","); sum_sc_1=0; for (x in sc_1) sum_sc_1 += sc_1[x]; split($5,sc_2,","); sum_sc_2=0; for (x in sc_2) sum_sc_2 += sc_2[x]; if ( NR==1 || ( (sum_ic_2 + sum_sc_2 > 0) && (sum_ic_1 + sum_sc_1 > 0) && (sum_ic_1 != 0 || sum_ic_2 != 0) && (sum_sc_1 != 0 || sum_sc_2 != 0) && $6!=0 && $7!=0 ) ) {print $0}}' JCEC.RNASeq.MXE.MATS.input.txt > MXE.JCEC.input.txt
/home/soft/rMATS.3.2.5/MATS/rMATS.sh -d MXE.JCEC.input.txt -o MXEout_JCEC/ -c 0.0001 -p 8 -t U
awk '{ split($2,ic_1,","); sum_ic_1=0; for (x in ic_1) sum_ic_1 += ic_1[x]; split($4,ic_2,","); sum_ic_2=0; for (x in ic_2) sum_ic_2 += ic_2[x]; split($3,sc_1,","); sum_sc_1=0; for (x in sc_1) sum_sc_1 += sc_1[x]; split($5,sc_2,","); sum_sc_2=0; for (x in sc_2) sum_sc_2 += sc_2[x]; if ( NR==1 || ( (sum_ic_2 + sum_sc_2 > 0) && (sum_ic_1 + sum_sc_1 > 0) && (sum_ic_1 != 0 || sum_ic_2 != 0) && (sum_sc_1 != 0 || sum_sc_2 != 0) && $6!=0 && $7!=0 ) ) {print $0}}' JCEC.RNASeq.A5SS.MATS.input.txt > A5SS.JCEC.input.txt
/home/soft/rMATS.3.2.5/MATS/rMATS.sh -d A5SS.JCEC.input.txt -o A5SSout_JCEC/ -c 0.0001 -p 8 -t U
awk '{ split($2,ic_1,","); sum_ic_1=0; for (x in ic_1) sum_ic_1 += ic_1[x]; split($4,ic_2,","); sum_ic_2=0; for (x in ic_2) sum_ic_2 += ic_2[x]; split($3,sc_1,","); sum_sc_1=0; for (x in sc_1) sum_sc_1 += sc_1[x]; split($5,sc_2,","); sum_sc_2=0; for (x in sc_2) sum_sc_2 += sc_2[x]; if ( NR==1 || ( (sum_ic_2 + sum_sc_2 > 0) && (sum_ic_1 + sum_sc_1 > 0) && (sum_ic_1 != 0 || sum_ic_2 != 0) && (sum_sc_1 != 0 || sum_sc_2 != 0) && $6!=0 && $7!=0 ) ) {print $0}}' JCEC.RNASeq.A3SS.MATS.input.txt > A3SS.JCEC.input.txt
/home/soft/rMATS.3.2.5/MATS/rMATS.sh -d A3SS.JCEC.input.txt -o A3SSout_JCEC/ -c 0.0001 -p 8 -t U
cd ..
rm -rf Root_H6piMinus_vs_H6piPlus
mkdir Root_H6piMinus_vs_H6piPlus
cd Root_H6piMinus_vs_H6piPlus
Rscript /home/pcd/merge_ASE_for_diff.r "H6piPlusR1_SRR1005308,H6piPlusR2_SRR1005309,H6piPlusR3_SRR1005310" "H6piMinusR1_SRR1005311,H6piMinusR2_SRR1005312,H6piMinusR3_SRR1005313"
awk '{ split($2,ic_1,","); sum_ic_1=0; for (x in ic_1) sum_ic_1 += ic_1[x]; split($4,ic_2,","); sum_ic_2=0; for (x in ic_2) sum_ic_2 += ic_2[x]; split($3,sc_1,","); sum_sc_1=0; for (x in sc_1) sum_sc_1 += sc_1[x]; split($5,sc_2,","); sum_sc_2=0; for (x in sc_2) sum_sc_2 += sc_2[x]; if ( NR==1 || ( (sum_ic_2 + sum_sc_2 > 0) && (sum_ic_1 + sum_sc_1 > 0) && (sum_ic_1 != 0 || sum_ic_2 != 0) && (sum_sc_1 != 0 || sum_sc_2 != 0) && $6!=0 && $7!=0 ) ) {print $0}}' JCEC.RNASeq.RI.MATS.input.txt > RI.JCEC.input.txt
/home/soft/rMATS.3.2.5/MATS/rMATS.sh -d RI.JCEC.input.txt -o RIout_JCEC/ -c 0.0001 -p 8 -t U
awk '{ split($2,ic_1,","); sum_ic_1=0; for (x in ic_1) sum_ic_1 += ic_1[x]; split($4,ic_2,","); sum_ic_2=0; for (x in ic_2) sum_ic_2 += ic_2[x]; split($3,sc_1,","); sum_sc_1=0; for (x in sc_1) sum_sc_1 += sc_1[x]; split($5,sc_2,","); sum_sc_2=0; for (x in sc_2) sum_sc_2 += sc_2[x]; if ( NR==1 || ( (sum_ic_2 + sum_sc_2 > 0) && (sum_ic_1 + sum_sc_1 > 0) && (sum_ic_1 != 0 || sum_ic_2 != 0) && (sum_sc_1 != 0 || sum_sc_2 != 0) && $6!=0 && $7!=0 ) ) {print $0}}' JCEC.RNASeq.SE.MATS.input.txt > SE.JCEC.input.txt
/home/soft/rMATS.3.2.5/MATS/rMATS.sh -d SE.JCEC.input.txt -o SEout_JCEC/ -c 0.0001 -p 8 -t U
awk '{ split($2,ic_1,","); sum_ic_1=0; for (x in ic_1) sum_ic_1 += ic_1[x]; split($4,ic_2,","); sum_ic_2=0; for (x in ic_2) sum_ic_2 += ic_2[x]; split($3,sc_1,","); sum_sc_1=0; for (x in sc_1) sum_sc_1 += sc_1[x]; split($5,sc_2,","); sum_sc_2=0; for (x in sc_2) sum_sc_2 += sc_2[x]; if ( NR==1 || ( (sum_ic_2 + sum_sc_2 > 0) && (sum_ic_1 + sum_sc_1 > 0) && (sum_ic_1 != 0 || sum_ic_2 != 0) && (sum_sc_1 != 0 || sum_sc_2 != 0) && $6!=0 && $7!=0 ) ) {print $0}}' JCEC.RNASeq.MXE.MATS.input.txt > MXE.JCEC.input.txt
/home/soft/rMATS.3.2.5/MATS/rMATS.sh -d MXE.JCEC.input.txt -o MXEout_JCEC/ -c 0.0001 -p 8 -t U
awk '{ split($2,ic_1,","); sum_ic_1=0; for (x in ic_1) sum_ic_1 += ic_1[x]; split($4,ic_2,","); sum_ic_2=0; for (x in ic_2) sum_ic_2 += ic_2[x]; split($3,sc_1,","); sum_sc_1=0; for (x in sc_1) sum_sc_1 += sc_1[x]; split($5,sc_2,","); sum_sc_2=0; for (x in sc_2) sum_sc_2 += sc_2[x]; if ( NR==1 || ( (sum_ic_2 + sum_sc_2 > 0) && (sum_ic_1 + sum_sc_1 > 0) && (sum_ic_1 != 0 || sum_ic_2 != 0) && (sum_sc_1 != 0 || sum_sc_2 != 0) && $6!=0 && $7!=0 ) ) {print $0}}' JCEC.RNASeq.A5SS.MATS.input.txt > A5SS.JCEC.input.txt
/home/soft/rMATS.3.2.5/MATS/rMATS.sh -d A5SS.JCEC.input.txt -o A5SSout_JCEC/ -c 0.0001 -p 8 -t U
awk '{ split($2,ic_1,","); sum_ic_1=0; for (x in ic_1) sum_ic_1 += ic_1[x]; split($4,ic_2,","); sum_ic_2=0; for (x in ic_2) sum_ic_2 += ic_2[x]; split($3,sc_1,","); sum_sc_1=0; for (x in sc_1) sum_sc_1 += sc_1[x]; split($5,sc_2,","); sum_sc_2=0; for (x in sc_2) sum_sc_2 += sc_2[x]; if ( NR==1 || ( (sum_ic_2 + sum_sc_2 > 0) && (sum_ic_1 + sum_sc_1 > 0) && (sum_ic_1 != 0 || sum_ic_2 != 0) && (sum_sc_1 != 0 || sum_sc_2 != 0) && $6!=0 && $7!=0 ) ) {print $0}}' JCEC.RNASeq.A3SS.MATS.input.txt > A3SS.JCEC.input.txt
/home/soft/rMATS.3.2.5/MATS/rMATS.sh -d A3SS.JCEC.input.txt -o A3SSout_JCEC/ -c 0.0001 -p 8 -t U
cd ..
rm -rf Root_H24piMinus_vs_H24piPlus
mkdir Root_H24piMinus_vs_H24piPlus
cd Root_H24piMinus_vs_H24piPlus
Rscript /home/pcd/merge_ASE_for_diff.r "H24piPlusR1_SRR1005296,H24piPlusR2_SRR1005297,H24piPlusR3_SRR1005298" "H24piMinusR1_SRR1005299,H24piMinusR2_SRR1005300,H24piMinusR3_SRR1005301"
awk '{ split($2,ic_1,","); sum_ic_1=0; for (x in ic_1) sum_ic_1 += ic_1[x]; split($4,ic_2,","); sum_ic_2=0; for (x in ic_2) sum_ic_2 += ic_2[x]; split($3,sc_1,","); sum_sc_1=0; for (x in sc_1) sum_sc_1 += sc_1[x]; split($5,sc_2,","); sum_sc_2=0; for (x in sc_2) sum_sc_2 += sc_2[x]; if ( NR==1 || ( (sum_ic_2 + sum_sc_2 > 0) && (sum_ic_1 + sum_sc_1 > 0) && (sum_ic_1 != 0 || sum_ic_2 != 0) && (sum_sc_1 != 0 || sum_sc_2 != 0) && $6!=0 && $7!=0 ) ) {print $0}}' JCEC.RNASeq.RI.MATS.input.txt > RI.JCEC.input.txt
/home/soft/rMATS.3.2.5/MATS/rMATS.sh -d RI.JCEC.input.txt -o RIout_JCEC/ -c 0.0001 -p 8 -t U
awk '{ split($2,ic_1,","); sum_ic_1=0; for (x in ic_1) sum_ic_1 += ic_1[x]; split($4,ic_2,","); sum_ic_2=0; for (x in ic_2) sum_ic_2 += ic_2[x]; split($3,sc_1,","); sum_sc_1=0; for (x in sc_1) sum_sc_1 += sc_1[x]; split($5,sc_2,","); sum_sc_2=0; for (x in sc_2) sum_sc_2 += sc_2[x]; if ( NR==1 || ( (sum_ic_2 + sum_sc_2 > 0) && (sum_ic_1 + sum_sc_1 > 0) && (sum_ic_1 != 0 || sum_ic_2 != 0) && (sum_sc_1 != 0 || sum_sc_2 != 0) && $6!=0 && $7!=0 ) ) {print $0}}' JCEC.RNASeq.SE.MATS.input.txt > SE.JCEC.input.txt
/home/soft/rMATS.3.2.5/MATS/rMATS.sh -d SE.JCEC.input.txt -o SEout_JCEC/ -c 0.0001 -p 8 -t U
awk '{ split($2,ic_1,","); sum_ic_1=0; for (x in ic_1) sum_ic_1 += ic_1[x]; split($4,ic_2,","); sum_ic_2=0; for (x in ic_2) sum_ic_2 += ic_2[x]; split($3,sc_1,","); sum_sc_1=0; for (x in sc_1) sum_sc_1 += sc_1[x]; split($5,sc_2,","); sum_sc_2=0; for (x in sc_2) sum_sc_2 += sc_2[x]; if ( NR==1 || ( (sum_ic_2 + sum_sc_2 > 0) && (sum_ic_1 + sum_sc_1 > 0) && (sum_ic_1 != 0 || sum_ic_2 != 0) && (sum_sc_1 != 0 || sum_sc_2 != 0) && $6!=0 && $7!=0 ) ) {print $0}}' JCEC.RNASeq.MXE.MATS.input.txt > MXE.JCEC.input.txt
/home/soft/rMATS.3.2.5/MATS/rMATS.sh -d MXE.JCEC.input.txt -o MXEout_JCEC/ -c 0.0001 -p 8 -t U
awk '{ split($2,ic_1,","); sum_ic_1=0; for (x in ic_1) sum_ic_1 += ic_1[x]; split($4,ic_2,","); sum_ic_2=0; for (x in ic_2) sum_ic_2 += ic_2[x]; split($3,sc_1,","); sum_sc_1=0; for (x in sc_1) sum_sc_1 += sc_1[x]; split($5,sc_2,","); sum_sc_2=0; for (x in sc_2) sum_sc_2 += sc_2[x]; if ( NR==1 || ( (sum_ic_2 + sum_sc_2 > 0) && (sum_ic_1 + sum_sc_1 > 0) && (sum_ic_1 != 0 || sum_ic_2 != 0) && (sum_sc_1 != 0 || sum_sc_2 != 0) && $6!=0 && $7!=0 ) ) {print $0}}' JCEC.RNASeq.A5SS.MATS.input.txt > A5SS.JCEC.input.txt
/home/soft/rMATS.3.2.5/MATS/rMATS.sh -d A5SS.JCEC.input.txt -o A5SSout_JCEC/ -c 0.0001 -p 8 -t U
awk '{ split($2,ic_1,","); sum_ic_1=0; for (x in ic_1) sum_ic_1 += ic_1[x]; split($4,ic_2,","); sum_ic_2=0; for (x in ic_2) sum_ic_2 += ic_2[x]; split($3,sc_1,","); sum_sc_1=0; for (x in sc_1) sum_sc_1 += sc_1[x]; split($5,sc_2,","); sum_sc_2=0; for (x in sc_2) sum_sc_2 += sc_2[x]; if ( NR==1 || ( (sum_ic_2 + sum_sc_2 > 0) && (sum_ic_1 + sum_sc_1 > 0) && (sum_ic_1 != 0 || sum_ic_2 != 0) && (sum_sc_1 != 0 || sum_sc_2 != 0) && $6!=0 && $7!=0 ) ) {print $0}}' JCEC.RNASeq.A3SS.MATS.input.txt > A3SS.JCEC.input.txt
/home/soft/rMATS.3.2.5/MATS/rMATS.sh -d A3SS.JCEC.input.txt -o A3SSout_JCEC/ -c 0.0001 -p 8 -t U
cd ..
rm -rf Root_D3piMinus_vs_D3piPlus
mkdir Root_D3piMinus_vs_D3piPlus
cd Root_D3piMinus_vs_D3piPlus
Rscript /home/pcd/merge_ASE_for_diff.r "D3piPlusR1_SRR1005302,D3piPlusR2_SRR1005303,D3piPlusR3_SRR1005304" "D3piMinusR1_SRR1005305,D3piMinusR2_SRR1005306,D3piMinusR3_SRR1005307"
awk '{ split($2,ic_1,","); sum_ic_1=0; for (x in ic_1) sum_ic_1 += ic_1[x]; split($4,ic_2,","); sum_ic_2=0; for (x in ic_2) sum_ic_2 += ic_2[x]; split($3,sc_1,","); sum_sc_1=0; for (x in sc_1) sum_sc_1 += sc_1[x]; split($5,sc_2,","); sum_sc_2=0; for (x in sc_2) sum_sc_2 += sc_2[x]; if ( NR==1 || ( (sum_ic_2 + sum_sc_2 > 0) && (sum_ic_1 + sum_sc_1 > 0) && (sum_ic_1 != 0 || sum_ic_2 != 0) && (sum_sc_1 != 0 || sum_sc_2 != 0) && $6!=0 && $7!=0 ) ) {print $0}}' JCEC.RNASeq.RI.MATS.input.txt > RI.JCEC.input.txt
/home/soft/rMATS.3.2.5/MATS/rMATS.sh -d RI.JCEC.input.txt -o RIout_JCEC/ -c 0.0001 -p 8 -t U
awk '{ split($2,ic_1,","); sum_ic_1=0; for (x in ic_1) sum_ic_1 += ic_1[x]; split($4,ic_2,","); sum_ic_2=0; for (x in ic_2) sum_ic_2 += ic_2[x]; split($3,sc_1,","); sum_sc_1=0; for (x in sc_1) sum_sc_1 += sc_1[x]; split($5,sc_2,","); sum_sc_2=0; for (x in sc_2) sum_sc_2 += sc_2[x]; if ( NR==1 || ( (sum_ic_2 + sum_sc_2 > 0) && (sum_ic_1 + sum_sc_1 > 0) && (sum_ic_1 != 0 || sum_ic_2 != 0) && (sum_sc_1 != 0 || sum_sc_2 != 0) && $6!=0 && $7!=0 ) ) {print $0}}' JCEC.RNASeq.SE.MATS.input.txt > SE.JCEC.input.txt
/home/soft/rMATS.3.2.5/MATS/rMATS.sh -d SE.JCEC.input.txt -o SEout_JCEC/ -c 0.0001 -p 8 -t U
awk '{ split($2,ic_1,","); sum_ic_1=0; for (x in ic_1) sum_ic_1 += ic_1[x]; split($4,ic_2,","); sum_ic_2=0; for (x in ic_2) sum_ic_2 += ic_2[x]; split($3,sc_1,","); sum_sc_1=0; for (x in sc_1) sum_sc_1 += sc_1[x]; split($5,sc_2,","); sum_sc_2=0; for (x in sc_2) sum_sc_2 += sc_2[x]; if ( NR==1 || ( (sum_ic_2 + sum_sc_2 > 0) && (sum_ic_1 + sum_sc_1 > 0) && (sum_ic_1 != 0 || sum_ic_2 != 0) && (sum_sc_1 != 0 || sum_sc_2 != 0) && $6!=0 && $7!=0 ) ) {print $0}}' JCEC.RNASeq.MXE.MATS.input.txt > MXE.JCEC.input.txt
/home/soft/rMATS.3.2.5/MATS/rMATS.sh -d MXE.JCEC.input.txt -o MXEout_JCEC/ -c 0.0001 -p 8 -t U
awk '{ split($2,ic_1,","); sum_ic_1=0; for (x in ic_1) sum_ic_1 += ic_1[x]; split($4,ic_2,","); sum_ic_2=0; for (x in ic_2) sum_ic_2 += ic_2[x]; split($3,sc_1,","); sum_sc_1=0; for (x in sc_1) sum_sc_1 += sc_1[x]; split($5,sc_2,","); sum_sc_2=0; for (x in sc_2) sum_sc_2 += sc_2[x]; if ( NR==1 || ( (sum_ic_2 + sum_sc_2 > 0) && (sum_ic_1 + sum_sc_1 > 0) && (sum_ic_1 != 0 || sum_ic_2 != 0) && (sum_sc_1 != 0 || sum_sc_2 != 0) && $6!=0 && $7!=0 ) ) {print $0}}' JCEC.RNASeq.A5SS.MATS.input.txt > A5SS.JCEC.input.txt
/home/soft/rMATS.3.2.5/MATS/rMATS.sh -d A5SS.JCEC.input.txt -o A5SSout_JCEC/ -c 0.0001 -p 8 -t U
awk '{ split($2,ic_1,","); sum_ic_1=0; for (x in ic_1) sum_ic_1 += ic_1[x]; split($4,ic_2,","); sum_ic_2=0; for (x in ic_2) sum_ic_2 += ic_2[x]; split($3,sc_1,","); sum_sc_1=0; for (x in sc_1) sum_sc_1 += sc_1[x]; split($5,sc_2,","); sum_sc_2=0; for (x in sc_2) sum_sc_2 += sc_2[x]; if ( NR==1 || ( (sum_ic_2 + sum_sc_2 > 0) && (sum_ic_1 + sum_sc_1 > 0) && (sum_ic_1 != 0 || sum_ic_2 != 0) && (sum_sc_1 != 0 || sum_sc_2 != 0) && $6!=0 && $7!=0 ) ) {print $0}}' JCEC.RNASeq.A3SS.MATS.input.txt > A3SS.JCEC.input.txt
/home/soft/rMATS.3.2.5/MATS/rMATS.sh -d A3SS.JCEC.input.txt -o A3SSout_JCEC/ -c 0.0001 -p 8 -t U
cd ..
rm -rf Root_D7piMinus_vs_D7piPlus
mkdir Root_D7piMinus_vs_D7piPlus
cd Root_D7piMinus_vs_D7piPlus
Rscript /home/pcd/merge_ASE_for_diff.r "D7piPlusR1_SRR1005314,D7piPlusR2_SRR1005315,D7piPlusR3_SRR1005316" "D7piMinusR1_SRR1005317,D7piMinusR2_SRR1005318,D7piMinusR3_SRR1005319"
awk '{ split($2,ic_1,","); sum_ic_1=0; for (x in ic_1) sum_ic_1 += ic_1[x]; split($4,ic_2,","); sum_ic_2=0; for (x in ic_2) sum_ic_2 += ic_2[x]; split($3,sc_1,","); sum_sc_1=0; for (x in sc_1) sum_sc_1 += sc_1[x]; split($5,sc_2,","); sum_sc_2=0; for (x in sc_2) sum_sc_2 += sc_2[x]; if ( NR==1 || ( (sum_ic_2 + sum_sc_2 > 0) && (sum_ic_1 + sum_sc_1 > 0) && (sum_ic_1 != 0 || sum_ic_2 != 0) && (sum_sc_1 != 0 || sum_sc_2 != 0) && $6!=0 && $7!=0 ) ) {print $0}}' JCEC.RNASeq.RI.MATS.input.txt > RI.JCEC.input.txt
/home/soft/rMATS.3.2.5/MATS/rMATS.sh -d RI.JCEC.input.txt -o RIout_JCEC/ -c 0.0001 -p 8 -t U
awk '{ split($2,ic_1,","); sum_ic_1=0; for (x in ic_1) sum_ic_1 += ic_1[x]; split($4,ic_2,","); sum_ic_2=0; for (x in ic_2) sum_ic_2 += ic_2[x]; split($3,sc_1,","); sum_sc_1=0; for (x in sc_1) sum_sc_1 += sc_1[x]; split($5,sc_2,","); sum_sc_2=0; for (x in sc_2) sum_sc_2 += sc_2[x]; if ( NR==1 || ( (sum_ic_2 + sum_sc_2 > 0) && (sum_ic_1 + sum_sc_1 > 0) && (sum_ic_1 != 0 || sum_ic_2 != 0) && (sum_sc_1 != 0 || sum_sc_2 != 0) && $6!=0 && $7!=0 ) ) {print $0}}' JCEC.RNASeq.SE.MATS.input.txt > SE.JCEC.input.txt
/home/soft/rMATS.3.2.5/MATS/rMATS.sh -d SE.JCEC.input.txt -o SEout_JCEC/ -c 0.0001 -p 8 -t U
awk '{ split($2,ic_1,","); sum_ic_1=0; for (x in ic_1) sum_ic_1 += ic_1[x]; split($4,ic_2,","); sum_ic_2=0; for (x in ic_2) sum_ic_2 += ic_2[x]; split($3,sc_1,","); sum_sc_1=0; for (x in sc_1) sum_sc_1 += sc_1[x]; split($5,sc_2,","); sum_sc_2=0; for (x in sc_2) sum_sc_2 += sc_2[x]; if ( NR==1 || ( (sum_ic_2 + sum_sc_2 > 0) && (sum_ic_1 + sum_sc_1 > 0) && (sum_ic_1 != 0 || sum_ic_2 != 0) && (sum_sc_1 != 0 || sum_sc_2 != 0) && $6!=0 && $7!=0 ) ) {print $0}}' JCEC.RNASeq.MXE.MATS.input.txt > MXE.JCEC.input.txt
/home/soft/rMATS.3.2.5/MATS/rMATS.sh -d MXE.JCEC.input.txt -o MXEout_JCEC/ -c 0.0001 -p 8 -t U
awk '{ split($2,ic_1,","); sum_ic_1=0; for (x in ic_1) sum_ic_1 += ic_1[x]; split($4,ic_2,","); sum_ic_2=0; for (x in ic_2) sum_ic_2 += ic_2[x]; split($3,sc_1,","); sum_sc_1=0; for (x in sc_1) sum_sc_1 += sc_1[x]; split($5,sc_2,","); sum_sc_2=0; for (x in sc_2) sum_sc_2 += sc_2[x]; if ( NR==1 || ( (sum_ic_2 + sum_sc_2 > 0) && (sum_ic_1 + sum_sc_1 > 0) && (sum_ic_1 != 0 || sum_ic_2 != 0) && (sum_sc_1 != 0 || sum_sc_2 != 0) && $6!=0 && $7!=0 ) ) {print $0}}' JCEC.RNASeq.A5SS.MATS.input.txt > A5SS.JCEC.input.txt
/home/soft/rMATS.3.2.5/MATS/rMATS.sh -d A5SS.JCEC.input.txt -o A5SSout_JCEC/ -c 0.0001 -p 8 -t U
awk '{ split($2,ic_1,","); sum_ic_1=0; for (x in ic_1) sum_ic_1 += ic_1[x]; split($4,ic_2,","); sum_ic_2=0; for (x in ic_2) sum_ic_2 += ic_2[x]; split($3,sc_1,","); sum_sc_1=0; for (x in sc_1) sum_sc_1 += sc_1[x]; split($5,sc_2,","); sum_sc_2=0; for (x in sc_2) sum_sc_2 += sc_2[x]; if ( NR==1 || ( (sum_ic_2 + sum_sc_2 > 0) && (sum_ic_1 + sum_sc_1 > 0) && (sum_ic_1 != 0 || sum_ic_2 != 0) && (sum_sc_1 != 0 || sum_sc_2 != 0) && $6!=0 && $7!=0 ) ) {print $0}}' JCEC.RNASeq.A3SS.MATS.input.txt > A3SS.JCEC.input.txt
/home/soft/rMATS.3.2.5/MATS/rMATS.sh -d A3SS.JCEC.input.txt -o A3SSout_JCEC/ -c 0.0001 -p 8 -t U
cd ..
rm -rf Root_D21piMinus_vs_D21piPlus
mkdir Root_D21piMinus_vs_D21piPlus
cd Root_D21piMinus_vs_D21piPlus
Rscript /home/pcd/merge_ASE_for_diff.r "D21piPlusR1_SRR1005290,D21piPlusR2_SRR1005291,D21piPlusR3_SRR1005292" "D21piMinusR1_SRR1005293,D21piMinusR2_SRR1005294,D21piMinusR3_SRR1005295"
awk '{ split($2,ic_1,","); sum_ic_1=0; for (x in ic_1) sum_ic_1 += ic_1[x]; split($4,ic_2,","); sum_ic_2=0; for (x in ic_2) sum_ic_2 += ic_2[x]; split($3,sc_1,","); sum_sc_1=0; for (x in sc_1) sum_sc_1 += sc_1[x]; split($5,sc_2,","); sum_sc_2=0; for (x in sc_2) sum_sc_2 += sc_2[x]; if ( NR==1 || ( (sum_ic_2 + sum_sc_2 > 0) && (sum_ic_1 + sum_sc_1 > 0) && (sum_ic_1 != 0 || sum_ic_2 != 0) && (sum_sc_1 != 0 || sum_sc_2 != 0) && $6!=0 && $7!=0 ) ) {print $0}}' JCEC.RNASeq.RI.MATS.input.txt > RI.JCEC.input.txt
/home/soft/rMATS.3.2.5/MATS/rMATS.sh -d RI.JCEC.input.txt -o RIout_JCEC/ -c 0.0001 -p 8 -t U
awk '{ split($2,ic_1,","); sum_ic_1=0; for (x in ic_1) sum_ic_1 += ic_1[x]; split($4,ic_2,","); sum_ic_2=0; for (x in ic_2) sum_ic_2 += ic_2[x]; split($3,sc_1,","); sum_sc_1=0; for (x in sc_1) sum_sc_1 += sc_1[x]; split($5,sc_2,","); sum_sc_2=0; for (x in sc_2) sum_sc_2 += sc_2[x]; if ( NR==1 || ( (sum_ic_2 + sum_sc_2 > 0) && (sum_ic_1 + sum_sc_1 > 0) && (sum_ic_1 != 0 || sum_ic_2 != 0) && (sum_sc_1 != 0 || sum_sc_2 != 0) && $6!=0 && $7!=0 ) ) {print $0}}' JCEC.RNASeq.SE.MATS.input.txt > SE.JCEC.input.txt
/home/soft/rMATS.3.2.5/MATS/rMATS.sh -d SE.JCEC.input.txt -o SEout_JCEC/ -c 0.0001 -p 8 -t U
awk '{ split($2,ic_1,","); sum_ic_1=0; for (x in ic_1) sum_ic_1 += ic_1[x]; split($4,ic_2,","); sum_ic_2=0; for (x in ic_2) sum_ic_2 += ic_2[x]; split($3,sc_1,","); sum_sc_1=0; for (x in sc_1) sum_sc_1 += sc_1[x]; split($5,sc_2,","); sum_sc_2=0; for (x in sc_2) sum_sc_2 += sc_2[x]; if ( NR==1 || ( (sum_ic_2 + sum_sc_2 > 0) && (sum_ic_1 + sum_sc_1 > 0) && (sum_ic_1 != 0 || sum_ic_2 != 0) && (sum_sc_1 != 0 || sum_sc_2 != 0) && $6!=0 && $7!=0 ) ) {print $0}}' JCEC.RNASeq.MXE.MATS.input.txt > MXE.JCEC.input.txt
/home/soft/rMATS.3.2.5/MATS/rMATS.sh -d MXE.JCEC.input.txt -o MXEout_JCEC/ -c 0.0001 -p 8 -t U
awk '{ split($2,ic_1,","); sum_ic_1=0; for (x in ic_1) sum_ic_1 += ic_1[x]; split($4,ic_2,","); sum_ic_2=0; for (x in ic_2) sum_ic_2 += ic_2[x]; split($3,sc_1,","); sum_sc_1=0; for (x in sc_1) sum_sc_1 += sc_1[x]; split($5,sc_2,","); sum_sc_2=0; for (x in sc_2) sum_sc_2 += sc_2[x]; if ( NR==1 || ( (sum_ic_2 + sum_sc_2 > 0) && (sum_ic_1 + sum_sc_1 > 0) && (sum_ic_1 != 0 || sum_ic_2 != 0) && (sum_sc_1 != 0 || sum_sc_2 != 0) && $6!=0 && $7!=0 ) ) {print $0}}' JCEC.RNASeq.A5SS.MATS.input.txt > A5SS.JCEC.input.txt
/home/soft/rMATS.3.2.5/MATS/rMATS.sh -d A5SS.JCEC.input.txt -o A5SSout_JCEC/ -c 0.0001 -p 8 -t U
awk '{ split($2,ic_1,","); sum_ic_1=0; for (x in ic_1) sum_ic_1 += ic_1[x]; split($4,ic_2,","); sum_ic_2=0; for (x in ic_2) sum_ic_2 += ic_2[x]; split($3,sc_1,","); sum_sc_1=0; for (x in sc_1) sum_sc_1 += sc_1[x]; split($5,sc_2,","); sum_sc_2=0; for (x in sc_2) sum_sc_2 += sc_2[x]; if ( NR==1 || ( (sum_ic_2 + sum_sc_2 > 0) && (sum_ic_1 + sum_sc_1 > 0) && (sum_ic_1 != 0 || sum_ic_2 != 0) && (sum_sc_1 != 0 || sum_sc_2 != 0) && $6!=0 && $7!=0 ) ) {print $0}}' JCEC.RNASeq.A3SS.MATS.input.txt > A3SS.JCEC.input.txt
/home/soft/rMATS.3.2.5/MATS/rMATS.sh -d A3SS.JCEC.input.txt -o A3SSout_JCEC/ -c 0.0001 -p 8 -t U
cd ..
rm -rf Root_D21H1piMinus_vs_D21H1piPlus
mkdir Root_D21H1piMinus_vs_D21H1piPlus
cd Root_D21H1piMinus_vs_D21H1piPlus
Rscript /home/pcd/merge_ASE_for_diff.r "D21H1piPlusR1_SRR1005263,D21H1piPlusR2_SRR1005264,D21H1piPlusR3_SRR1005265" "D21H1piMinusR1_SRR1005266,D21H1piMinusR2_SRR1005267,D21H1piMinusR3_SRR1005268"
awk '{ split($2,ic_1,","); sum_ic_1=0; for (x in ic_1) sum_ic_1 += ic_1[x]; split($4,ic_2,","); sum_ic_2=0; for (x in ic_2) sum_ic_2 += ic_2[x]; split($3,sc_1,","); sum_sc_1=0; for (x in sc_1) sum_sc_1 += sc_1[x]; split($5,sc_2,","); sum_sc_2=0; for (x in sc_2) sum_sc_2 += sc_2[x]; if ( NR==1 || ( (sum_ic_2 + sum_sc_2 > 0) && (sum_ic_1 + sum_sc_1 > 0) && (sum_ic_1 != 0 || sum_ic_2 != 0) && (sum_sc_1 != 0 || sum_sc_2 != 0) && $6!=0 && $7!=0 ) ) {print $0}}' JCEC.RNASeq.RI.MATS.input.txt > RI.JCEC.input.txt
/home/soft/rMATS.3.2.5/MATS/rMATS.sh -d RI.JCEC.input.txt -o RIout_JCEC/ -c 0.0001 -p 8 -t U
awk '{ split($2,ic_1,","); sum_ic_1=0; for (x in ic_1) sum_ic_1 += ic_1[x]; split($4,ic_2,","); sum_ic_2=0; for (x in ic_2) sum_ic_2 += ic_2[x]; split($3,sc_1,","); sum_sc_1=0; for (x in sc_1) sum_sc_1 += sc_1[x]; split($5,sc_2,","); sum_sc_2=0; for (x in sc_2) sum_sc_2 += sc_2[x]; if ( NR==1 || ( (sum_ic_2 + sum_sc_2 > 0) && (sum_ic_1 + sum_sc_1 > 0) && (sum_ic_1 != 0 || sum_ic_2 != 0) && (sum_sc_1 != 0 || sum_sc_2 != 0) && $6!=0 && $7!=0 ) ) {print $0}}' JCEC.RNASeq.SE.MATS.input.txt > SE.JCEC.input.txt
/home/soft/rMATS.3.2.5/MATS/rMATS.sh -d SE.JCEC.input.txt -o SEout_JCEC/ -c 0.0001 -p 8 -t U
awk '{ split($2,ic_1,","); sum_ic_1=0; for (x in ic_1) sum_ic_1 += ic_1[x]; split($4,ic_2,","); sum_ic_2=0; for (x in ic_2) sum_ic_2 += ic_2[x]; split($3,sc_1,","); sum_sc_1=0; for (x in sc_1) sum_sc_1 += sc_1[x]; split($5,sc_2,","); sum_sc_2=0; for (x in sc_2) sum_sc_2 += sc_2[x]; if ( NR==1 || ( (sum_ic_2 + sum_sc_2 > 0) && (sum_ic_1 + sum_sc_1 > 0) && (sum_ic_1 != 0 || sum_ic_2 != 0) && (sum_sc_1 != 0 || sum_sc_2 != 0) && $6!=0 && $7!=0 ) ) {print $0}}' JCEC.RNASeq.MXE.MATS.input.txt > MXE.JCEC.input.txt
/home/soft/rMATS.3.2.5/MATS/rMATS.sh -d MXE.JCEC.input.txt -o MXEout_JCEC/ -c 0.0001 -p 8 -t U
awk '{ split($2,ic_1,","); sum_ic_1=0; for (x in ic_1) sum_ic_1 += ic_1[x]; split($4,ic_2,","); sum_ic_2=0; for (x in ic_2) sum_ic_2 += ic_2[x]; split($3,sc_1,","); sum_sc_1=0; for (x in sc_1) sum_sc_1 += sc_1[x]; split($5,sc_2,","); sum_sc_2=0; for (x in sc_2) sum_sc_2 += sc_2[x]; if ( NR==1 || ( (sum_ic_2 + sum_sc_2 > 0) && (sum_ic_1 + sum_sc_1 > 0) && (sum_ic_1 != 0 || sum_ic_2 != 0) && (sum_sc_1 != 0 || sum_sc_2 != 0) && $6!=0 && $7!=0 ) ) {print $0}}' JCEC.RNASeq.A5SS.MATS.input.txt > A5SS.JCEC.input.txt
/home/soft/rMATS.3.2.5/MATS/rMATS.sh -d A5SS.JCEC.input.txt -o A5SSout_JCEC/ -c 0.0001 -p 8 -t U
awk '{ split($2,ic_1,","); sum_ic_1=0; for (x in ic_1) sum_ic_1 += ic_1[x]; split($4,ic_2,","); sum_ic_2=0; for (x in ic_2) sum_ic_2 += ic_2[x]; split($3,sc_1,","); sum_sc_1=0; for (x in sc_1) sum_sc_1 += sc_1[x]; split($5,sc_2,","); sum_sc_2=0; for (x in sc_2) sum_sc_2 += sc_2[x]; if ( NR==1 || ( (sum_ic_2 + sum_sc_2 > 0) && (sum_ic_1 + sum_sc_1 > 0) && (sum_ic_1 != 0 || sum_ic_2 != 0) && (sum_sc_1 != 0 || sum_sc_2 != 0) && $6!=0 && $7!=0 ) ) {print $0}}' JCEC.RNASeq.A3SS.MATS.input.txt > A3SS.JCEC.input.txt
/home/soft/rMATS.3.2.5/MATS/rMATS.sh -d A3SS.JCEC.input.txt -o A3SSout_JCEC/ -c 0.0001 -p 8 -t U
cd ..
rm -rf Root_D21H6piMinus_vs_D21H6piPlus
mkdir Root_D21H6piMinus_vs_D21H6piPlus
cd Root_D21H6piMinus_vs_D21H6piPlus
Rscript /home/pcd/merge_ASE_for_diff.r "D21H6piPlusR1_SRR1005281,D21H6piPlusR2_SRR1005282,D21H6piPlusR3_SRR1005283" "D21H6piMinusR1_SRR1005284,D21H6piMinusR2_SRR1005285,D21H6piMinusR3_SRR1005286"
awk '{ split($2,ic_1,","); sum_ic_1=0; for (x in ic_1) sum_ic_1 += ic_1[x]; split($4,ic_2,","); sum_ic_2=0; for (x in ic_2) sum_ic_2 += ic_2[x]; split($3,sc_1,","); sum_sc_1=0; for (x in sc_1) sum_sc_1 += sc_1[x]; split($5,sc_2,","); sum_sc_2=0; for (x in sc_2) sum_sc_2 += sc_2[x]; if ( NR==1 || ( (sum_ic_2 + sum_sc_2 > 0) && (sum_ic_1 + sum_sc_1 > 0) && (sum_ic_1 != 0 || sum_ic_2 != 0) && (sum_sc_1 != 0 || sum_sc_2 != 0) && $6!=0 && $7!=0 ) ) {print $0}}' JCEC.RNASeq.RI.MATS.input.txt > RI.JCEC.input.txt
/home/soft/rMATS.3.2.5/MATS/rMATS.sh -d RI.JCEC.input.txt -o RIout_JCEC/ -c 0.0001 -p 8 -t U
awk '{ split($2,ic_1,","); sum_ic_1=0; for (x in ic_1) sum_ic_1 += ic_1[x]; split($4,ic_2,","); sum_ic_2=0; for (x in ic_2) sum_ic_2 += ic_2[x]; split($3,sc_1,","); sum_sc_1=0; for (x in sc_1) sum_sc_1 += sc_1[x]; split($5,sc_2,","); sum_sc_2=0; for (x in sc_2) sum_sc_2 += sc_2[x]; if ( NR==1 || ( (sum_ic_2 + sum_sc_2 > 0) && (sum_ic_1 + sum_sc_1 > 0) && (sum_ic_1 != 0 || sum_ic_2 != 0) && (sum_sc_1 != 0 || sum_sc_2 != 0) && $6!=0 && $7!=0 ) ) {print $0}}' JCEC.RNASeq.SE.MATS.input.txt > SE.JCEC.input.txt
/home/soft/rMATS.3.2.5/MATS/rMATS.sh -d SE.JCEC.input.txt -o SEout_JCEC/ -c 0.0001 -p 8 -t U
awk '{ split($2,ic_1,","); sum_ic_1=0; for (x in ic_1) sum_ic_1 += ic_1[x]; split($4,ic_2,","); sum_ic_2=0; for (x in ic_2) sum_ic_2 += ic_2[x]; split($3,sc_1,","); sum_sc_1=0; for (x in sc_1) sum_sc_1 += sc_1[x]; split($5,sc_2,","); sum_sc_2=0; for (x in sc_2) sum_sc_2 += sc_2[x]; if ( NR==1 || ( (sum_ic_2 + sum_sc_2 > 0) && (sum_ic_1 + sum_sc_1 > 0) && (sum_ic_1 != 0 || sum_ic_2 != 0) && (sum_sc_1 != 0 || sum_sc_2 != 0) && $6!=0 && $7!=0 ) ) {print $0}}' JCEC.RNASeq.MXE.MATS.input.txt > MXE.JCEC.input.txt
/home/soft/rMATS.3.2.5/MATS/rMATS.sh -d MXE.JCEC.input.txt -o MXEout_JCEC/ -c 0.0001 -p 8 -t U
awk '{ split($2,ic_1,","); sum_ic_1=0; for (x in ic_1) sum_ic_1 += ic_1[x]; split($4,ic_2,","); sum_ic_2=0; for (x in ic_2) sum_ic_2 += ic_2[x]; split($3,sc_1,","); sum_sc_1=0; for (x in sc_1) sum_sc_1 += sc_1[x]; split($5,sc_2,","); sum_sc_2=0; for (x in sc_2) sum_sc_2 += sc_2[x]; if ( NR==1 || ( (sum_ic_2 + sum_sc_2 > 0) && (sum_ic_1 + sum_sc_1 > 0) && (sum_ic_1 != 0 || sum_ic_2 != 0) && (sum_sc_1 != 0 || sum_sc_2 != 0) && $6!=0 && $7!=0 ) ) {print $0}}' JCEC.RNASeq.A5SS.MATS.input.txt > A5SS.JCEC.input.txt
/home/soft/rMATS.3.2.5/MATS/rMATS.sh -d A5SS.JCEC.input.txt -o A5SSout_JCEC/ -c 0.0001 -p 8 -t U
awk '{ split($2,ic_1,","); sum_ic_1=0; for (x in ic_1) sum_ic_1 += ic_1[x]; split($4,ic_2,","); sum_ic_2=0; for (x in ic_2) sum_ic_2 += ic_2[x]; split($3,sc_1,","); sum_sc_1=0; for (x in sc_1) sum_sc_1 += sc_1[x]; split($5,sc_2,","); sum_sc_2=0; for (x in sc_2) sum_sc_2 += sc_2[x]; if ( NR==1 || ( (sum_ic_2 + sum_sc_2 > 0) && (sum_ic_1 + sum_sc_1 > 0) && (sum_ic_1 != 0 || sum_ic_2 != 0) && (sum_sc_1 != 0 || sum_sc_2 != 0) && $6!=0 && $7!=0 ) ) {print $0}}' JCEC.RNASeq.A3SS.MATS.input.txt > A3SS.JCEC.input.txt
/home/soft/rMATS.3.2.5/MATS/rMATS.sh -d A3SS.JCEC.input.txt -o A3SSout_JCEC/ -c 0.0001 -p 8 -t U
cd ..
rm -rf Root_D21H24piMinus_vs_D21H24piPlus
mkdir Root_D21H24piMinus_vs_D21H24piPlus
cd Root_D21H24piMinus_vs_D21H24piPlus
Rscript /home/pcd/merge_ASE_for_diff.r "D21H24piPlusR1_SRR1005272,D21H24piPlusR2_SRR1005273,D21H24piPlusR3_SRR1005274" "D21H24piMinusR1_SRR1005275,D21H24piMinusR2_SRR1005276,D21H24piMinusR3_SRR1005277"
awk '{ split($2,ic_1,","); sum_ic_1=0; for (x in ic_1) sum_ic_1 += ic_1[x]; split($4,ic_2,","); sum_ic_2=0; for (x in ic_2) sum_ic_2 += ic_2[x]; split($3,sc_1,","); sum_sc_1=0; for (x in sc_1) sum_sc_1 += sc_1[x]; split($5,sc_2,","); sum_sc_2=0; for (x in sc_2) sum_sc_2 += sc_2[x]; if ( NR==1 || ( (sum_ic_2 + sum_sc_2 > 0) && (sum_ic_1 + sum_sc_1 > 0) && (sum_ic_1 != 0 || sum_ic_2 != 0) && (sum_sc_1 != 0 || sum_sc_2 != 0) && $6!=0 && $7!=0 ) ) {print $0}}' JCEC.RNASeq.RI.MATS.input.txt > RI.JCEC.input.txt
/home/soft/rMATS.3.2.5/MATS/rMATS.sh -d RI.JCEC.input.txt -o RIout_JCEC/ -c 0.0001 -p 8 -t U
awk '{ split($2,ic_1,","); sum_ic_1=0; for (x in ic_1) sum_ic_1 += ic_1[x]; split($4,ic_2,","); sum_ic_2=0; for (x in ic_2) sum_ic_2 += ic_2[x]; split($3,sc_1,","); sum_sc_1=0; for (x in sc_1) sum_sc_1 += sc_1[x]; split($5,sc_2,","); sum_sc_2=0; for (x in sc_2) sum_sc_2 += sc_2[x]; if ( NR==1 || ( (sum_ic_2 + sum_sc_2 > 0) && (sum_ic_1 + sum_sc_1 > 0) && (sum_ic_1 != 0 || sum_ic_2 != 0) && (sum_sc_1 != 0 || sum_sc_2 != 0) && $6!=0 && $7!=0 ) ) {print $0}}' JCEC.RNASeq.SE.MATS.input.txt > SE.JCEC.input.txt
/home/soft/rMATS.3.2.5/MATS/rMATS.sh -d SE.JCEC.input.txt -o SEout_JCEC/ -c 0.0001 -p 8 -t U
awk '{ split($2,ic_1,","); sum_ic_1=0; for (x in ic_1) sum_ic_1 += ic_1[x]; split($4,ic_2,","); sum_ic_2=0; for (x in ic_2) sum_ic_2 += ic_2[x]; split($3,sc_1,","); sum_sc_1=0; for (x in sc_1) sum_sc_1 += sc_1[x]; split($5,sc_2,","); sum_sc_2=0; for (x in sc_2) sum_sc_2 += sc_2[x]; if ( NR==1 || ( (sum_ic_2 + sum_sc_2 > 0) && (sum_ic_1 + sum_sc_1 > 0) && (sum_ic_1 != 0 || sum_ic_2 != 0) && (sum_sc_1 != 0 || sum_sc_2 != 0) && $6!=0 && $7!=0 ) ) {print $0}}' JCEC.RNASeq.MXE.MATS.input.txt > MXE.JCEC.input.txt
/home/soft/rMATS.3.2.5/MATS/rMATS.sh -d MXE.JCEC.input.txt -o MXEout_JCEC/ -c 0.0001 -p 8 -t U
awk '{ split($2,ic_1,","); sum_ic_1=0; for (x in ic_1) sum_ic_1 += ic_1[x]; split($4,ic_2,","); sum_ic_2=0; for (x in ic_2) sum_ic_2 += ic_2[x]; split($3,sc_1,","); sum_sc_1=0; for (x in sc_1) sum_sc_1 += sc_1[x]; split($5,sc_2,","); sum_sc_2=0; for (x in sc_2) sum_sc_2 += sc_2[x]; if ( NR==1 || ( (sum_ic_2 + sum_sc_2 > 0) && (sum_ic_1 + sum_sc_1 > 0) && (sum_ic_1 != 0 || sum_ic_2 != 0) && (sum_sc_1 != 0 || sum_sc_2 != 0) && $6!=0 && $7!=0 ) ) {print $0}}' JCEC.RNASeq.A5SS.MATS.input.txt > A5SS.JCEC.input.txt
/home/soft/rMATS.3.2.5/MATS/rMATS.sh -d A5SS.JCEC.input.txt -o A5SSout_JCEC/ -c 0.0001 -p 8 -t U
awk '{ split($2,ic_1,","); sum_ic_1=0; for (x in ic_1) sum_ic_1 += ic_1[x]; split($4,ic_2,","); sum_ic_2=0; for (x in ic_2) sum_ic_2 += ic_2[x]; split($3,sc_1,","); sum_sc_1=0; for (x in sc_1) sum_sc_1 += sc_1[x]; split($5,sc_2,","); sum_sc_2=0; for (x in sc_2) sum_sc_2 += sc_2[x]; if ( NR==1 || ( (sum_ic_2 + sum_sc_2 > 0) && (sum_ic_1 + sum_sc_1 > 0) && (sum_ic_1 != 0 || sum_ic_2 != 0) && (sum_sc_1 != 0 || sum_sc_2 != 0) && $6!=0 && $7!=0 ) ) {print $0}}' JCEC.RNASeq.A3SS.MATS.input.txt > A3SS.JCEC.input.txt
/home/soft/rMATS.3.2.5/MATS/rMATS.sh -d A3SS.JCEC.input.txt -o A3SSout_JCEC/ -c 0.0001 -p 8 -t U
cd ..



# Identify differential alternative splicing events under Pi deficiency in shoots

rm -rf Shoot_H1piMinus_vs_H1piPlus
mkdir Shoot_H1piMinus_vs_H1piPlus
cd Shoot_H1piMinus_vs_H1piPlus
Rscript /home/pcd/merge_ASE_for_diff.r "H1piPlusR1_SRR1005320,H1piPlusR2_SRR1005321,H1piPlusR3_SRR1005322" "H1piMinusR1_SRR1005323,H1piMinusR2_SRR1005324,H1piMinusR3_SRR1005325"
awk '{ split($2,ic_1,","); sum_ic_1=0; for (x in ic_1) sum_ic_1 += ic_1[x]; split($4,ic_2,","); sum_ic_2=0; for (x in ic_2) sum_ic_2 += ic_2[x]; split($3,sc_1,","); sum_sc_1=0; for (x in sc_1) sum_sc_1 += sc_1[x]; split($5,sc_2,","); sum_sc_2=0; for (x in sc_2) sum_sc_2 += sc_2[x]; if ( NR==1 || ( (sum_ic_2 + sum_sc_2 > 0) && (sum_ic_1 + sum_sc_1 > 0) && (sum_ic_1 != 0 || sum_ic_2 != 0) && (sum_sc_1 != 0 || sum_sc_2 != 0) && $6!=0 && $7!=0 ) ) {print $0}}' JCEC.RNASeq.RI.MATS.input.txt > RI.JCEC.input.txt
/home/soft/rMATS.3.2.5/MATS/rMATS.sh -d RI.JCEC.input.txt -o RIout_JCEC/ -c 0.0001 -p 8 -t U
awk '{ split($2,ic_1,","); sum_ic_1=0; for (x in ic_1) sum_ic_1 += ic_1[x]; split($4,ic_2,","); sum_ic_2=0; for (x in ic_2) sum_ic_2 += ic_2[x]; split($3,sc_1,","); sum_sc_1=0; for (x in sc_1) sum_sc_1 += sc_1[x]; split($5,sc_2,","); sum_sc_2=0; for (x in sc_2) sum_sc_2 += sc_2[x]; if ( NR==1 || ( (sum_ic_2 + sum_sc_2 > 0) && (sum_ic_1 + sum_sc_1 > 0) && (sum_ic_1 != 0 || sum_ic_2 != 0) && (sum_sc_1 != 0 || sum_sc_2 != 0) && $6!=0 && $7!=0 ) ) {print $0}}' JCEC.RNASeq.SE.MATS.input.txt > SE.JCEC.input.txt
/home/soft/rMATS.3.2.5/MATS/rMATS.sh -d SE.JCEC.input.txt -o SEout_JCEC/ -c 0.0001 -p 8 -t U
awk '{ split($2,ic_1,","); sum_ic_1=0; for (x in ic_1) sum_ic_1 += ic_1[x]; split($4,ic_2,","); sum_ic_2=0; for (x in ic_2) sum_ic_2 += ic_2[x]; split($3,sc_1,","); sum_sc_1=0; for (x in sc_1) sum_sc_1 += sc_1[x]; split($5,sc_2,","); sum_sc_2=0; for (x in sc_2) sum_sc_2 += sc_2[x]; if ( NR==1 || ( (sum_ic_2 + sum_sc_2 > 0) && (sum_ic_1 + sum_sc_1 > 0) && (sum_ic_1 != 0 || sum_ic_2 != 0) && (sum_sc_1 != 0 || sum_sc_2 != 0) && $6!=0 && $7!=0 ) ) {print $0}}' JCEC.RNASeq.MXE.MATS.input.txt > MXE.JCEC.input.txt
/home/soft/rMATS.3.2.5/MATS/rMATS.sh -d MXE.JCEC.input.txt -o MXEout_JCEC/ -c 0.0001 -p 8 -t U
awk '{ split($2,ic_1,","); sum_ic_1=0; for (x in ic_1) sum_ic_1 += ic_1[x]; split($4,ic_2,","); sum_ic_2=0; for (x in ic_2) sum_ic_2 += ic_2[x]; split($3,sc_1,","); sum_sc_1=0; for (x in sc_1) sum_sc_1 += sc_1[x]; split($5,sc_2,","); sum_sc_2=0; for (x in sc_2) sum_sc_2 += sc_2[x]; if ( NR==1 || ( (sum_ic_2 + sum_sc_2 > 0) && (sum_ic_1 + sum_sc_1 > 0) && (sum_ic_1 != 0 || sum_ic_2 != 0) && (sum_sc_1 != 0 || sum_sc_2 != 0) && $6!=0 && $7!=0 ) ) {print $0}}' JCEC.RNASeq.A5SS.MATS.input.txt > A5SS.JCEC.input.txt
/home/soft/rMATS.3.2.5/MATS/rMATS.sh -d A5SS.JCEC.input.txt -o A5SSout_JCEC/ -c 0.0001 -p 8 -t U
awk '{ split($2,ic_1,","); sum_ic_1=0; for (x in ic_1) sum_ic_1 += ic_1[x]; split($4,ic_2,","); sum_ic_2=0; for (x in ic_2) sum_ic_2 += ic_2[x]; split($3,sc_1,","); sum_sc_1=0; for (x in sc_1) sum_sc_1 += sc_1[x]; split($5,sc_2,","); sum_sc_2=0; for (x in sc_2) sum_sc_2 += sc_2[x]; if ( NR==1 || ( (sum_ic_2 + sum_sc_2 > 0) && (sum_ic_1 + sum_sc_1 > 0) && (sum_ic_1 != 0 || sum_ic_2 != 0) && (sum_sc_1 != 0 || sum_sc_2 != 0) && $6!=0 && $7!=0 ) ) {print $0}}' JCEC.RNASeq.A3SS.MATS.input.txt > A3SS.JCEC.input.txt
/home/soft/rMATS.3.2.5/MATS/rMATS.sh -d A3SS.JCEC.input.txt -o A3SSout_JCEC/ -c 0.0001 -p 8 -t U
cd ..
rm -rf Shoot_H6piMinus_vs_H6piPlus
mkdir Shoot_H6piMinus_vs_H6piPlus
cd Shoot_H6piMinus_vs_H6piPlus
Rscript /home/pcd/merge_ASE_for_diff.r "H6piPlusR1_SRR1005371,H6piPlusR2_SRR1005372,H6piPlusR3_SRR1005373" "H6piMinusR1_SRR1005374,H6piMinusR2_SRR1005375,H6piMinusR3_SRR1005376"
awk '{ split($2,ic_1,","); sum_ic_1=0; for (x in ic_1) sum_ic_1 += ic_1[x]; split($4,ic_2,","); sum_ic_2=0; for (x in ic_2) sum_ic_2 += ic_2[x]; split($3,sc_1,","); sum_sc_1=0; for (x in sc_1) sum_sc_1 += sc_1[x]; split($5,sc_2,","); sum_sc_2=0; for (x in sc_2) sum_sc_2 += sc_2[x]; if ( NR==1 || ( (sum_ic_2 + sum_sc_2 > 0) && (sum_ic_1 + sum_sc_1 > 0) && (sum_ic_1 != 0 || sum_ic_2 != 0) && (sum_sc_1 != 0 || sum_sc_2 != 0) && $6!=0 && $7!=0 ) ) {print $0}}' JCEC.RNASeq.RI.MATS.input.txt > RI.JCEC.input.txt
/home/soft/rMATS.3.2.5/MATS/rMATS.sh -d RI.JCEC.input.txt -o RIout_JCEC/ -c 0.0001 -p 8 -t U
awk '{ split($2,ic_1,","); sum_ic_1=0; for (x in ic_1) sum_ic_1 += ic_1[x]; split($4,ic_2,","); sum_ic_2=0; for (x in ic_2) sum_ic_2 += ic_2[x]; split($3,sc_1,","); sum_sc_1=0; for (x in sc_1) sum_sc_1 += sc_1[x]; split($5,sc_2,","); sum_sc_2=0; for (x in sc_2) sum_sc_2 += sc_2[x]; if ( NR==1 || ( (sum_ic_2 + sum_sc_2 > 0) && (sum_ic_1 + sum_sc_1 > 0) && (sum_ic_1 != 0 || sum_ic_2 != 0) && (sum_sc_1 != 0 || sum_sc_2 != 0) && $6!=0 && $7!=0 ) ) {print $0}}' JCEC.RNASeq.SE.MATS.input.txt > SE.JCEC.input.txt
/home/soft/rMATS.3.2.5/MATS/rMATS.sh -d SE.JCEC.input.txt -o SEout_JCEC/ -c 0.0001 -p 8 -t U
awk '{ split($2,ic_1,","); sum_ic_1=0; for (x in ic_1) sum_ic_1 += ic_1[x]; split($4,ic_2,","); sum_ic_2=0; for (x in ic_2) sum_ic_2 += ic_2[x]; split($3,sc_1,","); sum_sc_1=0; for (x in sc_1) sum_sc_1 += sc_1[x]; split($5,sc_2,","); sum_sc_2=0; for (x in sc_2) sum_sc_2 += sc_2[x]; if ( NR==1 || ( (sum_ic_2 + sum_sc_2 > 0) && (sum_ic_1 + sum_sc_1 > 0) && (sum_ic_1 != 0 || sum_ic_2 != 0) && (sum_sc_1 != 0 || sum_sc_2 != 0) && $6!=0 && $7!=0 ) ) {print $0}}' JCEC.RNASeq.MXE.MATS.input.txt > MXE.JCEC.input.txt
/home/soft/rMATS.3.2.5/MATS/rMATS.sh -d MXE.JCEC.input.txt -o MXEout_JCEC/ -c 0.0001 -p 8 -t U
awk '{ split($2,ic_1,","); sum_ic_1=0; for (x in ic_1) sum_ic_1 += ic_1[x]; split($4,ic_2,","); sum_ic_2=0; for (x in ic_2) sum_ic_2 += ic_2[x]; split($3,sc_1,","); sum_sc_1=0; for (x in sc_1) sum_sc_1 += sc_1[x]; split($5,sc_2,","); sum_sc_2=0; for (x in sc_2) sum_sc_2 += sc_2[x]; if ( NR==1 || ( (sum_ic_2 + sum_sc_2 > 0) && (sum_ic_1 + sum_sc_1 > 0) && (sum_ic_1 != 0 || sum_ic_2 != 0) && (sum_sc_1 != 0 || sum_sc_2 != 0) && $6!=0 && $7!=0 ) ) {print $0}}' JCEC.RNASeq.A5SS.MATS.input.txt > A5SS.JCEC.input.txt
/home/soft/rMATS.3.2.5/MATS/rMATS.sh -d A5SS.JCEC.input.txt -o A5SSout_JCEC/ -c 0.0001 -p 8 -t U
awk '{ split($2,ic_1,","); sum_ic_1=0; for (x in ic_1) sum_ic_1 += ic_1[x]; split($4,ic_2,","); sum_ic_2=0; for (x in ic_2) sum_ic_2 += ic_2[x]; split($3,sc_1,","); sum_sc_1=0; for (x in sc_1) sum_sc_1 += sc_1[x]; split($5,sc_2,","); sum_sc_2=0; for (x in sc_2) sum_sc_2 += sc_2[x]; if ( NR==1 || ( (sum_ic_2 + sum_sc_2 > 0) && (sum_ic_1 + sum_sc_1 > 0) && (sum_ic_1 != 0 || sum_ic_2 != 0) && (sum_sc_1 != 0 || sum_sc_2 != 0) && $6!=0 && $7!=0 ) ) {print $0}}' JCEC.RNASeq.A3SS.MATS.input.txt > A3SS.JCEC.input.txt
/home/soft/rMATS.3.2.5/MATS/rMATS.sh -d A3SS.JCEC.input.txt -o A3SSout_JCEC/ -c 0.0001 -p 8 -t U
cd ..
rm -rf Shoot_H24piMinus_vs_H24piPlus
mkdir Shoot_H24piMinus_vs_H24piPlus
cd Shoot_H24piMinus_vs_H24piPlus
Rscript /home/pcd/merge_ASE_for_diff.r "H24piPlusR1_SRR1005359,H24piPlusR2_SRR1005360,H24piPlusR3_SRR1005361" "H24piMinusR1_SRR1005362,H24piMinusR2_SRR1005363,H24piMinusR3_SRR1005364"
awk '{ split($2,ic_1,","); sum_ic_1=0; for (x in ic_1) sum_ic_1 += ic_1[x]; split($4,ic_2,","); sum_ic_2=0; for (x in ic_2) sum_ic_2 += ic_2[x]; split($3,sc_1,","); sum_sc_1=0; for (x in sc_1) sum_sc_1 += sc_1[x]; split($5,sc_2,","); sum_sc_2=0; for (x in sc_2) sum_sc_2 += sc_2[x]; if ( NR==1 || ( (sum_ic_2 + sum_sc_2 > 0) && (sum_ic_1 + sum_sc_1 > 0) && (sum_ic_1 != 0 || sum_ic_2 != 0) && (sum_sc_1 != 0 || sum_sc_2 != 0) && $6!=0 && $7!=0 ) ) {print $0}}' JCEC.RNASeq.RI.MATS.input.txt > RI.JCEC.input.txt
/home/soft/rMATS.3.2.5/MATS/rMATS.sh -d RI.JCEC.input.txt -o RIout_JCEC/ -c 0.0001 -p 8 -t U
awk '{ split($2,ic_1,","); sum_ic_1=0; for (x in ic_1) sum_ic_1 += ic_1[x]; split($4,ic_2,","); sum_ic_2=0; for (x in ic_2) sum_ic_2 += ic_2[x]; split($3,sc_1,","); sum_sc_1=0; for (x in sc_1) sum_sc_1 += sc_1[x]; split($5,sc_2,","); sum_sc_2=0; for (x in sc_2) sum_sc_2 += sc_2[x]; if ( NR==1 || ( (sum_ic_2 + sum_sc_2 > 0) && (sum_ic_1 + sum_sc_1 > 0) && (sum_ic_1 != 0 || sum_ic_2 != 0) && (sum_sc_1 != 0 || sum_sc_2 != 0) && $6!=0 && $7!=0 ) ) {print $0}}' JCEC.RNASeq.SE.MATS.input.txt > SE.JCEC.input.txt
/home/soft/rMATS.3.2.5/MATS/rMATS.sh -d SE.JCEC.input.txt -o SEout_JCEC/ -c 0.0001 -p 8 -t U
awk '{ split($2,ic_1,","); sum_ic_1=0; for (x in ic_1) sum_ic_1 += ic_1[x]; split($4,ic_2,","); sum_ic_2=0; for (x in ic_2) sum_ic_2 += ic_2[x]; split($3,sc_1,","); sum_sc_1=0; for (x in sc_1) sum_sc_1 += sc_1[x]; split($5,sc_2,","); sum_sc_2=0; for (x in sc_2) sum_sc_2 += sc_2[x]; if ( NR==1 || ( (sum_ic_2 + sum_sc_2 > 0) && (sum_ic_1 + sum_sc_1 > 0) && (sum_ic_1 != 0 || sum_ic_2 != 0) && (sum_sc_1 != 0 || sum_sc_2 != 0) && $6!=0 && $7!=0 ) ) {print $0}}' JCEC.RNASeq.MXE.MATS.input.txt > MXE.JCEC.input.txt
/home/soft/rMATS.3.2.5/MATS/rMATS.sh -d MXE.JCEC.input.txt -o MXEout_JCEC/ -c 0.0001 -p 8 -t U
awk '{ split($2,ic_1,","); sum_ic_1=0; for (x in ic_1) sum_ic_1 += ic_1[x]; split($4,ic_2,","); sum_ic_2=0; for (x in ic_2) sum_ic_2 += ic_2[x]; split($3,sc_1,","); sum_sc_1=0; for (x in sc_1) sum_sc_1 += sc_1[x]; split($5,sc_2,","); sum_sc_2=0; for (x in sc_2) sum_sc_2 += sc_2[x]; if ( NR==1 || ( (sum_ic_2 + sum_sc_2 > 0) && (sum_ic_1 + sum_sc_1 > 0) && (sum_ic_1 != 0 || sum_ic_2 != 0) && (sum_sc_1 != 0 || sum_sc_2 != 0) && $6!=0 && $7!=0 ) ) {print $0}}' JCEC.RNASeq.A5SS.MATS.input.txt > A5SS.JCEC.input.txt
/home/soft/rMATS.3.2.5/MATS/rMATS.sh -d A5SS.JCEC.input.txt -o A5SSout_JCEC/ -c 0.0001 -p 8-t U
awk '{ split($2,ic_1,","); sum_ic_1=0; for (x in ic_1) sum_ic_1 += ic_1[x]; split($4,ic_2,","); sum_ic_2=0; for (x in ic_2) sum_ic_2 += ic_2[x]; split($3,sc_1,","); sum_sc_1=0; for (x in sc_1) sum_sc_1 += sc_1[x]; split($5,sc_2,","); sum_sc_2=0; for (x in sc_2) sum_sc_2 += sc_2[x]; if ( NR==1 || ( (sum_ic_2 + sum_sc_2 > 0) && (sum_ic_1 + sum_sc_1 > 0) && (sum_ic_1 != 0 || sum_ic_2 != 0) && (sum_sc_1 != 0 || sum_sc_2 != 0) && $6!=0 && $7!=0 ) ) {print $0}}' JCEC.RNASeq.A3SS.MATS.input.txt > A3SS.JCEC.input.txt
/home/soft/rMATS.3.2.5/MATS/rMATS.sh -d A3SS.JCEC.input.txt -o A3SSout_JCEC/ -c 0.0001 -p 8 -t U
cd ..
rm -rf Shoot_D3piMinus_vs_D3piPlus
mkdir Shoot_D3piMinus_vs_D3piPlus
cd Shoot_D3piMinus_vs_D3piPlus
Rscript /home/pcd/merge_ASE_for_diff.r "D3piPlusR1_SRR1005365,D3piPlusR2_SRR1005366,D3piPlusR3_SRR1005367" "D3piMinusR1_SRR1005368,D3piMinusR2_SRR1005369,D3piMinusR3_SRR1005370"
awk '{ split($2,ic_1,","); sum_ic_1=0; for (x in ic_1) sum_ic_1 += ic_1[x]; split($4,ic_2,","); sum_ic_2=0; for (x in ic_2) sum_ic_2 += ic_2[x]; split($3,sc_1,","); sum_sc_1=0; for (x in sc_1) sum_sc_1 += sc_1[x]; split($5,sc_2,","); sum_sc_2=0; for (x in sc_2) sum_sc_2 += sc_2[x]; if ( NR==1 || ( (sum_ic_2 + sum_sc_2 > 0) && (sum_ic_1 + sum_sc_1 > 0) && (sum_ic_1 != 0 || sum_ic_2 != 0) && (sum_sc_1 != 0 || sum_sc_2 != 0) && $6!=0 && $7!=0 ) ) {print $0}}' JCEC.RNASeq.RI.MATS.input.txt > RI.JCEC.input.txt
/home/soft/rMATS.3.2.5/MATS/rMATS.sh -d RI.JCEC.input.txt -o RIout_JCEC/ -c 0.0001 -p 8 -t U
awk '{ split($2,ic_1,","); sum_ic_1=0; for (x in ic_1) sum_ic_1 += ic_1[x]; split($4,ic_2,","); sum_ic_2=0; for (x in ic_2) sum_ic_2 += ic_2[x]; split($3,sc_1,","); sum_sc_1=0; for (x in sc_1) sum_sc_1 += sc_1[x]; split($5,sc_2,","); sum_sc_2=0; for (x in sc_2) sum_sc_2 += sc_2[x]; if ( NR==1 || ( (sum_ic_2 + sum_sc_2 > 0) && (sum_ic_1 + sum_sc_1 > 0) && (sum_ic_1 != 0 || sum_ic_2 != 0) && (sum_sc_1 != 0 || sum_sc_2 != 0) && $6!=0 && $7!=0 ) ) {print $0}}' JCEC.RNASeq.SE.MATS.input.txt > SE.JCEC.input.txt
/home/soft/rMATS.3.2.5/MATS/rMATS.sh -d SE.JCEC.input.txt -o SEout_JCEC/ -c 0.0001 -p 8 -t U
awk '{ split($2,ic_1,","); sum_ic_1=0; for (x in ic_1) sum_ic_1 += ic_1[x]; split($4,ic_2,","); sum_ic_2=0; for (x in ic_2) sum_ic_2 += ic_2[x]; split($3,sc_1,","); sum_sc_1=0; for (x in sc_1) sum_sc_1 += sc_1[x]; split($5,sc_2,","); sum_sc_2=0; for (x in sc_2) sum_sc_2 += sc_2[x]; if ( NR==1 || ( (sum_ic_2 + sum_sc_2 > 0) && (sum_ic_1 + sum_sc_1 > 0) && (sum_ic_1 != 0 || sum_ic_2 != 0) && (sum_sc_1 != 0 || sum_sc_2 != 0) && $6!=0 && $7!=0 ) ) {print $0}}' JCEC.RNASeq.MXE.MATS.input.txt > MXE.JCEC.input.txt
/home/soft/rMATS.3.2.5/MATS/rMATS.sh -d MXE.JCEC.input.txt -o MXEout_JCEC/ -c 0.0001 -p 8 -t U
awk '{ split($2,ic_1,","); sum_ic_1=0; for (x in ic_1) sum_ic_1 += ic_1[x]; split($4,ic_2,","); sum_ic_2=0; for (x in ic_2) sum_ic_2 += ic_2[x]; split($3,sc_1,","); sum_sc_1=0; for (x in sc_1) sum_sc_1 += sc_1[x]; split($5,sc_2,","); sum_sc_2=0; for (x in sc_2) sum_sc_2 += sc_2[x]; if ( NR==1 || ( (sum_ic_2 + sum_sc_2 > 0) && (sum_ic_1 + sum_sc_1 > 0) && (sum_ic_1 != 0 || sum_ic_2 != 0) && (sum_sc_1 != 0 || sum_sc_2 != 0) && $6!=0 && $7!=0 ) ) {print $0}}' JCEC.RNASeq.A5SS.MATS.input.txt > A5SS.JCEC.input.txt
/home/soft/rMATS.3.2.5/MATS/rMATS.sh -d A5SS.JCEC.input.txt -o A5SSout_JCEC/ -c 0.0001 -p 8 -t U
awk '{ split($2,ic_1,","); sum_ic_1=0; for (x in ic_1) sum_ic_1 += ic_1[x]; split($4,ic_2,","); sum_ic_2=0; for (x in ic_2) sum_ic_2 += ic_2[x]; split($3,sc_1,","); sum_sc_1=0; for (x in sc_1) sum_sc_1 += sc_1[x]; split($5,sc_2,","); sum_sc_2=0; for (x in sc_2) sum_sc_2 += sc_2[x]; if ( NR==1 || ( (sum_ic_2 + sum_sc_2 > 0) && (sum_ic_1 + sum_sc_1 > 0) && (sum_ic_1 != 0 || sum_ic_2 != 0) && (sum_sc_1 != 0 || sum_sc_2 != 0) && $6!=0 && $7!=0 ) ) {print $0}}' JCEC.RNASeq.A3SS.MATS.input.txt > A3SS.JCEC.input.txt
/home/soft/rMATS.3.2.5/MATS/rMATS.sh -d A3SS.JCEC.input.txt -o A3SSout_JCEC/ -c 0.0001 -p 8 -t U
cd ..
rm -rf Shoot_D7piMinus_vs_D7piPlus
mkdir Shoot_D7piMinus_vs_D7piPlus
cd Shoot_D7piMinus_vs_D7piPlus
Rscript /home/pcd/merge_ASE_for_diff.r "D7piPlusR1_SRR1005377,D7piPlusR2_SRR1005378,D7piPlusR3_SRR1005379" "D7piMinusR1_SRR1005380,D7piMinusR2_SRR1005381,D7piMinusR3_SRR1005382"
awk '{ split($2,ic_1,","); sum_ic_1=0; for (x in ic_1) sum_ic_1 += ic_1[x]; split($4,ic_2,","); sum_ic_2=0; for (x in ic_2) sum_ic_2 += ic_2[x]; split($3,sc_1,","); sum_sc_1=0; for (x in sc_1) sum_sc_1 += sc_1[x]; split($5,sc_2,","); sum_sc_2=0; for (x in sc_2) sum_sc_2 += sc_2[x]; if ( NR==1 || ( (sum_ic_2 + sum_sc_2 > 0) && (sum_ic_1 + sum_sc_1 > 0) && (sum_ic_1 != 0 || sum_ic_2 != 0) && (sum_sc_1 != 0 || sum_sc_2 != 0) && $6!=0 && $7!=0 ) ) {print $0}}' JCEC.RNASeq.RI.MATS.input.txt > RI.JCEC.input.txt
/home/soft/rMATS.3.2.5/MATS/rMATS.sh -d RI.JCEC.input.txt -o RIout_JCEC/ -c 0.0001 -p 8 -t U
awk '{ split($2,ic_1,","); sum_ic_1=0; for (x in ic_1) sum_ic_1 += ic_1[x]; split($4,ic_2,","); sum_ic_2=0; for (x in ic_2) sum_ic_2 += ic_2[x]; split($3,sc_1,","); sum_sc_1=0; for (x in sc_1) sum_sc_1 += sc_1[x]; split($5,sc_2,","); sum_sc_2=0; for (x in sc_2) sum_sc_2 += sc_2[x]; if ( NR==1 || ( (sum_ic_2 + sum_sc_2 > 0) && (sum_ic_1 + sum_sc_1 > 0) && (sum_ic_1 != 0 || sum_ic_2 != 0) && (sum_sc_1 != 0 || sum_sc_2 != 0) && $6!=0 && $7!=0 ) ) {print $0}}' JCEC.RNASeq.SE.MATS.input.txt > SE.JCEC.input.txt
/home/soft/rMATS.3.2.5/MATS/rMATS.sh -d SE.JCEC.input.txt -o SEout_JCEC/ -c 0.0001 -p 8 -t U
awk '{ split($2,ic_1,","); sum_ic_1=0; for (x in ic_1) sum_ic_1 += ic_1[x]; split($4,ic_2,","); sum_ic_2=0; for (x in ic_2) sum_ic_2 += ic_2[x]; split($3,sc_1,","); sum_sc_1=0; for (x in sc_1) sum_sc_1 += sc_1[x]; split($5,sc_2,","); sum_sc_2=0; for (x in sc_2) sum_sc_2 += sc_2[x]; if ( NR==1 || ( (sum_ic_2 + sum_sc_2 > 0) && (sum_ic_1 + sum_sc_1 > 0) && (sum_ic_1 != 0 || sum_ic_2 != 0) && (sum_sc_1 != 0 || sum_sc_2 != 0) && $6!=0 && $7!=0 ) ) {print $0}}' JCEC.RNASeq.MXE.MATS.input.txt > MXE.JCEC.input.txt
/home/soft/rMATS.3.2.5/MATS/rMATS.sh -d MXE.JCEC.input.txt -o MXEout_JCEC/ -c 0.0001 -p 8 -t U
awk '{ split($2,ic_1,","); sum_ic_1=0; for (x in ic_1) sum_ic_1 += ic_1[x]; split($4,ic_2,","); sum_ic_2=0; for (x in ic_2) sum_ic_2 += ic_2[x]; split($3,sc_1,","); sum_sc_1=0; for (x in sc_1) sum_sc_1 += sc_1[x]; split($5,sc_2,","); sum_sc_2=0; for (x in sc_2) sum_sc_2 += sc_2[x]; if ( NR==1 || ( (sum_ic_2 + sum_sc_2 > 0) && (sum_ic_1 + sum_sc_1 > 0) && (sum_ic_1 != 0 || sum_ic_2 != 0) && (sum_sc_1 != 0 || sum_sc_2 != 0) && $6!=0 && $7!=0 ) ) {print $0}}' JCEC.RNASeq.A5SS.MATS.input.txt > A5SS.JCEC.input.txt
/home/soft/rMATS.3.2.5/MATS/rMATS.sh -d A5SS.JCEC.input.txt -o A5SSout_JCEC/ -c 0.0001 -p 8 -t U
awk '{ split($2,ic_1,","); sum_ic_1=0; for (x in ic_1) sum_ic_1 += ic_1[x]; split($4,ic_2,","); sum_ic_2=0; for (x in ic_2) sum_ic_2 += ic_2[x]; split($3,sc_1,","); sum_sc_1=0; for (x in sc_1) sum_sc_1 += sc_1[x]; split($5,sc_2,","); sum_sc_2=0; for (x in sc_2) sum_sc_2 += sc_2[x]; if ( NR==1 || ( (sum_ic_2 + sum_sc_2 > 0) && (sum_ic_1 + sum_sc_1 > 0) && (sum_ic_1 != 0 || sum_ic_2 != 0) && (sum_sc_1 != 0 || sum_sc_2 != 0) && $6!=0 && $7!=0 ) ) {print $0}}' JCEC.RNASeq.A3SS.MATS.input.txt > A3SS.JCEC.input.txt
/home/soft/rMATS.3.2.5/MATS/rMATS.sh -d A3SS.JCEC.input.txt -o A3SSout_JCEC/ -c 0.0001 -p 8 -t U
cd ..

rm -rf Shoot_D21piMinus_vs_D21piPlus
mkdir Shoot_D21piMinus_vs_D21piPlus
cd Shoot_D21piMinus_vs_D21piPlus
Rscript /home/pcd/merge_ASE_for_diff.r "D21piPlusR1_SRR1005353,D21piPlusR2_SRR1005354,D21piPlusR3_SRR1005355" "D21piMinusR1_SRR1005356,D21piMinusR2_SRR1005357,D21piMinusR3_SRR1005358"
awk '{ split($2,ic_1,","); sum_ic_1=0; for (x in ic_1) sum_ic_1 += ic_1[x]; split($4,ic_2,","); sum_ic_2=0; for (x in ic_2) sum_ic_2 += ic_2[x]; split($3,sc_1,","); sum_sc_1=0; for (x in sc_1) sum_sc_1 += sc_1[x]; split($5,sc_2,","); sum_sc_2=0; for (x in sc_2) sum_sc_2 += sc_2[x]; if ( NR==1 || ( (sum_ic_2 + sum_sc_2 > 0) && (sum_ic_1 + sum_sc_1 > 0) && (sum_ic_1 != 0 || sum_ic_2 != 0) && (sum_sc_1 != 0 || sum_sc_2 != 0) && $6!=0 && $7!=0 ) ) {print $0}}' JCEC.RNASeq.RI.MATS.input.txt > RI.JCEC.input.txt
/home/soft/rMATS.3.2.5/MATS/rMATS.sh -d RI.JCEC.input.txt -o RIout_JCEC/ -c 0.0001 -p 8 -t U
awk '{ split($2,ic_1,","); sum_ic_1=0; for (x in ic_1) sum_ic_1 += ic_1[x]; split($4,ic_2,","); sum_ic_2=0; for (x in ic_2) sum_ic_2 += ic_2[x]; split($3,sc_1,","); sum_sc_1=0; for (x in sc_1) sum_sc_1 += sc_1[x]; split($5,sc_2,","); sum_sc_2=0; for (x in sc_2) sum_sc_2 += sc_2[x]; if ( NR==1 || ( (sum_ic_2 + sum_sc_2 > 0) && (sum_ic_1 + sum_sc_1 > 0) && (sum_ic_1 != 0 || sum_ic_2 != 0) && (sum_sc_1 != 0 || sum_sc_2 != 0) && $6!=0 && $7!=0 ) ) {print $0}}' JCEC.RNASeq.SE.MATS.input.txt > SE.JCEC.input.txt
/home/soft/rMATS.3.2.5/MATS/rMATS.sh -d SE.JCEC.input.txt -o SEout_JCEC/ -c 0.0001 -p 8 -t U
awk '{ split($2,ic_1,","); sum_ic_1=0; for (x in ic_1) sum_ic_1 += ic_1[x]; split($4,ic_2,","); sum_ic_2=0; for (x in ic_2) sum_ic_2 += ic_2[x]; split($3,sc_1,","); sum_sc_1=0; for (x in sc_1) sum_sc_1 += sc_1[x]; split($5,sc_2,","); sum_sc_2=0; for (x in sc_2) sum_sc_2 += sc_2[x]; if ( NR==1 || ( (sum_ic_2 + sum_sc_2 > 0) && (sum_ic_1 + sum_sc_1 > 0) && (sum_ic_1 != 0 || sum_ic_2 != 0) && (sum_sc_1 != 0 || sum_sc_2 != 0) && $6!=0 && $7!=0 ) ) {print $0}}' JCEC.RNASeq.MXE.MATS.input.txt > MXE.JCEC.input.txt
/home/soft/rMATS.3.2.5/MATS/rMATS.sh -d MXE.JCEC.input.txt -o MXEout_JCEC/ -c 0.0001 -p 8 -t U
awk '{ split($2,ic_1,","); sum_ic_1=0; for (x in ic_1) sum_ic_1 += ic_1[x]; split($4,ic_2,","); sum_ic_2=0; for (x in ic_2) sum_ic_2 += ic_2[x]; split($3,sc_1,","); sum_sc_1=0; for (x in sc_1) sum_sc_1 += sc_1[x]; split($5,sc_2,","); sum_sc_2=0; for (x in sc_2) sum_sc_2 += sc_2[x]; if ( NR==1 || ( (sum_ic_2 + sum_sc_2 > 0) && (sum_ic_1 + sum_sc_1 > 0) && (sum_ic_1 != 0 || sum_ic_2 != 0) && (sum_sc_1 != 0 || sum_sc_2 != 0) && $6!=0 && $7!=0 ) ) {print $0}}' JCEC.RNASeq.A5SS.MATS.input.txt > A5SS.JCEC.input.txt
/home/soft/rMATS.3.2.5/MATS/rMATS.sh -d A5SS.JCEC.input.txt -o A5SSout_JCEC/ -c 0.0001 -p 8 -t U
awk '{ split($2,ic_1,","); sum_ic_1=0; for (x in ic_1) sum_ic_1 += ic_1[x]; split($4,ic_2,","); sum_ic_2=0; for (x in ic_2) sum_ic_2 += ic_2[x]; split($3,sc_1,","); sum_sc_1=0; for (x in sc_1) sum_sc_1 += sc_1[x]; split($5,sc_2,","); sum_sc_2=0; for (x in sc_2) sum_sc_2 += sc_2[x]; if ( NR==1 || ( (sum_ic_2 + sum_sc_2 > 0) && (sum_ic_1 + sum_sc_1 > 0) && (sum_ic_1 != 0 || sum_ic_2 != 0) && (sum_sc_1 != 0 || sum_sc_2 != 0) && $6!=0 && $7!=0 ) ) {print $0}}' JCEC.RNASeq.A3SS.MATS.input.txt > A3SS.JCEC.input.txt
/home/soft/rMATS.3.2.5/MATS/rMATS.sh -d A3SS.JCEC.input.txt -o A3SSout_JCEC/ -c 0.0001 -p 8 -t U
cd ..
rm -rf Shoot_D21H1piMinus_vs_D21H1piPlus
mkdir Shoot_D21H1piMinus_vs_D21H1piPlus
cd Shoot_D21H1piMinus_vs_D21H1piPlus
Rscript /home/pcd/merge_ASE_for_diff.r "D21H1piPlusR1_SRR1005326,D21H1piPlusR2_SRR1005327,D21H1piPlusR3_SRR1005328" "D21H1piMinusR1_SRR1005329,D21H1piMinusR2_SRR1005330,D21H1piMinusR3_SRR1005331"
awk '{ split($2,ic_1,","); sum_ic_1=0; for (x in ic_1) sum_ic_1 += ic_1[x]; split($4,ic_2,","); sum_ic_2=0; for (x in ic_2) sum_ic_2 += ic_2[x]; split($3,sc_1,","); sum_sc_1=0; for (x in sc_1) sum_sc_1 += sc_1[x]; split($5,sc_2,","); sum_sc_2=0; for (x in sc_2) sum_sc_2 += sc_2[x]; if ( NR==1 || ( (sum_ic_2 + sum_sc_2 > 0) && (sum_ic_1 + sum_sc_1 > 0) && (sum_ic_1 != 0 || sum_ic_2 != 0) && (sum_sc_1 != 0 || sum_sc_2 != 0) && $6!=0 && $7!=0 ) ) {print $0}}' JCEC.RNASeq.RI.MATS.input.txt > RI.JCEC.input.txt
/home/soft/rMATS.3.2.5/MATS/rMATS.sh -d RI.JCEC.input.txt -o RIout_JCEC/ -c 0.0001 -p 8 -t U
awk '{ split($2,ic_1,","); sum_ic_1=0; for (x in ic_1) sum_ic_1 += ic_1[x]; split($4,ic_2,","); sum_ic_2=0; for (x in ic_2) sum_ic_2 += ic_2[x]; split($3,sc_1,","); sum_sc_1=0; for (x in sc_1) sum_sc_1 += sc_1[x]; split($5,sc_2,","); sum_sc_2=0; for (x in sc_2) sum_sc_2 += sc_2[x]; if ( NR==1 || ( (sum_ic_2 + sum_sc_2 > 0) && (sum_ic_1 + sum_sc_1 > 0) && (sum_ic_1 != 0 || sum_ic_2 != 0) && (sum_sc_1 != 0 || sum_sc_2 != 0) && $6!=0 && $7!=0 ) ) {print $0}}' JCEC.RNASeq.SE.MATS.input.txt > SE.JCEC.input.txt
/home/soft/rMATS.3.2.5/MATS/rMATS.sh -d SE.JCEC.input.txt -o SEout_JCEC/ -c 0.0001 -p 8 -t U
awk '{ split($2,ic_1,","); sum_ic_1=0; for (x in ic_1) sum_ic_1 += ic_1[x]; split($4,ic_2,","); sum_ic_2=0; for (x in ic_2) sum_ic_2 += ic_2[x]; split($3,sc_1,","); sum_sc_1=0; for (x in sc_1) sum_sc_1 += sc_1[x]; split($5,sc_2,","); sum_sc_2=0; for (x in sc_2) sum_sc_2 += sc_2[x]; if ( NR==1 || ( (sum_ic_2 + sum_sc_2 > 0) && (sum_ic_1 + sum_sc_1 > 0) && (sum_ic_1 != 0 || sum_ic_2 != 0) && (sum_sc_1 != 0 || sum_sc_2 != 0) && $6!=0 && $7!=0 ) ) {print $0}}' JCEC.RNASeq.MXE.MATS.input.txt > MXE.JCEC.input.txt
/home/soft/rMATS.3.2.5/MATS/rMATS.sh -d MXE.JCEC.input.txt -o MXEout_JCEC/ -c 0.0001 -p 8 -t U
awk '{ split($2,ic_1,","); sum_ic_1=0; for (x in ic_1) sum_ic_1 += ic_1[x]; split($4,ic_2,","); sum_ic_2=0; for (x in ic_2) sum_ic_2 += ic_2[x]; split($3,sc_1,","); sum_sc_1=0; for (x in sc_1) sum_sc_1 += sc_1[x]; split($5,sc_2,","); sum_sc_2=0; for (x in sc_2) sum_sc_2 += sc_2[x]; if ( NR==1 || ( (sum_ic_2 + sum_sc_2 > 0) && (sum_ic_1 + sum_sc_1 > 0) && (sum_ic_1 != 0 || sum_ic_2 != 0) && (sum_sc_1 != 0 || sum_sc_2 != 0) && $6!=0 && $7!=0 ) ) {print $0}}' JCEC.RNASeq.A5SS.MATS.input.txt > A5SS.JCEC.input.txt
/home/soft/rMATS.3.2.5/MATS/rMATS.sh -d A5SS.JCEC.input.txt -o A5SSout_JCEC/ -c 0.0001 -p 8 -t U
awk '{ split($2,ic_1,","); sum_ic_1=0; for (x in ic_1) sum_ic_1 += ic_1[x]; split($4,ic_2,","); sum_ic_2=0; for (x in ic_2) sum_ic_2 += ic_2[x]; split($3,sc_1,","); sum_sc_1=0; for (x in sc_1) sum_sc_1 += sc_1[x]; split($5,sc_2,","); sum_sc_2=0; for (x in sc_2) sum_sc_2 += sc_2[x]; if ( NR==1 || ( (sum_ic_2 + sum_sc_2 > 0) && (sum_ic_1 + sum_sc_1 > 0) && (sum_ic_1 != 0 || sum_ic_2 != 0) && (sum_sc_1 != 0 || sum_sc_2 != 0) && $6!=0 && $7!=0 ) ) {print $0}}' JCEC.RNASeq.A3SS.MATS.input.txt > A3SS.JCEC.input.txt
/home/soft/rMATS.3.2.5/MATS/rMATS.sh -d A3SS.JCEC.input.txt -o A3SSout_JCEC/ -c 0.0001 -p 8 -t U
cd ..
rm -rf Shoot_D21H6piMinus_vs_D21H6piPlus
mkdir Shoot_D21H6piMinus_vs_D21H6piPlus
cd Shoot_D21H6piMinus_vs_D21H6piPlus
Rscript /home/pcd/merge_ASE_for_diff.r "D21H6piPlusR1_SRR1005344,D21H6piPlusR2_SRR1005345,D21H6piPlusR3_SRR1005346" "D21H6piMinusR1_SRR1005347,D21H6piMinusR2_SRR1005348,D21H6piMinusR3_SRR1005349"
awk '{ split($2,ic_1,","); sum_ic_1=0; for (x in ic_1) sum_ic_1 += ic_1[x]; split($4,ic_2,","); sum_ic_2=0; for (x in ic_2) sum_ic_2 += ic_2[x]; split($3,sc_1,","); sum_sc_1=0; for (x in sc_1) sum_sc_1 += sc_1[x]; split($5,sc_2,","); sum_sc_2=0; for (x in sc_2) sum_sc_2 += sc_2[x]; if ( NR==1 || ( (sum_ic_2 + sum_sc_2 > 0) && (sum_ic_1 + sum_sc_1 > 0) && (sum_ic_1 != 0 || sum_ic_2 != 0) && (sum_sc_1 != 0 || sum_sc_2 != 0) && $6!=0 && $7!=0 ) ) {print $0}}' JCEC.RNASeq.RI.MATS.input.txt > RI.JCEC.input.txt
/home/soft/rMATS.3.2.5/MATS/rMATS.sh -d RI.JCEC.input.txt -o RIout_JCEC/ -c 0.0001 -p 8 -t U
awk '{ split($2,ic_1,","); sum_ic_1=0; for (x in ic_1) sum_ic_1 += ic_1[x]; split($4,ic_2,","); sum_ic_2=0; for (x in ic_2) sum_ic_2 += ic_2[x]; split($3,sc_1,","); sum_sc_1=0; for (x in sc_1) sum_sc_1 += sc_1[x]; split($5,sc_2,","); sum_sc_2=0; for (x in sc_2) sum_sc_2 += sc_2[x]; if ( NR==1 || ( (sum_ic_2 + sum_sc_2 > 0) && (sum_ic_1 + sum_sc_1 > 0) && (sum_ic_1 != 0 || sum_ic_2 != 0) && (sum_sc_1 != 0 || sum_sc_2 != 0) && $6!=0 && $7!=0 ) ) {print $0}}' JCEC.RNASeq.SE.MATS.input.txt > SE.JCEC.input.txt
/home/soft/rMATS.3.2.5/MATS/rMATS.sh -d SE.JCEC.input.txt -o SEout_JCEC/ -c 0.0001 -p 8 -t U
awk '{ split($2,ic_1,","); sum_ic_1=0; for (x in ic_1) sum_ic_1 += ic_1[x]; split($4,ic_2,","); sum_ic_2=0; for (x in ic_2) sum_ic_2 += ic_2[x]; split($3,sc_1,","); sum_sc_1=0; for (x in sc_1) sum_sc_1 += sc_1[x]; split($5,sc_2,","); sum_sc_2=0; for (x in sc_2) sum_sc_2 += sc_2[x]; if ( NR==1 || ( (sum_ic_2 + sum_sc_2 > 0) && (sum_ic_1 + sum_sc_1 > 0) && (sum_ic_1 != 0 || sum_ic_2 != 0) && (sum_sc_1 != 0 || sum_sc_2 != 0) && $6!=0 && $7!=0 ) ) {print $0}}' JCEC.RNASeq.MXE.MATS.input.txt > MXE.JCEC.input.txt
/home/soft/rMATS.3.2.5/MATS/rMATS.sh -d MXE.JCEC.input.txt -o MXEout_JCEC/ -c 0.0001 -p 8 -t U
awk '{ split($2,ic_1,","); sum_ic_1=0; for (x in ic_1) sum_ic_1 += ic_1[x]; split($4,ic_2,","); sum_ic_2=0; for (x in ic_2) sum_ic_2 += ic_2[x]; split($3,sc_1,","); sum_sc_1=0; for (x in sc_1) sum_sc_1 += sc_1[x]; split($5,sc_2,","); sum_sc_2=0; for (x in sc_2) sum_sc_2 += sc_2[x]; if ( NR==1 || ( (sum_ic_2 + sum_sc_2 > 0) && (sum_ic_1 + sum_sc_1 > 0) && (sum_ic_1 != 0 || sum_ic_2 != 0) && (sum_sc_1 != 0 || sum_sc_2 != 0) && $6!=0 && $7!=0 ) ) {print $0}}' JCEC.RNASeq.A5SS.MATS.input.txt > A5SS.JCEC.input.txt
/home/soft/rMATS.3.2.5/MATS/rMATS.sh -d A5SS.JCEC.input.txt -o A5SSout_JCEC/ -c 0.0001 -p 8 -t U
awk '{ split($2,ic_1,","); sum_ic_1=0; for (x in ic_1) sum_ic_1 += ic_1[x]; split($4,ic_2,","); sum_ic_2=0; for (x in ic_2) sum_ic_2 += ic_2[x]; split($3,sc_1,","); sum_sc_1=0; for (x in sc_1) sum_sc_1 += sc_1[x]; split($5,sc_2,","); sum_sc_2=0; for (x in sc_2) sum_sc_2 += sc_2[x]; if ( NR==1 || ( (sum_ic_2 + sum_sc_2 > 0) && (sum_ic_1 + sum_sc_1 > 0) && (sum_ic_1 != 0 || sum_ic_2 != 0) && (sum_sc_1 != 0 || sum_sc_2 != 0) && $6!=0 && $7!=0 ) ) {print $0}}' JCEC.RNASeq.A3SS.MATS.input.txt > A3SS.JCEC.input.txt
/home/soft/rMATS.3.2.5/MATS/rMATS.sh -d A3SS.JCEC.input.txt -o A3SSout_JCEC/ -c 0.0001 -p 8 -t U
cd ..
rm -rf Shoot_D21H24piMinus_vs_D21H24piPlus
mkdir Shoot_D21H24piMinus_vs_D21H24piPlus
cd Shoot_D21H24piMinus_vs_D21H24piPlus
Rscript /home/pcd/merge_ASE_for_diff.r "D21H24piPlusR1_SRR1005335,D21H24piPlusR2_SRR1005336,D21H24piPlusR3_SRR1005337" "D21H24piMinusR1_SRR1005338,D21H24piMinusR2_SRR1005339,D21H24piMinusR3_SRR1005340"
awk '{ split($2,ic_1,","); sum_ic_1=0; for (x in ic_1) sum_ic_1 += ic_1[x]; split($4,ic_2,","); sum_ic_2=0; for (x in ic_2) sum_ic_2 += ic_2[x]; split($3,sc_1,","); sum_sc_1=0; for (x in sc_1) sum_sc_1 += sc_1[x]; split($5,sc_2,","); sum_sc_2=0; for (x in sc_2) sum_sc_2 += sc_2[x]; if ( NR==1 || ( (sum_ic_2 + sum_sc_2 > 0) && (sum_ic_1 + sum_sc_1 > 0) && (sum_ic_1 != 0 || sum_ic_2 != 0) && (sum_sc_1 != 0 || sum_sc_2 != 0) && $6!=0 && $7!=0 ) ) {print $0}}' JCEC.RNASeq.RI.MATS.input.txt > RI.JCEC.input.txt
/home/soft/rMATS.3.2.5/MATS/rMATS.sh -d RI.JCEC.input.txt -o RIout_JCEC/ -c 0.0001 -p 8 -t U
awk '{ split($2,ic_1,","); sum_ic_1=0; for (x in ic_1) sum_ic_1 += ic_1[x]; split($4,ic_2,","); sum_ic_2=0; for (x in ic_2) sum_ic_2 += ic_2[x]; split($3,sc_1,","); sum_sc_1=0; for (x in sc_1) sum_sc_1 += sc_1[x]; split($5,sc_2,","); sum_sc_2=0; for (x in sc_2) sum_sc_2 += sc_2[x]; if ( NR==1 || ( (sum_ic_2 + sum_sc_2 > 0) && (sum_ic_1 + sum_sc_1 > 0) && (sum_ic_1 != 0 || sum_ic_2 != 0) && (sum_sc_1 != 0 || sum_sc_2 != 0) && $6!=0 && $7!=0 ) ) {print $0}}' JCEC.RNASeq.SE.MATS.input.txt > SE.JCEC.input.txt
/home/soft/rMATS.3.2.5/MATS/rMATS.sh -d SE.JCEC.input.txt -o SEout_JCEC/ -c 0.0001 -p 8 -t U
awk '{ split($2,ic_1,","); sum_ic_1=0; for (x in ic_1) sum_ic_1 += ic_1[x]; split($4,ic_2,","); sum_ic_2=0; for (x in ic_2) sum_ic_2 += ic_2[x]; split($3,sc_1,","); sum_sc_1=0; for (x in sc_1) sum_sc_1 += sc_1[x]; split($5,sc_2,","); sum_sc_2=0; for (x in sc_2) sum_sc_2 += sc_2[x]; if ( NR==1 || ( (sum_ic_2 + sum_sc_2 > 0) && (sum_ic_1 + sum_sc_1 > 0) && (sum_ic_1 != 0 || sum_ic_2 != 0) && (sum_sc_1 != 0 || sum_sc_2 != 0) && $6!=0 && $7!=0 ) ) {print $0}}' JCEC.RNASeq.MXE.MATS.input.txt > MXE.JCEC.input.txt
/home/soft/rMATS.3.2.5/MATS/rMATS.sh -d MXE.JCEC.input.txt -o MXEout_JCEC/ -c 0.0001 -p 8 -t U
awk '{ split($2,ic_1,","); sum_ic_1=0; for (x in ic_1) sum_ic_1 += ic_1[x]; split($4,ic_2,","); sum_ic_2=0; for (x in ic_2) sum_ic_2 += ic_2[x]; split($3,sc_1,","); sum_sc_1=0; for (x in sc_1) sum_sc_1 += sc_1[x]; split($5,sc_2,","); sum_sc_2=0; for (x in sc_2) sum_sc_2 += sc_2[x]; if ( NR==1 || ( (sum_ic_2 + sum_sc_2 > 0) && (sum_ic_1 + sum_sc_1 > 0) && (sum_ic_1 != 0 || sum_ic_2 != 0) && (sum_sc_1 != 0 || sum_sc_2 != 0) && $6!=0 && $7!=0 ) ) {print $0}}' JCEC.RNASeq.A5SS.MATS.input.txt > A5SS.JCEC.input.txt
/home/soft/rMATS.3.2.5/MATS/rMATS.sh -d A5SS.JCEC.input.txt -o A5SSout_JCEC/ -c 0.0001 -p 8 -t U
awk '{ split($2,ic_1,","); sum_ic_1=0; for (x in ic_1) sum_ic_1 += ic_1[x]; split($4,ic_2,","); sum_ic_2=0; for (x in ic_2) sum_ic_2 += ic_2[x]; split($3,sc_1,","); sum_sc_1=0; for (x in sc_1) sum_sc_1 += sc_1[x]; split($5,sc_2,","); sum_sc_2=0; for (x in sc_2) sum_sc_2 += sc_2[x]; if ( NR==1 || ( (sum_ic_2 + sum_sc_2 > 0) && (sum_ic_1 + sum_sc_1 > 0) && (sum_ic_1 != 0 || sum_ic_2 != 0) && (sum_sc_1 != 0 || sum_sc_2 != 0) && $6!=0 && $7!=0 ) ) {print $0}}' JCEC.RNASeq.A3SS.MATS.input.txt > A3SS.JCEC.input.txt
/home/soft/rMATS.3.2.5/MATS/rMATS.sh -d A3SS.JCEC.input.txt -o A3SSout_JCEC/ -c 0.0001 -p 8 -t U
cd ..




# Identify DASE (FDR <= 0.05 and |DiffPSI| >= 0.1) under -Fe, -Zn, -Cu, -Mn, and -Pi, and GO, KEGG enrichment analysis
ASEdir=/public/analysis/hefei/Project_RNAseqstudy/study2/PiFeZnCuMn_merge_cufflinksV3/Stringtie/
annotation=/public/analysis/hefei/Project_RNAseqstudy/study2/PiFeZnCuMn_merge_cufflinksV3/Stringtie/merged_filter_multiexon_junction40_TPM_annotation.xls
ASEexprs=/public/analysis/hefei/Project_RNAseqstudy/study2/PiFeZnCuMn_merge_cufflinksV3/DASG/ASE/PiFeZnCuMn_in_ex_filter.xls
bedfile=/public/analysis/hefei/Project_RNAseqstudy/study2/PiFeZnCuMn_merge_cufflinksV3/Stringtie/merged_filter_multiexon_junction40_TPM.bed
unset arr
arr=("Fe_vs_Control;ControlR1,ControlR2,ControlR3;FeR1,FeR2,FeR3")
arr[${#arr[@]}]="Zn_vs_Control;ControlR1,ControlR2,ControlR3;ZnR1,ZnR2,ZnR3"
arr[${#arr[@]}]="Cu_vs_Control;ControlR1,ControlR2,ControlR3;CuR1,CuR2,CuR3"
arr[${#arr[@]}]="Mn_vs_Control;ControlR1,ControlR2,ControlR3;MnR1,MnR2,MnR3"
arr[${#arr[@]}]="Root_H1piMinus_vs_H1piPlus;H1piPlusR1,H1piPlusR2,H1piPlusR3;H1piMinusR1,H1piMinusR2,H1piMinusR3"
arr[${#arr[@]}]="Root_H6piMinus_vs_H6piPlus;H6piPlusR1,H6piPlusR2,H6piPlusR3;H6piMinusR1,H6piMinusR2,H6piMinusR3"
arr[${#arr[@]}]="Root_H24piMinus_vs_H24piPlus;H24piPlusR1,H24piPlusR2,H24piPlusR3;H24piMinusR1,H24piMinusR2,H24piMinusR3"
arr[${#arr[@]}]="Root_D3piMinus_vs_D3piPlus;D3piPlusR1,D3piPlusR2,D3piPlusR3;D3piMinusR1,D3piMinusR2,D3piMinusR3"
arr[${#arr[@]}]="Root_D7piMinus_vs_D7piPlus;D7piPlusR1,D7piPlusR2,D7piPlusR3;D7piMinusR1,D7piMinusR2,D7piMinusR3"
arr[${#arr[@]}]="Root_D21piMinus_vs_D21piPlus;D21piPlusR1,D21piPlusR2,D21piPlusR3;D21piMinusR1,D21piMinusR2,D21piMinusR3"
arr[${#arr[@]}]="Root_D21H1piMinus_vs_D21H1piPlus;D21H1piPlusR1,D21H1piPlusR2,D21H1piPlusR3;D21H1piMinusR1,D21H1piMinusR2,D21H1piMinusR3"
arr[${#arr[@]}]="Root_D21H6piMinus_vs_D21H6piPlus;D21H6piPlusR1,D21H6piPlusR2,D21H6piPlusR3;D21H6piMinusR2,D21H6piMinusR2,D21H6piMinusR3"
arr[${#arr[@]}]="Root_D21H24piMinus_vs_D21H24piPlus;D21H24piPlusR1,D21H24piPlusR2,D21H24piPlusR3;D21H24piMinusR1,D21H24piMinusR2,D21H24piMinusR3"
arr[${#arr[@]}]="Shoot_H1piMinus_vs_H1piPlus;H1piPlusR1,H1piPlusR2,H1piPlusR3;H1piMinusR1,H1piMinusR2,H1piMinusR3"
arr[${#arr[@]}]="Shoot_H6piMinus_vs_H6piPlus;H6piPlusR1,H6piPlusR2,H6piPlusR3;H6piMinusR1,H6piMinusR2,H6piMinusR3"
arr[${#arr[@]}]="Shoot_H24piMinus_vs_H24piPlus;H24piPlusR1,H24piPlusR2,H24piPlusR3;H24piMinusR1,H24piMinusR2,H24piMinusR3"
arr[${#arr[@]}]="Shoot_D3piMinus_vs_D3piPlus;D3piPlusR1,D3piPlusR2,D3piPlusR3;D3piMinusR1,D3piMinusR2,D3piMinusR3"
arr[${#arr[@]}]="Shoot_D7piMinus_vs_D7piPlus;D7piPlusR1,D7piPlusR2,D7piPlusR3;D7piMinusR1,D7piMinusR2,D7piMinusR3"
arr[${#arr[@]}]="Shoot_D21piMinus_vs_D21piPlus;D21piPlusR1,D21piPlusR2,D21piPlusR3;D21piMinusR1,D21piMinusR2,D21piMinusR3"
arr[${#arr[@]}]="Shoot_D21H1piMinus_vs_D21H1piPlus;D21H1piPlusR1,D21H1piPlusR2,D21H1piPlusR3;D21H1piMinusR1,D21H1piMinusR2,D21H1piMinusR3"
arr[${#arr[@]}]="Shoot_D21H6piMinus_vs_D21H6piPlus;D21H6piPlusR1,D21H6piPlusR2,D21H6piPlusR3;D21H6piMinusR2,D21H6piMinusR3,D21H6piMinusR3"
arr[${#arr[@]}]="Shoot_D21H24piMinus_vs_D21H24piPlus;D21H24piPlusR1,D21H24piPlusR2,D21H24piPlusR3;D21H24piMinusR1,D21H24piMinusR2,D21H24piMinusR3"
for i in ${arr[@]}
do
comparname=`echo $i | awk '{split($0,a,";" ); print a[1]}'`
cd $comparname
cd SEout_JCEC
perl /home/pcd/add_AS_aseid_v4.pl -dir /home/project/ -AStype "SE"
cd ../RIout_JCEC
perl /home/pcd/add_AS_aseid_v4.pl -dir /home/project/ -AStype "RI"
cd ../A5SSout_JCEC
perl /home/pcd/add_AS_aseid_v4.pl -dir /home/project/ -AStype "A5SS"
cd ../A3SSout_JCEC
perl /home/pcd/add_AS_aseid_v4.pl -dir /home/project/ -AStype "A3SS"
cd ../MXEout_JCEC
perl /home/pcd/add_AS_aseid_v4.pl -dir /home/project/ -AStype "MXE"
cd ..
rm -f AS*
perl /home/pcd/MATS_AS_exprs_each_condition_JC.pl -mode "JCEC" -samplelist $i -cutoff 3
Rscript /home/pcd/AS_diff_by_pvalue_psi.r  
perl /home/pcd/filter_ASE_by_id.pl -target /home/project/DASG/ASE/PiFeZnCuMn_in_ex_filter.xls -AS AS_diff.xls > AS_diff_filter.xls 
Rscript /home/pcd/AS_diff_PSI.r
perl /home/pcd/add_ASE_tid.pl -AS AS_diff_filter.xls > AS_diff_filter_with_major.xls
perl /home/pcd/add_anno_AS.pl -annotation /home/project/merged_filter_multiexon_junction40_TPM_annotation.xls -degene AS_diff_filter_with_major.xls > AS_diff_filter_with_major_annotation.xls
perl /home/pcd/degene_for_DASG_forGO.pl -degene AS_diff_filter_with_major_annotation.xls > diffgene.txt
perl /home/pcd/GO_enrich.pl -ispvalue "Pvalue" -bl /home/database/allgene.txt -gl diffgene.txt -symbol2GO /home/database/symbol2GO.txt
perl /home/pcd/KEGG_enrich.pl -ispvalue "Pvalue" -bl /home/database/allgene.txt -gl diffgene.txt -symbol2KEGG /home/database/symbol2KEGG.txt
cd ..
done



# DASG of 22 SRs
unset arr
arr[${#arr[@]}]="G41233-Os06g0187900-OsRSZ21a"
arr[${#arr[@]}]="G75476-Os11g0704700-OsSCL57"
arr[${#arr[@]}]="G34292-Os05g0162600-OsRS2Z39"
arr[${#arr[@]}]="G52301-Os07g0623300-OsSC32"
arr[${#arr[@]}]="G26033-Os04g0118900-OsRS29"
arr[${#arr[@]}]="G793-Os01g0155600-OsRS2Z37"
arr[${#arr[@]}]="G20979-Os03g0374575-OsSCL26"
arr[${#arr[@]}]="G52353-Os07g0633200-OsSCL25"
arr[${#arr[@]}]="G14932-Os02g0610600-OsRSZ23"
arr[${#arr[@]}]="G80417-Os12g0572400-OsSCL30"
arr[${#arr[@]}]="G9853-Os02g0122800-OsRS33"
arr[${#arr[@]}]="G33733-Os05g0120100-OsRS2Z36"
arr[${#arr[@]}]="G37219-Os05g0364600-OsSR33a"
arr[${#arr[@]}]="G20575-Os03g0344100-OsSR32"
arr[${#arr[@]}]="G20791-Os03g0363800-OsSCL28"
arr[${#arr[@]}]="G2706-Os01g0316600-OsSR40"
arr[${#arr[@]}]="G21073-Os03g0388000-OsSC25"
arr[${#arr[@]}]="G11532-Os02g0252100-OsSCL30a"
arr[${#arr[@]}]="G19901-Os03g0285900-OsRS2Z38"
arr[${#arr[@]}]="G58169-Os08g0486200-OsSC34"
arr[${#arr[@]}]="G52973-Os07g0673500-OsSR33"
arr[${#arr[@]}]="G17071-Os02g0789400-OsRSZ21"
unset gnames
gnames=${arr[0]}
for i in ${arr[@]:1}  
do
   gnames=$gnames","${i}
done;
perl /home/pcd/filter_DASE_by_targetids.pl -targetids $gnames -DASEtable PiRoot_DASE_overlap.xls > SRs_PiRoot_DASG_overlap.xls
perl /home/pcd/filter_DASE_by_targetids.pl -targetids $gnames -DASEtable PiShoot_DASE_overlap.xls > SRs_PiShoot_DASG_overlap.xls
perl /home/pcd/filter_DASE_by_targetids.pl -targetids $gnames -DASEtable FeZnCuMn_DASE_overlap.xls > SRs_FeZnCuMn_DASG_overlap.xls


