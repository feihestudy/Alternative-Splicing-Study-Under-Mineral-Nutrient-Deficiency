# ASE quantification for all 141 samples using rMATS

FeZnCuMnSampleids=("ControlR1" "ControlR2" "ControlR3" "FeR1" "FeR2" "FeR3" "ZnR1" "ZnR2" "ZnR3" "CuR1" "CuR2" "CuR3" "MnR1" "MnR2" "MnR3")
PiRootSamples=("H1piPlusR1_SRR1005257" "H1piPlusR2_SRR1005258" "H1piPlusR3_SRR1005259" "H1piMinusR1_SRR1005260" "H1piMinusR2_SRR1005261" "H1piMinusR3_SRR1005262" "H6piPlusR1_SRR1005308" "H6piPlusR2_SRR1005309" "H6piPlusR3_SRR1005310" "H6piMinusR1_SRR1005311" "H6piMinusR2_SRR1005312" "H6piMinusR3_SRR1005313" "H24piPlusR1_SRR1005296" "H24piPlusR2_SRR1005297" "H24piPlusR3_SRR1005298" "H24piMinusR1_SRR1005299" "H24piMinusR2_SRR1005300" "H24piMinusR3_SRR1005301" "D3piPlusR1_SRR1005302" "D3piPlusR2_SRR1005303" "D3piPlusR3_SRR1005304" "D3piMinusR1_SRR1005305" "D3piMinusR2_SRR1005306" "D3piMinusR3_SRR1005307" "D7piPlusR1_SRR1005314" "D7piPlusR2_SRR1005315" "D7piPlusR3_SRR1005316" "D7piMinusR1_SRR1005317" "D7piMinusR2_SRR1005318" "D7piMinusR3_SRR1005319" "D21piPlusR1_SRR1005290" "D21piPlusR2_SRR1005291" "D21piPlusR3_SRR1005292" "D21piMinusR1_SRR1005293" "D21piMinusR2_SRR1005294" "D21piMinusR3_SRR1005295" "D21H1piPlusR1_SRR1005263" "D21H1piPlusR2_SRR1005264" "D21H1piPlusR3_SRR1005265" "D21H1piMinusR1_SRR1005266" "D21H1piMinusR2_SRR1005267" "D21H1piMinusR3_SRR1005268" "D21H1piPlusRecR1_SRR1005269" "D21H1piPlusRecR2_SRR1005270" "D21H1piPlusRecR3_SRR1005271" "D21H6piPlusR1_SRR1005281" "D21H6piPlusR2_SRR1005282" "D21H6piPlusR3_SRR1005283" "D21H6piMinusR1_SRR1005284" "D21H6piMinusR2_SRR1005285" "D21H6piMinusR3_SRR1005286" "D21H6piPlusRecR1_SRR1005287" "D21H6piPlusRecR2_SRR1005288" "D21H6piPlusRecR3_SRR1005289" "D21H24piPlusR1_SRR1005272" "D21H24piPlusR2_SRR1005273" "D21H24piPlusR3_SRR1005274" "D21H24piMinusR1_SRR1005275" "D21H24piMinusR2_SRR1005276" "D21H24piMinusR3_SRR1005277" "D21H24piPlusRecR1_SRR1005278" "D21H24piPlusRecR2_SRR1005279" "D21H24piPlusRecR3_SRR1005280")
PiShootSamples=("H1piPlusR1_SRR1005320" "H1piPlusR2_SRR1005321" "H1piPlusR3_SRR1005322" "H1piMinusR1_SRR1005323" "H1piMinusR2_SRR1005324" "H1piMinusR3_SRR1005325" "H6piPlusR1_SRR1005371" "H6piPlusR2_SRR1005372" "H6piPlusR3_SRR1005373" "H6piMinusR1_SRR1005374" "H6piMinusR2_SRR1005375" "H6piMinusR3_SRR1005376" "H24piPlusR1_SRR1005359" "H24piPlusR2_SRR1005360" "H24piPlusR3_SRR1005361" "H24piMinusR1_SRR1005362" "H24piMinusR2_SRR1005363" "H24piMinusR3_SRR1005364" "D3piPlusR1_SRR1005365" "D3piPlusR2_SRR1005366" "D3piPlusR3_SRR1005367" "D3piMinusR1_SRR1005368" "D3piMinusR2_SRR1005369" "D3piMinusR3_SRR1005370" "D7piPlusR1_SRR1005377" "D7piPlusR2_SRR1005378" "D7piPlusR3_SRR1005379" "D7piMinusR1_SRR1005380" "D7piMinusR2_SRR1005381" "D7piMinusR3_SRR1005382" "D21piPlusR1_SRR1005353" "D21piPlusR2_SRR1005354" "D21piPlusR3_SRR1005355" "D21piMinusR1_SRR1005356" "D21piMinusR2_SRR1005357" "D21piMinusR3_SRR1005358" "D21H1piPlusR1_SRR1005326" "D21H1piPlusR2_SRR1005327" "D21H1piPlusR3_SRR1005328" "D21H1piMinusR1_SRR1005329" "D21H1piMinusR2_SRR1005330" "D21H1piMinusR3_SRR1005331" "D21H1piPlusRecR1_SRR1005332" "D21H1piPlusRecR2_SRR1005333" "D21H1piPlusRecR3_SRR1005334" "D21H6piPlusR1_SRR1005344" "D21H6piPlusR2_SRR1005345" "D21H6piPlusR3_SRR1005346" "D21H6piMinusR1_SRR1005347" "D21H6piMinusR2_SRR1005348" "D21H6piMinusR3_SRR1005349" "D21H6piPlusRecR1_SRR1005350" "D21H6piPlusRecR2_SRR1005351" "D21H6piPlusRecR3_SRR1005352" "D21H24piPlusR1_SRR1005335" "D21H24piPlusR2_SRR1005336" "D21H24piPlusR3_SRR1005337" "D21H24piMinusR1_SRR1005338" "D21H24piMinusR2_SRR1005339" "D21H24piMinusR3_SRR1005340" "D21H24piPlusRecR1_SRR1005341" "D21H24piPlusRecR2_SRR1005342" "D21H24piPlusRecR3_SRR1005343")
AllSamples=(${FeZnCuMnSampleids[@]} ${PiRootSamples[@]} ${PiShootSamples[@]})


mkdir DASG
cd DASG
mkdir ASE
cd ASE
for((i=0; i<=140; i++));
do
mkdir ${AllSamples[$i]}
cd ${AllSamples[$i]}
samfile="/home/project/"${AllSamples[$i]}"/accepted_hits75.bam"
perl /home/pcd/MATS_config.pl -dir /home/project/ -readlen 75 -sam $samfile -dataType paired > config.txt
cp "/home/project/junctions.per.sample.pickle" ./
cp "/home/project/totalJunctions.pickle" ./
python /home/soft/rMATS.3.2.5/bin/MATS.processsUnique.bam.py config.txt
cd ..
done

rm -f PiFeZnCuMn_samplelist_condition
echo -e "ControlR1\0011ControlR1" > PiFeZnCuMn_samplelist_condition
echo -e "ControlR2\0011ControlR2" >> PiFeZnCuMn_samplelist_condition
echo -e "ControlR3\0011\ControlR3" >> PiFeZnCuMn_samplelist_condition
echo -e "FeR1\0011FeR1" >> PiFeZnCuMn_samplelist_condition
echo -e "FeR2\0011FeR2" >> PiFeZnCuMn_samplelist_condition
echo -e "FeR3\0011FeR3" >> PiFeZnCuMn_samplelist_condition
echo -e "ZnR1\0011ZnR1" >> PiFeZnCuMn_samplelist_condition
echo -e "ZnR2\0011ZnR2" >> PiFeZnCuMn_samplelist_condition
echo -e "ZnR3\0011ZnR3" >> PiFeZnCuMn_samplelist_condition
echo -e "CuR1\0011CuR1" >> PiFeZnCuMn_samplelist_condition
echo -e "CuR2\0011CuR2" >> PiFeZnCuMn_samplelist_condition
echo -e "CuR3\0011CuR3" >> PiFeZnCuMn_samplelist_condition
echo -e "MnR1\0011MnR1" >> PiFeZnCuMn_samplelist_condition
echo -e "MnR2\0011MnR2" >> PiFeZnCuMn_samplelist_condition
echo -e "MnR3\0011MnR3" >> PiFeZnCuMn_samplelist_condition
echo -e "MnR1\0011MnR1" >> PiFeZnCuMn_samplelist_condition
echo -e "MnR2\0011MnR2" >> PiFeZnCuMn_samplelist_condition
echo -e "H1piPlusR1_SRR1005257\0011H1piPlusR1_R" >> PiFeZnCuMn_samplelist_condition
echo -e "H1piPlusR2_SRR1005258\0011H1piPlusR2_R" >> PiFeZnCuMn_samplelist_condition
echo -e "H1piPlusR3_SRR1005259\0011H1piPlusR3_R" >> PiFeZnCuMn_samplelist_condition
echo -e "H1piMinusR1_SRR1005260\0011H1piMinusR1_R" >> PiFeZnCuMn_samplelist_condition
echo -e "H1piMinusR2_SRR1005261\0011H1piMinusR2_R" >> PiFeZnCuMn_samplelist_condition
echo -e "H1piMinusR3_SRR1005262\0011H1piMinusR3_R" >> PiFeZnCuMn_samplelist_condition
echo -e "H6piPlusR1_SRR1005308\0011H6piPlusR1_R" >> PiFeZnCuMn_samplelist_condition
echo -e "H6piPlusR2_SRR1005309\0011H6piPlusR2_R" >> PiFeZnCuMn_samplelist_condition
echo -e "H6piPlusR3_SRR1005310\0011H6piPlusR3_R" >> PiFeZnCuMn_samplelist_condition
echo -e "H6piMinusR1_SRR1005311\0011H6piMinusR1_R" >> PiFeZnCuMn_samplelist_condition
echo -e "H6piMinusR2_SRR1005312\0011H6piMinusR2_R" >> PiFeZnCuMn_samplelist_condition
echo -e "H6piMinusR3_SRR1005313\0011H6piMinusR3_R" >> PiFeZnCuMn_samplelist_condition
echo -e "H24piPlusR1_SRR1005296\0011H24piPlusR1_R" >> PiFeZnCuMn_samplelist_condition
echo -e "H24piPlusR2_SRR1005297\0011H24piPlusR2_R" >> PiFeZnCuMn_samplelist_condition
echo -e "H24piPlusR3_SRR1005298\0011H24piPlusR3_R" >> PiFeZnCuMn_samplelist_condition
echo -e "H24piMinusR1_SRR1005299\0011H24piMinusR1_R" >> PiFeZnCuMn_samplelist_condition
echo -e "H24piMinusR2_SRR1005300\0011H24piMinusR2_R" >> PiFeZnCuMn_samplelist_condition
echo -e "H24piMinusR3_SRR1005301\0011H24piMinusR3_R" >> PiFeZnCuMn_samplelist_condition
echo -e "D3piPlusR1_SRR1005302\0011D3piPlusR1_R" >> PiFeZnCuMn_samplelist_condition
echo -e "D3piPlusR2_SRR1005303\0011D3piPlusR2_R" >> PiFeZnCuMn_samplelist_condition
echo -e "D3piPlusR3_SRR1005304\0011D3piPlusR3_R" >> PiFeZnCuMn_samplelist_condition
echo -e "D3piMinusR1_SRR1005305\0011D3piMinusR1_R" >> PiFeZnCuMn_samplelist_condition
echo -e "D3piMinusR2_SRR1005306\0011D3piMinusR2_R" >> PiFeZnCuMn_samplelist_condition
echo -e "D3piMinusR3_SRR1005307\0011D3piMinusR3_R" >> PiFeZnCuMn_samplelist_condition
echo -e "D7piPlusR1_SRR1005314\0011D7piPlusR1_R" >> PiFeZnCuMn_samplelist_condition
echo -e "D7piPlusR2_SRR1005315\0011D7piPlusR2_R" >> PiFeZnCuMn_samplelist_condition
echo -e "D7piPlusR3_SRR1005316\0011D7piPlusR3_R" >> PiFeZnCuMn_samplelist_condition
echo -e "D7piMinusR1_SRR1005317\0011D7piMinusR1_R" >> PiFeZnCuMn_samplelist_condition
echo -e "D7piMinusR2_SRR1005318\0011D7piMinusR2_R" >> PiFeZnCuMn_samplelist_condition
echo -e "D7piMinusR3_SRR1005319\0011D7piMinusR3_R" >> PiFeZnCuMn_samplelist_condition
echo -e "D21piPlusR1_SRR1005290\0011D21piPlusR1_R" >> PiFeZnCuMn_samplelist_condition
echo -e "D21piPlusR2_SRR1005291\0011D21piPlusR2_R" >> PiFeZnCuMn_samplelist_condition
echo -e "D21piPlusR3_SRR1005292\0011D21piPlusR3_R" >> PiFeZnCuMn_samplelist_condition
echo -e "D21piMinusR1_SRR1005293\0011D21piMinusR1_R" >> PiFeZnCuMn_samplelist_condition
echo -e "D21piMinusR2_SRR1005294\0011D21piMinusR2_R" >> PiFeZnCuMn_samplelist_condition
echo -e "D21piMinusR3_SRR1005295\0011D21piMinusR3_R" >> PiFeZnCuMn_samplelist_condition
echo -e "D21H1piPlusR1_SRR1005263\0011D21H1piPlusR1_R" >> PiFeZnCuMn_samplelist_condition
echo -e "D21H1piPlusR2_SRR1005264\0011D21H1piPlusR2_R" >> PiFeZnCuMn_samplelist_condition
echo -e "D21H1piPlusR3_SRR1005265\0011D21H1piPlusR3_R" >> PiFeZnCuMn_samplelist_condition
echo -e "D21H1piMinusR1_SRR1005266\0011D21H1piMinusR1_R" >> PiFeZnCuMn_samplelist_condition
echo -e "D21H1piMinusR2_SRR1005267\0011D21H1piMinusR2_R" >> PiFeZnCuMn_samplelist_condition
echo -e "D21H1piMinusR3_SRR1005268\0011D21H1piMinusR3_R" >> PiFeZnCuMn_samplelist_condition
echo -e "D21H1piPlusRecR1_SRR1005269\0011D21H1piPlusRecR1_R" >> PiFeZnCuMn_samplelist_condition
echo -e "D21H1piPlusRecR2_SRR1005270\0011D21H1piPlusRecR2_R" >> PiFeZnCuMn_samplelist_condition
echo -e "D21H1piPlusRecR3_SRR1005271\0011D21H1piPlusRecR3_R" >> PiFeZnCuMn_samplelist_condition
echo -e "D21H6piPlusR1_SRR1005281\0011D21H6piPlusR1_R" >> PiFeZnCuMn_samplelist_condition
echo -e "D21H6piPlusR2_SRR1005282\0011D21H6piPlusR2_R" >> PiFeZnCuMn_samplelist_condition
echo -e "D21H6piPlusR3_SRR1005283\0011D21H6piPlusR3_R" >> PiFeZnCuMn_samplelist_condition
echo -e "D21H6piMinusR1_SRR1005284\0011D21H6piMinusR1_R" >> PiFeZnCuMn_samplelist_condition
echo -e "D21H6piMinusR2_SRR1005285\0011D21H6piMinusR2_R" >> PiFeZnCuMn_samplelist_condition
echo -e "D21H6piMinusR3_SRR1005286\0011D21H6piMinusR3_R" >> PiFeZnCuMn_samplelist_condition
echo -e "D21H6piPlusRecR1_SRR1005287\0011D21H6piPlusRecR1_R" >> PiFeZnCuMn_samplelist_condition
echo -e "D21H6piPlusRecR2_SRR1005288\0011D21H6piPlusRecR2_R" >> PiFeZnCuMn_samplelist_condition
echo -e "D21H6piPlusRecR3_SRR1005289\0011D21H6piPlusRecR3_R" >> PiFeZnCuMn_samplelist_condition
echo -e "D21H24piPlusR1_SRR1005272\0011D21H24piPlusR1_R" >> PiFeZnCuMn_samplelist_condition
echo -e "D21H24piPlusR2_SRR1005273\0011D21H24piPlusR2_R" >> PiFeZnCuMn_samplelist_condition
echo -e "D21H24piPlusR3_SRR1005274\0011D21H24piPlusR3_R" >> PiFeZnCuMn_samplelist_condition
echo -e "D21H24piMinusR1_SRR1005275\0011D21H24piMinusR1_R" >> PiFeZnCuMn_samplelist_condition
echo -e "D21H24piMinusR2_SRR1005276\0011D21H24piMinusR2_R" >> PiFeZnCuMn_samplelist_condition
echo -e "D21H24piMinusR3_SRR1005277\0011D21H24piMinusR3_R" >> PiFeZnCuMn_samplelist_condition
echo -e "D21H24piPlusRecR1_SRR1005278\0011D21H24piPlusRecR1_R" >> PiFeZnCuMn_samplelist_condition
echo -e "D21H24piPlusRecR2_SRR1005279\0011D21H24piPlusRecR2_R" >> PiFeZnCuMn_samplelist_condition
echo -e "D21H24piPlusRecR3_SRR1005280\0011D21H24piPlusRecR3_R" >> PiFeZnCuMn_samplelist_condition
echo -e "H1piPlusR1_SRR1005320\0011H1piPlusR1_S" >> PiFeZnCuMn_samplelist_condition
echo -e "H1piPlusR2_SRR1005321\0011H1piPlusR2_S" >> PiFeZnCuMn_samplelist_condition
echo -e "H1piPlusR3_SRR1005322\0011H1piPlusR3_S" >> PiFeZnCuMn_samplelist_condition
echo -e "H1piMinusR1_SRR1005323\0011H1piMinusR1_S" >> PiFeZnCuMn_samplelist_condition
echo -e "H1piMinusR2_SRR1005324\0011H1piMinusR2_S" >> PiFeZnCuMn_samplelist_condition
echo -e "H1piMinusR3_SRR1005325\0011H1piMinusR3_S" >> PiFeZnCuMn_samplelist_condition
echo -e "H6piPlusR1_SRR1005371\0011H6piPlusR1_S" >> PiFeZnCuMn_samplelist_condition
echo -e "H6piPlusR2_SRR1005372\0011H6piPlusR2_S" >> PiFeZnCuMn_samplelist_condition
echo -e "H6piPlusR3_SRR1005373\0011H6piPlusR3_S" >> PiFeZnCuMn_samplelist_condition
echo -e "H6piMinusR1_SRR1005374\0011H6piMinusR1_S" >> PiFeZnCuMn_samplelist_condition
echo -e "H6piMinusR2_SRR1005375\0011H6piMinusR2_S" >> PiFeZnCuMn_samplelist_condition
echo -e "H6piMinusR3_SRR1005376\0011H6piMinusR3_S" >> PiFeZnCuMn_samplelist_condition
echo -e "H24piPlusR1_SRR1005359\0011H24piPlusR1_S" >> PiFeZnCuMn_samplelist_condition
echo -e "H24piPlusR2_SRR1005360\0011H24piPlusR2_S" >> PiFeZnCuMn_samplelist_condition
echo -e "H24piPlusR3_SRR1005361\0011H24piPlusR3_S" >> PiFeZnCuMn_samplelist_condition
echo -e "H24piMinusR1_SRR1005362\0011H24piMinusR1_S" >> PiFeZnCuMn_samplelist_condition
echo -e "H24piMinusR2_SRR1005363\0011H24piMinusR2_S" >> PiFeZnCuMn_samplelist_condition
echo -e "H24piMinusR3_SRR1005364\0011H24piMinusR3_S" >> PiFeZnCuMn_samplelist_condition
echo -e "D3piPlusR1_SRR1005365\0011D3piPlusR1_S" >> PiFeZnCuMn_samplelist_condition
echo -e "D3piPlusR2_SRR1005366\0011D3piPlusR2_S" >> PiFeZnCuMn_samplelist_condition
echo -e "D3piPlusR3_SRR1005367\0011D3piPlusR3_S" >> PiFeZnCuMn_samplelist_condition
echo -e "D3piMinusR1_SRR1005368\0011D3piMinusR1_S" >> PiFeZnCuMn_samplelist_condition
echo -e "D3piMinusR2_SRR1005369\0011D3piMinusR2_S" >> PiFeZnCuMn_samplelist_condition
echo -e "D3piMinusR3_SRR1005370\0011D3piMinusR3_S" >> PiFeZnCuMn_samplelist_condition
echo -e "D7piPlusR1_SRR1005377\0011D7piPlusR1_S" >> PiFeZnCuMn_samplelist_condition
echo -e "D7piPlusR2_SRR1005378\0011D7piPlusR2_S" >> PiFeZnCuMn_samplelist_condition
echo -e "D7piPlusR3_SRR1005379\0011D7piPlusR3_S" >> PiFeZnCuMn_samplelist_condition
echo -e "D7piMinusR1_SRR1005380\0011D7piMinusR1_S" >> PiFeZnCuMn_samplelist_condition
echo -e "D7piMinusR2_SRR1005381\0011D7piMinusR2_S" >> PiFeZnCuMn_samplelist_condition
echo -e "D7piMinusR3_SRR1005382\0011D7piMinusR3_S" >> PiFeZnCuMn_samplelist_condition
echo -e "D21piPlusR1_SRR1005353\0011D21piPlusR1_S" >> PiFeZnCuMn_samplelist_condition
echo -e "D21piPlusR2_SRR1005354\0011D21piPlusR2_S" >> PiFeZnCuMn_samplelist_condition
echo -e "D21piPlusR3_SRR1005355\0011D21piPlusR3_S" >> PiFeZnCuMn_samplelist_condition
echo -e "D21piMinusR1_SRR1005356\0011D21piMinusR1_S" >> PiFeZnCuMn_samplelist_condition
echo -e "D21piMinusR2_SRR1005357\0011D21piMinusR2_S" >> PiFeZnCuMn_samplelist_condition
echo -e "D21piMinusR3_SRR1005358\0011D21piMinusR3_S" >> PiFeZnCuMn_samplelist_condition
echo -e "D21H1piPlusR1_SRR1005326\0011D21H1piPlusR1_S" >> PiFeZnCuMn_samplelist_condition
echo -e "D21H1piPlusR2_SRR1005327\0011D21H1piPlusR2_S" >> PiFeZnCuMn_samplelist_condition
echo -e "D21H1piPlusR3_SRR1005328\0011D21H1piPlusR3_S" >> PiFeZnCuMn_samplelist_condition
echo -e "D21H1piMinusR1_SRR1005329\0011D21H1piMinusR1_S" >> PiFeZnCuMn_samplelist_condition
echo -e "D21H1piMinusR2_SRR1005330\0011D21H1piMinusR2_S" >> PiFeZnCuMn_samplelist_condition
echo -e "D21H1piMinusR3_SRR1005331\0011D21H1piMinusR3_S" >> PiFeZnCuMn_samplelist_condition
echo -e "D21H1piPlusRecR1_SRR1005332\0011D21H1piPlusRecR1_S" >> PiFeZnCuMn_samplelist_condition
echo -e "D21H1piPlusRecR2_SRR1005333\0011D21H1piPlusRecR2_S" >> PiFeZnCuMn_samplelist_condition
echo -e "D21H1piPlusRecR3_SRR1005334\0011D21H1piPlusRecR3_S" >> PiFeZnCuMn_samplelist_condition
echo -e "D21H6piPlusR1_SRR1005344\0011D21H6piPlusR1_S" >> PiFeZnCuMn_samplelist_condition
echo -e "D21H6piPlusR2_SRR1005345\0011D21H6piPlusR2_S" >> PiFeZnCuMn_samplelist_condition
echo -e "D21H6piPlusR3_SRR1005346\0011D21H6piPlusR3_S" >> PiFeZnCuMn_samplelist_condition
echo -e "D21H6piMinusR1_SRR1005347\0011D21H6piMinusR1_S" >> PiFeZnCuMn_samplelist_condition
echo -e "D21H6piMinusR2_SRR1005348\0011D21H6piMinusR2_S" >> PiFeZnCuMn_samplelist_condition
echo -e "D21H6piMinusR3_SRR1005349\0011D21H6piMinusR3_S" >> PiFeZnCuMn_samplelist_condition
echo -e "D21H6piPlusRecR1_SRR1005350\0011D21H6piPlusRecR1_S" >> PiFeZnCuMn_samplelist_condition
echo -e "D21H6piPlusRecR2_SRR1005351\0011D21H6piPlusRecR2_S" >> PiFeZnCuMn_samplelist_condition
echo -e "D21H6piPlusRecR3_SRR1005352\0011D21H6piPlusRecR3_S" >> PiFeZnCuMn_samplelist_condition
echo -e "D21H24piPlusR1_SRR1005335\0011D21H24piPlusR1_S" >> PiFeZnCuMn_samplelist_condition
echo -e "D21H24piPlusR2_SRR1005336\0011D21H24piPlusR2_S" >> PiFeZnCuMn_samplelist_condition
echo -e "D21H24piPlusR3_SRR1005337\0011D21H24piPlusR3_S" >> PiFeZnCuMn_samplelist_condition
echo -e "D21H24piMinusR1_SRR1005338\0011D21H24piMinusR1_S" >> PiFeZnCuMn_samplelist_condition
echo -e "D21H24piMinusR2_SRR1005339\0011D21H24piMinusR2_S" >> PiFeZnCuMn_samplelist_condition
echo -e "D21H24piMinusR3_SRR1005340\0011D21H24piMinusR3_S" >> PiFeZnCuMn_samplelist_condition
echo -e "D21H24piPlusRecR1_SRR1005341\0011D21H24piPlusRecR1_S" >> PiFeZnCuMn_samplelist_condition
echo -e "D21H24piPlusRecR2_SRR1005342\0011D21H24piPlusRecR2_S" >> PiFeZnCuMn_samplelist_condition
echo -e "D21H24piPlusRecR3_SRR1005343\0011D21H24piPlusRecR3_S" >> PiFeZnCuMn_samplelist_condition

#IJC and SJC quantification model (only spanned splicing junctions)
perl /home/pcd/inclusive_incl_excl_all_sample.pl -ASEdir /home/project/ -samplelist PiFeZnCuMn_samplelist_condition > PiFeZnCuMn_in_ex.xls





#IC and SC quantification model (including spanned the splicing junction or occurred in the alternative exon or intron target)
perl /home/pcd/inclusive_incl_excl_all_sample_jcec.pl -ASEdir /home/project/ -samplelist PiFeZnCuMn_samplelist_condition > PiFeZnCuMn_in_ex_JCEC.xls

perl /home/pcd/inclusive_ratio_all_sample.pl -ASEdir /home/project/ -samplelist PiFeZnCuMn_samplelist_condition > PiFeZnCuMn_psi.xls
Rscript /home/pcd/transform_PSI.r
perl /home/pcd/ASE_call_all_sample.pl -ASEdir /home/project/ -samplelist PiFeZnCuMn_samplelist_condition -cutoff 30 > PiFeZnCuMn_call.xls
perl /home/pcd/filter_ASE.pl -inputfile PiFeZnCuMn_in_ex.xls -totalcutoff 30 > PiFeZnCuMn_in_ex_filter.xls


# Detected ASEs, in total or by each condition(-Fe, -Zn, -Cu, -Mn, -Pi)  
perl /home/pcd/detected_ASE.pl -inputfile PiFeZnCuMn_in_ex.xls -totalcutoff 0 > AS_event_list_sample.txt
perl /home/pcd/detected_ASE_all.pl -inputfile PiFeZnCuMn_in_ex.xls -totalcutoff 0 > AS_event_list_all.txt
cut -f 1,2,3,4 PiFeZnCuMn_in_ex.xls > one.xls
perl /home/pcd/detected_ASE_all.pl -inputfile one.xls -totalcutoff 0 > AS_event_list_Control.txt
cut -f 1,5,6,7 PiFeZnCuMn_in_ex.xls > one.xls
perl /home/pcd/detected_ASE_all.pl -inputfile one.xls -totalcutoff 0 > AS_event_list_Fe.txt
cut -f 1,8,9,10 PiFeZnCuMn_in_ex.xls > one.xls
perl /home/pcd/detected_ASE_all.pl -inputfile one.xls -totalcutoff 0 > AS_event_list_Zn.txt
cut -f 1,11,12,13 PiFeZnCuMn_in_ex.xls > one.xls
perl /home/pcd/detected_ASE_all.pl -inputfile one.xls -totalcutoff 0 > AS_event_list_Cu.txt
cut -f 1,14,15,16 PiFeZnCuMn_in_ex.xls > one.xls
perl /home/pcd/detected_ASE_all.pl -inputfile one.xls -totalcutoff 0 > AS_event_list_Mn.txt
cut -f 1,5,6,7,8,9,10,11,12,13,14,15,16 PiFeZnCuMn_in_ex.xls > one.xls
perl /home/pcd/detected_ASE_all.pl -inputfile one.xls -totalcutoff 0 > AS_event_list_FeZnCuMn.txt
cut -f 1,17,18,19,23,24,25,29,30,31,35,36,37,41,42,43,47,48,49,53,54,55,62,63,64,71,72,73 PiFeZnCuMn_in_ex.xls > one.xls
perl /home/pcd/detected_ASE_all.pl -inputfile one.xls -totalcutoff 0 > AS_event_list_RootPiPlus.txt
cut -f 1,20,21,22,26,27,28,32,33,34,38,39,40,44,45,46,50,51,52,56,57,58,65,66,67,74,75,76 PiFeZnCuMn_in_ex.xls > one.xls
perl /home/pcd/detected_ASE_all.pl -inputfile one.xls -totalcutoff 0 > AS_event_list_RootPiMinus.txt
cut -f 1,80,81,82,86,87,88,92,93,94,98,99,100,104,105,106,110,111,112,116,117,118,125,126,127,134,135,136 PiFeZnCuMn_in_ex.xls > one.xls
perl /home/pcd/detected_ASE_all.pl -inputfile one.xls -totalcutoff 0 > AS_event_list_ShootPiPlus.txt
cut -f 1,83,84,85,89,90,91,95,96,97,101,102,103,107,108,109,113,114,115,119,120,121,128,129,130,137,138,139 PiFeZnCuMn_in_ex.xls > one.xls
perl /home/pcd/detected_ASE_all.pl -inputfile one.xls -totalcutoff 0 > AS_event_list_ShootPiMinus.txt
cut -f 1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,38,39,40,41,42,43,44,45,46,47,48,49,50,51,52,53,54,55,56,57,58,59,60,61,62,63,64,65,66,67,68,69,70,71,72,73,74,75,76,77,78,79 PiFeZnCuMn_in_ex.xls > one.xls
perl /home/pcd/detected_ASE_all.pl -inputfile one.xls -totalcutoff 0 > AS_event_list_Root.txt
cut -f 1,80,81,82,83,84,85,86,87,88,89,90,91,92,93,94,95,96,97,98,99,100,101,102,103,104,105,106,107,108,109,110,111,112,113,114,115,116,117,118,119,120,121,122,123,124,125,126,127,128,129,130,131,132,133,134,135,136,137,138,139,140,141,142 PiFeZnCuMn_in_ex.xls > one.xls
perl /home/pcd/detected_ASE_all.pl -inputfile one.xls -totalcutoff 0 > AS_event_list_Shoot.txt


perl /home/pcd/AS_event_stat.pl -ASE AS_event_list_all.txt > PiFeZnCuMn_ASE_type.xls

Rscript /home/pcd/Figure1A_ASE_structure_and_ratio.r


#Splice Site type distribution of ASE (known or novel)
perl /home/pcd/get_splicesite_SE.pl -ASE PiFeZnCuMn_in_ex.xls -genome /home/database/genome.fa > splice_site_SE.txt
perl /home/pcd/get_splicesite_RI.pl -ASE PiFeZnCuMn_in_ex.xls -genome /home/database/genome.fa > splice_site_RI.txt
perl /home/pcd/get_splicesite_A5SS.pl -ASE PiFeZnCuMn_in_ex.xls -genome /home/database/genome.fa > splice_site_A5SS.txt
perl /home/pcd/get_splicesite_A3SS.pl -ASE PiFeZnCuMn_in_ex.xls -genome /home/database/genome.fa > splice_site_A3SS.txt
perl /home/pcd/get_splicesite_MXE.pl -ASE PiFeZnCuMn_in_ex.xls -genome /home/database/genome.fa > splice_site_MXE.txt
perl /home/pcd/event_is_known_novel.pl -ASE PiFeZnCuMn_in_ex.xls
cat RI_known.txt A3SS_known.txt SE_known.txt A5SS_known.txt MXE_known.txt > ASE_known.txt
cat RI_novel.txt A3SS_novel.txt SE_novel.txt A5SS_novel.txt MXE_novel.txt > ASE_novel.txt
n1=`awk 'END{print NR}' ASE_known.txt`
n2=`awk 'END{print NR}' ASE_novel.txt`
echo $n1
echo $n2
let n3=$n2*1000/$n1
echo $n3
perl /home/pcd/add_ASE_known_novel.pl -AStype splice_site_RI.txt > splice_site_RI_known_novel.txt
perl /home/pcd/add_ASE_known_novel.pl -AStype splice_site_A5SS.txt > splice_site_A5SS_known_novel.txt
perl /home/pcd/add_ASE_known_novel.pl -AStype splice_site_A3SS.txt > splice_site_A3SS_known_novel.txt
perl /home/pcd/add_ASE_known_novel.pl -AStype splice_site_SE.txt > splice_site_SE_known_novel.txt
perl /home/pcd/add_ASE_known_novel.pl -AStype splice_site_MXE.txt > splice_site_MXE_known_novel.txt
perl /home/pcd/splice_site_stat.pl > splice_site_stat.xls

Rscript /home/pcd/SupplementalFigure2C_SpliceSite_proportion_known_novel.r

Rscript /home/pcd/SupplementalFigure3_ASG_stat.r /home/project/intron_list.txt

perl /home/pcd/AS_gene_stat_by_all_v3.pl -intron /home/project/intron_list.txt -inputfile PiFeZnCuMn_in_ex.xls -readcutoff 0 > gene2ASE_num.txt

Rscript /home/pcd/hdl_AS_gene_all_cutoff.r "gene2ASE_num.txt" /home/project/intron_list.txt

# The correlation between ASE frequency, GC percentage, expression level (TPM), exon number, intron number, exon length, intron length
perl /home/pcd/add_GC_exprs_exonnum.pl -genefpkm /home/project/DEG/FeZnCuMnPiRootShoot/PiFeZnCuMn_gene_sample_TPM_filter.xls -genefeature /home/project/gene_feature.xls -ASEnum gene2ASE_num.txt > gene2ASE_num_gene_structure.xls
perl /home/pcd/add_expression_ASE.pl -genefpkm /home/project/DEG/FeZnCuMnPiRootShoot/PiFeZnCuMn_gene_sample_TPM_filter.xls -diffgene high_AS_genes.txt > high_AS_genes_exprs.xls
perl /home/pcd/add_expression_ASE.pl -genefpkm /home/project/DEG/FeZnCuMnPiRootShoot/PiFeZnCuMn_gene_sample_TPM_filter.xls -diffgene medium_AS_genes.txt > medium_AS_genes_exprs.xls
perl /home/pcd/add_expression_ASE.pl -genefpkm /home/project/DEG/FeZnCuMnPiRootShoot/PiFeZnCuMn_gene_sample_TPM_filter.xls -diffgene low_AS_genes.txt > low_AS_genes_exprs.xls
perl /home/pcd/add_expression_ASE.pl -genefpkm /home/project/DEG/FeZnCuMnPiRootShoot/PiFeZnCuMn_gene_sample_TPM_filter.xls -diffgene no_AS_genes.txt > no_AS_genes_exprs.xls
perl /home/pcd/add_seq_feature_ASE.pl -genefeature /home/project/gene_feature.xls -diffgene high_AS_genes.txt > high_AS_genes_seq_feature.xls
perl /home/pcd/add_seq_feature_ASE.pl -genefeature /home/project/gene_feature.xls -diffgene medium_AS_genes.txt > medium_AS_genes_seq_feature.xls
perl /home/pcd/add_seq_feature_ASE.pl -genefeature /home/project/gene_feature.xls -diffgene low_AS_genes.txt > low_AS_genes_seq_feature.xls
perl /home/pcd/add_seq_feature_ASE.pl -genefeature /home/project/gene_feature.xls -diffgene no_AS_genes.txt > no_AS_genes_seq_feature.xls
perl /home/pcd/add_GC_ASE.pl -genefeature /home/project/gene_feature.xls -diffgene high_AS_genes.txt > high_AS_genes_GC.xls
perl /home/pcd/add_GC_ASE.pl -genefeature /home/project/gene_feature.xls -diffgene medium_AS_genes.txt > medium_AS_genes_GC.xls
perl /home/pcd/add_GC_ASE.pl -genefeature /home/project/gene_feature.xls -diffgene low_AS_genes.txt > low_AS_genes_GC.xls
perl /home/pcd/add_GC_ASE.pl -genefeature /home/project/gene_feature.xls -diffgene no_AS_genes.txt > no_AS_genes_GC.xls

 
Rscript /home/pcd/SupplementalFigure2D_ASE_feature_boxplot.r
 


# Domain enrichment analysis of high AS frequency genes
Rscript /home/pcd/hdl_AS_gene_all_cutoff.r "gene2ASE_num.txt" /home/project/intron_list.txt
perl /home/pcd/degene2known.pl -annotation /home/project/merged_filter_multiexon_junction40_TPM_annotation.xls -degene high_AS_genes.txt  > high_AS_genes_known.txt
perl /home/pcd/domain_enrich.pl -ispvalue "Pvalue" -bl /home/database/allgene.txt -gl high_AS_genes_known.txt -symbol2domain /home/database/symbol2interpro_pfam.txt
mv domain_enrichment_significant.xls high_AS_genes_domain_enrichment_significant.xls
mv domain_enrichment.xls high_AS_genes_domain_enrichment.xls
Rscript /home/pcd/Figure5A_high_ASG_domain_barplot_high.r


#Two representative inclusive and exclusive transcripts of ASEs
perl /home/pcd/find_transcript_by_SE_event.pl -gtf /home/database/merged_filter_multiexon_junction40_TPM.gtf -SE PiFeZnCuMn_in_ex_filter.xls > SE2transcript.xls
perl /home/pcd/find_transcript_by_RI_event.pl -gtf /home/database/merged_filter_multiexon_junction40_TPM.gtf -RI PiFeZnCuMn_in_ex_filter.xls > RI2transcript.xls
perl /home/pcd/find_transcript_by_A3SS_event.pl -gtf /home/database/merged_filter_multiexon_junction40_TPM.gtf -A3SS PiFeZnCuMn_in_ex_filter.xls > A3SS2transcript.xls
perl /home/pcd/find_transcript_by_A5SS_event.pl -gtf /home/database/merged_filter_multiexon_junction40_TPM.gtf -A5SS PiFeZnCuMn_in_ex_filter.xls > A5SS2transcript.xls
perl /home/pcd/find_transcript_by_MXE_event.pl -gtf /home/database/merged_filter_multiexon_junction40_TPM.gtf -MXE PiFeZnCuMn_in_ex_filter.xls > MXE2transcript.xls
perl /home/pcd/add_ASE_tid_v4.pl -FPKM /home/project/DEG/FeZnCuMnPiRootShoot/PiFeZnCuMn_gene_sample_TPM_filter.xls

