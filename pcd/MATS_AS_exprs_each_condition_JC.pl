use strict;
use Getopt::Long;
use List::Uniq ':all';
use vars qw($gtffile $cutoff $pvaluecutoff $mode $diffcutoff $ASEexprsfile $samplelist);
Getopt::Long::GetOptions(
    'mode=s' => \$mode,
	'samplelist=s' => \$samplelist,
	'cutoff=s' => \$cutoff,
);

my @tmp = split(/;/, $samplelist);
my $comname = $tmp[0];
my @controllist = split(/\,/, $tmp[1]);
my @caselist = split(/\,/, $tmp[2]);
my @controllist_in_ex;
for(my $i = 0; $i <= $#controllist; $i++)
{
    push(@controllist_in_ex, $controllist[$i]."_inclusion_skipping_count");
}
my @controllist_psi;
for(my $i = 0; $i <= $#controllist; $i++)
{
    push(@controllist_psi, $controllist[$i]."_PSI");
}
my @caselist_in_ex;
for(my $i = 0; $i <= $#caselist; $i++)
{
    push(@caselist_in_ex, $caselist[$i]."_inclusion_skipping_count");
}
my @caselist_psi;
for(my $i = 0; $i <= $#caselist; $i++)
{
    push(@caselist_psi, $caselist[$i]."_PSI");
}



my %exon_skip;
open(IN, "./SEout_".$mode."/SE.xls");
my $tt = <IN>;
$tt=~s/\n$//;
my @tmp = split(/\t/, $tt);
my $n = scalar(@tmp);
while(<IN>)
{
    chomp;
	s/\s+$//;
	my @arr = split(/\t/, $_);
	#print $_."\n";
	my $key = "SE&".$arr[0]."&".$arr[2]."&".$arr[3]."&".$arr[4]."&".$arr[5]."&".$arr[6]."&".$arr[7]."&".$arr[8]."&".$arr[9];
	my @tmp = split(/\,/, $arr[($n-11)]);
	my $s2_in = $tmp[0];
	my $s3_in = $tmp[1];	
	my $s4_in = $tmp[2];
	my @tmp = split(/\,/, $arr[($n-10)]);
	my $s2_ex = $tmp[0];
	my $s3_ex = $tmp[1];
	my $s4_ex = $tmp[2];	
	my @tmp = split(/\,/, $arr[($n-9)]);
	my $s8_in = $tmp[0];
	my $s9_in = $tmp[1];	
	my $s10_in = $tmp[2];
	my @tmp = split(/\,/, $arr[($n-8)]);
	my $s8_ex = $tmp[0];
	my $s9_ex = $tmp[1];	
	my $s10_ex = $tmp[2];	
	my @tmp = split(/\,/, $arr[($n-3)]);
	my $s2_ratio = $tmp[0];
	my $s3_ratio = $tmp[1];	
	my $s4_ratio = $tmp[2];
	my @tmp = split(/\,/, $arr[($n-2)]);
	my $s8_ratio = $tmp[0];
	my $s9_ratio = $tmp[1];	
	my $s10_ratio = $tmp[2];	
	#$exon_skip{$key}->{"C1"}=$s2_in."\t".$s2_ex."\t".($s2_in+$s2_ex)."\t".$s2_ratio;
	$exon_skip{$key}->{"C1"}->{"IN"}=$s2_in;
	$exon_skip{$key}->{"C1"}->{"EX"}=$s2_ex;
	$exon_skip{$key}->{"C1"}->{"TT"}=$s2_in+$s2_ex;
	$exon_skip{$key}->{"C1"}->{"PSI"}=$s2_ratio;
	$exon_skip{$key}->{"C2"}->{"IN"}=$s3_in;
	$exon_skip{$key}->{"C2"}->{"EX"}=$s3_ex;
	$exon_skip{$key}->{"C2"}->{"TT"}=$s3_in+$s3_ex;
	$exon_skip{$key}->{"C2"}->{"PSI"}=$s3_ratio;
	$exon_skip{$key}->{"C3"}->{"IN"}=$s4_in;
	$exon_skip{$key}->{"C3"}->{"EX"}=$s4_ex;
	$exon_skip{$key}->{"C3"}->{"TT"}=$s4_in+$s4_ex;
	$exon_skip{$key}->{"C3"}->{"PSI"}=$s4_ratio;
	$exon_skip{$key}->{"T1"}->{"IN"}=$s8_in;
	$exon_skip{$key}->{"T1"}->{"EX"}=$s8_ex;
	$exon_skip{$key}->{"T1"}->{"TT"}=$s8_in+$s8_ex;
	$exon_skip{$key}->{"T1"}->{"PSI"}=$s8_ratio;
	$exon_skip{$key}->{"T2"}->{"IN"}=$s9_in;
	$exon_skip{$key}->{"T2"}->{"EX"}=$s9_ex;
	$exon_skip{$key}->{"T2"}->{"TT"}=$s9_in+$s9_ex;
	$exon_skip{$key}->{"T2"}->{"PSI"}=$s9_ratio;
	$exon_skip{$key}->{"T3"}->{"IN"}=$s10_in;
	$exon_skip{$key}->{"T3"}->{"EX"}=$s10_ex;
	$exon_skip{$key}->{"T3"}->{"TT"}=$s10_in+$s10_ex;
	$exon_skip{$key}->{"T3"}->{"PSI"}=$s10_ratio;	
	$exon_skip{$key}->{"T_vs_C_diff"}=(-1)*$arr[($n-1)];
	$exon_skip{$key}->{"T_vs_C_Pvalue"}=$arr[($n-5)];
	$exon_skip{$key}->{"T_vs_C_FDR"}=$arr[($n-4)];	
}
close(IN);



my $file = "./RIout_".$mode."/RI.xls";
open(IN, $file);
my $tt = <IN>;
$tt=~s/\n$//;
my @tmp = split(/\t/, $tt);
my $n = scalar(@tmp);
#print $tmp[($n-5)]."\n";
while(<IN>)
{
    chomp;
	s/\s+$//;
	my @arr = split(/\t/, $_);
	my $key = "RI&".$arr[0]."&".$arr[2]."&".$arr[3]."&".$arr[4]."&".$arr[5]."&".$arr[6]."&".$arr[7]."&".$arr[8]."&".$arr[9];
	my @tmp = split(/\,/, $arr[($n-11)]);
	my $s2_in = $tmp[0];
	my $s3_in = $tmp[1];	
	my $s4_in = $tmp[2];
	my @tmp = split(/\,/, $arr[($n-10)]);
	my $s2_ex = $tmp[0];
	my $s3_ex = $tmp[1];
	my $s4_ex = $tmp[2];	
	my @tmp = split(/\,/, $arr[($n-9)]);
	my $s8_in = $tmp[0];
	my $s9_in = $tmp[1];	
	my $s10_in = $tmp[2];
	my @tmp = split(/\,/, $arr[($n-8)]);
	my $s8_ex = $tmp[0];
	my $s9_ex = $tmp[1];	
	my $s10_ex = $tmp[2];	
	my @tmp = split(/\,/, $arr[($n-3)]);
	my $s2_ratio = $tmp[0];
	my $s3_ratio = $tmp[1];	
	my $s4_ratio = $tmp[2];
	my @tmp = split(/\,/, $arr[($n-2)]);
	my $s8_ratio = $tmp[0];
	my $s9_ratio = $tmp[1];	
	my $s10_ratio = $tmp[2];	
	$exon_skip{$key}->{"C1"}->{"IN"}=$s2_in;
	$exon_skip{$key}->{"C1"}->{"EX"}=$s2_ex;
	$exon_skip{$key}->{"C1"}->{"TT"}=$s2_in+$s2_ex;
	$exon_skip{$key}->{"C1"}->{"PSI"}=$s2_ratio;
	$exon_skip{$key}->{"C2"}->{"IN"}=$s3_in;
	$exon_skip{$key}->{"C2"}->{"EX"}=$s3_ex;
	$exon_skip{$key}->{"C2"}->{"TT"}=$s3_in+$s3_ex;
	$exon_skip{$key}->{"C2"}->{"PSI"}=$s3_ratio;
	$exon_skip{$key}->{"C3"}->{"IN"}=$s4_in;
	$exon_skip{$key}->{"C3"}->{"EX"}=$s4_ex;
	$exon_skip{$key}->{"C3"}->{"TT"}=$s4_in+$s4_ex;
	$exon_skip{$key}->{"C3"}->{"PSI"}=$s4_ratio;
	$exon_skip{$key}->{"T1"}->{"IN"}=$s8_in;
	$exon_skip{$key}->{"T1"}->{"EX"}=$s8_ex;
	$exon_skip{$key}->{"T1"}->{"TT"}=$s8_in+$s8_ex;
	$exon_skip{$key}->{"T1"}->{"PSI"}=$s8_ratio;
	$exon_skip{$key}->{"T2"}->{"IN"}=$s9_in;
	$exon_skip{$key}->{"T2"}->{"EX"}=$s9_ex;
	$exon_skip{$key}->{"T2"}->{"TT"}=$s9_in+$s9_ex;
	$exon_skip{$key}->{"T2"}->{"PSI"}=$s9_ratio;
	$exon_skip{$key}->{"T3"}->{"IN"}=$s10_in;
	$exon_skip{$key}->{"T3"}->{"EX"}=$s10_ex;
	$exon_skip{$key}->{"T3"}->{"TT"}=$s10_in+$s10_ex;
	$exon_skip{$key}->{"T3"}->{"PSI"}=$s10_ratio;
	$exon_skip{$key}->{"T_vs_C_diff"}=(-1)*$arr[($n-1)];
	$exon_skip{$key}->{"T_vs_C_Pvalue"}=$arr[($n-5)];	
	$exon_skip{$key}->{"T_vs_C_FDR"}=$arr[($n-4)];	
}
close(IN);



my $file = "./A5SSout_".$mode."/A5SS.xls";
open(IN, $file);
my $tt = <IN>;
$tt=~s/\n$//;
my @tmp = split(/\t/, $tt);
my $n = scalar(@tmp);
#print $tmp[($n-5)]."\n";
while(<IN>)
{
    chomp;
	s/\s+$//;
	my @arr = split(/\t/, $_);
	my $key = "A5SS&".$arr[0]."&".$arr[2]."&".$arr[3]."&".$arr[4]."&".$arr[5]."&".$arr[6]."&".$arr[7]."&".$arr[8]."&".$arr[9];
	my @tmp = split(/\,/, $arr[($n-11)]);
	my $s2_in = $tmp[0];
	my $s3_in = $tmp[1];	
	my $s4_in = $tmp[2];
	my @tmp = split(/\,/, $arr[($n-10)]);
	my $s2_ex = $tmp[0];
	my $s3_ex = $tmp[1];
	my $s4_ex = $tmp[2];	
	my @tmp = split(/\,/, $arr[($n-9)]);
	my $s8_in = $tmp[0];
	my $s9_in = $tmp[1];	
	my $s10_in = $tmp[2];
	my @tmp = split(/\,/, $arr[($n-8)]);
	my $s8_ex = $tmp[0];
	my $s9_ex = $tmp[1];	
	my $s10_ex = $tmp[2];	
	my @tmp = split(/\,/, $arr[($n-3)]);
	my $s2_ratio = $tmp[0];
	my $s3_ratio = $tmp[1];	
	my $s4_ratio = $tmp[2];
	my @tmp = split(/\,/, $arr[($n-2)]);
	my $s8_ratio = $tmp[0];
	my $s9_ratio = $tmp[1];	
	my $s10_ratio = $tmp[2];	
	$exon_skip{$key}->{"C1"}->{"IN"}=$s2_in;
	$exon_skip{$key}->{"C1"}->{"EX"}=$s2_ex;
	$exon_skip{$key}->{"C1"}->{"TT"}=$s2_in+$s2_ex;
	$exon_skip{$key}->{"C1"}->{"PSI"}=$s2_ratio;
	$exon_skip{$key}->{"C2"}->{"IN"}=$s3_in;
	$exon_skip{$key}->{"C2"}->{"EX"}=$s3_ex;
	$exon_skip{$key}->{"C2"}->{"TT"}=$s3_in+$s3_ex;
	$exon_skip{$key}->{"C2"}->{"PSI"}=$s3_ratio;
	$exon_skip{$key}->{"C3"}->{"IN"}=$s4_in;
	$exon_skip{$key}->{"C3"}->{"EX"}=$s4_ex;
	$exon_skip{$key}->{"C3"}->{"TT"}=$s4_in+$s4_ex;
	$exon_skip{$key}->{"C3"}->{"PSI"}=$s4_ratio;
	$exon_skip{$key}->{"T1"}->{"IN"}=$s8_in;
	$exon_skip{$key}->{"T1"}->{"EX"}=$s8_ex;
	$exon_skip{$key}->{"T1"}->{"TT"}=$s8_in+$s8_ex;
	$exon_skip{$key}->{"T1"}->{"PSI"}=$s8_ratio;
	$exon_skip{$key}->{"T2"}->{"IN"}=$s9_in;
	$exon_skip{$key}->{"T2"}->{"EX"}=$s9_ex;
	$exon_skip{$key}->{"T2"}->{"TT"}=$s9_in+$s9_ex;
	$exon_skip{$key}->{"T2"}->{"PSI"}=$s9_ratio;
	$exon_skip{$key}->{"T3"}->{"IN"}=$s10_in;
	$exon_skip{$key}->{"T3"}->{"EX"}=$s10_ex;
	$exon_skip{$key}->{"T3"}->{"TT"}=$s10_in+$s10_ex;
	$exon_skip{$key}->{"T3"}->{"PSI"}=$s10_ratio;	
	$exon_skip{$key}->{"T_vs_C_diff"}=(-1)*$arr[($n-1)];
	$exon_skip{$key}->{"T_vs_C_Pvalue"}=$arr[($n-5)];
	$exon_skip{$key}->{"T_vs_C_FDR"}=$arr[($n-4)];		
;
}
close(IN);



my $file = "./A3SSout_".$mode."/A3SS.xls";
open(IN, $file);
my $tt = <IN>;
$tt=~s/\n$//;
my @tmp = split(/\t/, $tt);
my $n = scalar(@tmp);
#print $tmp[($n-5)]."\n";
while(<IN>)
{
    chomp;
	s/\s+$//;
	my @arr = split(/\t/, $_);
	my $key = "A3SS&".$arr[0]."&".$arr[2]."&".$arr[3]."&".$arr[4]."&".$arr[5]."&".$arr[6]."&".$arr[7]."&".$arr[8]."&".$arr[9];
	my @tmp = split(/\,/, $arr[($n-11)]);
	my $s2_in = $tmp[0];
	my $s3_in = $tmp[1];	
	my $s4_in = $tmp[2];
	my @tmp = split(/\,/, $arr[($n-10)]);
	my $s2_ex = $tmp[0];
	my $s3_ex = $tmp[1];
	my $s4_ex = $tmp[2];	
	my @tmp = split(/\,/, $arr[($n-9)]);
	my $s8_in = $tmp[0];
	my $s9_in = $tmp[1];	
	my $s10_in = $tmp[2];
	my @tmp = split(/\,/, $arr[($n-8)]);
	my $s8_ex = $tmp[0];
	my $s9_ex = $tmp[1];	
	my $s10_ex = $tmp[2];	
	my @tmp = split(/\,/, $arr[($n-3)]);
	my $s2_ratio = $tmp[0];
	my $s3_ratio = $tmp[1];	
	my $s4_ratio = $tmp[2];
	my @tmp = split(/\,/, $arr[($n-2)]);
	my $s8_ratio = $tmp[0];
	my $s9_ratio = $tmp[1];	
	my $s10_ratio = $tmp[2];	
	$exon_skip{$key}->{"C1"}->{"IN"}=$s2_in;
	$exon_skip{$key}->{"C1"}->{"EX"}=$s2_ex;
	$exon_skip{$key}->{"C1"}->{"TT"}=$s2_in+$s2_ex;
	$exon_skip{$key}->{"C1"}->{"PSI"}=$s2_ratio;
	$exon_skip{$key}->{"C2"}->{"IN"}=$s3_in;
	$exon_skip{$key}->{"C2"}->{"EX"}=$s3_ex;
	$exon_skip{$key}->{"C2"}->{"TT"}=$s3_in+$s3_ex;
	$exon_skip{$key}->{"C2"}->{"PSI"}=$s3_ratio;
	$exon_skip{$key}->{"C3"}->{"IN"}=$s4_in;
	$exon_skip{$key}->{"C3"}->{"EX"}=$s4_ex;
	$exon_skip{$key}->{"C3"}->{"TT"}=$s4_in+$s4_ex;
	$exon_skip{$key}->{"C3"}->{"PSI"}=$s4_ratio;
	$exon_skip{$key}->{"T1"}->{"IN"}=$s8_in;
	$exon_skip{$key}->{"T1"}->{"EX"}=$s8_ex;
	$exon_skip{$key}->{"T1"}->{"TT"}=$s8_in+$s8_ex;
	$exon_skip{$key}->{"T1"}->{"PSI"}=$s8_ratio;
	$exon_skip{$key}->{"T2"}->{"IN"}=$s9_in;
	$exon_skip{$key}->{"T2"}->{"EX"}=$s9_ex;
	$exon_skip{$key}->{"T2"}->{"TT"}=$s9_in+$s9_ex;
	$exon_skip{$key}->{"T2"}->{"PSI"}=$s9_ratio;
	$exon_skip{$key}->{"T3"}->{"IN"}=$s10_in;
	$exon_skip{$key}->{"T3"}->{"EX"}=$s10_ex;
	$exon_skip{$key}->{"T3"}->{"TT"}=$s10_in+$s10_ex;
	$exon_skip{$key}->{"T3"}->{"PSI"}=$s10_ratio;	
	$exon_skip{$key}->{"T_vs_C_diff"}=(-1)*$arr[($n-1)];
	$exon_skip{$key}->{"T_vs_C_Pvalue"}=$arr[($n-5)];
	$exon_skip{$key}->{"T_vs_C_FDR"}=$arr[($n-4)];		

}
close(IN);




my $file = "./MXEout_".$mode."/MXE.xls";
open(IN, $file);
my $tt = <IN>;
$tt=~s/\n$//;
my @tmp = split(/\t/, $tt);
my $n = scalar(@tmp);
#print $tmp[($n-5)]."\n";
while(<IN>)
{
    chomp;
	s/\s+$//;
	my @arr = split(/\t/, $_);
	my $key = "MXE&".$arr[0]."&".$arr[2]."&".$arr[3]."&".$arr[4]."&".$arr[5]."&".$arr[6]."&".$arr[7]."&".$arr[8]."&".$arr[9]."&".$arr[10]."&".$arr[11];
	my @tmp = split(/\,/, $arr[($n-11)]);
	my $s2_in = $tmp[0];
	my $s3_in = $tmp[1];	
	my $s4_in = $tmp[2];
	my @tmp = split(/\,/, $arr[($n-10)]);
	my $s2_ex = $tmp[0];
	my $s3_ex = $tmp[1];
	my $s4_ex = $tmp[2];	
	my @tmp = split(/\,/, $arr[($n-9)]);
	my $s8_in = $tmp[0];
	my $s9_in = $tmp[1];	
	my $s10_in = $tmp[2];
	my @tmp = split(/\,/, $arr[($n-8)]);
	my $s8_ex = $tmp[0];
	my $s9_ex = $tmp[1];	
	my $s10_ex = $tmp[2];	
	my @tmp = split(/\,/, $arr[($n-3)]);
	my $s2_ratio = $tmp[0];
	my $s3_ratio = $tmp[1];	
	my $s4_ratio = $tmp[2];
	my @tmp = split(/\,/, $arr[($n-2)]);
	my $s8_ratio = $tmp[0];
	my $s9_ratio = $tmp[1];	
	my $s10_ratio = $tmp[2];	
	$exon_skip{$key}->{"C1"}->{"IN"}=$s2_in;
	$exon_skip{$key}->{"C1"}->{"EX"}=$s2_ex;
	$exon_skip{$key}->{"C1"}->{"TT"}=$s2_in+$s2_ex;
	$exon_skip{$key}->{"C1"}->{"PSI"}=$s2_ratio;
	$exon_skip{$key}->{"C2"}->{"IN"}=$s3_in;
	$exon_skip{$key}->{"C2"}->{"EX"}=$s3_ex;
	$exon_skip{$key}->{"C2"}->{"TT"}=$s3_in+$s3_ex;
	$exon_skip{$key}->{"C2"}->{"PSI"}=$s3_ratio;
	$exon_skip{$key}->{"C3"}->{"IN"}=$s4_in;
	$exon_skip{$key}->{"C3"}->{"EX"}=$s4_ex;
	$exon_skip{$key}->{"C3"}->{"TT"}=$s4_in+$s4_ex;
	$exon_skip{$key}->{"C3"}->{"PSI"}=$s4_ratio;
	$exon_skip{$key}->{"T1"}->{"IN"}=$s8_in;
	$exon_skip{$key}->{"T1"}->{"EX"}=$s8_ex;
	$exon_skip{$key}->{"T1"}->{"TT"}=$s8_in+$s8_ex;
	$exon_skip{$key}->{"T1"}->{"PSI"}=$s8_ratio;
	$exon_skip{$key}->{"T2"}->{"IN"}=$s9_in;
	$exon_skip{$key}->{"T2"}->{"EX"}=$s9_ex;
	$exon_skip{$key}->{"T2"}->{"TT"}=$s9_in+$s9_ex;
	$exon_skip{$key}->{"T2"}->{"PSI"}=$s9_ratio;
	$exon_skip{$key}->{"T3"}->{"IN"}=$s10_in;
	$exon_skip{$key}->{"T3"}->{"EX"}=$s10_ex;
	$exon_skip{$key}->{"T3"}->{"TT"}=$s10_in+$s10_ex;
	$exon_skip{$key}->{"T3"}->{"PSI"}=$s10_ratio;	
	$exon_skip{$key}->{"T_vs_C_diff"}=(-1)*$arr[($n-1)];
	$exon_skip{$key}->{"T_vs_C_Pvalue"}=$arr[($n-5)];
	$exon_skip{$key}->{"T_vs_C_FDR"}=$arr[($n-4)];		

}
close(IN);




open(OUT, ">AS.xls");
print OUT "ASE_id"."\t".join("\t",@controllist_in_ex)."\t".join("\t",@caselist_in_ex)."\t";
print OUT join("\t",@controllist_psi)."\t".join("\t",@caselist_psi)."\t";
print OUT $comname."_diffPSI"."\t".$comname."_Pvalue"."\n";
my @samplearr=("C1", "C2", "C3", "T1", "T2", "T3");
foreach my $key (keys %exon_skip)
{
    #next if(not defined $targetaseid{$key});
	my $tt = 0;
	for(my $i = 0; $i <= $#samplearr; $i++)
	{
	if(defined $exon_skip{$key}->{$samplearr[$i]}->{"TT"})
	{
	    #$tt=$tt+$exon_skip{$key}->{$samplearr[$i]}->{"TT"};
		$tt++ if($exon_skip{$key}->{$samplearr[$i]}->{"TT"} >= 30);
	}
	}
	next if($tt<=$cutoff);	
	
	print OUT $key."\t";
	my @tmp=();
	for(my $i = 0; $i <= $#samplearr; $i++)
	{
	if(defined $exon_skip{$key}->{$samplearr[$i]}->{"IN"} and defined $exon_skip{$key}->{$samplearr[$i]}->{"EX"})
	{
	    push(@tmp, $exon_skip{$key}->{$samplearr[$i]}->{"IN"}.",".$exon_skip{$key}->{$samplearr[$i]}->{"EX"});
	}
	else
	{
		push(@tmp, "---");
	}
	}
	print OUT join("\t", @tmp)."\t";
	my @tmp=();
	for(my $i = 0; $i <= $#samplearr; $i++)
	{
	if(defined $exon_skip{$key}->{$samplearr[$i]}->{"PSI"})
	{
	    push(@tmp, $exon_skip{$key}->{$samplearr[$i]}->{"PSI"});
	}
	else
	{
		push(@tmp, "---");
	}
	}
	print OUT join("\t", @tmp)."\t";
	if(defined $exon_skip{$key}->{"T_vs_C_diff"})
	{
	    print OUT $exon_skip{$key}->{"T_vs_C_diff"}."\t";
	}
	else
	{
	    print OUT "---"."\t";
	}
	if(defined $exon_skip{$key}->{"T_vs_C_Pvalue"})
	{
	    print OUT $exon_skip{$key}->{"T_vs_C_Pvalue"}."\n";
	}
	else
	{
	    print OUT "---"."\n";
	}
}
close(OUT);







open(OUT, ">AS_raw.xls");
print OUT "ASE_id"."\t".join("\t",@controllist_in_ex)."\t".join("\t",@caselist_in_ex)."\t";
print OUT join("\t",@controllist_psi)."\t".join("\t",@caselist_psi)."\t";
print OUT $comname."_diffPSI"."\t".$comname."_Pvalue"."\n";
my @samplearr=("C1", "C2", "C3", "T1", "T2", "T3");
foreach my $key (keys %exon_skip)
{
    #next if(not defined $targetaseid{$key});
	my $in_tt = 0;
	my $ex_tt = 0;	
	for(my $i = 0; $i <= $#samplearr; $i++)
	{
	if(defined $exon_skip{$key}->{$samplearr[$i]}->{"IN"})
	{
	    $in_tt=$in_tt+$exon_skip{$key}->{$samplearr[$i]}->{"IN"};
	}
	}
	for(my $i = 0; $i <= $#samplearr; $i++)
	{
	if(defined $exon_skip{$key}->{$samplearr[$i]}->{"EX"})
	{
	    $ex_tt=$ex_tt+$exon_skip{$key}->{$samplearr[$i]}->{"EX"};
	}
	}
	next if($in_tt == 0 or $ex_tt == 0);	
	print OUT $key."\t";
	my @tmp=();
	for(my $i = 0; $i <= $#samplearr; $i++)
	{
	if(defined $exon_skip{$key}->{$samplearr[$i]}->{"IN"} and defined $exon_skip{$key}->{$samplearr[$i]}->{"EX"})
	{
	    push(@tmp, $exon_skip{$key}->{$samplearr[$i]}->{"IN"}.",".$exon_skip{$key}->{$samplearr[$i]}->{"EX"});
	}
	else
	{
		push(@tmp, "---");
	}
	}
	print OUT join("\t", @tmp)."\t";
	my @tmp=();
	for(my $i = 0; $i <= $#samplearr; $i++)
	{
	if(defined $exon_skip{$key}->{$samplearr[$i]}->{"PSI"})
	{
	    push(@tmp, $exon_skip{$key}->{$samplearr[$i]}->{"PSI"});
	}
	else
	{
		push(@tmp, "---");
	}
	}
	print OUT join("\t", @tmp)."\t";
	if(defined $exon_skip{$key}->{"T_vs_C_diff"})
	{
	    print OUT $exon_skip{$key}->{"T_vs_C_diff"}."\t";
	}
	else
	{
	    print OUT "---"."\t";
	}
	if(defined $exon_skip{$key}->{"T_vs_C_Pvalue"})
	{
	    print OUT $exon_skip{$key}->{"T_vs_C_Pvalue"}."\n";
	}
	else
	{
	    print OUT "---"."\n";
	}
}
close(OUT);
