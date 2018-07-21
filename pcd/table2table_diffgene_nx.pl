use strict;
use Getopt::Long;
use vars qw($samplepairfile $targetfile $rowcountfile $normcountfile $fpkmcountfile);
Getopt::Long::GetOptions(
    'samplepair=s' => \$samplepairfile,
	'target=s' => \$targetfile,
	'rawcount=s' => \$rowcountfile,
	'normcount=s' => \$normcountfile,
	'TPM=s' => \$fpkmcountfile,
);
my %targetids;
open(IN, $targetfile);
<IN>;
while(<IN>)
{
    chomp;
	s/\s+$//;
	my @arr = split(/\t/, $_);
	$arr[0]=~s/\s+$//;
	$targetids{$arr[0]}=1;
}
close(IN);


my %gene2exist;
my %rawcount_table;
my %normcount_table;
my %FPKM_table;
open(IN, $rowcountfile);
my $tt = <IN>;
$tt=~s/\n$//;
$tt=~s/\s+$//;
my @tmp = split(/\t/, $tt);
shift(@tmp);
my @recindex;
for(my $i = 0; $i <= $#tmp; $i++)
{
   $tmp[$i]=~s/\_\S+$//;
   if($tmp[$i]=~m/RecR/)
   {
        
   }
   else
   {
		push(@recindex, $i);
   }
   $tmp[$i]=$tmp[$i]."_raw_count";
   
}
my $rawcount_tt = join("\t", @tmp[@recindex]);
while(<IN>)
{
    chomp;
	s/\s+$//;
	s/\"//g;
	#print $_."\n";
	my @arr = split(/\t/, $_);
	my $id = shift(@arr);
	$id=~s/\s+$//;
	if(defined $targetids{$id})
	{
		$rawcount_table{$id}->{"rawcount"}=join("\t",@arr[@recindex]);
	#print $_."\n";
	}
}
close(IN);


foreach my $key (keys %rawcount_table)
{
   #print $key."\t".
}


open(IN, $normcountfile);
my $tt = <IN>;
$tt=~s/\n$//;
$tt=~s/\s+$//;
my @tmp = split(/\t/, $tt);
shift(@tmp);
my @recindex;
for(my $i = 0; $i <= $#tmp; $i++)
{
   $tmp[$i]=~s/\_\S+$//;
   if($tmp[$i]=~m/RecR/)
   {
        
   }
   else
   {
		push(@recindex, $i);
   }   
   $tmp[$i]=$tmp[$i]."_normalized_count";
}
my $normcount_tt = join("\t", @tmp[@recindex]);
while(<IN>)
{
    chomp;
	s/\s+$//;
	s/\"//g;
	my @arr = split(/\t/, $_);
	my $id = shift(@arr);
	$id=~s/\s+$//;
	if(defined $targetids{$id})
	{
	$normcount_table{$id}->{"normcount"}=join("\t", @arr[@recindex]);
	}
}
close(IN);




open(IN, $fpkmcountfile);
my $tt = <IN>;
$tt=~s/\n$//;
$tt=~s/\s+$//;
my @tmp = split(/\t/, $tt);
shift(@tmp);
my @recindex;
for(my $i = 0; $i <= $#tmp; $i++)
{
   $tmp[$i]=~s/\_\S+$//;
    if($tmp[$i]=~m/RecR/)
   {
        
   }
   else
   {
		push(@recindex, $i);
   }   
   $tmp[$i]=$tmp[$i]."_TPM";
}
my $fpkm_tt = join("\t", @tmp[@recindex]);
while(<IN>)
{
    chomp;
	s/\s+$//;
	s/\"//g;
	my @arr = split(/\t/, $_);
	my $id = shift(@arr);
	$id=~s/\s+$//;
	if(defined $targetids{$id})
	{
	    $FPKM_table{$id}->{"fpkm"}=join("\t", @arr[@recindex]);
		#print $_."\n";
	}
}
close(IN);



my @samplepairs = split(/\,/, $samplepairfile);
#print join("\n", @samplepairs)."\n";

for(my $i = 0; $i <= $#samplepairs; $i++)
{
     my $file = $samplepairs[$i]."_gene_exp.xls";
	 #print $file."\n";
	open(IN, $file);
	<IN>;
	while(<IN>)
	{
	    chomp;
		s/\s+$//;
		s/\"//g;
		my @arr = split(/\t/, $_);
		#print $arr[0]."\n";
	$arr[0]=~s/\s+$//;
		my $fc = $arr[($#arr-2)];
		if(defined $targetids{$arr[0]})
		{
		$gene2exist{$arr[0]}->{$samplepairs[$i]}=$arr[($#arr-2)]."\t".$arr[($#arr-1)]."\t".$arr[($#arr)];
		}
   }
	close(IN);	 
}
#print join("\n", keys %gene2exist)."\n";
#print "Gene_locus_id"."\t".join("\t", @samplepairs)."\n";
print "Gene_locus_id"."\t";
print $rawcount_tt."\t";
print $normcount_tt."\t";
print $fpkm_tt."\t";
print "H1piMinus_vs_H1piPlus_normalized_count_log2FoldChange"."\t"."H1piMinus_vs_H1piPlus_Pvalue"."\t"."H1piMinus_vs_H1piPlus_FDR"."\t";
print "H6piMinus_vs_H6piPlus_normalized_count_log2FoldChange"."\t"."H6piMinus_vs_H6piPlus_Pvalue"."\t"."H6piMinus_vs_H6piPlus_FDR"."\t";
print "H24piMinus_vs_H24piPlus_normalized_count_log2FoldChange"."\t"."H24piMinus_vs_H24piPlus_Pvalue"."\t"."H24piMinus_vs_H24piPlus_FDR"."\t";
print "D3piMinus_vs_D3piPlus_normalized_count_log2FoldChange"."\t"."D3piMinus_vs_D3piPlus_Pvalue"."\t"."D3piMinus_vs_D3piPlus_FDR"."\t";
print "D7piMinus_vs_D7piPlus_normalized_count_log2FoldChange"."\t"."D7piMinus_vs_D7piPlus_Pvalue"."\t"."D7piMinus_vs_D7piPlus_FDR"."\t";
print "D21piMinus_vs_D21piPlus_normalized_count_log2FoldChange"."\t"."D21piMinus_vs_D21piPlus_Pvalue"."\t"."D21piMinus_vs_D7piPlus_FDR"."\t";
print "D21H1piMinus_vs_D21H1piPlus_normalized_count_log2FoldChange"."\t"."D21H1piMinus_vs_D21H1piPlus_Pvalue"."\t"."D21H1piMinus_vs_D21H1piPlus_FDR"."\t";
print "D21H6piMinus_vs_D21H6piPlus_normalized_count_log2FoldChange"."\t"."D21H6piMinus_vs_D21H6piPlus_Pvalue"."\t"."D21H6piMinus_vs_D21H6piPlus_FDR"."\t";
print "D21H24piMinus_vs_D21H24piPlus_normalized_count_log2FoldChange"."\t"."D21H24piMinus_vs_D21H24piPlus_Pvalue"."\t"."D21H24piMinus_vs_D21H24piPlus_FDR"."\n";
foreach my $key (keys %gene2exist)
{
    print $key."\t";
	print $rawcount_table{$key}->{"rawcount"}."\t";
	print $normcount_table{$key}->{"normcount"}."\t";
	print $FPKM_table{$key}->{"fpkm"}."\t";	
    my @values;  
    for(my $i = 0; $i <= $#samplepairs; $i++)
	{
	     if(defined $gene2exist{$key}->{$samplepairs[$i]})
		 {
		    push(@values, $gene2exist{$key}->{$samplepairs[$i]});
		 }
		 else
		 {
		    push(@values, "---"."\t"."---"."\t"."---");		 
		 }
	}
	print join("\t", @values)."\n";
}



