use strict;
use Getopt::Long;
use List::Util qw(first max maxstr min minstr reduce shuffle sum);
use List::Uniq ':all';
use vars qw($ASE $inexinputfile $psiinputfile $samplelist);
Getopt::Long::GetOptions(
    'ASE=s' => \$ASE, 
	'inex=s' => \$inexinputfile,
	'PSI=s' => \$psiinputfile,
	'samplelist=s' => \$samplelist,
);
my %sample2info;
open(IN, $inexinputfile);
my $tt = <IN>;
$tt=~s/\n$//;
$tt=~s/\s+$//;
my @tmp = split(/\t/, $tt);
shift(@tmp);
my @samples = split(/\,/, $samplelist);
my @index=();
for(my $i = 0; $i <= $#samples; $i++)
{
    for(my $j = 0; $j <= $#tmp; $j++)
	{
	    push(@index, $j) if($samples[$i] eq $tmp[$j]);	
	}
}
while(<IN>)
{
   chomp;
   s/\s+$//;
   my @arr = split(/\t/, $_);
   my $aseid = shift(@arr);
   if($aseid eq $ASE)
   {
       #print $_."\n";
	   for(my $i = 0; $i <= $#samples; $i++)
	   {
           	$sample2info{$samples[$i]}->{"in_ex"}=$arr[$index[$i]];
	   }
  
   }
   
}
close(IN);



open(IN, $psiinputfile);
my $tt = <IN>;
$tt=~s/\n$//;
$tt=~s/\s+$//;
my @tmp = split(/\t/, $tt);
shift(@tmp);
my @samples = split(/\,/, $samplelist);
my @index=();
for(my $i = 0; $i <= $#samples; $i++)
{
    for(my $j = 0; $j <= $#tmp; $j++)
	{
	    push(@index, $j) if($samples[$i] eq $tmp[$j]);	
	}
}
while(<IN>)
{
   chomp;
   s/\s+$//;
   my @arr = split(/\t/, $_);
   my $aseid = shift(@arr);
   if($aseid eq $ASE)
   {
	   for(my $i = 0; $i <= $#samples; $i++)
	   {
           	$sample2info{$samples[$i]}->{"PSI"}=$arr[$index[$i]];
	   }
  
   }
   
}
close(IN);


for(my $i = 0; $i <= $#samples; $i++)
{
    if(defined $sample2info{$samples[$i]}->{"in_ex"} and defined $sample2info{$samples[$i]}->{"PSI"})
	{
	    my $ssname = $samples[$i];
		$ssname=~s/\_(R|S)$//;
	    open(OUT, ">".$ssname."_junction.txt");
	    print OUT $sample2info{$samples[$i]}->{"in_ex"}.",".$sample2info{$samples[$i]}->{"PSI"}."\n";
		close(OUT);
	
	}


}


