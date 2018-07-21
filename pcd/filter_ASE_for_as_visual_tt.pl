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
my @ases = split(/\;/, $ASE);
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
   for(my $j = 0; $j <= $#ases; $j++)
   {
   if($aseid eq $ases[$j])
   {
       #print $_."\n";
	   for(my $i = 0; $i <= $#samples; $i++)
	   {
		   my @tmp = split(/\,/, $arr[$index[$i]]);
           $sample2info{$ases[$j]}->{"in"}=$sample2info{$ases[$j]}->{"in"}+$tmp[0];
           $sample2info{$ases[$j]}->{"ex"}=$sample2info{$ases[$j]}->{"ex"}+$tmp[1];
	   }  
   }
   }   
}
close(IN);


foreach my $key (keys %sample2info)
{
}



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
   #print $aseid."\t".$arr[0]."\n";
   for(my $j = 0; $j <= $#ases; $j++)
   {
   if($aseid eq $ases[$j])
   {
       #print $_."\n";
	   for(my $i = 0; $i <= $#samples; $i++)
	   {
           $sample2info{$ases[$j]}->{"PSI"}=$sample2info{$ases[$j]}->{"PSI"}+$arr[$index[$i]];
	   }  
   }
   }
}
close(IN);


foreach my $key (keys %sample2info)
{
    print $key."\t".$sample2info{$key}->{"in"}."\t".$sample2info{$key}->{"ex"}."\t".$sample2info{$key}->{"PSI"}/scalar(@index)."\n";
	
}

