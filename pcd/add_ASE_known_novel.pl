use strict;
use List::Uniq ':all';
use List::Compare;
use Getopt::Long;
use vars qw($asefile $astype);
Getopt::Long::GetOptions(
    'AStype=s'    => \$astype,
);

my %known;
open(IN, "ASE_known.txt");
while(<IN>)
{
    chomp;
	s/\s+$//;
	my @arr = split(/\t/, $_);
	$known{$arr[0]}=1;
}
close(IN);

my %novel;
open(IN, "ASE_novel.txt");
while(<IN>)
{
    chomp;
	s/\s+$//;
	my @arr = split(/\t/, $_);
	$novel{$arr[0]}=1;
}
close(IN);




open(IN, $astype);
while(<IN>)
{
    chomp;
	s/\s+$//;
	my @arr = split(/\t/, $_);
	if(defined $known{$arr[0]})
	{
	    print $_."\t"."known"."\n";
	}
	if(defined $novel{$arr[0]})
	{
	    print $_."\t"."novel"."\n";
	}
}
close(IN);


