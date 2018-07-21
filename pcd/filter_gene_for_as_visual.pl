use strict;
use Getopt::Long;
use List::Util qw(first max maxstr min minstr reduce shuffle sum);
use List::Uniq ':all';
use vars qw($gene $bedfile);
Getopt::Long::GetOptions(
    'gene=s' => \$gene,
	'bed=s' => \$bedfile,
);
open(IN, $bedfile);
while(<IN>)
{
    chomp;
	s/\s+$//;
	my @arr = split(/\t/, $_);
	if($arr[3] eq $gene)
	{
       print $_."\n";	   
	}

}
close(IN);



