use strict;
use List::Util qw(first max maxstr min minstr reduce shuffle sum);
use List::Uniq ':all';
use strict;
use Getopt::Long;
use vars qw($targetfile $ASfile);
Getopt::Long::GetOptions(
    'target=s' => \$targetfile,
	'AS=s' => \$ASfile,
);

my %ASElist;
open(IN, $targetfile);
<IN>;
while(<IN>)
{
    chomp;;
	s/\s+$//;
	my @arr = split(/\t/, $_);
	$ASElist{$arr[0]}=1;
}
close(IN);


open(IN, $ASfile);
my $tt = <IN>;
print $tt;
while(<IN>)
{
    chomp;
	s/\s+$//;
	my @arr = split(/\t/, $_);
	my $aseid = shift(@arr);
    #my @arr2 = split(/\&/, $aseid);
	#print $aseid."\n";
	if(defined $ASElist{$aseid})
	{
	   print $_."\n";
	}
}
close(IN);
