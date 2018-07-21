use strict;
use Getopt::Long;
use List::Uniq ':all';
use vars qw($tid $transfile);
Getopt::Long::GetOptions(
	'tid=s' =>\$tid,
	'trans=s' =>\$transfile,
);


$/=">";
open(IN, $transfile);
<IN>;
my $n;
while(<IN>)
{
    chomp;
	s/\s+$//;
	my @arr = split(/\n/, $_);
	my $idone = shift(@arr);
	my $seq = join("", @arr);
	$seq=~s/\>$//;
	if($idone eq $tid)
	{
		print ">".$idone."\n".$seq."\n";
	}
}
close(IN);





