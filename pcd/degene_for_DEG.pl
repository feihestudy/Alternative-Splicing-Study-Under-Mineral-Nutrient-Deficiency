use strict;
use List::Uniq ':all';
use Getopt::Long;
use List::Util qw(first max maxstr min minstr reduce shuffle sum);
use Getopt::Long;
use Statistics::R;
use vars qw($degenfile);
Getopt::Long::GetOptions(
    'degene=s'    =>\$degenfile,
);



open(IN, $degenfile);
#<IN>;
my $tt = <IN>;
$tt=~s/\s+$//;
my %exist;
while(<IN>)
{
    chomp;
	s/\s+$//;
	s/\"//g;
	s/\'//g;	
	my @arr = split(/\t/, $_);
	my $gid = $arr[($#arr-5)];
	next if($gid eq "NULL");
	my @gids = split(/\/\//, $gid);
	for(my $i = 0; $i <= $#gids; $i++)
	{
	$exist{$gids[$i]}=1;
	}
}
close(IN);

print join("\n", keys %exist)."\n";
