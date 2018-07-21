use strict;
use Getopt::Long;
use vars qw($AStype $dir);
Getopt::Long::GetOptions(
	'AStype=s' => \$AStype,
	'dir=s' => \$dir,
);

my $otfile = $AStype.".xls";
open(OUT1, ">$otfile");
my %se;
my $file = $dir."/"."AS.".$AStype.".txt";
open(IN, $file);
my $tt = <IN>;
$tt=~s/\n$//;
$tt=~s/\s+$//;
my @tmp = split(/\t/, $tt);
print OUT1 join("\t", @tmp[1..$#tmp])."\t";
while(<IN>)
{
    chomp;
	s/\s+$//;
	s/\"//g;
	my @arr = split(/\t/, $_);
	$se{$arr[0]}=join("\t", @arr[1..$#arr]);
}
close(IN);



open(IN, "rMATS_Result.txt");
my $tt = <IN>;
print OUT1 $tt;
while(<IN>)
{
    chomp;
	s/\s+$//;
	my @arr = split(/\t/, $_);
	if(defined $se{$arr[0]})
	{
	   print OUT1 $se{$arr[0]}."\t".$_."\n";
	}
}
close(IN);
close(OUT1);

