use strict;
use List::Util qw(first max maxstr min minstr reduce shuffle sum);
use List::Uniq ':all';
use strict;
use Getopt::Long;
use vars qw($inputfile $incutoff $excutoff $ttcutoff);
Getopt::Long::GetOptions(
	'inputfile=s' =>\$inputfile,
	'totalcutoff=s' => \$ttcutoff,
);
open(IN, $inputfile);
my $tt = <IN>;
#print $tt;
while(<IN>)
{
    chomp;
	s/\s+$//;
	my @arr = split(/\t/, $_);
	my $aseid = shift(@arr);
	my $int=0;
	my $ext=0;
	for(my $i = 0; $i <= $#arr; $i++)
	{
	     my @arr3 = split(/\,/, $arr[$i]);
		 $int=$int+$arr3[0];
		 $ext=$ext+$arr3[1];		 
	}
	if($int > 0 and $ext > 0 and ($int+$ext)>=$ttcutoff)
    {
	    print $aseid."\n";
	}
}
close(IN);
