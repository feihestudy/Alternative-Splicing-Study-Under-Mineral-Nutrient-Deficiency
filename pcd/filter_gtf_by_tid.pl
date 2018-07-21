use strict;
use Getopt::Long;
use List::Uniq ':all';
use vars qw($targettid $allgtf);
Getopt::Long::GetOptions(
    'targettid=s' => \$targettid,
	'allgtf=s'    => \$allgtf,
);
my %transcript;
open(IN, $allgtf);
while(<IN>)
{
	chomp;
	my @arr = split(/\t/, $_);
	my $chr = $arr[0];
	my @tmp = split(/;\s/, $arr[8]);
    my $tid;
	my $gid;
    for(my $i = 0; $i <= $#tmp; $i++)
    {
       if($tmp[$i]=~m/^gene\_id/)
       {
           $gid = $tmp[$i];
           $gid =~s/^gene\_id//;
		   $gid=~s/^\s+//;
		   $gid=~s/\s+$//;
		   $gid=~s/\"//g;
       }
       if($tmp[$i]=~m/^transcript\_id/)
       {
           $tid = $tmp[$i];
           $tid =~s/^transcript\_id//;
		   $tid=~s/^\s+//;
		   $tid=~s/\;$//;
		   $tid=~s/\s+$//;
		   $tid=~s/\"//g;
       }
    }
	if($tid eq $targettid)
	{
		print $_."\n";
	}
}
close(IN);







