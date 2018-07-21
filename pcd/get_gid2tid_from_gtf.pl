use strict;
use Getopt::Long;
use List::Uniq ':all';
use vars qw($gtffile);
Getopt::Long::GetOptions(
    'gtf=s'    =>\$gtffile,
);
my %tids;
my %gid2tids;
open(IN, $gtffile);
while(<IN>)
{
    chomp;
	#print $_."\n";
	my @arr = split(/\t/, $_);
	$arr[8]=~s/\;$//;
	#print $arr[8]."\n";
	#next if ($arr[2] ne "mRNA");
	next if(m/^\#/);
	next if($arr[2] ne "transcript");
	my @tmp = split(/;\s+/, $arr[8]);
	my $tid = 'NULL';
	my $gid = 'NULL';
	for(my $i = 0; $i <= $#tmp; $i++)
	{
	    $tmp[$i]=~s/\"//g;
	    if($tmp[$i]=~m/^transcript\_id/)
		{
			$tid = $tmp[$i];
			$tid=~s/transcript\_id\s+//; 
			$tid=~s/\;$//g;
			$tid=~s/\"//g;
			
		}
	    if($tmp[$i]=~m/^gene\_id/)
		{
			$gid = $tmp[$i];
			$gid=~s/gene\_id\s+//; 
			$gid=~s/\;$//g;
			$gid=~s/\"//g;
		}
	}
	#print $tid."\t".$gid."\n";
	my $key = $gid."\t".$tid;
	$gid2tids{$key}=1;
}
close(IN);

print join("\n", keys %gid2tids)."\n";

