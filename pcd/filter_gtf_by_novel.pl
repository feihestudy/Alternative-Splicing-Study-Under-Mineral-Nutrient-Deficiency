use strict;
use Getopt::Long;
use List::Uniq ':all';
use List::Util qw(first max maxstr min minstr reduce shuffle sum);
use vars qw($gtffile $lencutoff);
Getopt::Long::GetOptions(
    'GTF=s' => \$gtffile,
);
my %exon;
open(IN, $gtffile);
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
		   $gid=~s/\;$//;
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
    #print $tid."\t".$arr[2]."\n";
	if($tid=~m/^TU\d+/)
	{
	   print $_."\n";
	}
}
close(IN);

