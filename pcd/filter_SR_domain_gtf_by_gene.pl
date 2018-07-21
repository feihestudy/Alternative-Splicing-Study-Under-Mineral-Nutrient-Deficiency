use strict;
use Getopt::Long;
use List::Util qw(first max maxstr min minstr reduce shuffle sum);
use vars qw($SRdomaingtf $geneid);
Getopt::Long::GetOptions(
    'SRdomain=s' =>\$SRdomaingtf,
	'gene=s' =>\$geneid,
);
my $targetgeneid=$geneid;

open(IN, $SRdomaingtf);
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
	if($targetgeneid eq $gid)
	{
		print join("\t", @arr[0..7])."\t"."gene_id "."\"".$gid."\"; "."transcript_id "."\"".$tid."\";"."\n";
	}
}
close(IN);
