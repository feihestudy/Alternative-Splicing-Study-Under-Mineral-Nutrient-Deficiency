use strict;
use Getopt::Long;
use List::Uniq ':all';
my %targetid;
open(IN, "uniq.txt");
while(<IN>)
{
	chomp;
	my @arr = split(/\t/, $_);
	if($arr[0]=~m/^LOC/ and $arr[1]=~m/^Os/)
	{
	   $targetid{$arr[0]}=1;
	}
	if($arr[0]=~m/^LOC/ and $arr[1]=~m/^LOC/)
	{
	   $targetid{$arr[0]}=1;
	}
}
close(IN);



open(IN, "RAPDBMSU_mRNA_nx.gtf");
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
	if(defined $targetid{$tid})
	{
	}
	else
	{
		print join("\t", @arr[0..7])."\t"."gene_id "."\"".$gid."\"; "."transcript_id "."\"".$tid."\";"."\n";
	}
}
close(IN);

