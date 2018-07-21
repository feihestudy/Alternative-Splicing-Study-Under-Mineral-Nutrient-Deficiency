use strict;
use Getopt::Long;
open(IN, "all_MSU.gff3");
my %tid2gid;
while(<IN>)
{
    chomp;
    my @arr = split(/\t/, $_);
    my $refid = $arr[0];
	my $type = $arr[2];
    my $st = $arr[3];
    my $en = $arr[4];
	my $strand = $arr[6];
    my $trnaid = $arr[8];
	my $tid;
	my $gid;
	if ($type eq "mRNA")
	{
		$trnaid=~m/ID\=(\S+)\;Name\=(\S+)\;Parent\=(\S+)$/;
		$tid=$1;
		$gid=$3;
		$tid2gid{$tid}=$gid;
	}
}
close(IN);



open(IN, "all_MSU.gff3");
while(<IN>)
{
    chomp;
    my @arr = split(/\t/, $_);
    my $refid = $arr[0];
	my $type = $arr[2];
    my $st = $arr[3];
    my $en = $arr[4];
	my $strand = $arr[6];
    my $trnaid = $arr[8];
	my $tid;
	my $gid;
	if ($type eq "mRNA")
	{
		$trnaid=~m/ID\=(\S+)\;Name\=(\S+)\;Parent\=(\S+)$/;
		$tid=$1;
		$gid=$3;
        print join("\t", @arr[0..7])."\t"."gene_id \"".$gid."\";"." "."transcript_id \"".$tid."\";"."\n";
	}
	if ($type eq "exon")
	{
		$trnaid=~m/\;Parent\=(\S+)$/;
		$tid=$1;
        print join("\t", @arr[0..7])."\t"."gene_id \"".$tid2gid{$tid}."\";"." "."transcript_id \"".$tid."\";"."\n";
	}
}
close(IN);



