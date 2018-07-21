use strict;
use Getopt::Long;
open(IN, "predicted_transcripts.gff");
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
		my @tmp = split(/\;/, $trnaid);
		for(my $i = 0; $i <= $#tmp; $i++)
		{
		     if($tmp[$i]=~m/^ID\=/)
			 {
			    $tid=$tmp[$i];
				$tid=~s/^ID\=//;
			 }
		     if($tmp[$i]=~m/^Locus\_id\=/)
			 {
			    $gid=$tmp[$i];
				$gid=~s/^Locus\_id\=//;
			 }
		}
		$tid2gid{$tid}=$gid;
	}
}
close(IN);



open(IN, "predicted_transcripts.gff");
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
		my @tmp = split(/\;/, $trnaid);
		for(my $i = 0; $i <= $#tmp; $i++)
		{
		     if($tmp[$i]=~m/^ID\=/)
			 {
			    $tid=$tmp[$i];
				$tid=~s/^ID\=//;
			 }
		     if($tmp[$i]=~m/^Locus\_id\=/)
			 {
			    $gid=$tmp[$i];
				$gid=~s/^Locus\_id\=//;
			 }
		}
        print join("\t", @arr[0..7])."\t"."gene_id \"".$gid."\";"." "."transcript_id \"".$tid."\";"."\n";
	}
	if ($type eq "CDS")
	{
		$trnaid=~m/Parent\=(\S+)$/;
		$tid=$1;
		$arr[2]="exon";
        print join("\t", @arr[0..7])."\t"."gene_id \"".$tid2gid{$tid}."\";"." "."transcript_id \"".$tid."\";"."\n";
	}
}
close(IN);



