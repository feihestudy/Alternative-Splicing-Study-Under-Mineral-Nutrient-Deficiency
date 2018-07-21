use strict;
use Getopt::Long;
$/=">";
my %targetid;
open(IN, "repre_protein.fasta");
<IN>;
while(<IN>)
{
    my @arr = split(/\n/, $_);
    my $id = shift(@arr);
    $id =~s/^>//;
	#next if(not defined $chrs{$id});
	#$id = $chrs{$id};
	my @tmp = split(/\s+/, $id);
	$id=$tmp[0];
    my $seq = join("", @arr);
    $seq=~s/>$//;
	$seq=~s/\*$//;
    #print ">".$id."\n".$seq."\n";
	$targetid{$id}=1;
}
close(IN);
$/="\n";
open(IN, "repre_transcripts_exon.gff");
my %tid2gid;
while(<IN>)
{
    chomp;
    my @arr = split(/\t/, $_);
    my $refid = $arr[0];
	#next if(not defined $chrs{$refid});
	#$arr[0] = $chrs{$refid};
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



open(IN, "repre_transcripts_exon.gff");
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
		if(defined $targetid{$tid})
		{
        print join("\t", @arr[0..7])."\t"."gene_id \"".$gid."\";"." "."transcript_id \"".$tid."\";"."\n";
		}
	}
	if ($type eq "exon")
	{
		$trnaid=~m/Parent\=(\S+)$/;
		$tid=$1;
		#$arr[2]="exon";
		if(defined $targetid{$tid})
		{
        print join("\t", @arr[0..7])."\t"."gene_id \"".$tid2gid{$tid}."\";"." "."transcript_id \"".$tid."\";"."\n";
		}
	}
}
close(IN);



