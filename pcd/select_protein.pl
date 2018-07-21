use strict;
my %cpcids;
open(IN, "novel.fa.cpc2");
<IN>;
while(<IN>)
{
    chomp;
	s/\s+$//;
	my @arr = split(/\t/, $_);
	if($arr[$#arr] eq "coding")
	{
	   $cpcids{$arr[0]}=1;
	}
}
close(IN);

my %pid2tid;
$/ = ">";
open(IN, "novel.fa.transdecoder.pep.aa");
<IN>;
while(<IN>)
{
    my @arr = split(/\n/, $_);
    my $id = shift(@arr);
	my $seq = join("", @arr);
	$seq=~s/\>$//;
	$seq=~s/\*//g;
	my @tmp = split(/\:\:/, $id);
	if(defined $cpcids{$tmp[1]})
	{
	$pid2tid{$id}=$tmp[1];
	}
}
close(IN);


$/="\n";


open(OUT, ">novel.fa.transdecoder_CDS.gtf");
open(IN, "novel.fa.transdecoder.gff3");
while(<IN>)
{
    chomp;
	s/\s+$//;
	my @arr = split(/\t/, $_);
	my $type = $arr[2];
	my $id = $arr[0];
	if($type eq "CDS")
	{
		 my @tmp = split(/\;/, $arr[8]);
		 $tmp[1]=~s/^Parent\=//;
		 if(defined $pid2tid{$tmp[1]})
		 {
				print OUT join("\t", @arr[0..7])."\t"."gene_id "."\"".$pid2tid{$tmp[1]}."\"; "."transcript_id "."\"".$pid2tid{$tmp[1]}."\";"."\n";
		 }
	}
	if($type eq "five_prime_UTR")
	{
		 my @tmp = split(/\;/, $arr[8]);
		 $tmp[1]=~s/^Parent\=//;
		 if(defined $pid2tid{$tmp[1]})
		 {
				print OUT join("\t", @arr[0..7])."\t"."gene_id "."\"".$pid2tid{$tmp[1]}."\"; "."transcript_id "."\"".$pid2tid{$tmp[1]}."\";"."\n";
		 }
	}
	if($type eq "three_prime_UTR")
	{
		 my @tmp = split(/\;/, $arr[8]);
		 $tmp[1]=~s/^Parent\=//;
		 if(defined $pid2tid{$tmp[1]})
		 {
				print OUT join("\t", @arr[0..7])."\t"."gene_id "."\"".$pid2tid{$tmp[1]}."\"; "."transcript_id "."\"".$pid2tid{$tmp[1]}."\";"."\n";
		 }
	}
}
close(IN);
close(OUT);





