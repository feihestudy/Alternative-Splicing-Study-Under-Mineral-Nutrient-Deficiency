use strict;
use Getopt::Long;
use List::Util qw(first max maxstr min minstr reduce shuffle sum);
use vars qw($tmapfile $gtffile);
Getopt::Long::GetOptions(
	'tmap=s'    => \$tmapfile,
	'gtf=s' => \$gtffile,
);
$/="\n";
open(IN, $tmapfile);
my %leftid;
my %id2id;
while(<IN>)
{
	#print $_;
	chomp;
	my @arr = split(/\t/, $_);
	my $type = $arr[2];
	#print $type."\n";
    if($type eq "=")
	{
	   $leftid{$arr[1]} = $arr[3];
	}
}
close(IN);



open(IN, $gtffile);
my %reftranscript;
while(<IN>)
{
    chomp;
	my @arr = split(/\t/, $_);
	$arr[8]=~s/\;$//;
	#next if(not defined $chr2seq{$arr[0]});
	#next if($arr[2] ne "exon");
	if($arr[2] eq "mRNA")
	{
	    $arr[2]="transcript";
	}
	my @tmp = split(/;\s+/, $arr[8]);
	$tmp[0]=~s/gene\_id\s+//;
	$tmp[0]=~s/\"//g;	
	$tmp[1]=~s/transcript\_id\s+//;
	$tmp[1]=~s/\"//g;
	#my $tid = $tmp[0]."//".$tmp[1];
	my $tid = $tmp[1];
	my $gid = $tmp[0];
	
	if(defined $leftid{$tid})
	{
		print join("\t", @arr[0..7])."\t"."gene_id \"".$leftid{$tid}."\";"." "."transcript_id \"".$tid."\";"."\n";
	}
}
close(IN);


