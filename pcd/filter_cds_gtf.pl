use strict;
use Getopt::Long;
use List::Uniq ':all';
use vars qw($gtffile);
Getopt::Long::GetOptions(
    'GTF=s' => \$gtffile,
);
my %targetid;
open(IN, $gtffile);
while(<IN>)
{
	chomp;
	s/\s+$//;
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
	$targetid{$tid}=$gid;
}
close(IN);


open(IN, "/home/database/all_MSU.gff3");
while(<IN>)
{
	chomp;
	s/\s+$//;
	my @arr = split(/\t/, $_);
	my $chr = $arr[0];
	my @tmp = split(/\;/, $arr[8]);
	$arr[8]=$tmp[1];
    my $tid;
	my $gid;
	#print $arr[8]."\n";
    if($arr[2] eq "five_prime_UTR" or $arr[2] eq "CDS" or $arr[2] eq "three_prime_UTR")
	{
	    #print $_."\n";
		$arr[8]=~s/^Parent\=//;
		#print $arr[8]."\n";
		if(defined $targetid{$arr[8]})
		{
			print join("\t", @arr[0..7])."\t"."gene_id "."\"".$targetid{$arr[8]}."\"; "."transcript_id "."\"".$arr[8]."\";"."\n";
		}
	}
}
close(IN);



open(IN, "/home/database/predicted_transcripts.gff");
while(<IN>)
{
	chomp;
	my @arr = split(/\t/, $_);
	my $chr = $arr[0];
	my @tmp = split(/;\s/, $arr[8]);
    my $tid;
	my $gid;
    if($arr[2] eq "five_prime_UTR" or $arr[2] eq "CDS" or $arr[2] eq "three_prime_UTR")
	{
	    #print $_."\n";
		$arr[8]=~s/^Parent\=//;
		#print $arr[8]."\n";
		if(defined $targetid{$arr[8]})
		{
			print join("\t", @arr[0..7])."\t"."gene_id "."\"".$targetid{$arr[8]}."\"; "."transcript_id "."\"".$arr[8]."\";"."\n";
		}
	}
}
close(IN);



open(IN, "/home/database/repre_transcripts.gff");
while(<IN>)
{
	chomp;
	my @arr = split(/\t/, $_);
	my $chr = $arr[0];
	my @tmp = split(/;\s/, $arr[8]);
    my $tid;
	my $gid;
    if($arr[2] eq "five_prime_UTR" or $arr[2] eq "CDS" or $arr[2] eq "three_prime_UTR")
	{
	    #print $_."\n";
		$arr[8]=~s/^Parent\=//;
		#print $arr[8]."\n";
		if(defined $targetid{$arr[8]})
		{
			print join("\t", @arr[0..7])."\t"."gene_id "."\"".$targetid{$arr[8]}."\"; "."transcript_id "."\"".$arr[8]."\";"."\n";
		}
	}
}
close(IN);

open(IN, "novel_AS_CDS.gtf");
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
    if($arr[2] eq "five_prime_UTR" or $arr[2] eq "CDS" or $arr[2] eq "three_prime_UTR")
	{
		if(defined $targetid{$tid})
		{
			print join("\t", @arr[0..7])."\t"."gene_id "."\"".$targetid{$tid}."\"; "."transcript_id "."\"".$tid."\";"."\n";
		}
	}
}
