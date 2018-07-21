use strict;
use Getopt::Long;
use List::Util qw(first max maxstr min minstr reduce shuffle sum);
use vars qw($inGTFfile $exGTFfile);
Getopt::Long::GetOptions(
	'inGTF=s'    => \$inGTFfile,
	'exGTF=s'    => \$exGTFfile,
);
my $strand="+";
my $inlen=0;
my @inSts;
my @inEns;
open(IN, $inGTFfile);
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
	if($arr[2] eq "CDS")
	{
	    #print $tid."\t".$arr[3]."\t".$arr[4]."\n";
		$inlen = $inlen + $arr[4] - $arr[3] + 1;
		push(@inSts, $arr[3]);
		push(@inEns, $arr[4]);
		
	}
	$strand=$arr[6];
}
close(IN);



my $exlen=0;
my @exSts;
my @exEns;
open(IN, $exGTFfile);
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
	if($arr[2] eq "CDS")
	{
	    #print $tid."\t".$arr[3]."\t".$arr[4]."\n";
		$exlen = $exlen + $arr[4] - $arr[3] + 1;
		push(@exSts, $arr[3]);
		push(@exEns, $arr[4]);
	}
}
close(IN);


my $inCDS_st=min(@inSts);
my $inCDS_en=max(@inEns);


my $exCDS_st=min(@exSts);
my $exCDS_en=max(@exEns);



my $start_codon=0;
if($strand eq "+" and $inlen >= $exlen)
{
    $start_codon=$inCDS_st;
}
if($strand eq "+" and $inlen <= $exlen)
{
    $start_codon=$exCDS_st;
}

if($strand eq "-" and $inlen >= $exlen)
{
    $start_codon=$inCDS_en;
}
if($strand eq "-" and $inlen <= $exlen)
{
    $start_codon=$exCDS_en;
}


print $start_codon."\n";







