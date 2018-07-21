use strict;
use Getopt::Long;
use List::Uniq ':all';
use List::Util qw(first max maxstr min minstr reduce shuffle sum);
use vars qw($gtffile $lencutoff);
Getopt::Long::GetOptions(
    'GTF=s' => \$gtffile,
	'LEN=s' => \$lencutoff,
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
	if($arr[2] eq "exon")
	{
		$exon{$tid} = $exon{$tid} + $arr[4] - $arr[3] + 1;
	}
}
close(IN);

my %targetids;
foreach my $key (keys %exon)
{
   if($exon{$key} >= 200)
   {
	   $targetids{$key}=1;
   }
}


open(IN, $gtffile);
while(<IN>)
{
	chomp;
	my @arr = split(/\t/, $_);
	my $chr = $arr[0];
	my @tmp = split(/;\s/, $arr[8]);
    my $tid;
	my $gid;
	my $fpkm;
    for(my $i = 0; $i <= $#tmp; $i++)
    {
       if($tmp[$i]=~m/^transcript\_id/)
       {
           $tid= $tmp[$i];		   
           $tid=~s/^transcript\_id//;
		   $tid=~s/^\s+//;
		   $tid=~s/\"//g;
		   $tid=~s/\;$//;
       }
    }
	if(defined $targetids{$tid})
	{
	   print $_."\n";
	}
}
close(IN);










