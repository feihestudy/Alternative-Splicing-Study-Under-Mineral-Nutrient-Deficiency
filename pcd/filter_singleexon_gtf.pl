use strict;
use Getopt::Long;
use List::Uniq ':all';
use vars qw($gtffile);
Getopt::Long::GetOptions(
    'gtf=s'    =>\$gtffile,
);
use List::Uniq ':all';
my %tid2gid;
my %gids;
my %transcript;
open(IN, $gtffile);
while(<IN>)
{
    chomp;
	my @arr = split(/\t/, $_);
	$arr[8]=~s/\;$//;
	my @tmp = split(/;\s+/, $arr[8]);
	my $tid = '';
	my $gid = '';
	for(my $i = 0; $i <= $#tmp; $i++)
	{
	    $tmp[$i]=~s/\"//g;
	    if($tmp[$i]=~m/^transcript\_id/)
		{
			$tid = $tmp[$i];
			$tid=~s/transcript\_id\s+//; 
			$tid=~s/\;$//g;
			$tid=~s/\"//g;
			
		}
	    if($tmp[$i]=~m/^gene\_id/)
		{
			$gid = $tmp[$i];
			$gid=~s/gene\_id\s+//; 
			$gid=~s/\;$//g;
			$gid=~s/\"//g;
		}
		#print $tid."\t".$gid."\n";
	}
    my $strand = $arr[6];  	
	unless(defined $transcript{$tid})
	{
	    $transcript{$tid} = {
			chr => $arr[0],
			strand => $strand,
			geneid=>$gid,
		};
	}
	if($arr[2] eq "exon")
	{
		push(@{$transcript{$tid}->{exst}}, $arr[3]);
		push(@{$transcript{$tid}->{exen}}, $arr[4]);
		if($arr[3] >= $arr[4])
		{
		   # print $_."\n";
		}
	}
}
close(IN);


my %keytids;
foreach my $key (keys %transcript)
{
  if(scalar(@{$transcript{$key}->{exst}}) > 1)
  {
    #$exonm2++;
	$keytids{$key}=1;
  }
}

open(IN, $gtffile);
while(<IN>)
{
    chomp;
    my @arr = split(/\t/, $_);
    my $refid = $arr[0];
	#next if(not defined $chrs{$refid});
	#$arr[0] = $chrs{$refid};
	my @tmp = split(/;\s/, $arr[8]);
    my $tid;
	my $gid;
    for(my $i = 0; $i <= $#tmp; $i++)
    {
       if($tmp[$i]=~m/^transcript\_id/)
       {
           $tid = $tmp[$i];
           $tid =~s/^transcript\_id//;
		   $tid=~s/^\s+//;
		   $tid=~s/\s+$//;
		   $tid=~s/\"//g;
		   $tid=~s/\;//g;
       }
       if($tmp[$i]=~m/^gene\_id/)
       {
           $gid = $tmp[$i];
           $gid =~s/^gene\_id//;
		   $gid=~s/^\s+//;
		   $gid=~s/\s+$//;
		   $gid=~s/\"//g;
		   $gid=~s/\;//g;
       }
    }
	if(defined $keytids{$tid})
	{
	   print $_."\n";
	}
}
close(IN);




