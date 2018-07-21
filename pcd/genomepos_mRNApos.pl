#use Set::Infinite;
use List::Util qw(first max maxstr min minstr reduce shuffle sum);
use Getopt::Long;
use vars qw($targetpos $gtffile $mRNAfile);
Getopt::Long::GetOptions(
	'targetpos=s'    => \$targetpos,
	'gtf=s' => \$gtffile,
	'mRNA=s' => \$mRNAfile,
);


my @exonSts;
my @exonEns;

my $strand;
open(IN, $gtffile);
while(<IN>)
{
    chomp;
    my @arr = split(/\t/, $_);
    my $refid = $arr[0];
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
	    if($arr[6] eq ".")
		{
		    $arr[6]="+";
		}	
    $strand = $arr[6];
	if($arr[2] eq "exon")
	{
	    $arr[3]=~s/\s+$//;
	    $arr[3]=~s/^\s+//;
	    $arr[4]=~s/\s+$//;
	    $arr[4]=~s/^\s+//;		
		push(@exonSts, $arr[3]);
		push(@exonEns, $arr[4]);
	}
}
close(IN);

#print $strand."\n";

my @sorted_exsts = sort {$exonSts[$a] <=> $exonSts[$b]} 0..$#exonSts;
@exonSts=@exonSts[@sorted_exsts];
@exonEns=@exonEns[@sorted_exsts];

my $len = 0;
for(my $i = 0; $i <= $#exonSts; $i++)
{
    my $st = $exonSts[$i];
	my $en = $exonEns[$i];
    $len=$len+$en-$st+1;
}
my $mRNApos=0;

if($strand eq "+")
{

for(my $i = 0; $i <= $#exonSts; $i++)
{
    my $st = $exonSts[$i];
	my $en = $exonEns[$i];
    if($targetpos >= $st and $targetpos <= $en)
	{
	   $mRNApos=$mRNApos+$targetpos-$st+1;
	   last;
	}
    if($targetpos < $st)
	{
	   #$mRNApos=$mRNApos+$en-$st+1;
	   last;
	}
    if($targetpos > $en)
	{
	   $mRNApos=$mRNApos+$en-$st+1;
	   #last;
	}	
	
}
$mRNApos=$mRNApos-1;
}


if($strand eq "-")
{

for(my $i = 0; $i <= $#exonSts; $i++)
{
    my $st = $exonSts[$i];
	my $en = $exonEns[$i];
    if($targetpos >= $st and $targetpos <= $en)
	{
	   $mRNApos=$mRNApos+$targetpos-$st+1;
	   last;
	}
    if($targetpos < $st)
	{
	   #$mRNApos=$mRNApos+$en-$st+1;
	   last;
	}
    if($targetpos > $en)
	{
	   $mRNApos=$mRNApos+$en-$st+1;
	   #last;
	}	
	
}

$mRNApos=$len-$mRNApos;
}

open(IN, $mRNAfile);
<IN>;
my $seq = <IN>;
$seq=~s/\n$//;
$seq=~s/\s+$//;


print $mRNApos."\n";
#print substr($seq, $mRNApos, 3)."\n";
