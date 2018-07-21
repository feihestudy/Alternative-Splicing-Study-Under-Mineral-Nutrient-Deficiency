use strict;
use List::Util qw(first max maxstr min minstr reduce shuffle sum);
use strict;
use Getopt::Long;
use vars qw($gtffile $sjfile $sumcutoff);
Getopt::Long::GetOptions(
    'GTF=s' => \$gtffile,
	'SJ=s' => \$sjfile,
	'sumcutoff=s' => \$sumcutoff,
);

use List::Util qw(first max maxstr min minstr reduce shuffle sum);
use List::Uniq ':all';

my %sj2c;
open(IN, $sjfile);
<IN>;
while(<IN>)
{
    chomp;
	s/\s+$//;
	my @arr = split(/\t/, $_);
	my $id =shift(@arr);
	my $sum=0;
	for(my $i = 0; $i <= $#arr; $i++)
	{
	   $sum=$sum+$arr[$i];
	}
	if($sum >= $sumcutoff)
	{
	    $sj2c{$id}=1;
	}
}
close(IN);


my %transcript;
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


my %gene2ok;
my %singleexon;

my %ASgene;
my %introngene;



foreach my $tid (keys %transcript)
{
	#my $gid = "Os05g0481100";
	#print join(",", @{$transcript{"Os05g0481100"} -> {exst}})."\n";
	#print $tid."\n";
	my %exon;
	my $chr = $transcript{$tid}->{chr};
	my $strand = $transcript{$tid}->{strand};
	my $gid = $transcript{$tid}->{geneid};
	my @exsts = @{$transcript{$tid}->{exst}};
	#print join(",",@exsts)."\n";
	my @exens = @{$transcript{$tid} -> {exen}};		
    my @sorted_exsts = sort {$exsts[$a] <=> $exsts[$b]} 0..$#exsts;	
	@exsts = @exsts[@sorted_exsts];
	@exens = @exens[@sorted_exsts];
	my $nexon = scalar(@exens);
	if($nexon == 1)
	{
	$singleexon{$tid}=1;
	next;
	}
	#print $tid."\t".join(",", @exsts)."\t".join(",",@exens)."\n";
	for(my $i = 0; $i < $#exsts; $i++)
	{
	   my $in_st = $exens[$i]+1;
	   my $in_en = $exsts[($i+1)]-1;
	   #print $tid."\t".$strand."\t".($i+1)."\t".$in_st."\t".$in_en."\n";
	   my $key = $chr."&".$in_st."&".$in_en."&".$strand;
	   if(defined $sj2c{$key})
	   {
		   push(@{$gene2ok{$tid}}, $sj2c{$key});
	   }
	   else
	   {
		  push(@{$gene2ok{$tid}}, "0");
	   }
	}
}


my %targetids;

foreach my $key (keys %gene2ok)
{
   #print $key."\t".join(",", @{$gene2ok{$key}})."\n";
   my @tmp = @{$gene2ok{$key}};
   my $sum = 0;
   for(my $i = 0; $i <= $#tmp; $i++)
   {
      $sum++ if($tmp[$i] == 0);
   }
   if($sum == 0)
   {
       #print $key."\n";
	   $targetids{$key}=1;
   }
}

foreach my $key (keys %singleexon)
{
   #print $key."\n";
   $targetids{$key}=1;
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
	if(defined $targetids{$tid})
	{
	   print $_."\n";
	}
}
close(IN);



