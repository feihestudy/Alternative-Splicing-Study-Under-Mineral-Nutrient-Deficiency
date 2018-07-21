use List::Uniq ':all';
use List::Util qw(first max maxstr min minstr reduce shuffle sum);
use Getopt::Long;
use vars qw($gtf);
Getopt::Long::GetOptions(
    'gtf=s'    =>\$gtf,
);
my %transcript;
open(IN, $gtf);
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
		push(@{$transcript{$tid} -> {exst}}, $arr[3]);
		push(@{$transcript{$tid} -> {exen}}, $arr[4]);
	}
}
close(IN);

my %transcript2;
foreach my $key (keys %transcript)
{
   my $chr = $transcript{$key} -> {chr};
   my $strand = $transcript{$key} -> {strand};
   my @sts = @{$transcript{$key} -> {exst}};
   my @ens = @{$transcript{$key} -> {exen}};
   my @sorted_exst = sort {$sts[$a] <=> $sts[$b]} 0..$#sts;
   @sts = @sts[@sorted_exst];
   @ens = @ens[@sorted_exst];
   @{$transcript{$key} -> {exst}} = @sts;
   @{$transcript{$key} -> {exen}} = @ens;
   my $vv = $chr."//".$strand."//".join("//",@sts)."//".join("//",@ens); 
   push(@{$transcript2{$vv}}, $key);
}


foreach my $key (keys %transcript2)
{
    #print $key."\t".scalar(@{$transcript2{$key}})."\t".join(",", sort @{$transcript2{$key}})."\n";
	my $num = scalar(@{$transcript2{$key}});
	if($num >= 2)
	{
	    print join("\t", sort @{$transcript2{$key}})."\n";
	}
}


