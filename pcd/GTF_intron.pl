use strict;
use List::Util qw(first max maxstr min minstr reduce shuffle sum);
use List::Uniq ':all';
use strict;
use Getopt::Long;
use vars qw($gtffile);
Getopt::Long::GetOptions(
    'GTF=s' => \$gtffile,
);

my %transcript;
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
		}
	}
}
close(IN);


my %gene2tr;
my %ASgene;
my %introngene;
foreach my $tid (keys %transcript)
{
	my %exon;
	my $chr = $transcript{$tid}->{chr};
	my $strand = $transcript{$tid}->{strand};
	my $gid = $transcript{$tid}->{geneid};
	my @exsts = @{$transcript{$tid}->{exst}};
	my @exens = @{$transcript{$tid} -> {exen}};		
    my @sorted_exsts = sort {$exsts[$a] <=> $exsts[$b]} 0..$#exsts;	
	@exsts = @exsts[@sorted_exsts];
	@exens = @exens[@sorted_exsts];
	my $nexon = scalar(@exens);
	if($nexon > 1)
	{
	    $introngene{$gid} = 1;
	}
	push(@{$gene2tr{$gid}}, $tid);
}
foreach my $key (keys %gene2tr)
{
     @{$gene2tr{$key}}=uniq(@{$gene2tr{$key}});
	 if(defined $introngene{$key})
	 {
		 print $key."\n";
		 
	 }
}
