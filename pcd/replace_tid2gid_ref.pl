#use Set::Infinite;
use List::Uniq ':all';
use List::Util qw(first max maxstr min minstr reduce shuffle sum);
use Getopt::Long;
use Statistics::R;
my $index=1;
#63698
my %tid2gid;
open(IN, "RAPDB_MSU_mRNA_merge.bed");
while(<IN>)
{
    chomp;
	s/\s+$//;
	#print $_."\n";
	my @arr = split(/\t/, $_);
	#print $arr[0]."\t".$arr[4]."\n";
	$arr[3]=~s/\"//g;
	my @tmp =split(/\;/, $arr[3]);
	for(my $i = 0; $i <= $#tmp; $i++)
	{
		$tid2gid{$tmp[$i]}=$index;
	}
	$index=$index+1;
}
close(IN);



foreach my $key (keys %tid2gid)
{
   my $v = $tid2gid{$key};
   my $len = length($tid2gid{$key});
   if($len == 1)
   {
		$v="M0000".$v;
   }
   if($len == 2)
   {
		$v="M000".$v;
   }
   if($len == 3)
   {
		$v="M00".$v;
   }
   if($len == 4)
   {
		$v="M0".$v;
   }
   if($len == 5)
   {
		$v="M".$v;
   }
   $tid2gid{$key}=$v;
   #print $key."\t".$v."\n";
}


#=cut;
open(IN, "RAPDBMSU_mRNA.gtf");
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
	if(defined $tid2gid{$tid})
	{
        print join("\t", @arr[0..7])."\t"."gene_id \"".$tid2gid{$tid}."\";"." "."transcript_id \"".$tid."\";"."\n";
	    
	}
}
close(IN);
