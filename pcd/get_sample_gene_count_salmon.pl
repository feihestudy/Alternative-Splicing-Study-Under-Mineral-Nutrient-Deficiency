use strict;
use Getopt::Long;
use vars qw($samplelist $dir);
Getopt::Long::GetOptions(
    'sample=s'    => \$samplelist,
	'dir=s' => \$dir,
);

my %tid2gid;
open(IN, "gid2tid.txt");
while(<IN>)
{
   chomp;
   s/\s+$//;
   my @arr = split(/\t/, $_);
   $tid2gid{$arr[1]}=$arr[0];
}
close(IN);



my %refid2count;

my @samples = split(/\,/, $samplelist);
#$/="\n";
for(my $i = 0; $i <= $#samples; $i++)
{
	#my $samfile = "C3.sam";
	#print $samples[$i]."\n";
	my $samfile = $dir."/".$samples[$i]."/"."quant.sf";
	#print $samfile."\n";
	my %samplecount;
	open(IN, $samfile);
	<IN>;
	while(<IN>)
	{
		chomp;
		s/\s+$//;
		my @arr = split(/\t/, $_);
		$arr[4]=~s/\s+$//;
		if(defined $tid2gid{$arr[0]})
		{
		    push(@{$samplecount{$tid2gid{$arr[0]}}}, $arr[4]);
		}
	}
	close(IN);
	foreach my $key (keys %samplecount)
	{
	   my @tmp = @{$samplecount{$key}};
	   #print $key."\t".join(",",@tmp)."\n";
	   my $ttcount = 0;
	   for(my $j=0; $j<=$#tmp; $j++)
	   {
	       #$tmp[$j]=ord($tmp[$j]);
	       $ttcount=$ttcount+$tmp[$j];
	   }	   
	   #my $ttcount = sum(@tmp);
	   $ttcount=round($ttcount,0);
	   #print $key."\t".$ttcount."\n";
	   $refid2count{$key}->{$samples[$i]}=$ttcount;
	}
}
print "gene_ID"."\t".join("\t", @samples)."\n";


#=cut;
foreach my $key (keys %refid2count)
{
    #my %tmp = %{$refid2count{$key}};
	my @tmpcount;
	for(my $j = 0; $j <= $#samples; $j++)
	{
	    if(defined $refid2count{$key}->{$samples[$j]})
		{
			push(@tmpcount, $refid2count{$key}->{$samples[$j]});
		}
		else
		{
		    push(@tmpcount, 0);
		}
	}
	print $key."\t".join("\t", @tmpcount)."\n";
}



sub round {
    my $val = shift;
    my $col = shift;
    my $r = 10 ** $col;
    my $a = ($val > 0) ? 0.5 : -0.5;
    return int($val * $r + $a) / $r;
}
       