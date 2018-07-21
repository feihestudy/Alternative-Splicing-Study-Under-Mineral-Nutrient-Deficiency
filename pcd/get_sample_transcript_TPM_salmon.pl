use strict;
use Getopt::Long;
use vars qw($samplelist $dir);
Getopt::Long::GetOptions(
    'sample=s'    => \$samplelist,
	'dir=s' => \$dir,
);

my %refid2count;

my @samples = split(/\,/, $samplelist);
for(my $i = 0; $i <= $#samples; $i++)
{
	my $samfile = $dir."/".$samples[$i]."/"."quant.sf";
	open(IN, $samfile);
	<IN>;
	while(<IN>)
	{
		chomp;
		s/\s+$//;
		my @arr = split(/\t/, $_);
		$arr[4]=~s/\s+$//;
	    $refid2count{$arr[0]}->{$samples[$i]}=$arr[3];
	}
	close(IN);
}
print "transcript_ID"."\t".join("\t", @samples)."\n";
foreach my $key (keys %refid2count)
{
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
       