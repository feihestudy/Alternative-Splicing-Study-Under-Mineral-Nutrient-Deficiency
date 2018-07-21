use strict;
use List::Util qw(first max maxstr min minstr reduce shuffle sum);
use Getopt::Long;
use vars qw($dir $samplelist $samplecutoff);
Getopt::Long::GetOptions(
    'dir=s' => \$dir,
	'samplelist=s' => \$samplelist,
	'samplecutoff=s' => \$samplecutoff,
);
my %junction;
my @samples = split(/\,/, $samplelist);

for(my $j = 0; $j <= $#samples; $j++)
{
	my $file = $dir."/".$samples[$j]."/SJ.out.tab";
	open(IN, $file);
	while(<IN>)
	{
		chomp;
		s/\s+$//;
		#print $_."\n";
		my @arr = split(/\t/, $_);
		if($arr[3] eq "1")
		{
		   $arr[3]="+";
		}
		if($arr[3] eq "2")
		{
		   $arr[3]="-";
		}
		if($arr[3] eq "0")
		{
		   $arr[3]=".";
		}
		my $key = $arr[0]."&".$arr[1]."&".$arr[2]."&".$arr[3];
		#print $key."\n";
		$junction{$key}->{$samples[$j]}=$arr[6];
	}
	close(IN);
}
#=cut;
print "Jun"."\t".join("\t", @samples)."\n";
#print "Jun"."\t"."Total_reads"."\n";
foreach my $key (keys %junction)
{
my @arr;  
for(my $i = 0; $i <= $#samples; $i++)
{
	 if(defined $junction{$key}->{$samples[$i]})
	 {
	    push(@arr, $junction{$key}->{$samples[$i]});
	 }
	 else
	 {
	    push(@arr, 0);
	 }
}
print $key."\t".join("\t", @arr)."\n";
}

sub round {
    my $val = shift;
    my $col = shift;
    my $r = 10 ** $col;
    my $a = ($val > 0) ? 0.5 : -0.5;
    return int($val * $r + $a) / $r;
}