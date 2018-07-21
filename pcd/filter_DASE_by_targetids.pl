use strict;
use Getopt::Long;
use vars qw($targetids $DASEtable);
Getopt::Long::GetOptions(
	'DASEtable=s' =>\$DASEtable,
	'targetids=s' =>\$targetids,	
);

my %target;
my @syms;
my @arr = split(/\,/, $targetids);
for(my $i = 0; $i <= $#arr; $i++)
{
    my @tmp2 = split(/\-/, $arr[$i]);
	$target{$tmp2[0]}->{"gene"}=$tmp2[1];
	$target{$tmp2[0]}->{"symbol"}=$tmp2[2];	
	push(@syms, $tmp2[2]);
}

my %ASE;
open(IN, $DASEtable);
my $tt = <IN>;
$tt=~s/\n$//;
$tt=~s/\s+$//;
#print $tt;
my @tt_tmp = split(/\t/, $tt);
shift(@tt_tmp);
while(<IN>)
{
    chomp;
	s/\s+$//;
	my @arr = split(/\t/, $_);
	my $aseid=shift(@arr);
	my @arr2 = split(/\&/,$aseid);
	if(defined $target{$arr2[1]})
	{
		my $sym = $target{$arr2[1]}->{"symbol"};
		$arr2[1]=$target{$arr2[1]}->{"gene"};
		for(my $j = 0; $j <= $#arr; $j++)
		{
		    push(@{$ASE{$sym}->{$tt_tmp[$j]}}, $arr[$j]);
		}
	}
}
close(IN);

print "Gene"."\t".join("\t", @tt_tmp)."\n";
for(my $k = 0; $k <= $#syms; $k++)
{
    my $key = $syms[$k];
    my @value;
	for(my $i = 0; $i <= $#tt_tmp; $i++)
	{
	    if(defined $ASE{$key}->{$tt_tmp[$i]})
		{
		    my @zeroone=@{$ASE{$key}->{$tt_tmp[$i]}};
			my $hit = 0;
			for(my $l=0; $l <= $#zeroone; $l++)
			{
			   $hit=1 if($zeroone[$l] != 0);
			}
		    push(@value, $hit);
		}
		else
		{
		    push(@value, 0);		
		}
	}
	print $key."\t".join("\t", @value)."\n";
}