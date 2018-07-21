use strict;
use Getopt::Long;
use vars qw($samplepairfile);
Getopt::Long::GetOptions(
    'samplepair=s' => \$samplepairfile,
);
my @samplepairs = split(/\,/, $samplepairfile);
my %gene2exist;
for(my $i = 0; $i <= $#samplepairs; $i++)
{
     my $file = $samplepairs[$i]."_gene_exp_significant.xls";
	 #print $file."\n";
	open(IN, $file);
	<IN>;
	while(<IN>)
	{
	    chomp;
		s/\s+$//;
		s/\"//g;
		my @arr = split(/\t/, $_);
		#print $arr[0]."\n";
		my $fc = $arr[($#arr-2)];
		if($fc >= 1)
		{
			$gene2exist{$arr[0]}->{$samplepairs[$i]}=1;
		}
		if($fc <= -1)
		{
			$gene2exist{$arr[0]}->{$samplepairs[$i]}=-1;
		}		
	}
	close(IN);	 
}
#print join("\n", keys %gene2exist)."\n";
print "Gene_id"."\t".join("\t", @samplepairs)."\n";
foreach my $key (keys %gene2exist)
{
    my @values;  
    for(my $i = 0; $i <= $#samplepairs; $i++)
	{
	     if(defined $gene2exist{$key}->{$samplepairs[$i]})
		 {
		    push(@values, $gene2exist{$key}->{$samplepairs[$i]});
		 }
		 else
		 {
		    push(@values, 0);		 
		 }
	}
	print $key."\t".join("\t", @values)."\n";
}



