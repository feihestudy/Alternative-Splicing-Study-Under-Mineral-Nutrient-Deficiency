use strict;
use Getopt::Long;
use vars qw($genefeature $diffgene);
Getopt::Long::GetOptions(
    'genefeature=s' => \$genefeature,
	'diffgene=s'    => \$diffgene,
);
my %id2ex3;
my $file = $genefeature;
open(IN, $file);
my $tt2 = <IN>;
while(<IN>)
{
    chomp;
	s/\s+$//;
	my @arr = split(/\t/, $_);
    my $id = shift(@arr);
	$id=~s/\s+$//;
	$id=~s/^\s+//;
	$id2ex3{$id}=join("\t", @arr[0..($#arr-2)]);
}
close(IN);


my $file = $diffgene;
open(IN, $file);
print "geneid"."\t"."exon_num"."\t"."intron_num"."\t"."exon_len"."\t"."intron_len"."\n";
while(<IN>)
{
   chomp;
   s/\s+$//;
   my @arr = split(/\t/, $_);
   $arr[0]=~s/\s+$//;
   $arr[0]=~s/^\s+//;
   #print $arr[0]."\t"."---"."\n";
   if(defined $id2ex3{$arr[0]})
   {
       print $_."\t".$id2ex3{$arr[0]}."\n";
   }
   else
   {
      print $_."\t"."---"."\t"."---"."\t"."---"."\t"."---"."\n";	   
   }
}
close(IN);	
