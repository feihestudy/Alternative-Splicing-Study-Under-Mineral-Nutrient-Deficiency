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
	$id2ex3{$id}=join("\t", $arr[($#arr-1)]);
}
close(IN);

#print scalar(keys %id2ex3)."\n";



my $file = $diffgene;
open(IN, $file);
my $tt = <IN>;
$tt=~s/\n$//;
$tt=~s/\s+$//;
print "geneid"."\t"."GC"."\n";
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
   print $_."\t"."---"."\n";	   
   }
}
close(IN);	
