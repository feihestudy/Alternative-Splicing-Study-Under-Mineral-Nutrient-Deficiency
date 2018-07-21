use strict;
use Getopt::Long;
use List::Util qw(first max maxstr min minstr reduce shuffle sum);
use vars qw($genefpkm $diffgene);
Getopt::Long::GetOptions(
    'genefpkm=s' => \$genefpkm,
	'diffgene=s'    => \$diffgene,
);
my %id2ex3;
my $file = $genefpkm;
open(IN, $file);
my $tt2 = <IN>;
$tt2=~s/\n$//;
$tt2=~s/\s+$//;
my @tmp = split(/\t/, $tt2);
shift(@tmp);
$tt2=join("\t", @tmp);
while(<IN>)
{
    chomp;
	s/\s+$//;
	my @arr = split(/\t/, $_);
    my $id = shift(@arr);
	$id=~s/\s+$//;
	$id=~s/^\s+//;
	for(my $i = 0; $i <= $#arr; $i++)
	{
	    $arr[$i] = log($arr[$i]+1)/log(2)
	}
	$id2ex3{$id}=sum(@arr)/scalar(@arr);
}
close(IN);

#print scalar(keys %id2ex3)."\n";



my $file = $diffgene;
open(IN, $file);
my $tt = <IN>;
$tt=~s/\n$//;
$tt=~s/\s+$//;
print "geneid"."\t"."expression"."\n";
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
