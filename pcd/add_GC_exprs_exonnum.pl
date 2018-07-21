use strict;
use Getopt::Long;
use List::Util qw(first max maxstr min minstr reduce shuffle sum);
use vars qw($genefpkm $genefeature $ASEnum);
Getopt::Long::GetOptions(
    'genefpkm=s' => \$genefpkm,
	'genefeature=s' => \$genefeature,
	'ASEnum=s' =>\$ASEnum,
);
my %exprs;
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
	#@arr=@arr[0..72];
	for(my $i = 0; $i <= $#arr; $i++)
	{
	    $arr[$i] = log($arr[$i]+1)/log(2)
	}
	$exprs{$id}=sum(@arr)/scalar(@arr);
}
close(IN);

my %gfs;
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
	$gfs{$id}=join("\t", @arr[0..($#arr-1)]);
}
close(IN);


my %gene2ASE;
my $file = $ASEnum;
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
	$gene2ASE{$id}=$arr[0];
}
close(IN);



my %targetid;
open(IN, "high_AS_genes.txt");
while(<IN>)
{
   chomp;
   s/\s+$//;
   $targetid{$_}=1;

}
close(IN);


open(IN, "low_AS_genes.txt");
while(<IN>)
{
   chomp;
   s/\s+$//;
   $targetid{$_}=1;

}
close(IN);


open(IN, "no_AS_genes.txt");
while(<IN>)
{
   chomp;
   s/\s+$//;
   $targetid{$_}=1;

}
close(IN);

print "gene_id\tFPKM\texon_number\tintron_number\texon_length\tintron_length\tGC_percentage\tASEnum"."\n";
foreach my $key (keys %targetid)
{
   print $key."\t";
   if(defined $exprs{$key})
   {
       print $exprs{$key}."\t";
   }
   else
   {
	   print "---"."\t";	   
   }
    if(defined $gfs{$key})
   {
       print $gfs{$key}."\t";
   }
   else
   {
	   print "---"."\t"."---"."\t"."---"."\t"."---"."\t"."----"."\t";	   
   }  
    if(defined $gene2ASE{$key})
   {
       print $gene2ASE{$key}."\n";
   }
   else
   {
	   print "0"."\n";	   
   } 
}
