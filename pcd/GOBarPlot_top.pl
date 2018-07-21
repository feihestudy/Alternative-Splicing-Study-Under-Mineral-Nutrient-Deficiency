use strict;
use Getopt::Long;
use vars qw($tissue);
Getopt::Long::GetOptions(
    'tissue=s'    => \$tissue,
);

sub log_base { my ($base, $value) = @_; return log($value)/log($base); } 


my $cutoff=log_base(10, 0.05)*(-1);
 
my %is_diff;

my %diff_per;

my %pvs;

my %targetids;
open(IN, $tissue."_DEG_all_GO_BP_simple.txt");
<IN>;
my $i = 0;
while(<IN>)
{
   chomp;
   s/\s+$//;
   my @arr = split(/\t/, $_);
   $i++;
   if($i > 10)
   {
       last;   
   }
   $targetids{$arr[1]}=1;
}
close(IN);

open(IN, $tissue."_DASG_all_GO_BP_simple.txt");
<IN>;
$i = 0;
while(<IN>)
{
   chomp;
   s/\s+$//;
   my @arr = split(/\t/, $_);
   $i++;
   if($i > 10)
   {
       last;   
   }
   $targetids{$arr[1]}=1;
}
close(IN);

open(IN, $tissue."_DEG_all_GO_BP.txt");
<IN>;
while(<IN>)
{
   chomp;
   s/\s+$//;
   my @arr = split(/\t/, $_);
   if(defined $targetids{$arr[1]})
   {
       $pvs{$arr[1]}->{"DEG"}=log_base(10, $arr[($#arr-1)])*(-1);
   }
}
close(IN);
open(IN, $tissue."_DASG_all_GO_BP.txt");
<IN>;
while(<IN>)
{
   chomp;
   s/\s+$//;
   my @arr = split(/\t/, $_);
   if(defined $targetids{$arr[1]})
   {
       $pvs{$arr[1]}->{"DASG"}=log_base(10, $arr[($#arr-1)])*(-1);
   }
}
close(IN);
open(OUT1, ">GO_PVs_tmp.xls");
print OUT1 "GO"."\t"."DEG"."\t"."DASG"."\n";


foreach my $key (keys %pvs)
{
   my @tmp = ();
   
   
   if(defined $pvs{$key}->{"DEG"})
   {
       push(@tmp, $pvs{$key}->{"DEG"});      
   }
   else
   {
       push(@tmp, 0);      
   
   }
   
    if(defined $pvs{$key}->{"DASG"})
   {
       push(@tmp, $pvs{$key}->{"DASG"});      
   }
   else
   {
       push(@tmp, 0);      
   
   }
  
   my $yes=0;
   
   for(my $i = 0; $i <= $#tmp; $i++)
   {
       $yes=1 if($tmp[$i] >= $cutoff);
   
   }
   if($yes eq 1)
   {
   print OUT1 $key."\t".join("\t", @tmp)."\n";
   }
}
close(OUT1);


