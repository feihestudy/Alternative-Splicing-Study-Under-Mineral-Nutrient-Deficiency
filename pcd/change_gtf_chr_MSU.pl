####
use strict;
my %chrs;
$chrs{"Chr1"} = "chr01";
$chrs{"Chr10"} = "chr10";
$chrs{"Chr11"} = "chr11";
$chrs{"Chr12"} = "chr12";
$chrs{"Chr2"} = "chr02";
$chrs{"Chr3"} = "chr03";
$chrs{"Chr4"} = "chr04";
$chrs{"Chr5"} = "chr05";
$chrs{"Chr6"} = "chr06";
$chrs{"Chr7"} = "chr07";
$chrs{"Chr8"} = "chr08";
$chrs{"Chr9"} = "chr09";
open(IN, "all.gff3");
while(<IN>)
{
   chomp;
   my @arr = split(/\t/, $_);
   if(defined $chrs{$arr[0]})
   {
       $arr[0] = $chrs{$arr[0]};
	   print join("\t", @arr)."\n";
   }

}
close(IN);

