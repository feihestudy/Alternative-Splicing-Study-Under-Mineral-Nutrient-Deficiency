use strict;
use Getopt::Long;
use vars qw($tmapfile);
Getopt::Long::GetOptions(
	'tmap=s' => \$tmapfile,
);
use List::Util qw(first max maxstr min minstr reduce shuffle sum);
$/="\n";
my %classcode;
my $otfile = "class_code_stat.xls";
open(OUT, ">$otfile"); 
print OUT "Class code"."\t"."Descrition of class code"."\t"."Number of transcripts within each class code"."\n";
open(IN, $tmapfile);
<IN>;
while(<IN>)
{
	#print $_;
	chomp;
	my @arr = split(/\t/, $_);
	my $type = $arr[2];
	$classcode{$type}++;
}
close(IN);
if(defined $classcode{"c"})
{
    print OUT "c"."\t"."Contained"."\t".$classcode{"c"}."\n";
}
else
{
    print OUT "c"."\t"."Contained"."\t"."0"."\n";
}
if(defined $classcode{"="})
{
    print OUT "="."\t"."Complete match of intron chain"."\t".$classcode{"="}."\n";
}
else
{
    print OUT "="."\t"."Complete match of intron chain"."\t"."0"."\n";
}
if(defined $classcode{"u"})
{
    print OUT "u"."\t"."Unknown, intergenic transcript"."\t".$classcode{"u"}."\n";
}
else
{
    print OUT "u"."\t"."Unknown, intergenic transcript"."\t"."0"."\n";
}
if(defined $classcode{"j"})
{
    print OUT "j"."\t"."Potentially novel isoform (fragment): at least one splice junction is shared with a reference transcript"."\t".$classcode{"j"}."\n";
}
else
{
    print OUT "j"."\t"."Potentially novel isoform (fragment): at least one splice junction is shared with a reference transcript"."\t"."0"."\n";
}
if(defined $classcode{"i"})
{
    print OUT "i"."\t"."A transfrag falling entirely within a reference intron"."\t".$classcode{"i"}."\n";
}
else
{
    print OUT "i"."\t"."A transfrag falling entirely within a reference intron"."\t"."0"."\n";
}
if(defined $classcode{"o"})
{
    print OUT "o"."\t"."Generic exonic overlap with a reference transcript"."\t".$classcode{"o"}."\n";
}
else
{
    print OUT "o"."\t"."Generic exonic overlap with a reference transcript"."\t"."0"."\n";
}
if(defined $classcode{"e"})
{
    print OUT "e"."\t"."Single exon transfrag overlapping a reference exon and at least 10 bp of a reference intron, indicating a possible pre-mRNA fragment."."\t".$classcode{"e"}."\n";
}
else
{
    print OUT "e"."\t"."Single exon transfrag overlapping a reference exon and at least 10 bp of a reference intron, indicating a possible pre-mRNA fragment"."\t"."0"."\n";
}
if(defined $classcode{"x"})
{
    print OUT "x"."\t"."Exonic overlap with reference on the opposite strand"."\t".$classcode{"x"}."\n";
}
else
{
    print OUT "x"."\t"."Exonic overlap with reference on the opposite strand"."\t"."0"."\n";
}
if(defined $classcode{"s"})
{
    print OUT "s"."\t"."An intron of the transfrag overlaps a reference intron on the opposite strand (likely due to read mapping errors)"."\t".$classcode{"s"}."\n";
}
else
{
    print OUT "s"."\t"."An intron of the transfrag overlaps a reference intron on the opposite strand (likely due to read mapping errors)"."\t"."0"."\n";
}
if(defined $classcode{"p"})
{
    print OUT "p"."\t"."Possible polymerase run-on fragment (within 2Kbases of a reference transcript)"."\t".$classcode{"p"}."\n";
}
else
{
    print OUT "p"."\t"."Possible polymerase run-on fragment (within 2Kbases of a reference transcript)"."\t"."0"."\n";
}
close(OUT);
