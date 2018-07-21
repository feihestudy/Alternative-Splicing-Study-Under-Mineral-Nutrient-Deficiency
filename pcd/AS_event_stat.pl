use strict;
use Getopt::Long;
use vars qw($ASE);
Getopt::Long::GetOptions(
    'ASE=s'    => \$ASE,
);
use List::Uniq ':all';
#use List::Compare;
my $ri=0;
my $se=0;
my $a3ss=0;
my $a5ss=0;
my $mxe=0;
open(IN, $ASE);
#<IN>;
while(<IN>)
{
    chomp;
	s/\s+$//;
	my @arr = split(/\t/, $_);
	my @arr2 = split(/\&/, $arr[0]);
	$ri++ if($arr2[0] eq "RI");
	$a3ss++ if($arr2[0] eq "A3SS");
	$se++ if($arr2[0] eq "SE");
	$a5ss++ if($arr2[0] eq "A5SS");
	$mxe++ if($arr2[0] eq "MXE");
}
close(IN);
my $tt = ($ri+$a3ss+$se+$a5ss+$mxe);
print "Event Type"."\t"."Percentage(%)"."\n";
$ri=$ri/$tt;
$ri=round($ri*100, 2)."%";
print "RI\t".$ri."\n";
$a3ss=$a3ss/$tt;
$a3ss=round($a3ss*100, 2)."%";
print "A3SS\t".$a3ss."\n";
$a5ss=$a5ss/$tt;
$a5ss=round($a5ss*100, 2)."%";
print "A5SS\t".$a5ss."\n";
$se=$se/$tt;
$se=round($se*100, 2)."%";
print "SE\t".$se."\n";
$mxe=$mxe/$tt;
$mxe=round($mxe*100, 2)."%";
print "MXE\t".$mxe."\n";
#print "Total\t".($ri+$a3ss+$se+$a5ss+$mxe)."\n";

sub round 
{
    my $val = shift;
    my $col = shift;
    my $r = 10 ** $col;
    my $a = ($val > 0) ? 0.5 : -0.5;
    return int($val * $r + $a) / $r;
}

	


