use strict;
use Getopt::Long;
use vars qw($fq1file $fq2file $readlen);
Getopt::Long::GetOptions(
    'fq1=s' => \$fq1file,
    'fq2=s' => \$fq2file,
	'readlen=s' => \$readlen,
);
open(IN1, $fq1file);
open(IN2, $fq2file);
my $ot1fqfile = $fq1file.".".$readlen;
my $ot2fqfile = $fq2file.".".$readlen;
open(OUT1, ">$ot1fqfile");
open(OUT2, ">$ot2fqfile");
while(!eof(IN1))
{
    my $id1 = <IN1>;
    my $fa1 = <IN1>;
    $fa1 =~s/\n$//;
    my $plus1 = <IN1>;
    my $qv_ascii1 = <IN1>;
    $qv_ascii1 =~s/\n$//;
    my $n1 = length($fa1);
	my $subfa1 = $fa1;
	my $subqv_ascii1 = $qv_ascii1;	
    my $id2 = <IN2>;
    my $fa2 = <IN2>;
    $fa2 =~s/\n$//;
    my $plus2 = <IN2>;
    my $qv_ascii2 = <IN2>;
    $qv_ascii2 =~s/\n$//;
    my $n2 = length($fa2);
	my $subfa2 = $fa2;
	my $subqv_ascii2 = $qv_ascii2;	
	if($n1 >= $readlen and $n2 >= $readlen)
	{
	    $subfa1 = substr($fa1, 0, $readlen);
	    $subqv_ascii1 = substr($qv_ascii1, 0, $readlen);
	    $subfa2 = substr($fa2, 0, $readlen);
	    $subqv_ascii2 = substr($qv_ascii2, 0, $readlen);
        print OUT1 $id1.$subfa1."\n".$plus1.$subqv_ascii1."\n";
        print OUT2 $id2.$subfa2."\n".$plus2.$subqv_ascii2."\n";
    }
}
close(IN);
close(OUT1);
close(OUT2);

