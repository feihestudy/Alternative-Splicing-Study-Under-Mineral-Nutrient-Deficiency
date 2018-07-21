use strict;
use List::Util qw(first max maxstr min minstr reduce shuffle sum);
use List::Uniq ':all';
use strict;
use Getopt::Long;
use vars qw($input1 $input2 $input3);
Getopt::Long::GetOptions(
    'input1=s' => \$input1,
	'input2=s' => \$input2,
	'input3=s' => \$input3,	
);

open(IN, $input1);
my $tt1 = <IN>;
close(IN);
$tt1=~s/\n$//;
$tt1=~s/\s+$//;
my @tmp1 = split(/\,/, $tt1);


open(IN, $input2);
my $tt2 = <IN>;
close(IN);
$tt2=~s/\n$//;
$tt2=~s/\s+$//;
my @tmp2 = split(/\,/, $tt2);

open(IN, $input3);
my $tt3 = <IN>;
close(IN);
$tt3=~s/\n$//;
$tt3=~s/\s+$//;
my @tmp3 = split(/\,/, $tt3);



my $sum1 = $tmp1[0]+$tmp2[0]+$tmp3[0];
my $sum2 = $tmp1[1]+$tmp2[1]+$tmp3[1];
my $psi = $tmp1[2]+$tmp2[2]+$tmp3[2];
$psi = $psi/3;

print $sum1.",".$sum2.",".$psi."\n";






