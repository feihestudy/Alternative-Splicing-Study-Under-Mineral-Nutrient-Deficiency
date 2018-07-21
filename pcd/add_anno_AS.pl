use strict;
use List::Uniq ':all';
use List::Util qw(first max maxstr min minstr reduce shuffle sum);
use Getopt::Long;
use Statistics::R;
use vars qw($degenfile $annotationfile);
Getopt::Long::GetOptions(
    'degene=s'    =>\$degenfile,
	'annotation=s' =>\$annotationfile,
);
my %rice2ann;
my $file = $annotationfile;
open(IN, $file);
my $tt1 = <IN>;
$tt1=~s/\n$//;
$tt1=~s/\s+$//;
my @tttmp = split(/\t/, $tt1);
$tt1=join("\t", @tttmp[1..$#tttmp]);
while(<IN>)
{
   chomp;
   s/\s+$//;
   my @arr = split(/\t/, $_);
   my $id = shift(@arr);
   $rice2ann{$id}= join("\t", @arr);
}
close(IN);	

open(IN, $degenfile);
my $tt = <IN>;
$tt=~s/\n$//;
$tt=~s/\s+$//;
print $tt."\t".$tt1."\n";
while(<IN>)
{
    chomp;
	s/\s+$//;
	my @arr = split(/\t/, $_);
	my @arr2 = split(/\&/, $arr[0]);
	if(defined $rice2ann{$arr2[1]})
	{
	     print $_."\t".$rice2ann{$arr2[1]}."\n";
	}
}
close(IN);






