use strict;
use List::Uniq ':all';
use Getopt::Long;
use List::Util qw(first max maxstr min minstr reduce shuffle sum);
use Getopt::Long;
use Statistics::R;
use vars qw($degenfile $annotationfile);
Getopt::Long::GetOptions(
    'degene=s'    =>\$degenfile,
	'annotation=s' => \$annotationfile,
);


my %exist;
open(IN, $degenfile);
while(<IN>)
{
    chomp;
	s/\s+$//;
	s/\"//g;
	s/\'//g;	
	$exist{$_}=1;
}
close(IN);

my %exist2;

my %rice2ann;
open(IN,$annotationfile);
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
   next if ($arr[2] eq "NULL");
   my @gids = split(/\/\//, $arr[2]);
   next if(not defined $exist{$arr[0]});
   for(my $i = 0; $i <= $#gids; $i++)
   {
		#print $gids[$i]."\n";
		$exist2{$gids[$i]}=1;
   }
  
}
close(IN);	


print join("\n",keys %exist2)."\n";

