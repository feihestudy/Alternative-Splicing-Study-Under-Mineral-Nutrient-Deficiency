use strict;
use Getopt::Long;
use List::Util qw(first max maxstr min minstr reduce shuffle sum);
use vars qw($ASE $dir);
Getopt::Long::GetOptions(
	'dir=s' => \$dir,
	'ASE=s' => \$ASE,
);

my %ASE2trans;
open(IN, $dir."/RI2transcript_major.xls");
<IN>;
while(<IN>)
{
    chomp;
	s/\s+$//;
	my @arr = split(/\t/, $_);
	$ASE2trans{$arr[0]}->{"include"}=$arr[1];
	$ASE2trans{$arr[0]}->{"exclude"}=$arr[2];	
}
close(IN);
open(IN, $dir."/SE2transcript_major.xls");
<IN>;
while(<IN>)
{
    chomp;
	s/\s+$//;
	my @arr = split(/\t/, $_);
	$ASE2trans{$arr[0]}->{"include"}=$arr[1];
	$ASE2trans{$arr[0]}->{"exclude"}=$arr[2];	
}
close(IN);
open(IN, $dir."/A3SS2transcript_major.xls");
<IN>;
while(<IN>)
{
    chomp;
	s/\s+$//;
	my @arr = split(/\t/, $_);
	$ASE2trans{$arr[0]}->{"include"}=$arr[1];
	$ASE2trans{$arr[0]}->{"exclude"}=$arr[2];	
}
close(IN);
open(IN, $dir."/A5SS2transcript_major.xls");
<IN>;
while(<IN>)
{
    chomp;
	s/\s+$//;
	my @arr = split(/\t/, $_);
	$ASE2trans{$arr[0]}->{"include"}=$arr[1];
	$ASE2trans{$arr[0]}->{"exclude"}=$arr[2];	
}
close(IN);
open(IN, $dir."/MXE2transcript_major.xls");
<IN>;
while(<IN>)
{
    chomp;
	s/\s+$//;
	my @arr = split(/\t/, $_);
	$ASE2trans{$arr[0]}->{"include"}=$arr[1];
	$ASE2trans{$arr[0]}->{"exclude"}=$arr[2];	
}
close(IN);

foreach my $key (keys %ASE2trans)
{
   if($key eq $ASE)
   {
      print $key."\t".$ASE2trans{$key}->{"include"}."\t".$ASE2trans{$key}->{"exclude"}."\n";
   }
}




