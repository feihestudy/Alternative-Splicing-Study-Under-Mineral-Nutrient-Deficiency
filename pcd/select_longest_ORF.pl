use strict;
use Getopt::Long;
use List::Util qw(first max maxstr min minstr reduce shuffle sum);
use vars qw($fafile);
Getopt::Long::GetOptions(
    'fa=s'    => \$fafile,
);

open(IN, $fafile);
$/=">";
my %transid2score;
my %id2len;
<IN>;
while(<IN>)
{
  my @arr = split(/\n/, $_);
  my $id = shift(@arr);
  my @tmp = split(/\s+/, $id);
  $id=$tmp[0];
  my $seq = join("", @arr);
  $seq=~s/\>$//;
  $seq=~s/\*$//;
  my @tmp = split(/\:\:/, $id);
  my $id2 = $tmp[0]."::".$tmp[1];
  push(@{$transid2score{$id2}}, $id);
  $id2len{$id}=length($seq);
}
close(IN);


my %targetid;
foreach my $key (keys %transid2score)
{
	my @tmp = @{$transid2score{$key}};
	my $max = 0;
	my $maxindex = 0;
	for(my $i = 0; $i <= $#tmp; $i++)
	{
		my $len = $id2len{$tmp[$i]};
		if($len > $max)
		{
		    $max = $len;
			$maxindex = $i;
		}
	}
	$targetid{$tmp[$maxindex]} = 1;
}


open(IN, $fafile);
$/=">";
<IN>;
while(<IN>)
{
  my @arr = split(/\n/, $_);
  my $id = shift(@arr);
  my @tmp = split(/\s+/, $id);
  $id=$tmp[0];
  my $seq = join("", @arr);
  $seq=~s/\>$//;
  $seq=~s/\*$//;  
  if(defined $targetid{$id})
  {
      print ">".$id."\n".$seq."\n";
  }
}
close(IN);

