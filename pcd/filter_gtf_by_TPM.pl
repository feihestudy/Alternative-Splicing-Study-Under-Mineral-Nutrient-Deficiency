use strict;
use Getopt::Long;
use List::Uniq ':all';
use List::Util qw(first max maxstr min minstr reduce shuffle sum);
use vars qw($gtffile $tpmcutoff $ncutoff $tpmfile);
Getopt::Long::GetOptions(
    'GTF=s' => \$gtffile,
	'TPM=s' => \$tpmcutoff,
	'N=s' => \$ncutoff,
	'TPMfile=s' => \$tpmfile,
);
my %keytid;
open(IN, $tpmfile);
<IN>;
while(<IN>)
{
   chomp;
   s/\s+$//;
   my @arr = split(/\t/, $_);
   my $tid = shift(@arr);
   my $sum = 0;
   for(my $i = 0; $i <= $#arr; $i++)
   {
       $sum++ if($arr[$i] >= $tpmcutoff);
   }
   if($sum >= $ncutoff)
   {
       $keytid{$tid}=1;
   }
}
close(IN);


my %uniq;
open(IN, $gtffile);
while(<IN>)
{
	chomp;
	my @arr = split(/\t/, $_);
	my $chr = $arr[0];
	my @tmp = split(/;\s/, $arr[8]);
    my $tid;
	my $gid;
	my $fpkm;
    for(my $i = 0; $i <= $#tmp; $i++)
    {
       if($tmp[$i]=~m/^transcript\_id/)
       {
           $tid= $tmp[$i];		   
           $tid=~s/^transcript\_id//;
		   $tid=~s/^\s+//;
		   $tid=~s/\"//g;
		   $tid=~s/\;$//;
       }
    }
	if(defined $keytid{$tid})
	{
	   print $_."\n";
	}
}
close(IN);



