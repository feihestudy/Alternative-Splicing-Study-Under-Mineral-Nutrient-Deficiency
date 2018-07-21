use strict;
use Getopt::Long;
use vars qw($GTFfile $inputfile);
Getopt::Long::GetOptions(
    'input=s'    => \$inputfile,
	'GTF=s' => \$GTFfile,
);
my %transcript;
open(IN, "$GTFfile");
while(<IN>)
{
   chomp;
   my @arr = split(/\t/, $_);
   #print $arr[8]."\n";
   my @tmp = split(/;\s+/, $arr[8]);
   my $tid;
   my $gid;
   my $exon = $arr[2];
   next if($exon ne "exon");
   for(my $i = 0; $i <= $#tmp; $i++)
    {
       if($tmp[$i]=~m/^transcript\_id/)
       {
           $tid = $tmp[$i];
           $tid =~s/^transcript\_id//;
		   $tid=~s/^\s+//;
		   $tid=~s/\s+$//;
		   $tid=~s/\"//g;
		   $tid=~s/\;//g;
       }
       if($tmp[$i]=~m/^gene\_id/)
       {
           $gid = $tmp[$i];
           $gid =~s/^gene\_id//;
		   $gid=~s/^\s+//;
		   $gid=~s/\s+$//;
		   $gid=~s/\"//g;
		   $gid=~s/\;//g;
       }
	   
    }
	$transcript{$gid}=1;

}
close(IN);

#print join("\n", keys %transcript)."\n";

open(IN, $inputfile);
my $tt = <IN>;
$tt=~s/\n$//;
print $tt."\n";
while(<IN>)
{
    chomp;
	s/\s+$//;
	my @arr = split(/\t/, $_);
	#print $arr[0]."\n";
    if(defined $transcript{$arr[0]})
	{
	    print $_."\n";
	}
	else
	{
	     #print $_."\n";
	}
}
close(IN);

