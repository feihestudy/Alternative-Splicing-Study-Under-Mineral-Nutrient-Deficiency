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
	$transcript{$tid}=1;

}
close(IN);

$/ = ">";
open(IN, $inputfile);
<IN>;
while(<IN>)
{
    my @arr = split(/\n/, $_);
    my $id = shift(@arr);
	my $seq = join("", @arr);
	$seq=~s/\>$//;
	$seq=~s/\*//g;
	if(defined $transcript{$id})
	{
	print ">".$id."\n".$seq."\n";
	}
}
close(IN);

