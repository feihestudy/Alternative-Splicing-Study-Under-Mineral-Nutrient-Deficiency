use strict;
use Getopt::Long;
use List::Util qw(first max maxstr min minstr reduce shuffle sum);
use vars qw($tmapfile $gtffile);
Getopt::Long::GetOptions(
	'tmap=s'    => \$tmapfile,
	'gtf=s' => \$gtffile,
);
$/="\n";
open(IN, $tmapfile);
my %leftid;
while(<IN>)
{
	#print $_;
	chomp;
	my @arr = split(/\t/, $_);
	my $type = $arr[2];
    if($type eq "j" or $type eq "o" or $type eq "u")
	{
       $leftid{$arr[4]} = 1;
	}
}
close(IN);

open(IN, $gtffile);
my %transcript;
while(<IN>)
{
    chomp;
	my @arr = split(/\t/, $_);
	next if($arr[6] eq ".");
	my @tmp = split(/;\s+/, $arr[8]);
    my $tid;
	my $gid;
	my $fpkm;
    for(my $i = 0; $i <= $#tmp; $i++)
    {
	   $tmp[$i]=~s/^\s+//;
       if($tmp[$i]=~m/^transcript\_id/)
       {
           $tid= $tmp[$i];
		   $tid=~s/\;$//;
		   $tid=~s/^\s+//;
           $tid=~s/^transcript\_id//;
		   $tid=~s/^\s+//;
		   $tid=~s/\s+$//;
		   $tid=~s/\"//g;
       }
       if($tmp[$i]=~m/^gene\_id/)
       {
           $gid= $tmp[$i];
		   $gid=~s/\;$//;
		   $gid=~s/^\s+//;		   
           $gid=~s/^gene\_id//;
		   $gid=~s/^\s+//;
		   $gid=~s/\s+$//;
		   $gid=~s/\"//g;
       }   
    }
	
	if(defined $leftid{$tid})
	{
		print join("\t", @arr[0..7])."\t"."gene_id \"".$gid."\";"." "."transcript_id \"".$tid."\";"."\n";
	}
}
close(IN);
