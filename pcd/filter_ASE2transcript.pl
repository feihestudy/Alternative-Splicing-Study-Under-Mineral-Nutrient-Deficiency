use strict;
use Getopt::Long;
use List::Util qw(first max maxstr min minstr reduce shuffle sum);
use vars qw($CDSgtf $exongtf);
Getopt::Long::GetOptions(
    'CDSgtf=s' =>\$CDSgtf,
	'exongtf=s' =>\$exongtf,
);
my $intr;
my $extr;
open(IN, "ASE2transcript_major.txt");
my $line = <IN>;
$line=~s/\n$//;
$line=~s/\s+$//;
my @arr = split(/\t/, $line);
$intr=$arr[1];
$extr=$arr[2];


my $intr_isok = 0;
my $extr_isok = 0;
open(IN, $CDSgtf);
while(<IN>)
{
	chomp;
	my @arr = split(/\t/, $_);
	my $chr = $arr[0];
	my @tmp = split(/;\s/, $arr[8]);
    my $tid;
	my $gid;
    for(my $i = 0; $i <= $#tmp; $i++)
    {
       if($tmp[$i]=~m/^gene\_id/)
       {
           $gid = $tmp[$i];
           $gid =~s/^gene\_id//;
		   $gid=~s/^\s+//;
		   $gid=~s/\s+$//;
		   $gid=~s/\"//g;
		   $gid=~s/\;$//;
		   
       }
       if($tmp[$i]=~m/^transcript\_id/)
       {
           $tid = $tmp[$i];
           $tid =~s/^transcript\_id//;
		   $tid=~s/^\s+//;
		   $tid=~s/\;$//;
		   $tid=~s/\s+$//;
		   $tid=~s/\"//g;
       }
    }
	if($tid eq $intr)
	{
	$intr_isok=1;	
	}
	if($tid eq $extr)
	{
	$extr_isok=1;
	}
	
	
	
	
}
close(IN);


open(INCLUDE, ">include.gtf");
open(EXCLUDE, ">exclude.gtf");

if($intr_isok  == 1)
{
open(IN, $CDSgtf);
while(<IN>)
{
	chomp;
	my @arr = split(/\t/, $_);
	my $chr = $arr[0];
	my @tmp = split(/;\s/, $arr[8]);
    my $tid;
	my $gid;
    for(my $i = 0; $i <= $#tmp; $i++)
    {
       if($tmp[$i]=~m/^gene\_id/)
       {
           $gid = $tmp[$i];
           $gid =~s/^gene\_id//;
		   $gid=~s/^\s+//;
		   $gid=~s/\s+$//;
		   $gid=~s/\"//g;
		   $gid=~s/\;$//;
		   
       }
       if($tmp[$i]=~m/^transcript\_id/)
       {
           $tid = $tmp[$i];
           $tid =~s/^transcript\_id//;
		   $tid=~s/^\s+//;
		   $tid=~s/\;$//;
		   $tid=~s/\s+$//;
		   $tid=~s/\"//g;
       }
    }
	if($tid eq $intr)
	{
		print INCLUDE join("\t", @arr[0..7])."\t"."gene_id "."\"".$gid."\"; "."transcript_id "."\"".$tid."\";"."\n";
	}
}
close(IN);
close(INCLUDE);
}
else
{


open(IN, $exongtf);
while(<IN>)
{
	chomp;
	my @arr = split(/\t/, $_);
	my $chr = $arr[0];
	my @tmp = split(/;\s/, $arr[8]);
    my $tid;
	my $gid;
    for(my $i = 0; $i <= $#tmp; $i++)
    {
       if($tmp[$i]=~m/^gene\_id/)
       {
           $gid = $tmp[$i];
           $gid =~s/^gene\_id//;
		   $gid=~s/^\s+//;
		   $gid=~s/\s+$//;
		   $gid=~s/\"//g;
		   $gid=~s/\;$//;
		   
       }
       if($tmp[$i]=~m/^transcript\_id/)
       {
           $tid = $tmp[$i];
           $tid =~s/^transcript\_id//;
		   $tid=~s/^\s+//;
		   $tid=~s/\;$//;
		   $tid=~s/\s+$//;
		   $tid=~s/\"//g;
       }
    }
	if($tid eq $intr)
	{
		print INCLUDE join("\t", @arr[0..7])."\t"."gene_id "."\"".$gid."\"; "."transcript_id "."\"".$tid."\";"."\n";
	}
}
close(IN);
close(INCLUDE);

}




if($extr_isok  == 1)
{
open(IN, $CDSgtf);
while(<IN>)
{
	chomp;
	my @arr = split(/\t/, $_);
	my $chr = $arr[0];
	my @tmp = split(/;\s/, $arr[8]);
    my $tid;
	my $gid;
    for(my $i = 0; $i <= $#tmp; $i++)
    {
       if($tmp[$i]=~m/^gene\_id/)
       {
           $gid = $tmp[$i];
           $gid =~s/^gene\_id//;
		   $gid=~s/^\s+//;
		   $gid=~s/\s+$//;
		   $gid=~s/\"//g;
		   $gid=~s/\;$//;
		   
       }
       if($tmp[$i]=~m/^transcript\_id/)
       {
           $tid = $tmp[$i];
           $tid =~s/^transcript\_id//;
		   $tid=~s/^\s+//;
		   $tid=~s/\;$//;
		   $tid=~s/\s+$//;
		   $tid=~s/\"//g;
       }
    }
	if($tid eq $extr)
	{
		print EXCLUDE join("\t", @arr[0..7])."\t"."gene_id "."\"".$gid."\"; "."transcript_id "."\"".$tid."\";"."\n";
	}
}
close(IN);
close(EXCLUDE);
}
else
{


open(IN, $exongtf);
while(<IN>)
{
	chomp;
	my @arr = split(/\t/, $_);
	my $chr = $arr[0];
	my @tmp = split(/;\s/, $arr[8]);
    my $tid;
	my $gid;
    for(my $i = 0; $i <= $#tmp; $i++)
    {
       if($tmp[$i]=~m/^gene\_id/)
       {
           $gid = $tmp[$i];
           $gid =~s/^gene\_id//;
		   $gid=~s/^\s+//;
		   $gid=~s/\s+$//;
		   $gid=~s/\"//g;
		   $gid=~s/\;$//;
		   
       }
       if($tmp[$i]=~m/^transcript\_id/)
       {
           $tid = $tmp[$i];
           $tid =~s/^transcript\_id//;
		   $tid=~s/^\s+//;
		   $tid=~s/\;$//;
		   $tid=~s/\s+$//;
		   $tid=~s/\"//g;
       }
    }
	if($tid eq $extr)
	{
		print EXCLUDE join("\t", @arr[0..7])."\t"."gene_id "."\"".$gid."\"; "."transcript_id "."\"".$tid."\";"."\n";
	}
}
close(IN);
close(EXCLUDE);

}









