use strict;
use Getopt::Long;
use vars qw($ASEfile);
Getopt::Long::GetOptions(
    'ASE=s'    => \$ASEfile,
);
use List::Uniq ':all';
use List::Compare;
my %database_se_list;
open(IN, "/home/database/AS.RI.txt");
<IN>;
while(<IN>)
{
    chomp;
	s/\s+$//;
	s/\"//g;
	my @arr = split(/\t/, $_);
	my $key = $arr[1]."&".$arr[3]."&".$arr[4]."&".$arr[5]."&".$arr[6]."&".$arr[7]."&".$arr[8]."&".$arr[9]."&".$arr[10];
	$database_se_list{$key}=1;
}
close(IN);
my @Control1_se_list;
open(IN, $ASEfile);
<IN>;
while(<IN>)
{
    chomp;
	s/\s+$//;
	my @tmp = split(/\t/, $_);
	my @arr = split(/\&/, $tmp[0]);
	next if($arr[0] ne "RI");
	#push(@Control1_se_list, join("&", @arr[1..$#arr]));
	push(@Control1_se_list, join("&", @arr[0..$#arr]));	
}
close(IN);
#print join("\n", @Control1_se_list)."\n";
my %known;
for(my $i = 0; $i <= $#Control1_se_list; $i++)
{
	#print $Control1_se_list[$i]."\n";
	my @elements = split(/\&/, $Control1_se_list[$i]);
	my $target_gid = $elements[1];
	my $target_chr = $elements[2];
	my $target_strand =	$elements[3];	
	my $target_up_exst=$elements[6];
	$target_up_exst++;
	my $target_up_exen=$elements[7];	
	my $target_down_exst=$elements[8];
	$target_down_exst++;
	my $target_down_exen=$elements[9];	
	my $target_exst=$target_up_exen+1;
	my $target_exen=$target_down_exst-1;
	foreach my $key (keys %database_se_list)
	{
		#print $database_se_list[$j]."\n";
		my @elements2 = split(/\&/, $key);
		my $target2_gid = $elements2[0];
		my $target2_chr = $elements2[1];
		my $target2_strand = $elements2[2];	
		my $target2_up_exst=$elements2[5];
		$target2_up_exst++;
		my $target2_up_exen=$elements2[6];	
		my $target2_down_exst=$elements2[7];
		$target2_down_exst++;
		my $target2_down_exen=$elements2[8];	
		my $target2_exst=$target2_up_exen+1;
		my $target2_exen=$target2_down_exst-1;
		next if($target_chr ne $target2_chr);
		next if($target_strand ne $target2_strand);
		#next if($target_gid ne $target2_gid);
		my $left_overlap = overlap($target_up_exst, $target_up_exen, $target2_up_exst, $target2_up_exen);
		my $right_overlap = overlap($target_down_exst, $target_down_exen, $target2_down_exst, $target2_down_exen);
		if($left_overlap != 12 and $right_overlap != 12 and $target_exst eq $target2_exst and $target_exen eq $target2_exen)
		{
		      #print $Control1_se_list[$i]."\n";
			  $known{$Control1_se_list[$i]}=1;
			  last;
		}		
	}
}
#print scalar(keys %know)."\n";
my $otfile1 = $ASEfile;
$otfile1="RI_known.txt";
my $otfile2 = $ASEfile;
$otfile2="RI_novel.txt";
open(OUT1, ">$otfile1");
open(OUT2, ">$otfile2");
for(my $i = 0; $i<=$#Control1_se_list; $i++)
{
     if(defined $known{$Control1_se_list[$i]})
	 {
		print OUT1 $Control1_se_list[$i]."\n";
	 }
	 else
	 {
		print OUT2 $Control1_se_list[$i]."\n";
	 }	 
}
close(OUT1);
close(OUT2);





#### A3SS
my %database_se_list;
open(IN, "/home/database/AS.A3SS.txt");
<IN>;
while(<IN>)
{
    chomp;
	s/\s+$//;
	s/\"//g;
	my @arr = split(/\t/, $_);
	my $key = $arr[1]."&".$arr[3]."&".$arr[4]."&".$arr[5]."&".$arr[6]."&".$arr[7]."&".$arr[8]."&".$arr[9]."&".$arr[10];
	$database_se_list{$key}=1;
}
close(IN);
my @Control1_se_list;
open(IN, $ASEfile);
<IN>;
while(<IN>)
{
    chomp;
	s/\s+$//;
	my @tmp = split(/\t/, $_);
	my @arr = split(/\&/, $tmp[0]);
	next if($arr[0] ne "A3SS");
	#push(@Control1_se_list, join("&", @arr[1..$#arr]));
	push(@Control1_se_list, join("&", @arr[0..$#arr]));	
}
close(IN);
#print join("\n", @Control1_se_list)."\n";
my %known;
for(my $i = 0; $i <= $#Control1_se_list; $i++)
{
	#print $Control1_se_list[$i]."\n";
	my @elements = split(/\&/, $Control1_se_list[$i]);
	my $target_gid = $elements[1];
	my $target_chr = $elements[2];
	my $target_strand =	$elements[3];
    if($target_strand eq "-")
	{	
		my $target_up_exst=$elements[6];
		$target_up_exst++;
		my $target_up_exen=$elements[7];		
		my $target_down_exst=$elements[8];
		$target_down_exst++;
		my $target_down_exen=$elements[9];	
		my $target_exst=$elements[7]+1;
		my $target_exen=$elements[5];
		foreach my $key (keys %database_se_list)
		{
			#print $database_se_list[$j]."\n";
			my @elements2 = split(/\&/, $key);
			my $target2_gid = $elements2[0];
			my $target2_chr = $elements2[1];
			my $target2_strand =$elements2[2];
			my $target2_up_exst=$elements2[5];
			$target2_up_exst++;
			my $target2_up_exen=$elements2[6];	
			my $target2_down_exst=$elements2[7];
			$target2_down_exst++;
			my $target2_down_exen=$elements2[8];	
			my $target2_exst=$elements2[6]+1;
			my $target2_exen=$elements2[4];			
			next if($target_chr ne $target2_chr);
			next if($target_strand ne $target2_strand);
			my $left_overlap = overlap($target_up_exst, $target_up_exen, $target2_up_exst, $target2_up_exen);
			my $right_overlap = overlap($target_down_exst, $target_down_exen, $target2_down_exst, $target2_down_exen);
			if($left_overlap != 12 and $right_overlap != 12 and $target_exst eq $target2_exst and $target_exen eq $target2_exen)
			{
		      #print $Control1_se_list[$i]."\n";
			  $known{$Control1_se_list[$i]}=1;
			  last;
			}
        }			
	}
    if($target_strand eq "+")
	{	
		my $target_up_exst=$elements[8];
		$target_up_exst++;
		my $target_up_exen=$elements[9];		
		my $target_down_exst=$elements[6];
		$target_down_exst++;
		my $target_down_exen=$elements[7];		
		my $target_exst=$elements[4];
		$target_exst++;
		my $target_exen=$target_down_exst - 1;
		foreach my $key (keys %database_se_list)
		{
			#print $database_se_list[$j]."\n";

			my @elements2 = split(/\&/, $key);
			my $target2_gid = $elements2[0];
			my $target2_chr = $elements2[1];
			my $target2_strand =$elements2[2];
			
			my $target2_up_exst=$elements2[7];
			$target2_up_exst++;
			my $target2_up_exen=$elements2[8];
		
			my $target2_down_exst=$elements2[5];
			$target2_down_exst++;
			my $target2_down_exen=$elements2[6];
		
			my $target2_exst=$elements2[3];
			$target2_exst++;
			my $target2_exen=$target2_down_exst - 1;			
			next if($target_chr ne $target2_chr);
			next if($target_strand ne $target2_strand);
			my $left_overlap = overlap($target_up_exst, $target_up_exen, $target2_up_exst, $target2_up_exen);
			my $right_overlap = overlap($target_down_exst, $target_down_exen, $target2_down_exst, $target2_down_exen);
			if($left_overlap != 12 and $right_overlap != 12 and $target_exst eq $target2_exst and $target_exen eq $target2_exen)
			{
		      #print $Control1_se_list[$i]."\n";
			  $known{$Control1_se_list[$i]}=1;
			  last;
			}
        }			
	}
}

#print scalar(keys %know)."\n";
my $otfile1 = $ASEfile;
$otfile1="A3SS_known.txt";
my $otfile2 = $ASEfile;
$otfile2="A3SS_novel.txt";
open(OUT1, ">$otfile1");
open(OUT2, ">$otfile2");
for(my $i = 0; $i<=$#Control1_se_list; $i++)
{
     if(defined $known{$Control1_se_list[$i]})
	 {
		print OUT1 $Control1_se_list[$i]."\n";
	 }
	 else
	 {
		print OUT2 $Control1_se_list[$i]."\n";
	 }	 
}
close(OUT1);
close(OUT2);








my %database_se_list;
open(IN, "/home/database/AS.SE.txt");
<IN>;
while(<IN>)
{
    chomp;
	s/\s+$//;
	s/\"//g;
	my @arr = split(/\t/, $_);
	my $key = $arr[1]."&".$arr[3]."&".$arr[4]."&".$arr[5]."&".$arr[6]."&".$arr[7]."&".$arr[8]."&".$arr[9]."&".$arr[10];
	$database_se_list{$key}=1;
}
close(IN);
my @Control1_se_list;
open(IN, $ASEfile);
<IN>;
while(<IN>)
{
    chomp;
	s/\s+$//;
	my @tmp = split(/\t/, $_);
	my @arr = split(/\&/, $tmp[0]);
	next if($arr[0] ne "SE");
	#push(@Control1_se_list, join("&", @arr[1..$#arr]));
	push(@Control1_se_list, join("&", @arr[0..$#arr]));	
}
close(IN);
#print join("\n", @Control1_se_list)."\n";


my %known;
for(my $i = 0; $i <= $#Control1_se_list; $i++)
{
	#print $Control1_se_list[$i]."\n";
	my @elements = split(/\&/, $Control1_se_list[$i]);
	my $target_gid = $elements[1];
	my $target_chr = $elements[2];
	my $target_strand =	$elements[3];	
	my $target_up_exst=$elements[6];
	$target_up_exst++;
	my $target_up_exen=$elements[7];	
	my $target_down_exst=$elements[8];
	$target_down_exst++;
	my $target_down_exen=$elements[9];	
	my $target_exst=$elements[4];
	$target_exst++;
	my $target_exen=$elements[5];
	foreach my $key (keys %database_se_list)
	{
		#print $database_se_list[$j]."\n";
		my @elements2 = split(/\&/, $key);
		my $target2_gid = $elements2[0];
		my $target2_chr = $elements2[1];
		my $target2_strand = $elements2[2];	
		my $target2_up_exst=$elements2[5];
		$target2_up_exst++;
		my $target2_up_exen=$elements2[6];	
		my $target2_down_exst=$elements2[7];
		$target2_down_exst++;
		my $target2_down_exen=$elements2[8];	
		my $target2_exst=$elements2[3];
		$target2_exst++;
		my $target2_exen=$elements2[4];
		next if($target_chr ne $target2_chr);
		next if($target_strand ne $target2_strand);
		my $left_overlap = overlap($target_up_exst, $target_up_exen, $target2_up_exst, $target2_up_exen);
		my $right_overlap = overlap($target_down_exst, $target_down_exen, $target2_down_exst, $target2_down_exen);
		if($left_overlap != 12 and $right_overlap != 12 and $target_exst eq $target2_exst and $target_exen eq $target2_exen)
		{
		      #print $Control1_se_list[$i]."\n";
			  $known{$Control1_se_list[$i]}=1;
			  last;
		}		
	}
}


my $otfile1 = $ASEfile;
$otfile1="SE_known.txt";
my $otfile2 = $ASEfile;
$otfile2="SE_novel.txt";
open(OUT1, ">$otfile1");
open(OUT2, ">$otfile2");
for(my $i = 0; $i<=$#Control1_se_list; $i++)
{
     if(defined $known{$Control1_se_list[$i]})
	 {
		print OUT1 $Control1_se_list[$i]."\n";
	 }
	 else
	 {
		print OUT2 $Control1_se_list[$i]."\n";
	 }	 
}
close(OUT1);
close(OUT2);




my %database_se_list;
open(IN, "/home/database/AS.A5SS.txt");
<IN>;
while(<IN>)
{
    chomp;
	s/\s+$//;
	s/\"//g;
	my @arr = split(/\t/, $_);
	my $key = $arr[1]."&".$arr[3]."&".$arr[4]."&".$arr[5]."&".$arr[6]."&".$arr[7]."&".$arr[8]."&".$arr[9]."&".$arr[10];
	$database_se_list{$key}=1;
}
close(IN);
my @Control1_se_list;
open(IN, $ASEfile);
<IN>;
while(<IN>)
{
    chomp;
	s/\s+$//;
	my @tmp = split(/\t/, $_);
	my @arr = split(/\&/, $tmp[0]);
	next if($arr[0] ne "A5SS");
	#push(@Control1_se_list, join("&", @arr[1..$#arr]));
	push(@Control1_se_list, join("&", @arr[0..$#arr]));	
}
close(IN);
my %known;
for(my $i = 0; $i <= $#Control1_se_list; $i++)
{
	#print $Control1_se_list[$i]."\n";
	my @elements = split(/\&/, $Control1_se_list[$i]);
	my $target_gid = $elements[1];
	my $target_chr = $elements[2];
	my $target_strand =	$elements[3];
    if($target_strand eq "+")
	{	
		my $target_up_exst=$elements[6];
		$target_up_exst++;
		my $target_up_exen=$elements[7];	
		my $target_down_exst=$elements[8];
		$target_down_exst++;
		my $target_down_exen=$elements[9];	
		my $target_exst=$elements[7]+1;
		$target_exst++;
		my $target_exen=$elements[5];
		foreach my $key (keys %database_se_list)
		{
			#print $database_se_list[$j]."\n";

			my @elements2 = split(/\&/, $key);
			my $target2_gid = $elements2[0];
			my $target2_chr = $elements2[1];
			my $target2_strand =$elements2[2];
			my $target2_up_exst=$elements2[5];
			$target2_up_exst++;
			my $target2_up_exen=$elements2[6];	
			my $target2_down_exst=$elements2[7];
			$target2_down_exst++;
			my $target2_down_exen=$elements2[8];	
			my $target2_exst=$elements2[6]+1;
			$target2_exst++;
			my $target2_exen=$elements2[4];
			
			next if($target_chr ne $target2_chr);
			next if($target_strand ne $target2_strand);
			my $left_overlap = overlap($target_up_exst, $target_up_exen, $target2_up_exst, $target2_up_exen);
			my $right_overlap = overlap($target_down_exst, $target_down_exen, $target2_down_exst, $target2_down_exen);
			if($left_overlap != 12 and $right_overlap != 12 and $target_exst eq $target2_exst and $target_exen eq $target2_exen)
			{
		      #print $Control1_se_list[$i]."\n";
			  $known{$Control1_se_list[$i]}=1;
			  last;
			}
        }			
	}
    if($target_strand eq "-")
	{	
		my $target_up_exst=$elements[8];
		$target_up_exst++;
		my $target_up_exen=$elements[9];
		
		my $target_down_exst=$elements[6];
		$target_down_exst++;
		my $target_down_exen=$elements[7];
		
		my $target_exst=$elements[4]+1;
		$target_exst++;
		my $target_exen=$target_down_exst - 1;
		foreach my $key (keys %database_se_list)
		{
			#print $database_se_list[$j]."\n";

			my @elements2 = split(/\&/, $key);
			my $target2_gid = $elements2[0];
			my $target2_chr = $elements2[1];
			my $target2_strand =$elements2[2];
			
			my $target2_up_exst=$elements2[7];
			$target2_up_exst++;
			my $target2_up_exen=$elements2[8];
		
			my $target2_down_exst=$elements2[5];
			$target2_down_exst++;
			my $target2_down_exen=$elements2[6];
		
			my $target2_exst=$elements2[3]+1;
			$target2_exst++;
			my $target2_exen=$target2_down_exst - 1;
			
			next if($target_chr ne $target2_chr);
			next if($target_strand ne $target2_strand);
			my $left_overlap = overlap($target_up_exst, $target_up_exen, $target2_up_exst, $target2_up_exen);
			my $right_overlap = overlap($target_down_exst, $target_down_exen, $target2_down_exst, $target2_down_exen);
			if($left_overlap != 12 and $right_overlap != 12 and $target_exst eq $target2_exst and $target_exen eq $target2_exen)
			{
		      #print $Control1_se_list[$i]."\n";
			  $known{$Control1_se_list[$i]}=1;
			  last;
			}
        }			
	}
}





my $otfile1 = $ASEfile;
$otfile1="A5SS_known.txt";
my $otfile2 = $ASEfile;
$otfile2="A5SS_novel.txt";
open(OUT1, ">$otfile1");
open(OUT2, ">$otfile2");
for(my $i = 0; $i<=$#Control1_se_list; $i++)
{
     if(defined $known{$Control1_se_list[$i]})
	 {
		print OUT1 $Control1_se_list[$i]."\n";
	 }
	 else
	 {
		print OUT2 $Control1_se_list[$i]."\n";
	 }	 
}
close(OUT1);
close(OUT2);





my %database_se_list;
open(IN, "/home/database/AS.MXE.txt");
<IN>;
while(<IN>)
{
    chomp;
	s/\s+$//;
	s/\"//g;
	my @arr = split(/\t/, $_);
	my $key = $arr[1]."&".$arr[3]."&".$arr[4]."&".$arr[5]."&".$arr[6]."&".$arr[7]."&".$arr[8]."&".$arr[9]."&".$arr[10]."&".$arr[11]."&".$arr[12];
	$database_se_list{$key}=1;
}
close(IN);
my @Control1_se_list;
open(IN, $ASEfile);
<IN>;
while(<IN>)
{
    chomp;
	s/\s+$//;
	my @tmp = split(/\t/, $_);
	my @arr = split(/\&/, $tmp[0]);
	next if($arr[0] ne "MXE");
	#push(@Control1_se_list, join("&", @arr[1..$#arr]));
	push(@Control1_se_list, join("&", @arr[0..$#arr]));	
}
close(IN);

my %known;
for(my $i = 0; $i <= $#Control1_se_list; $i++)
{
	#print $Control1_se_list[$i]."\n";
	my @elements = split(/\&/, $Control1_se_list[$i]);
	my $target_gid = $elements[1];
	my $target_chr = $elements[2];
	my $target_strand =	$elements[3];	
	my $target_up_exst=$elements[8];
	$target_up_exst++;
	my $target_up_exen=$elements[9];	
	my $target_down_exst=$elements[10];
	$target_down_exst++;
	my $target_down_exen=$elements[11];	
	my $target_first_exst=$elements[4];
	$target_first_exst++;
	my $target_first_exen=$elements[5];
	my $target_sec_exst=$elements[6];
	$target_sec_exst++;
	my $target_sec_exen=$elements[7];
	foreach my $key (keys %database_se_list)
	{
		#print $database_se_list[$j]."\n";
		my @elements2 = split(/\&/, $key);
		my $target2_gid = $elements2[0];
		my $target2_chr = $elements2[1];
		my $target2_strand = $elements2[2];	
		my $target2_up_exst=$elements2[7];
		$target2_up_exst++;
		my $target2_up_exen=$elements2[8];	
		my $target2_down_exst=$elements2[9];
		$target2_down_exst++;
		my $target2_down_exen=$elements2[10];
		
		
		my $target2_first_exst=$elements2[3];
		$target2_first_exst++;
		my $target2_first_exen=$elements2[4];
		my $target2_sec_exst=$elements2[5];
		$target2_sec_exst++;
		my $target2_sec_exen=$elements2[6];
		
		
		next if($target_chr ne $target2_chr);
		next if($target_strand ne $target2_strand);
		my $left_overlap = overlap($target_up_exst, $target_up_exen, $target2_up_exst, $target2_up_exen);
		my $right_overlap = overlap($target_down_exst, $target_down_exen, $target2_down_exst, $target2_down_exen);
		if($left_overlap != 12 and $right_overlap != 12 and $target_first_exst eq $target2_first_exst and $target_first_exen eq $target2_first_exen and $target_sec_exst eq $target2_sec_exst and $target_sec_exen eq $target2_sec_exen)
		{
		      #print $Control1_se_list[$i]."\n";
			  $known{$Control1_se_list[$i]}=1;
			  last;
		}		
	}
}




my $otfile1 = $ASEfile;
$otfile1="MXE_known.txt";
my $otfile2 = $ASEfile;
$otfile2="MXE_novel.txt";
open(OUT1, ">$otfile1");
open(OUT2, ">$otfile2");
for(my $i = 0; $i<=$#Control1_se_list; $i++)
{
     if(defined $known{$Control1_se_list[$i]})
	 {
		print OUT1 $Control1_se_list[$i]."\n";
	 }
	 else
	 {
		print OUT2 $Control1_se_list[$i]."\n";
	 }	 
}
close(OUT1);
close(OUT2);


###### Sub function overlap ########
sub overlap{
	my @temp=@_;
	my @target=($temp[0], $temp[1]);
	my @query=($temp[2], $temp[3]);
	my $hit;
	if($target[1] < $target[0]){
		die "Target interval is wrong!";
	}
	if($query[1] < $query[0]){
		die "Query interval is wrong!";
	}
	if($query[0] < $target[0] && $query[1] == $target[0]){
        	$hit=1;
	}elsif($query[0] < $target[0] && $query[1] > $target[0] && $query[1] < $target[1]){
        	$hit=2;
	}elsif($query[0] < $target[0] && $query[1] == $target[1]){
        	$hit=3;
	}elsif($query[0] < $target[0] && $query[1] > $target[1]){
        	$hit=4;
	}elsif($query[0] == $target[0] && $query[1] > $target[0] && $query[1] < $target[1]){
        	$hit=5;
	}elsif($query[0] == $target[0] && $query[1] == $target[1]){
        	$hit=6;
	}elsif($query[0] == $target[0] && $query[1] > $target[1]){
        	$hit=7;
	}elsif($query[0] > $target[0] && $query[1] < $target[1]){
        	$hit=8;
	}elsif($query[0] > $target[0] && $query[1] == $target[1]){
        	$hit=9;
	}elsif($query[0] > $target[0] && $query[0] < $target[1] && $query[1] > $target[1]){
        	$hit=10;
	}elsif($query[0] == $target[1]){
        	$hit=11;
	}else{
        	$hit=12;
	}
	return $hit;
}
