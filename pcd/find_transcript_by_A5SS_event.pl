use List::Util qw(first max maxstr min minstr reduce shuffle sum);
use Getopt::Long;
use vars qw($gtffile $A5SSeventfile $bedfile);
Getopt::Long::GetOptions(
	'gtf=s'    => \$gtffile,
	'A5SS=s'     => \$A5SSeventfile,
);

my %transcript;
open(IN, "$gtffile");
while(<IN>)
{
    chomp;
    my @arr = split(/\t/, $_);
    my $refid = $arr[0];
	my @tmp = split(/;\s+/, $arr[8]);
    my $tid;
	my $gid;
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
    my $strand = $arr[6];
	$arr[3]=~s/\s+$//;
	$arr[3]=~s/^\s+//;
	$arr[4]=~s/\s+$//;
	$arr[4]=~s/^\s+//; 	
	unless(defined $transcript{$tid})
	{
	    $transcript{$tid} = {
			chr => $arr[0],
			strand => $strand,
			geneid=>$gid,
		};
	}

	if($arr[2] eq "exon")
	{
		push(@{$transcript{$tid}->{exonst}}, $arr[3]);
		push(@{$transcript{$tid}->{exonen}}, $arr[4]);
	}
}
close(IN);


foreach my $key (keys %transcript)
{

}

foreach my $key (keys %transcript)
{

}


open(IN, "$A5SSeventfile");
my $tt = <IN>;
print "Event"."\t"."Inclusive_trascript_ids"."\t"."Exclusive_trascript_ids"."\n";
while(<IN>)
{
    chomp;
	s/\s+$//;
	s/\"//g;
	my $line = $_;
	my @arr = split(/\t/, $line);	
	my @elements = split(/\&/, $arr[0]);
    next if($elements[0] ne "A5SS");
	my $target_strand =	$elements[3];
	if($target_strand eq "+")
	{
	
		my $target_gid = $elements[1];
		$target_gid=~s/\"//g;
		my $target_chr = $elements[2];	
		
		my $target_exst=$elements[7]+1;
		my $target_exen=$elements[5];
	
		my $target_up_exst=$elements[4];
		$target_up_exst++;
		my $target_up_exen=$elements[7];
	
		my $target_down_exst=$elements[8];
		$target_down_exst++;
		my $target_down_exen=$elements[9];
		#print $target_exst."\t".$target_exen."\n";
		my @transcripts_in=();
		my @transcripts_ex=();	
		foreach my $key (keys %transcript)
		{
			my $tid = $key;
			my $chr = $transcript{$tid}->{chr};
			my $strand = $transcript{$tid}->{strand};
			my $gid = $transcript{$tid}->{geneid};
			next if($gid ne $target_gid);
			next if($strand ne $target_strand);
			next if($chr ne $target_chr);
			my @exsst = @{$transcript{$tid}->{exonst}};
			my @exsen = @{$transcript{$tid}->{exonen}};
			my @sorted_exsts = sort {$exsst[$a] <=> $exsst[$b]} 0..$#exsst;
			@exsst=@exsst[@sorted_exsts];
			@exsen=@exsen[@sorted_exsts];
			if(scalar(@exsst) >= 2)
			{
				for(my $pos1 = 0; $pos1 <= ($#exsst - 1); $pos1++)
				{
					my $add1 = $pos1 + 1;
					#my $add2 = $pos1 + 2;				
					my @two1st = ($exsst[$pos1], $exsst[$add1]);
					my @two1en = ($exsen[$pos1], $exsen[$add1]);				
					my $left_overlap = overlap($target_up_exst, $target_up_exen, $two1st[0], $two1en[0]);
					my $right_overlap = overlap($target_down_exst, $target_down_exen, $two1st[1], $two1en[1]);
					my $middle_overlap = overlap($target_exst, $target_exen, $two1st[0], $two1en[0]);
					if($left_overlap != 12 and $right_overlap != 12 and $middle_overlap != 12 and $two1en[0] eq $target_exen and $two1st[1] eq $target_down_exst)
					{
						#print $tid."\n";
						push(@transcripts_in, $tid);
					}
				}
			}
			if(scalar(@exsst) >= 2)
			{
				for(my $pos1 = 0; $pos1 <= ($#exsst - 1); $pos1++)
				{
					my $add1 = $pos1 + 1;
					#my $add2 = $pos1 + 2;				
					my @two1st = ($exsst[$pos1], $exsst[$add1]);
					my @two1en = ($exsen[$pos1], $exsen[$add1]);			
					my $left_overlap = overlap($target_up_exst, $target_up_exen, $two1st[0], $two1en[0]);
					my $right_overlap = overlap($target_down_exst, $target_down_exen, $two1st[1], $two1en[1]);
					my $middle_overlap = overlap($target_exst, $target_exen, $two1st[0], $two1en[0]);
					if($left_overlap != 12 and $right_overlap != 12 and $middle_overlap == 12 and $two1en[0] eq $target_up_exen and $two1st[1] eq $target_down_exst)
					
					{
						#print $tid."\n";
						push(@transcripts_ex, $tid);
					}
				}
			}        		
		}
		print $arr[0];
		if(scalar(@transcripts_in) > 0)
		{
			print "\t".join(",", @transcripts_in);
		}
		else
		{
			print "\t"."---";	
		}
		if(scalar(@transcripts_ex) > 0)
		{
			print "\t".join(",", @transcripts_ex)."\n";
		}
		else
		{
			print "\t"."---"."\n";	
		}
	}
	if($target_strand eq "-")
	{
	
		my $target_gid = $elements[1];
		$target_gid=~s/\"//g;
	#if(not defined $gid2gid{$target_gid})
	#{
	     #print $_."\n";
	#	 next;
	#}
	#$target_gid=$gid2gid{$target_gid};
		my $target_chr = $elements[2];	
		
		my $target_exst=$elements[4]+1;
		my $target_exen=$elements[6];
	
		my $target_up_exst=$elements[8];
		$target_up_exst++;
		my $target_up_exen=$elements[9];
	
		my $target_down_exst=$elements[6];
		$target_down_exst++;
		my $target_down_exen=$elements[5];
		#print $target_exst."\t".$target_exen."\n";
		my @transcripts_in=();
		my @transcripts_ex=();	
		foreach my $key (keys %transcript)
		{
			my $tid = $key;
			my $chr = $transcript{$tid}->{chr};
			my $strand = $transcript{$tid}->{strand};
			my $gid = $transcript{$tid}->{geneid};
			next if($gid ne $target_gid);
			next if($strand ne $target_strand);
			next if($chr ne $target_chr);
			my @exsst = @{$transcript{$tid}->{exonst}};
			my @exsen = @{$transcript{$tid}->{exonen}};
			my @sorted_exsts = sort {$exsst[$a] <=> $exsst[$b]} 0..$#exsst;
			@exsst=@exsst[@sorted_exsts];
			@exsen=@exsen[@sorted_exsts];
			if(scalar(@exsst) >= 2)
			{
				for(my $pos1 = 0; $pos1 <= ($#exsst - 1); $pos1++)
				{
					my $add1 = $pos1 + 1;
					#my $add2 = $pos1 + 2;				
					my @two1st = ($exsst[$pos1], $exsst[$add1]);
					my @two1en = ($exsen[$pos1], $exsen[$add1]);				
					my $left_overlap = overlap($target_up_exst, $target_up_exen, $two1st[0], $two1en[0]);
					my $right_overlap = overlap($target_down_exst, $target_down_exen, $two1st[1], $two1en[1]);
					my $middle_overlap = overlap($target_exst, $target_exen, $two1st[1], $two1en[1]);
					if($left_overlap != 12 and $right_overlap != 12 and $middle_overlap != 12 and $two1st[1] eq $target_exst and $two1en[0] eq $target_up_exen)
					{
						#print $tid."\n";
						push(@transcripts_in, $tid);
					}
				}
			}
			if(scalar(@exsst) >= 2)
			{
				for(my $pos1 = 0; $pos1 <= ($#exsst - 1); $pos1++)
				{
					my $add1 = $pos1 + 1;
					#my $add2 = $pos1 + 2;				
					my @two1st = ($exsst[$pos1], $exsst[$add1]);
					my @two1en = ($exsen[$pos1], $exsen[$add1]);			
					my $left_overlap = overlap($target_up_exst, $target_up_exen, $two1st[0], $two1en[0]);
					my $right_overlap = overlap($target_down_exst, $target_down_exen, $two1st[1], $two1en[1]);
					my $middle_overlap = overlap($target_exst, $target_exen, $two1st[1], $two1en[1]);
					if($left_overlap != 12 and $right_overlap != 12 and $middle_overlap == 12 and $two1st[1] eq ($target_exen+1) and $two1en[0] eq $target_up_exen)
					
					{
						#print $tid."\n";
						push(@transcripts_ex, $tid);
					}
				}
			}        		
		}
		print $arr[0];
		if(scalar(@transcripts_in) > 0)
		{
			print "\t".join(",", @transcripts_in);
		}
		else
		{
			print "\t"."---";	
		}
		if(scalar(@transcripts_ex) > 0)
		{
			print "\t".join(",", @transcripts_ex)."\n";
		}
		else
		{
			print "\t"."---"."\n";	
		}
	}
}	
#print scalar(keys %hash)."\n";




######################################
sub merge{
	#my @set = ("2-10", , "12-17", "19-30", "31-40");
	my @set = @_;
	my @set2;
	for(my $i = 0; $i <= $#set; $i++)
	{
		if(scalar(@set2) == 0)
		{
			push(@set2, $set[$i]);
			next;
		}
	    my @tmpy = split(/\-/, $set2[$#set2]);
	    my @tmpx = split(/\-/, $set[$i]);
		if(($tmpx[0] - $tmpy[1]) == 1)
		{
		    my $ele = $tmpy[0]."-".$tmpx[1];
			$set2[$#set2] = $ele;
		}
		else
		{
		    push(@set2, $set[$i])
		}
	}
	return(@set2);
}






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


