use List::Util qw(first max maxstr min minstr reduce shuffle sum);
use Getopt::Long;
my %transcript;
open(IN, "novel.gtf");
while(<IN>)
{
    chomp;
    my @arr = split(/\t/, $_);
    my $refid = $arr[0];
	my @tmp = split(/;\s/, $arr[8]);
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
	    if($arr[6] eq ".")
		{
		    $arr[6]="+";
		}	
    my $strand = $arr[6];

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
	    $arr[3]=~s/\s+$//;
	    $arr[3]=~s/^\s+//;
	    $arr[4]=~s/\s+$//;
	    $arr[4]=~s/^\s+//;		
		push(@{$transcript{$tid}->{exonst}}, $arr[3]);
		push(@{$transcript{$tid}->{exonen}}, $arr[4]);
		if($arr[3] >= $arr[4])
		{
		}
	}
}
close(IN);

foreach my $key (keys %transcript)
{
	my $tid = $key;
	my $chr = $transcript{$tid}->{chr};
	my $strand = $transcript{$tid}->{strand};
	my $gid = $transcript{$tid}->{geneid};
	my @exsst = @{$transcript{$tid}->{exonst}};
	my @exsen = @{$transcript{$tid}->{exonen}};
	my @sorted_exsts = sort {$exsst[$a] <=> $exsst[$b]} 0..$#exsst;
	@exsst=@exsst[@sorted_exsts];
	@exsen=@exsen[@sorted_exsts];
	#print $tid."\t".$chr."\t".$gid."\t".$strand."\t".join(",", @exsst)."\t".join(",", @exsen)."\n";
}

my %CDS_pos;
open(IN, "novel.fa.transdecoder_CDS.gtf");
while(<IN>)
{
    chomp;
	s/\s+$//;
	my @arr = split(/\t/, $_);
	my @tmp = split(/;\s/, $arr[8]);
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
	my @exsst = @{$transcript{$tid}->{exonst}};
	my @exsen = @{$transcript{$tid}->{exonen}};
	my $strand = $transcript{$tid}->{strand};
	if($arr[2] eq "CDS")
	{
	    my $st_codon = 0;
	    my $en_codon = 0;		
		if($strand eq "+")
		{
	         #print $_."\n";
			 my $sum1 = 0;
			 my $sum2 = 0;
             for(my $i = 0; $i <= $#exsst; $i++)
             {
                for(my $j = $exsst[$i]; $j <= $exsen[$i]; $j++)
                {
                	$sum1++;$sum2++;					
                    if($arr[3] eq $sum1)
                    {
                        #print $arr[0]."\t".$tid."\t".$j."\n";
						$st_codon = $j
                    } 
                    if($arr[4] eq $sum1)
                    {
                        #print $arr[0]."\t".$tid."\t".$j."\n";
						$en_codon = $j
                    } 					
                }
             }
			$CDS_pos{$tid}->{"CDSst"}=$st_codon;
			$CDS_pos{$tid}->{"CDSen"}=$en_codon;			
		}
		if($strand eq "-")
		{
	         #print $_."\n";
			 my $sum1 = 0;
			 my $sum2 = 0;
             for(my $i = $#exsst; $i >= 0; $i--)
             {
                for(my $j = $exsen[$i]; $j >= $exsst[$i]; $j--)
                {
                	$sum1++;$sum2++;					
                    if($arr[3] eq $sum1)
                    {
                        #print $arr[0]."\t".$tid."\t".$j."\n";
						$st_codon = $j
                    } 
                    if($arr[4] eq $sum1)
                    {
                        #print $arr[0]."\t".$tid."\t".$j."\n";
						$en_codon = $j
                    } 					
                }
             }
			$CDS_pos{$tid}->{"CDSst"}=$en_codon;
			$CDS_pos{$tid}->{"CDSen"}=$st_codon;			 
		}
		
		#print $arr[0]."\t".$strand."\t".$tid."\t".$st_codon."\t".$en_codon."\n";
	}
	
}
close(IN);
foreach my $key (keys %transcript)
{
	my $tid = $key;
	#next if($key=~m/^LOC\_/);
	#next if($key=~m/^Os/);	
	my @exsst = @{$transcript{$tid}->{exonst}};
	my @exsen = @{$transcript{$tid}->{exonen}};
	next if (not defined $CDS_pos{$key});
    my $gid = $transcript{$tid}->{"geneid"};
    my $chr = $transcript{$tid}->{"chr"};	
	my $cdsst = $CDS_pos{$key}->{"CDSst"};
	my $cdsen = $CDS_pos{$key}->{"CDSen"};
	my $strand = $transcript{$tid}->{strand};
    for(my $i = 0; $i <= $#exsst; $i++)
    {
			my $exst = $exsst[$i];my $exen = $exsen[$i];
			my $ishit = overlap($cdsst, $cdsen, $exst, $exen);
			#print $ishit."\n";
			if($ishit == 1)
			{
			    if($strand eq "+") {print $chr."\t"."taco"."\t"."five_prime_UTR"."\t".$exst."\t".($exen-1)."\t"."."."\t".$strand."\t"."."."\t"."gene_id "."\"".$gid."\"; "."transcript_id "."\"".$tid."\";"."\n";}
				if($strand eq "-") {print $chr."\t"."taco"."\t"."three_prime_UTR"."\t".$exst."\t".($exen-1)."\t"."."."\t".$strand."\t"."."."\t"."gene_id "."\"".$gid."\"; "."transcript_id "."\"".$tid."\";"."\n";}
			    print $chr."\t"."taco"."\t"."CDS"."\t".$exen."\t".$exen."\t"."."."\t".$strand."\t"."."."\t"."gene_id "."\"".$gid."\"; "."transcript_id "."\"".$tid."\";"."\n";

			}
			if($ishit == 2 or $ishit == 3)
			{
			    if($strand eq "+") {print $chr."\t"."taco"."\t"."five_prime_UTR"."\t".$exst."\t".($cdsst-1)."\t"."."."\t".$strand."\t"."."."\t"."gene_id "."\"".$gid."\"; "."transcript_id "."\"".$tid."\";"."\n";}
			    if($strand eq "-") {print $chr."\t"."taco"."\t"."three_prime_UTR"."\t".$exst."\t".($cdsst-1)."\t"."."."\t".$strand."\t"."."."\t"."gene_id "."\"".$gid."\"; "."transcript_id "."\"".$tid."\";"."\n";}
				print $chr."\t"."taco"."\t"."CDS"."\t".$cdsst."\t".$exen."\t"."."."\t".$strand."\t"."."."\t"."gene_id "."\"".$gid."\"; "."transcript_id "."\"".$tid."\";"."\n";
			}
			if($ishit == 4)
			{
			    if($strand eq "+") {print $chr."\t"."taco"."\t"."five_prime_UTR"."\t".$exst."\t".($cdsst-1)."\t"."."."\t".$strand."\t"."."."\t"."gene_id "."\"".$gid."\"; "."transcript_id "."\"".$tid."\";"."\n";}
			    if($strand eq "-") {print $chr."\t"."taco"."\t"."three_prime_UTR"."\t".$exst."\t".($cdsst-1)."\t"."."."\t".$strand."\t"."."."\t"."gene_id "."\"".$gid."\"; "."transcript_id "."\"".$tid."\";"."\n";}
			    print $chr."\t"."taco"."\t"."CDS"."\t".$cdsst."\t".$cdsen."\t"."."."\t".$strand."\t"."."."\t"."gene_id "."\"".$gid."\"; "."transcript_id "."\"".$tid."\";"."\n";
			    if($strand eq "+") {print $chr."\t"."taco"."\t"."three_prime_UTR"."\t".($cdsen+1)."\t".$exen."\t"."."."\t".$strand."\t"."."."\t"."gene_id "."\"".$gid."\"; "."transcript_id "."\"".$tid."\";"."\n";}
			    if($strand eq "-") {print $chr."\t"."taco"."\t"."five_prime_UTR"."\t".($cdsen+1)."\t".$exen."\t"."."."\t".$strand."\t"."."."\t"."gene_id "."\"".$gid."\"; "."transcript_id "."\"".$tid."\";"."\n";}
			}
			if($ishit == 5 or $ishit == 6)
			{
			    print $chr."\t"."taco"."\t"."CDS"."\t".$cdsst."\t".$exen."\t"."."."\t".$strand."\t"."."."\t"."gene_id "."\"".$gid."\"; "."transcript_id "."\"".$tid."\";"."\n";
			}
			if($ishit == 7)
			{
			    print $chr."\t"."taco"."\t"."CDS"."\t".$cdsst."\t".$cdsen."\t"."."."\t".$strand."\t"."."."\t"."gene_id "."\"".$gid."\"; "."transcript_id "."\"".$tid."\";"."\n";
			    if($strand eq "+")  {print $chr."\t"."taco"."\t"."three_prime_UTR"."\t".($cdsen+1)."\t".$exen."\t"."."."\t".$strand."\t"."."."\t"."gene_id "."\"".$gid."\"; "."transcript_id "."\"".$tid."\";"."\n";}
			    if($strand eq "-")  {print $chr."\t"."taco"."\t"."five_prime_UTR"."\t".($cdsen+1)."\t".$exen."\t"."."."\t".$strand."\t"."."."\t"."gene_id "."\"".$gid."\"; "."transcript_id "."\"".$tid."\";"."\n";}
			}			
			if($ishit == 8 or $ishit == 9)
			{
			    print $chr."\t"."taco"."\t"."CDS"."\t".$exst."\t".$exen."\t"."."."\t".$strand."\t"."."."\t"."gene_id "."\"".$gid."\"; "."transcript_id "."\"".$tid."\";"."\n";
			}
			if($ishit == 10 or $ishit == 11)
			{
			    print $chr."\t"."taco"."\t"."CDS"."\t".$exst."\t".$cdsen."\t"."."."\t".$strand."\t"."."."\t"."gene_id "."\"".$gid."\"; "."transcript_id "."\"".$tid."\";"."\n";
			    if($strand eq "+") {print $chr."\t"."taco"."\t"."three_prime_UTR"."\t".($cdsen+1)."\t".$exen."\t"."."."\t".$strand."\t"."."."\t"."gene_id "."\"".$gid."\"; "."transcript_id "."\"".$tid."\";"."\n";}
			    if($strand eq "-") {print $chr."\t"."taco"."\t"."five_prime_UTR"."\t".($cdsen+1)."\t".$exen."\t"."."."\t".$strand."\t"."."."\t"."gene_id "."\"".$gid."\"; "."transcript_id "."\"".$tid."\";"."\n";}
			}
			if($ishit == 12)
			{
			    if($strand eq "+") {print $chr."\t"."taco"."\t"."five_prime_UTR"."\t".$exst."\t".$exen."\t"."."."\t".$strand."\t"."."."\t"."gene_id "."\"".$gid."\"; "."transcript_id "."\"".$tid."\";"."\n";}
			    if($strand eq "-") {print $chr."\t"."taco"."\t"."three_prime_UTR"."\t".$exst."\t".$exen."\t"."."."\t".$strand."\t"."."."\t"."gene_id "."\"".$gid."\"; "."transcript_id "."\"".$tid."\";"."\n";}
				
			}
			if($ishit == 13)
			{
			    if($strand eq "+") {print $chr."\t"."taco"."\t"."three_prime_UTR"."\t".$exst."\t".$exen."\t"."."."\t".$strand."\t"."."."\t"."gene_id "."\"".$gid."\"; "."transcript_id "."\"".$tid."\";"."\n";}
			    if($strand eq "-") {print $chr."\t"."taco"."\t"."five_prime_UTR"."\t".$exst."\t".$exen."\t"."."."\t".$strand."\t"."."."\t"."gene_id "."\"".$gid."\"; "."transcript_id "."\"".$tid."\";"."\n";}
			}

	}

	
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
	}elsif($query[1] < $target[0]){
        	$hit=12;
	}elsif($query[0] > $target[1]){
            $hit=13;	
	}
	return $hit;
}


	


