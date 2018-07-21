use List::Uniq ':all';
use List::Util qw(first max maxstr min minstr reduce shuffle sum);
use Getopt::Long;
use Statistics::R;
my %transcript;
open(IN, "all_MSU.gtf");
while(<IN>)
{
    chomp;
    my @arr = split(/\t/, $_);
    my $refid = $arr[0];
	#next if(not defined $chrs{$refid});
	#$arr[0] = $chrs{$refid};
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
		push(@{$transcript{$tid} -> {exst}}, $arr[3]);
		push(@{$transcript{$tid} -> {exen}}, $arr[4]);
	}
}
close(IN);


foreach my $key (keys %transcript)
{
    my $chr = $transcript{$key}->{chr};
    my $strand = $transcript{$key}->{strand};
    my @exsts = @{$transcript{$key}->{exst}};
    my @exens = @{$transcript{$key}->{exen}};
   	my $gst = min(@exsts);
	my $gen = max(@exens);
		print $chr."\t".($gst-1)."\t".$gen."\t".$key."\t"."."."\t".$strand."\n";
}
close(OUT);


my %transcript;
open(IN, "repre_transcripts_exon.gtf");
while(<IN>)
{
    chomp;
    my @arr = split(/\t/, $_);
    my $refid = $arr[0];
	#next if(not defined $chrs{$refid});
	#$arr[0] = $chrs{$refid};
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
		push(@{$transcript{$tid} -> {exst}}, $arr[3]);
		push(@{$transcript{$tid} -> {exen}}, $arr[4]);
	}
}
close(IN);


foreach my $key (keys %transcript)
{
    my $chr = $transcript{$key}->{chr};
    my $strand = $transcript{$key}->{strand};
    my @exsts = @{$transcript{$key}->{exst}};
    my @exens = @{$transcript{$key}->{exen}};
   	my $gst = min(@exsts);
	my $gen = max(@exens);
		print $chr."\t".($gst-1)."\t".$gen."\t".$key."\t"."."."\t".$strand."\n";
}
close(OUT);






my %transcript;
open(IN, "predicted_transcripts.gtf");
while(<IN>)
{
    chomp;
    my @arr = split(/\t/, $_);
    my $refid = $arr[0];
	#next if(not defined $chrs{$refid});
	#$arr[0] = $chrs{$refid};
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
		push(@{$transcript{$tid} -> {exst}}, $arr[3]);
		push(@{$transcript{$tid} -> {exen}}, $arr[4]);
	}
}
close(IN);


foreach my $key (keys %transcript)
{
    my $chr = $transcript{$key}->{chr};
    my $strand = $transcript{$key}->{strand};
    my @exsts = @{$transcript{$key}->{exst}};
    my @exens = @{$transcript{$key}->{exen}};
   	my $gst = min(@exsts);
	my $gen = max(@exens);
		print $chr."\t".($gst-1)."\t".$gen."\t".$key."\t"."."."\t".$strand."\n";
}
close(OUT);






