my $known_GTAG = 0;
my $known_GCAG = 0;
my $known_other = 0;
my $novel_GTAG = 0;
my $novel_GCAG = 0;
my $novel_other = 0;
open(IN, "splice_site_RI_known_novel.txt");
while(<IN>)
{
    chomp;
	s/\s+$//;
	my @arr = split(/\t/, $_);
	#print $arr[1]."\t".$arr[2]."\t".$arr[3]."\n";
	my $a = $arr[1]."->".$arr[2];
	my $b = $arr[3];
	$known_GTAG++ if($a eq "GT->AG" and $b eq "known");
	$known_GCAG++ if($a eq "GC->AG" and $b eq "known");
	$known_other++ if($a ne "GT->AG" and $a ne "GC->AG" and $b eq "known");

	$novel_GTAG++ if($a eq "GT->AG" and $b eq "novel");
	$novel_GCAG++ if($a eq "GC->AG" and $b eq "novel");
	$novel_other++ if($a ne "GT->AG" and $a ne "GC->AG" and $b eq "novel");
	
	#print $_."\n";
}
close(IN);


open(IN, "splice_site_SE_known_novel.txt");
while(<IN>)
{
    chomp;
	s/\s+$//;
	my @arr = split(/\t/, $_);
	#print $arr[1]."\t".$arr[2]."\t".$arr[3]."\n";
	my $a = $arr[1]."->".$arr[2];
	my $b = $arr[3];
	$known_GTAG++ if($a eq "GT->AG" and $b eq "known");
	$known_GCAG++ if($a eq "GC->AG" and $b eq "known");
	$known_other++ if($a ne "GT->AG" and $a ne "GC->AG" and $b eq "known");

	$novel_GTAG++ if($a eq "GT->AG" and $b eq "novel");
	$novel_GCAG++ if($a eq "GC->AG" and $b eq "novel");
	$novel_other++ if($a ne "GT->AG" and $a ne "GC->AG" and $b eq "novel");
	
	#print $_."\n";
}
close(IN);



open(IN, "splice_site_A5SS_known_novel.txt");
while(<IN>)
{
    chomp;
	s/\s+$//;
	my @arr = split(/\t/, $_);
	#print $arr[1]."\t".$arr[2]."\t".$arr[3]."\n";
	my $a = $arr[1]."->".$arr[2];
	my $b = $arr[3];
	$known_GTAG++ if($a eq "GT->AG" and $b eq "known");
	$known_GCAG++ if($a eq "GC->AG" and $b eq "known");
	$known_other++ if($a ne "GT->AG" and $a ne "GC->AG" and $b eq "known");

	$novel_GTAG++ if($a eq "GT->AG" and $b eq "novel");
	$novel_GCAG++ if($a eq "GC->AG" and $b eq "novel");
	$novel_other++ if($a ne "GT->AG" and $a ne "GC->AG" and $b eq "novel");
	
	#print $_."\n";
}
close(IN);



open(IN, "splice_site_A3SS_known_novel.txt");
while(<IN>)
{
    chomp;
	s/\s+$//;
	my @arr = split(/\t/, $_);
	#print $arr[1]."\t".$arr[2]."\t".$arr[3]."\n";
	my $a = $arr[1]."->".$arr[2];
	my $b = $arr[3];
	$known_GTAG++ if($a eq "GT->AG" and $b eq "known");
	$known_GCAG++ if($a eq "GC->AG" and $b eq "known");
	$known_other++ if($a ne "GT->AG" and $a ne "GC->AG" and $b eq "known");

	$novel_GTAG++ if($a eq "GT->AG" and $b eq "novel");
	$novel_GCAG++ if($a eq "GC->AG" and $b eq "novel");
	$novel_other++ if($a ne "GT->AG" and $a ne "GC->AG" and $b eq "novel");
	
	#print $_."\n";
}
close(IN);


open(IN, "splice_site_MXE_known_novel.txt");
while(<IN>)
{
    chomp;
	s/\s+$//;
	my @arr = split(/\t/, $_);
	#print $arr[1]."\t".$arr[2]."\t".$arr[3]."\n";
	my $a = $arr[1]."->".$arr[2];
	my $b = $arr[3];
	$known_GTAG++ if($a eq "GT->AG" and $b eq "known");
	$known_GCAG++ if($a eq "GC->AG" and $b eq "known");
	$known_other++ if($a ne "GT->AG" and $a ne "GC->AG" and $b eq "known");

	$novel_GTAG++ if($a eq "GT->AG" and $b eq "novel");
	$novel_GCAG++ if($a eq "GC->AG" and $b eq "novel");
	$novel_other++ if($a ne "GT->AG" and $a ne "GC->AG" and $b eq "novel");
	
	#print $_."\n";
}
close(IN);

#print $known_GTAG."\t".$known_GCAG."\t".$known_other."\t".$novel_GTAG."\t".$novel_GCAG."\t".$novel_other."\n";
print "canonical GT->AG"."\t".$known_GTAG."\n"."non-canonical GC->AG"."\t".$known_GCAG."\n"."non-canonical others"."\t".$known_other."\n"."canonical GT->AG"."\t".$novel_GTAG."\n"."non-canonical GC->AG"."\t".$novel_GCAG."\n"."non-canonical others"."\t".$novel_other."\n";

