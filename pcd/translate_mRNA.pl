use strict;
use Getopt::Long;
use vars qw($fafile $pos $gtffile);
Getopt::Long::GetOptions(
	'fa=s'     =>\$fafile,
	'pos=s'    => \$pos,
	'gtf=s'    => \$gtffile,
);
$/ = ">";
my $cdsst;
my $cdsen;

open(IN, $fafile);
<IN>;
while(<IN>)
{
    my @arr = split(/\n/, $_);
    my $id = shift(@arr);
	my $seq = join("", @arr);
	#$arr[$#arr]=~s/\>$//;
	$seq=~s/\>$//;
	#print ">".$id."\n".$seq."\n";
	my $protein = "";
	my $codon;
	$seq=substr($seq, $pos);
	for(my $i=0;$i<(length($seq)-2);$i+=3)
	{
		$codon=substr($seq,$i,3);
		my $symbol=&codon2aa($codon);
		$protein=$protein.$symbol;
		if($symbol eq "_")		{
		    
			#print $symbol."\n";
			last;
		}
	}
	$protein=~s/\_$//;
	my $cdslen=length($protein)*3;
	$cdsst = $pos+1;
	$cdsen=$cdslen+$cdsst-1+3;
	#print $cdsst."\t".$cdsen."\n";
	my $cdsseq = substr($seq, ($cdsst-1), ($cdsen-$cdsst+1));
}
close(IN);
$/="\n";
my $tid;
my $gid;
my @exsst;
my @exsen;	
my $strand="+";
my $chr;
open(IN, $gtffile);
while(<IN>)
{
    chomp;
    my @arr = split(/\t/, $_);
    my $refid = $arr[0];
	my @tmp = split(/;\s/, $arr[8]);
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
    $strand = $arr[6];
	$chr=$arr[0];
	if($arr[2] eq "exon")
	{
	    $arr[3]=~s/\s+$//;
	    $arr[3]=~s/^\s+//;
	    $arr[4]=~s/\s+$//;
	    $arr[4]=~s/^\s+//;		
		#print $tid."\t".$gid."\t".$arr[3]."\t".$arr[4]."\n";
		push(@exsst, $arr[3]);
		push(@exsen, $arr[4]);
    }
}
close(IN);
#print join(",", @exsst)."\n";
#print join(",", @exsen)."\n";
my @sorted_exsts = sort {$exsst[$a] <=> $exsst[$b]} 0..$#exsst;
@exsst=@exsst[@sorted_exsts];
@exsen=@exsen[@sorted_exsts];
my $st_codon;
my $en_codon;
if($strand eq "+")
{
	my $sum1 = 0;
	my $sum2 = 0;
    for(my $i = 0; $i <= $#exsst; $i++)
    {
        for(my $j = $exsst[$i]; $j <= $exsen[$i]; $j++)
        {
          	$sum1++;$sum2++;
            if($cdsst eq $sum1)
            {
				$st_codon = $j
            } 
            if($cdsen eq $sum1)
            {
				$en_codon = $j
            } 					
       }
   }
}
if($strand eq "-")
{
	my $sum1 = 0;
	my $sum2 = 0;
    for(my $i = $#exsst; $i >= 0; $i--)
    {
        for(my $j = $exsen[$i]; $j >= $exsst[$i]; $j--)
        {
          	$sum1++;$sum2++;
           if($cdsst eq $sum1)
           {
    		$en_codon = $j
            } 
            if($cdsen eq $sum1)
            {
 				$st_codon = $j
            } 					
        }
    }
}
for(my $i = 0; $i <= $#exsst; $i++)
{
	my $exst = $exsst[$i];my $exen = $exsen[$i];
	my $ishit = overlap($st_codon, $en_codon, $exst, $exen);
	if($ishit == 1)
	{
	    if($strand eq "+") {print $chr."\t"."taco"."\t"."five_prime_UTR"."\t".$exst."\t".($exen-1)."\t"."."."\t".$strand."\t"."."."\t"."gene_id "."\"".$gid."\"; "."transcript_id "."\"".$tid."\";"."\n";}
		if($strand eq "-") {print $chr."\t"."taco"."\t"."three_prime_UTR"."\t".$exst."\t".($exen-1)."\t"."."."\t".$strand."\t"."."."\t"."gene_id "."\"".$gid."\"; "."transcript_id "."\"".$tid."\";"."\n";}
		print $chr."\t"."taco"."\t"."CDS"."\t".$exen."\t".$exen."\t"."."."\t".$strand."\t"."."."\t"."gene_id "."\"".$gid."\"; "."transcript_id "."\"".$tid."\";"."\n";
	}
	if($ishit == 2 or $ishit == 3)
	{
	    if($strand eq "+") {print $chr."\t"."taco"."\t"."five_prime_UTR"."\t".$exst."\t".($st_codon-1)."\t"."."."\t".$strand."\t"."."."\t"."gene_id "."\"".$gid."\"; "."transcript_id "."\"".$tid."\";"."\n";}
		if($strand eq "-") {print $chr."\t"."taco"."\t"."three_prime_UTR"."\t".$exst."\t".($st_codon-1)."\t"."."."\t".$strand."\t"."."."\t"."gene_id "."\"".$gid."\"; "."transcript_id "."\"".$tid."\";"."\n";}
		print $chr."\t"."taco"."\t"."CDS"."\t".$st_codon."\t".$exen."\t"."."."\t".$strand."\t"."."."\t"."gene_id "."\"".$gid."\"; "."transcript_id "."\"".$tid."\";"."\n";
	}
	if($ishit == 4)
	{
	    if($strand eq "+") {print $chr."\t"."taco"."\t"."five_prime_UTR"."\t".$exst."\t".($st_codon-1)."\t"."."."\t".$strand."\t"."."."\t"."gene_id "."\"".$gid."\"; "."transcript_id "."\"".$tid."\";"."\n";}
	    if($strand eq "-") {print $chr."\t"."taco"."\t"."three_prime_UTR"."\t".$exst."\t".($st_codon-1)."\t"."."."\t".$strand."\t"."."."\t"."gene_id "."\"".$gid."\"; "."transcript_id "."\"".$tid."\";"."\n";}
	    print $chr."\t"."taco"."\t"."CDS"."\t".$st_codon."\t".$en_codon."\t"."."."\t".$strand."\t"."."."\t"."gene_id "."\"".$gid."\"; "."transcript_id "."\"".$tid."\";"."\n";
	    if($strand eq "+") {print $chr."\t"."taco"."\t"."three_prime_UTR"."\t".($en_codon+1)."\t".$exen."\t"."."."\t".$strand."\t"."."."\t"."gene_id "."\"".$gid."\"; "."transcript_id "."\"".$tid."\";"."\n";}
	    if($strand eq "-") {print $chr."\t"."taco"."\t"."five_prime_UTR"."\t".($en_codon+1)."\t".$exen."\t"."."."\t".$strand."\t"."."."\t"."gene_id "."\"".$gid."\"; "."transcript_id "."\"".$tid."\";"."\n";}
	}
			if($ishit == 5 or $ishit == 6)
			{
			    print $chr."\t"."taco"."\t"."CDS"."\t".$st_codon."\t".$exen."\t"."."."\t".$strand."\t"."."."\t"."gene_id "."\"".$gid."\"; "."transcript_id "."\"".$tid."\";"."\n";
			}
			if($ishit == 7)
			{
			    print $chr."\t"."taco"."\t"."CDS"."\t".$st_codon."\t".$en_codon."\t"."."."\t".$strand."\t"."."."\t"."gene_id "."\"".$gid."\"; "."transcript_id "."\"".$tid."\";"."\n";
			    if($strand eq "+")  {print $chr."\t"."taco"."\t"."three_prime_UTR"."\t".($en_codon+1)."\t".$exen."\t"."."."\t".$strand."\t"."."."\t"."gene_id "."\"".$gid."\"; "."transcript_id "."\"".$tid."\";"."\n";}
			    if($strand eq "-")  {print $chr."\t"."taco"."\t"."five_prime_UTR"."\t".($en_codon+1)."\t".$exen."\t"."."."\t".$strand."\t"."."."\t"."gene_id "."\"".$gid."\"; "."transcript_id "."\"".$tid."\";"."\n";}
			}			
			if($ishit == 8 or $ishit == 9)
			{
			    print $chr."\t"."taco"."\t"."CDS"."\t".$exst."\t".$exen."\t"."."."\t".$strand."\t"."."."\t"."gene_id "."\"".$gid."\"; "."transcript_id "."\"".$tid."\";"."\n";
			}
			if($ishit == 10 or $ishit == 11)
			{
			    print $chr."\t"."taco"."\t"."CDS"."\t".$exst."\t".$en_codon."\t"."."."\t".$strand."\t"."."."\t"."gene_id "."\"".$gid."\"; "."transcript_id "."\"".$tid."\";"."\n";
			    if($strand eq "+") {print $chr."\t"."taco"."\t"."three_prime_UTR"."\t".($en_codon+1)."\t".$exen."\t"."."."\t".$strand."\t"."."."\t"."gene_id "."\"".$gid."\"; "."transcript_id "."\"".$tid."\";"."\n";}
			    if($strand eq "-") {print $chr."\t"."taco"."\t"."five_prime_UTR"."\t".($en_codon+1)."\t".$exen."\t"."."."\t".$strand."\t"."."."\t"."gene_id "."\"".$gid."\"; "."transcript_id "."\"".$tid."\";"."\n";}
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

	
	
	

sub codon2aa{
	my($codon)=@_;
	$codon=uc($codon);
	my(%g)=('TCA'=>'S','TCC'=>'S','TCG'=>'S','TCT'=>'S','TTC'=>'F','TTT'=>'F','TTA'=>'L','TTG'=>'L','TAC'=>'Y','TAT'=>'Y','TAA'=>'_','TAG'=>'_','TGC'=>'C','TGT'=>'C','TGA'=>'_','TGG'=>'W','CTA'=>'L','CTC'=>'L','CTG'=>'L','CTT'=>'L','CCA'=>'P','CCC'=>'P','CCG'=>'P','CCT'=>'P','CAC'=>'H','CAT'=>'H','CAA'=>'Q','CAG'=>'Q','CGA'=>'R','CGC'=>'R','CGG'=>'R','CGT'=>'R','ATA'=>'I','ATC'=>'I','ATT'=>'I','ATG'=>'M','ACA'=>'T','ACC'=>'T','ACG'=>'T','ACT'=>'T','AAC'=>'N','AAT'=>'N','AAA'=>'K','AAG'=>'K','AGC'=>'S','AGT'=>'S','AGA'=>'R','AGG'=>'R','GTA'=>'V','GTC'=>'V','GTG'=>'V','GTT'=>'V','GCA'=>'A','GCC'=>'A','GCG'=>'A','GCT'=>'A','GAC'=>'D','GAT'=>'D','GAA'=>'E','GAG'=>'E','GGA'=>'G','GGC'=>'G','GGG'=>'G','GGT'=>'G');
	if(exists $g{$codon})
	{
		return $g{$codon};
	}
	else
	{
		#print STDERR "Bad codon \"$codon\"!!\n";
		#exit;
		return "*";
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

