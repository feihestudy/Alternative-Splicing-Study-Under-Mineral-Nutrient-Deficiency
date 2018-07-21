use strict;
use List::Util qw(first max maxstr min minstr reduce shuffle sum);
use List::Uniq ':all';
use strict;
use Getopt::Long;
use vars qw($gtffile);
Getopt::Long::GetOptions(
    'GTF=s' => \$gtffile,
);
$/="\n";
my %transcript;
open(IN, $gtffile);
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
		push(@{$transcript{$tid}->{exst}}, $arr[3]);
		push(@{$transcript{$tid}->{exen}}, $arr[4]);
		if($arr[3] >= $arr[4])
		{
		}
	}
}
close(IN);


my %gene2tr;
my %tid2nexon;
foreach my $tid (keys %transcript)
{
	my %exon;
	my $chr = $transcript{$tid}->{chr};
	my $strand = $transcript{$tid}->{strand};
	my $gid = $transcript{$tid}->{geneid};
	my @exsts = @{$transcript{$tid}->{exst}};
	my @exens = @{$transcript{$tid} -> {exen}};		
    my @sorted_exsts = sort {$exsts[$a] <=> $exsts[$b]} 0..$#exsts;	
	@exsts = @exsts[@sorted_exsts];
	@exens = @exens[@sorted_exsts];
	my $nexon = scalar(@exens);
	$tid2nexon{$tid}->{"exon_number"}=$nexon;
	$tid2nexon{$tid}->{"intron_number"}=$nexon-1;
	my @introns=();
	if($nexon > 1)
	{
		for(my $i = 0; $i <= ($#exsts-1); $i++)
		{
			my $intron_length = $exsts[($i+1)] - 1 - ($exens[$i] + 1) + 1;
			push(@introns, $intron_length);
		}
	}
	else
	{
			push(@introns, 0);
	}
	$tid2nexon{$tid}->{"intron_length"}=join(",", @introns);
	push(@{$gene2tr{$gid}->{"tids"}}, $tid);
	push(@{$gene2tr{$gid}->{"exsts"}}, @exsts);
	push(@{$gene2tr{$gid}->{"exens"}}, @exens);
	
	$gene2tr{$gid}->{"chr"}=$chr;
}

open(OUT, ">gene_intron_distribution1.xls");
print OUT "gene_id\tintron_number\n";
foreach my $key (keys %gene2tr)
{
     my @tids = @{$gene2tr{$key}->{"tids"}};
	 my @g_exsts = @{$gene2tr{$key}->{"exsts"}};
	 my @g_exens = @{$gene2tr{$key}->{"exens"}};
	 my $g_st = min(@g_exsts);
	 my $g_en = max(@g_exens);
	 my $chr = $gene2tr{$key}->{"chr"};
	 my @vs_ex_num;my @vs_in_num;my @vs_ex_len;my @vs_in_len;
	 for(my $i = 0; $i <= $#tids; $i++)
	 {
	    push(@vs_ex_num, $tid2nexon{$tids[$i]}->{"exon_number"});
	    push(@vs_in_num, $tid2nexon{$tids[$i]}->{"intron_number"});		
	    push(@vs_in_len, $tid2nexon{$tids[$i]}->{"intron_length"});		
		 
	 }
	 my $maxindex;
	 my $max = $vs_ex_num[0];
	 for(my $i = 0; $i <= $#vs_ex_num; $i++)
	 {
	     if($vs_ex_num[$i] > $max)
		 {
		    $maxindex=$i;		 
		 }	 
	 }
	 print OUT $key."\t".$vs_in_num[$maxindex]."\n";
}
close(OUT);


open(OUT, ">gene_intron_distribution2.xls");
print OUT "gene_id\tintron_length\n";
foreach my $key (keys %gene2tr)
{
     my @tids = @{$gene2tr{$key}->{"tids"}};
	 my @g_exsts = @{$gene2tr{$key}->{"exsts"}};
	 my @g_exens = @{$gene2tr{$key}->{"exens"}};
	 my $g_st = min(@g_exsts);
	 my $g_en = max(@g_exens);
	 my $chr = $gene2tr{$key}->{"chr"};
	 my @vs_ex_num;my @vs_in_num;my @vs_ex_len;my @vs_in_len;
	 for(my $i = 0; $i <= $#tids; $i++)
	 {
	    push(@vs_ex_num, $tid2nexon{$tids[$i]}->{"exon_number"});
	    push(@vs_in_num, $tid2nexon{$tids[$i]}->{"intron_number"});		
	    #push(@vs_ex_len, $tid2nexon{$tids[$i]}->{"exon_length"});		
	    push(@vs_in_len, $tid2nexon{$tids[$i]}->{"intron_length"});		
		 
	 }
	 my $maxindex;
	 my $max = $vs_ex_num[0];
	 for(my $i = 0; $i <= $#vs_ex_num; $i++)
	 {
	     if($vs_ex_num[$i] > $max)
		 {
		    $maxindex=$i;		 
		 }	 
	 }
	 my $len = $vs_in_len[$maxindex];
	 my @tmp = split(/\,/, $len);
	 for(my $i = 0; $i <= $#tmp; $i++)
	 {
	    print OUT $key."\t".$tmp[$i]."\n";
	 }
}
close(OUT);

