use strict;
use Getopt::Long;
use vars qw($database $gtffile);
Getopt::Long::GetOptions(
    'GTF=s' => \$gtffile,
	'database=s'    => \$database,
);
$/=">";
open(IN, $database);
<IN>;
my %chr2seq;
while(<IN>)
{
    my @arr = split(/\n/, $_);
    my $id = shift(@arr);
    $id =~s/^>//;
	my @tmp = split(/\s+/, $id);
	$id=$tmp[0];
    my $seq = join("", @arr);
    $seq =~s/>$//;
	$seq=uc($seq);
    #print $id."\n".length($seq)."\n";
    $chr2seq{$id} = $seq;
}
close(IN);

$/="\n";
open(IN, $gtffile);
my %transcript;
while(<IN>)
{
   chomp;
   my @arr = split(/\t/, $_);
   #print $arr[8]."\n";
   my @tmp = split(/;\s+/, $arr[8]);
   my $gid = $tmp[0];
   $gid=~s/gene\_id\s+//;
   $gid=~s/\"//g;
   $gid=~s/\;//g;
   my $tid = $tmp[1];
   $tid=~s/transcript\_id\s+//;
   $tid=~s/\s+$//;
   $tid=~s/\;$//;   
   $tid=~s/\"//g;
   my $exon = $arr[2];
   $transcript{$tid} -> {strand} = $arr[6];
   $transcript{$tid} -> {chr} = $arr[0];
   if($exon eq "exon")
   {
       push(@{$transcript{$tid} -> {exst}}, $arr[3]);
       push(@{$transcript{$tid} -> {exen}}, $arr[4]);
   }	   
}
close(IN);

foreach my $key (keys %transcript)
{
	my $strand = $transcript{$key}->{strand};
	my @sts = @{$transcript{$key}->{exst}};
	#@sts = sort {$a <=> $b} @sts;
	my @ens = @{$transcript{$key}->{exen}};
	#@ens = sort {$a <=> $b} @ens;	
	my @sorted_exst = sort {$sts[$a] <=> $sts[$b]} 0..$#sts;
	@sts = @sts[@sorted_exst];
	@ens = @ens[@sorted_exst];
	my $chr = $transcript{$key}->{chr};
    if($strand eq "+")
    {
        my $tseq = '';
        for(my $i = 0; $i <= $#sts; $i++)
        {
		    my $len = $ens[$i] - $sts[$i] + 1;
			my $st = $sts[$i] - 1;
            my $seq = substr($chr2seq{$chr}, $st, $len);
            $tseq = $tseq.$seq;
        }
		next if(length($tseq) eq 0);		
        print ">".$key."\n".$tseq."\n";
    }
    else
    {
        if($strand eq "-")
        {
            my $tseq = '';
            for(my $i = 0; $i <= $#sts; $i++)
            {
                my $st = $sts[$i] - 1;
                my $len = $ens[$i] - $sts[$i] + 1;
                my $seq = substr($chr2seq{$chr}, $st, $len);
                $tseq = $tseq.$seq;
            }
			next if(length($tseq) eq 0);
            my $revcomp = reverse($tseq);
            $revcomp =~ tr/ACGTacgt/TGCAtgca/;
            $tseq = $revcomp;
            print ">".$key."\n".$tseq."\n";
        }
        else
        {
              next;
        }
    }

}
close(OUT);
