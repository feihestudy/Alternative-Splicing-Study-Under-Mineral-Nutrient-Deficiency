use strict;
use Getopt::Long;
use vars qw($ase $genomefile);
Getopt::Long::GetOptions(
    'ASE=s' => \$ase,
    'genome=s' => \$genomefile,
);
$/=">";
open(IN, $genomefile);
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
    $chr2seq{$id} = $seq;
}
close(IN);


$/="\n";
my $tt_len = 0;
my $au_len = 0;
open(IN, $ase);
<IN>;
while(<IN>)
{
   chomp;
   s/\s+$//;
   my @arr3 = split(/\t/, $_);
   my @tmp = split(/\&/, $arr3[0]);
   if($tmp[0] eq "SE")
   {
	  my $st0 = $tmp[6] + 1;
	  my $en0 = $tmp[7];
	  my $st2 = $tmp[8] + 1;
	  my $en2 = $tmp[9];	  
	  my $st1 = $tmp[4] + 1;
	  my $en1 = $tmp[5];
	  my $chr = $tmp[2];
	  my $strand = $tmp[3];
	  if($strand eq "+")
	  {
			my $st = $en0 + 1 - 1;
			my $len = 2;
			my $seq5 = substr($chr2seq{$chr}, $st, $len);
			$seq5=uc($seq5);
			my $st = $st1-2-1;
			my $len = 2;
			my $seq3 = substr($chr2seq{$chr}, $st, $len);
			$seq3=uc($seq3);
            print $arr3[0]."\t".$seq5."\t".$seq3."\n";
			
			my $st = $en1 + 1 - 1;
			my $len = 2;
			my $seq5 = substr($chr2seq{$chr}, $st, $len);
			$seq5=uc($seq5);
			my $st = $st2-2-1;
			my $len = 2;
			my $seq3 = substr($chr2seq{$chr}, $st, $len);
			$seq3=uc($seq3);
            print $arr3[0]."\t".$seq5."\t".$seq3."\n";
			
	  }
	  if($strand eq "-")
	  {
			my $st5 = $st2 - 2 - 1;
			my $len = 2;
			my $seq5 = substr($chr2seq{$chr}, $st5, $len);
            my $revcomp = reverse($seq5);
            $revcomp =~ tr/ACGTacgt/TGCAtgca/;
            $seq5 = $revcomp;
			$seq5=uc($seq5);
			my $st3 = $en1+1-1;
			my $len = 2;
			my $seq3 = substr($chr2seq{$chr}, $st3, $len);
            my $revcomp = reverse($seq3);
            $revcomp =~ tr/ACGTacgt/TGCAtgca/;
            $seq3 = $revcomp;
			$seq3=uc($seq3);
            print $arr3[0]."\t".$seq5."\t".$seq3."\n";

			my $st5 = $st1 - 2 - 1;
			my $len = 2;
			my $seq5 = substr($chr2seq{$chr}, $st5, $len);
            my $revcomp = reverse($seq5);
            $revcomp =~ tr/ACGTacgt/TGCAtgca/;
            $seq5 = $revcomp;
			$seq5=uc($seq5);
			my $st3 = $en0+1-1;
			my $len = 2;
			my $seq3 = substr($chr2seq{$chr}, $st3, $len);
            my $revcomp = reverse($seq3);
            $revcomp =~ tr/ACGTacgt/TGCAtgca/;
            $seq3 = $revcomp;
			$seq3=uc($seq3);
            print $arr3[0]."\t".$seq5."\t".$seq3."\n";
	  }	  

    }

 }
 close(IN);
