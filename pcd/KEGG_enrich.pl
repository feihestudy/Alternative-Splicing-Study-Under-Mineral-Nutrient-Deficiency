use strict;
use List::Uniq ':all';
use Getopt::Long;
use Statistics::R;
use vars qw($blFile $glFile $symbol2KEGGfile $ispv);
Getopt::Long::GetOptions(
    'bl=s'    =>\$blFile,
    'gl=s'    => \$glFile,
    'symbol2KEGG=s'    => \$symbol2KEGGfile,
    'ispvalue=s' => \$ispv,
);
my %blsymbol;
my %allsymbol2KEGG;
my %desymbol2KEGG;
my %KEGG2allsymbol;
my %KEGG2desymbol;
open(IN,  "$blFile");
while(<IN>)
{
    chomp;
	s/\s+$//;
	$blsymbol{$_} =1;
}
close(IN);
##### All_Symol2KEGG
my $KEGGfile=$symbol2KEGGfile;
open(IN, $KEGGfile);
while(<IN>)
{
	chomp;
	my @arr = split(/\t/, $_);
	#next if(not defined $blsymbol{$arr[0]});
	my $keggid = shift(@arr);
	my @xtmp = split(/\/\//, $keggid);
	next if ($xtmp[0] eq "hsa01100");
	$keggid=~s/\/\//\t/g;
	for(my $i = 0; $i <= $#arr; $i++)
	{
	if(defined $blsymbol{$arr[$i]})
	{
        push(@{$allsymbol2KEGG{$arr[$i]}}, $keggid);
	    push(@{$KEGG2allsymbol{$keggid}}, $arr[$i]);
	}
	}
}
close(IN);
foreach my $key (keys %KEGG2allsymbol)
{
    @{$KEGG2allsymbol{$key}} = uniq(@{$KEGG2allsymbol{$key}});
	#print $key."\t".join(",", @{$KEGG2allsymbol{$key}})."\n";
}
#### My_symbol2KEGG
open(IN, $glFile);
while(<IN>)
{
    chomp;
	s/\s+$//;
    my $id = $_;
	print $id."\n";
    if(defined $allsymbol2KEGG{$id})
    {
        $desymbol2KEGG{$id} = 1;
        my @tmp = @{$allsymbol2KEGG{$id}};
        for(my $i = 0; $i <= $#tmp; $i++)
        {
            push(@{$KEGG2desymbol{$tmp[$i]}}, $id);
        }
    }
}
foreach my $key (keys %KEGG2desymbol)
{
    @{$KEGG2desymbol{$key}} = uniq(@{$KEGG2desymbol{$key}});
	print $key."\t".join(",", @{$KEGG2desymbol{$key}})."\n";
}


#=cut;
##### Run GO enrichment analysis #######
my $otFile = "KEGG_enrichment.xls";
open(OUT, ">$otFile");
print OUT "Pathway_id"."\t"."Pathway_name"."\t"."Pathway_class"."\t"."Annotated"."\t"."Significant"."\t"."Expected"."\t"."Geneids"."\t"."Pvalue_by_Hypergeometric_test"."\t"."FDR"."\n";



my $oddratio;
print "Run KEGG enrichment analysis"."\n";
my $num_symbols_in_blist = scalar(keys %allsymbol2KEGG);
my $num_symbols_in_list = scalar(keys %desymbol2KEGG);

my $R = Statistics::R->new();
$R->run(q'pvalues <- NULL');
my @tmp2;
foreach my $key (keys %KEGG2desymbol)
{
    #print $key."\n";
    my @x1 = @{$KEGG2desymbol{$key}};
	my $num_symbols_in_list_in_KEGG = scalar(@x1);
	my $num_symbols_in_blist_in_KEGG = scalar(@{$KEGG2allsymbol{$key}});
	my $num_symbols_in_blist_not_KEGG =  $num_symbols_in_blist - $num_symbols_in_blist_in_KEGG;
	my $countratio = $num_symbols_in_list_in_KEGG / $num_symbols_in_list;
    $R->set('x', $num_symbols_in_list_in_KEGG);
    $R->set('k', $num_symbols_in_list);
    $R->set('m', $num_symbols_in_blist_in_KEGG);
    $R->set('n', $num_symbols_in_blist_not_KEGG);
    $R->run(q'y <- phyper(x-1, m, n, k, lower.tail=F)');
    $R->run(q'pvalues <- c(pvalues, y)');
    #my $oddratio = ($num_symbols_in_list_in_KEGG / $num_symbols_in_list)/($num_symbols_in_blist_in_KEGG / $num_symbols_in_blist);
	my $expectnum = $num_symbols_in_blist_in_KEGG*($num_symbols_in_list / $num_symbols_in_blist);
	my @syms = @x1;
	if(scalar(@syms) >= 500)
	{
	   @syms = @syms[0..500];
	}
    push(@tmp2, $key."\t".$num_symbols_in_blist_in_KEGG."\t".$num_symbols_in_list_in_KEGG."\t".$expectnum."\t".join(",", @syms));

}
close(TM);
my $one = "q <- p.adjust(pvalues, method='BH')";
$R->run($one);
my $pv_ref = $R->get('pvalues');
my $qv_ref = $R->get('q');

my @pvs;
my @qvs;
if(scalar(keys %KEGG2desymbol) == 1)
{
   push(@pvs, $pv_ref);
   push(@qvs, $qv_ref);
}
else
{
	@pvs = @{$pv_ref};
	@qvs = @{$qv_ref};
}
for(my $i = 0; $i <= $#tmp2; $i++)
{
	$tmp2[$i] = $tmp2[$i]."\t".$pvs[$i]."\t".$qvs[$i];
}
my $nindex=8;
@tmp2 = sort { (split "\t",$a)[8] <=> (split "\t",$b)[8]} @tmp2;
print OUT join("\n", @tmp2)."\n";
close(OUT);



if($ispv eq "FDR")
{
open(IN, "KEGG_enrichment.xls");
my $otfile = "KEGG_enrichment_significant.xls";
open(OUT, ">$otfile");
my $tt = <IN>;
print OUT $tt;
while(<IN>)
{
    chomp;
	my @arr = split(/\t/, $_);
	if($arr[$#arr] <= 0.05)
	{
	    print OUT $_."\n";
	}
}
close(OUT);
}

if($ispv eq "Pvalue")
{
open(IN, "KEGG_enrichment.xls");
my $otfile = "KEGG_enrichment_significant.xls";
open(OUT, ">$otfile");
my $tt = <IN>;
print OUT $tt;
while(<IN>)
{
    chomp;
	my @arr = split(/\t/, $_);
	if($arr[($#arr-1)] <= 0.05)
	{
	    print OUT $_."\n";
	}
}
close(OUT);
}

	






