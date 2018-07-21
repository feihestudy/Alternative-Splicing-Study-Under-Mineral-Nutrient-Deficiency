use strict;
use List::Uniq ':all';
use Getopt::Long;
use Statistics::R;
use vars qw($blFile $glFile $symbol2GOfile $ispv);
Getopt::Long::GetOptions(
    'bl=s'    =>\$blFile,
    'gl=s'    => \$glFile,
	'symbol2GO=s' => \$symbol2GOfile,
    'ispvalue=s' => \$ispv,
);
my %blsymbol;
my %allsymbol2go;
my %desymbol2go;
my %go2allsymbol;
my %bp2allsymbol;
my %mf2allsymbol;
my %cc2allsymbol;
my %go2desymbol;
my %bp2desymbol;
my %mf2desymbol;
my %cc2desymbol;
open(IN,  "$blFile");
while(<IN>)
{
    chomp;
	s/\s+$//;
	s/^\s+//;
	$blsymbol{$_} =1;
	#print $_."\n";
}
close(IN);
##### All_Symol2GO
#my $GOfile = "F:/database/Oryza_sativa/symbol2GO.txt";
my $GOfile = $symbol2GOfile;
open(IN, $GOfile);
while(<IN>)
{
	chomp;
	#print $_."\n";
	my @arr = split(/\t/, $_);
	my $goid = shift(@arr);
	$goid=~s/\/\//\t/g;
	for(my $i = 0; $i <= $#arr; $i++)
	{
	$arr[$i]=~s/\s+$//;
	if(defined $blsymbol{$arr[$i]})
	{
        push(@{$allsymbol2go{$arr[$i]}}, $goid);
		#print $goid."\n";
	    push(@{$go2allsymbol{$goid}}, $arr[$i]);
	}
	}
}
close(IN);
foreach my $key (keys %go2allsymbol)
{
    @{$go2allsymbol{$key}} = uniq(@{$go2allsymbol{$key}});
	#print $key."\t".join(",", @{$go2allsymbol{$key}})."\n";
}
foreach my $key (keys %go2allsymbol)
{
   my @tmp = @{$go2allsymbol{$key}};
   my @tmp2 = split(/\t/, $key);
   if($tmp2[2] eq "biological_process")
   {
       push(@{$bp2allsymbol{$tmp2[2]}}, @tmp);
   }
   if($tmp2[2] eq "molecular_function")
   {
       push(@{$mf2allsymbol{$tmp2[2]}}, @tmp);
   }
   if($tmp2[2] eq "cellular_component")
   {
       push(@{$cc2allsymbol{$tmp2[2]}}, @tmp);
   }
}
@{$bp2allsymbol{biological_process}} = uniq(@{$bp2allsymbol{biological_process}});
@{$mf2allsymbol{molecular_function}} = uniq(@{$mf2allsymbol{molecular_function}});
@{$cc2allsymbol{cellular_component}} = uniq(@{$cc2allsymbol{cellular_component}});

#### My_symbol2GO
open(IN, $glFile);
#<IN>;
while(<IN>)
{
    chomp;
	s/\s+$//;
    my $id = $_;
    if(defined $allsymbol2go{$id})
    {
        $desymbol2go{$id} = 1;
        my @tmp = @{$allsymbol2go{$id}};
        for(my $i = 0; $i <= $#tmp; $i++)
        {
            push(@{$go2desymbol{$tmp[$i]}}, $id);
        }
    }
}
foreach my $key (keys %go2desymbol)
{
    @{$go2desymbol{$key}} = uniq(@{$go2desymbol{$key}});
}
foreach my $key (keys %go2desymbol)
{
   my @tmp = @{$go2desymbol{$key}};
   my @tmp2 = split(/\t/, $key);
   if($tmp2[2] eq "biological_process")
   {
       push(@{$bp2desymbol{$tmp2[2]}}, @tmp);
   }
   if($tmp2[2] eq "molecular_function")
   {
       push(@{$mf2desymbol{$tmp2[2]}}, @tmp);
   }
   if($tmp2[2] eq "cellular_component")
   {
       push(@{$cc2desymbol{$tmp2[2]}}, @tmp);
   }
}
@{$bp2desymbol{biological_process}} = uniq(@{$bp2desymbol{biological_process}});
@{$mf2desymbol{molecular_function}} = uniq(@{$mf2desymbol{molecular_function}});
@{$cc2desymbol{cellular_component}} = uniq(@{$cc2desymbol{cellular_component}});

##### Run GO enrichment analysis #######
my $otFile = "GO_enrichment.xls";
open(OUT, ">$otFile");
print OUT "GO_id"."\t"."GO_term"."\t"."GO_category"."\t"."GO_rank"."\t"."Annotated"."\t"."Significant"."\t"."Expected"."\t"."Geneids"."\t"."Pvalue_by_Hypergeometric_test"."\t"."FDR"."\n";


my $oddratio;
print "Run GO enrichment analysis"."\n";
my $num_symbols_in_bp_blist = scalar(@{$bp2allsymbol{biological_process}});
my $num_symbols_in_mf_blist = scalar(@{$mf2allsymbol{molecular_function}});
my $num_symbols_in_cc_blist = scalar(@{$cc2allsymbol{cellular_component}});
#my $num_symbols_in_list = scalar(keys %desymbol2go);
my $num_symbols_in_bp_list = scalar(@{$bp2desymbol{biological_process}});
my $num_symbols_in_mf_list = scalar(@{$mf2desymbol{molecular_function}});
my $num_symbols_in_cc_list = scalar(@{$cc2desymbol{cellular_component}});

print $num_symbols_in_bp_blist."\n";
print $num_symbols_in_mf_blist."\n";
print $num_symbols_in_cc_blist."\n";
print $num_symbols_in_bp_list."\n";
print $num_symbols_in_mf_list."\n";
print $num_symbols_in_cc_list."\n";

#=cut;
my $R = Statistics::R->new();
$R->run(q'pvalues <- NULL');
my @tmp2;
foreach my $key (keys %go2desymbol)
{
    #print $key."\n";
	my @tmpx = split(/\t/, $key);
	my $cat = $tmpx[2];
	my $num_symbols_in_blist2;
	my $num_symbols_in_list2;
	my $num_symbols_in_blist_in_GO;
	my $num_symbols_in_blist_not_GO;
	if($cat eq "biological_process")
	{
	    $num_symbols_in_list2 = $num_symbols_in_bp_list;
	    $num_symbols_in_blist2 = $num_symbols_in_bp_blist;
		$num_symbols_in_blist_in_GO = scalar(@{$go2allsymbol{$key}});
        $num_symbols_in_blist_not_GO = $num_symbols_in_blist2 - $num_symbols_in_blist_in_GO;

	}
	if($cat eq "molecular_function")
	{
	    $num_symbols_in_list2 = $num_symbols_in_mf_list;
	    $num_symbols_in_blist2 = $num_symbols_in_mf_blist;
		$num_symbols_in_blist_in_GO = @{$go2allsymbol{$key}};
        $num_symbols_in_blist_not_GO = $num_symbols_in_blist2 - $num_symbols_in_blist_in_GO;
}
	if($cat eq "cellular_component")
	{
	    $num_symbols_in_list2 = $num_symbols_in_cc_list;
	    $num_symbols_in_blist2 = $num_symbols_in_cc_blist;
		$num_symbols_in_blist_in_GO = @{$go2allsymbol{$key}};
        $num_symbols_in_blist_not_GO = $num_symbols_in_blist2 - $num_symbols_in_blist_in_GO;
	}
    my @x1 = @{$go2desymbol{$key}};
	my $num_symbols_in_list_in_GO = scalar(@x1);
	#next 
	my $countratio = $num_symbols_in_list_in_GO / $num_symbols_in_list2;
    $R->set('x', $num_symbols_in_list_in_GO);
    $R->set('k', $num_symbols_in_list2);
    $R->set('m', $num_symbols_in_blist_in_GO);
    $R->set('n', $num_symbols_in_blist_not_GO);
    $R->run(q'y <- phyper(x-1, m, n, k, lower.tail=F)');
    $R->run(q'pvalues <- c(pvalues, y)');
    my $oddratio = ($num_symbols_in_list_in_GO / $num_symbols_in_list2)/($num_symbols_in_blist_in_GO / $num_symbols_in_blist2);
	my $expectnum = $num_symbols_in_blist_in_GO*($num_symbols_in_list2 / $num_symbols_in_blist2);
	#$goidname=~s/\/\//\t/g;
	my @syms = @x1;
	if(scalar(@syms) >= 500)
	{
	   @syms = @syms[0..500];
	}
    #push(@tmp2, $key."\t".$num_symbols_in_list_in_GO."\t".$countratio."\t".$num_symbols_in_blist_in_GO."\t".$num_symbols_in_list2."\t".$num_symbols_in_blist2."\t".join(",", @syms)."\t".$oddratio);
    #push(@tmp2, $key."\t".$g_all_num."\t".$ref_all_num."\t".$g_num_in_path."\t".join(",",@syms)."\t".$ref_num_in_path."\t".$oddratio);
    push(@tmp2, $key."\t".$num_symbols_in_blist_in_GO."\t".$num_symbols_in_list_in_GO."\t".$expectnum."\t".join(",", @syms));

}
close(TM);

my $one = "q <- p.adjust(pvalues, method='BH')";
$R->run($one);
my $pv_ref = $R->get('pvalues');
my $qv_ref = $R->get('q');

my @pvs;
my @qvs;
if(scalar(keys %go2desymbol) == 1)
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
open(IN, "GO_enrichment.xls");
my $otfile = "GO_enrichment_significant.xls";
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
open(IN, "GO_enrichment.xls");
my $otfile = "GO_enrichment_significant.xls";
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






