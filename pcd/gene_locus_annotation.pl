use List::Uniq ':all';
use List::Util qw(first max maxstr min minstr reduce shuffle sum);
use Getopt::Long;
use vars qw($gtffile $bedfile $tmapfile);
Getopt::Long::GetOptions(
    'gtf=s'    =>\$gtffile,
	'tmap=s'   =>\$tmapfile,
);
my %knowngid2tid;
open(IN, "/home/database/RAPDBMSU_mRNA_nx.gtf");
while(<IN>)
{
    chomp;
	s/\s+$//;
	my @arr = split(/\t/, $_);
	my @tmp = split(/\;\s+/,$arr[8]);
	my $gid = $tmp[0];
	my $tid = $tmp[1];
	$gid=~s/\;$//;
	$gid=~s/^gene_id\s+//;
	$gid=~s/\"//g;
	$tid=~s/\;$//;
	$tid=~s/^transcript_id\s+//;
	$tid=~s/\"//g;		
	push(@{$knowngid2tid{$gid}}, $tid);
}
close(IN);

my %knowntid2gid;
open(IN, "/home/database/RAPDBMSU_mRNA.gtf");
while(<IN>)
{
    chomp;
	s/\s+$//;
	my @arr = split(/\t/, $_);
	my @tmp = split(/\;\s+/,$arr[8]);
	my $gid = $tmp[0];
	my $tid = $tmp[1];
	$gid=~s/\;$//;
	$gid=~s/^gene_id\s+//;
	$gid=~s/\"//g;
	$tid=~s/\;$//;
	$tid=~s/^transcript_id\s+//;
	$tid=~s/\"//g;		
    $knowntid2gid{$tid}=$gid;
}
close(IN);

my $gid2gid;
open(IN,$tmapfile);
while(<IN>)
{
    chomp;
	s/\s+$//;
	my @arr = split(/\t/, $_);
	if($arr[2] eq "j" or $arr[2] eq "o" or $arr[2] eq "=")
	{
		push(@{$gid2gid{$arr[3]}}, $arr[0]);
	}
}
close(IN);

foreach my $key (keys %knowngid2tid)
{
	@{$knowngid2tid{$key}}=uniq(@{$knowngid2tid{$key}});
}

foreach my $key (keys %gid2gid)
{
	@{$gid2gid{$key}}=uniq(@{$gid2gid{$key}});
}
my %tid2gid;
my %gtf2exon;
open(IN, $gtffile);
while(<IN>)
{
    chomp;
	s/\s+$//;
	my @arr = split(/\t/, $_);
	my @tmp = split(/\;\s+/,$arr[8]);
	my $gid = $tmp[0];
	my $tid = $tmp[1];
	$gid=~s/\;$//;
	$gid=~s/^gene_id\s+//;
	$gid=~s/\"//g;
	$tid=~s/\;$//;
	$tid=~s/^transcript_id\s+//;
	$tid=~s/\"//g;		
	if(defined $gid2gid{$gid})
	{
		my @knowgids = @{$gid2gid{$gid}};
        for(my $i = 0; $i <= $#knowgids; $i++)
        {
		    if(defined $knowngid2tid{$knowgids[$i]})
			{
				push(@{$tid2gid{$gid}}, @{$knowngid2tid{$knowgids[$i]}});
			}
		}	
	}
	else
	{
	
		push(@{$tid2gid{$gid}}, $tid);
	}
	if($arr[2] eq "exon")
	{
		push(@{$gtf2exon{$gid} -> {exst}}, $arr[3]);
		push(@{$gtf2exon{$gid} -> {exen}}, $arr[4]);
		$gtf2exon{$gid}->{chr}=$arr[0];
		$gtf2exon{$gid}->{strand}=$arr[6];
	}	
}
close(IN);



foreach my $key (keys %gtf2exon)
{
	my $exst_min = min(@{$gtf2exon{$key} -> {exst}});
	my $exen_max = max(@{$gtf2exon{$key} -> {exen}});
	my $chr =$gtf2exon{$key}->{chr};
	my $strand = $gtf2exon{$key}->{strand};
	$gtf2exon{$key}->{locus}=$chr.":".$exst_min."-".$exen_max."(".$strand.")";
}
	



#### MSU Gene information
my %gid2infor_msu;
open(IN, "/home/database/all.locus_brief_info.7.0");
<IN>;
while(<IN>)
{
     chomp;
	 s/\s+$//;
	 s/\'//g;
	 my @arr = split(/\t/, $_);
	 $gid2infor_msu{$arr[1]}=$arr[9];
}
close(IN);




#### RAPDB Gene information
my %gid2infor1;
my $file = "/home/database/repre_transcripts.gff";
open(IN, $file);
#my %go2osaid;
while(<IN>)
{
     chomp;
	 s/\'//g;
	 my @arr = split(/\t/, $_);
	 if($arr[8]=~m/^ID\=/)
	 {
	     #print $arr[8]."\n";
		 my @tmp = split(/\;/, $arr[8]);
		 my $rapdbid;my $symbol;my $goid;my $note;my $uniprotid;my $keggid;my $InterProid;my $genename;
		 for(my $i = 0; $i <= $#tmp; $i++)
		 {
			if($tmp[$i]=~m/^Locus\_id\=/)
			{
				$rapdbid = $tmp[$i];
			}
            if($tmp[$i]=~m/^Note\=/)
		    {
			    $note = $tmp[$i];
			}
            if($tmp[$i]=~m/^GO\=/)
		    {
			    $goid = $tmp[$i];
			}
            if($tmp[$i]=~m/^CGSNL\s+Gene\s+Symbol\=/)
		    {
			    $symbol = $tmp[$i];

			}
            if($tmp[$i]=~m/^KEGG\=/)
		    {
			    $keggid = $tmp[$i];
			}
            if($tmp[$i]=~m/^InterPro\=/)
		    {
			    $InterProid = $tmp[$i];
			}
	        if($tmp[$i]=~m/^ORF_evidence\=/)
		    {
			    $uniprotid = $tmp[$i];
			}
	        if($tmp[$i]=~m/^CGSNL\s+Gene\s+Name\=/)
		    {
			    $genename = $tmp[$i];

			}			
        }
		$rapdbid=~s/^Locus\_id\=//;
		$symbol=~s/^CGSNL\s+Gene\s+Symbol\=//;
		$genename=~s/^CGSNL\s+Gene\s+Name\=//;
		$note=~s/^Note\=//;
		
		next if($rapdbid eq "");
		if($symbol ne "")
		{
				if($symbol eq "_" or $symbol eq "-")
				{
				    $symbol="NULL";
				}
		$gid2infor1{$rapdbid}->{symbol}=$symbol;

		}
		else
		{
		     $gid2infor1{$rapdbid}->{symbol}="NULL";		
		}
		if($note ne "")
		{
		     $gid2infor1{$rapdbid}->{note}=$note;
		}
		else
		{
		     $gid2infor1{$rapdbid}->{note}="NULL";		
		}
		if($goid ne "")
		{
		     $gid2infor1{$rapdbid}->{GO}=$goid;
		}
		else
		{
		     $gid2infor1{$rapdbid}->{GO}="NULL";		
		}
		if($keggid ne "")
		{
		     $gid2infor1{$rapdbid}->{KEGG}=$keggid;
		}
		else
		{
		     $gid2infor1{$rapdbid}->{KEGG}="NULL";		
		}
		if($InterProid ne "")
		{
		     $gid2infor1{$rapdbid}->{InterPro}=$InterProid;
		}
		else
		{
		     $gid2infor1{$rapdbid}->{InterPro}="NULL";		
		}
		if($uniprotid ne "")
		{
		     $gid2infor1{$rapdbid}->{uniprot}=$uniprotid;
		}
		else
		{
		     $gid2infor1{$rapdbid}->{uniprot}="NULL";		
		}
		if($genename ne "")
		{
				if($genename eq "_" or $genename eq "-")
				{
				    $genename="NULL";
				}
		     $gid2infor1{$rapdbid}->{genename}=$genename;
		}
		else
		{
		     $gid2infor1{$rapdbid}->{genename}="NULL";		
		}
	    #print $rapdbid."\t".$symbol."\t".$note."\t".$goid."\n";
    }		
}
close(IN);








#my %gid2infor2;
my $file = "/home/database/predicted_transcripts.gff";
open(IN, $file);
#my %go2osaid;
while(<IN>)
{
     chomp;
	 s/\'//g;
	 my @arr = split(/\t/, $_);
	 if($arr[8]=~m/^ID\=/)
	 {
	     #print $arr[8]."\n";
		 my @tmp = split(/\;/, $arr[8]);
		 my $rapdbid;my $symbol;my $goid;my $note;my $uniprotid;my $keggid;my $InterProid;my $genename;
		 for(my $i = 0; $i <= $#tmp; $i++)
		 {
			if($tmp[$i]=~m/^Locus\_id\=/)
			{
				$rapdbid = $tmp[$i];
			}
            if($tmp[$i]=~m/^Note\=/)
		    {
			    $note = $tmp[$i];
			}
            if($tmp[$i]=~m/^GO\=/)
		    {
			    $goid = $tmp[$i];
			}
            if($tmp[$i]=~m/^CGSNL\s+Gene\s+Symbol\=/)
		    {
			    $symbol = $tmp[$i];
			}
            if($tmp[$i]=~m/^KEGG\=/)
		    {
			    $keggid = $tmp[$i];
			}
            if($tmp[$i]=~m/^InterPro\=/)
		    {
			    $InterProid = $tmp[$i];
			}
	        if($tmp[$i]=~m/^ORF_evidence\=/)
		    {
			    $uniprotid = $tmp[$i];
			}
	        if($tmp[$i]=~m/^CGSNL\s+Gene\s+Name\=/)
		    {
			    $genename = $tmp[$i];

			}			
        }
		$rapdbid=~s/^Locus\_id\=//;
		$symbol=~s/^CGSNL\s+Gene\s+Symbol\=//;
		$genename=~s/^CGSNL\s+Gene\s+Name\=//;
		$note=~s/^Note\=//;
		
		next if($rapdbid eq "");
		if($symbol ne "")
		{
				if($symbol eq "_" or $symbol eq "-")
				{
				    $symbol="NULL";
				}
		     $gid2infor1{$rapdbid}->{symbol}=$symbol;
		}
		else
		{
		     $gid2infor1{$rapdbid}->{symbol}="NULL";		
		}
		if($note ne "")
		{
		     $gid2infor1{$rapdbid}->{note}=$note;
		}
		else
		{
		     $gid2infor1{$rapdbid}->{note}="NULL";		
		}
		if($goid ne "")
		{
		     $gid2infor1{$rapdbid}->{GO}=$goid;
		}
		else
		{
		     $gid2infor1{$rapdbid}->{GO}="NULL";		
		}
		if($keggid ne "")
		{
		     $gid2infor1{$rapdbid}->{KEGG}=$keggid;
		}
		else
		{
		     $gid2infor1{$rapdbid}->{KEGG}="NULL";		
		}
		if($InterProid ne "")
		{
		     $gid2infor1{$rapdbid}->{InterPro}=$InterProid;
		}
		else
		{
		     $gid2infor1{$rapdbid}->{InterPro}="NULL";		
		}
		if($uniprotid ne "")
		{
		     $gid2infor1{$rapdbid}->{uniprot}=$uniprotid;
		}
		else
		{
		     $gid2infor1{$rapdbid}->{uniprot}="NULL";		
		}
		if($genename ne "")
		{
				if($genename eq "_" or $genename eq "-")
				{
				    $genename="NULL";
				}
		     $gid2infor1{$rapdbid}->{genename}=$genename;
		}
		else
		{
		     $gid2infor1{$rapdbid}->{genename}="NULL";		
		}
	    #print $rapdbid."\t".$symbol."\t".$note."\t".$goid."\n";
    }		
}
close(IN);



print "Gene_Locus\tLocus_infor\tRAPDBid\tMSUid\tRAPDB_symbol\tRAPDB_genename\tRAPDB_function_description\tMSU_function_description"."\n";
foreach my $key (keys %tid2gid)
{
    
	@{$tid2gid{$key}}=uniq(@{$tid2gid{$key}});
	#print $key."\t".join(",", @{$tid2gid{$key}})."\n";
	my @tmp = @{$tid2gid{$key}};
	my @msu;
	my @rapdb;
	for(my $i = 0; $i <= $#tmp; $i++)
	{
	    if(defined $knowntid2gid{$tmp[$i]})
		{
			if($knowntid2gid{$tmp[$i]}=~m/^LOC/)
			{
			   push(@msu, $knowntid2gid{$tmp[$i]});
			}
			if($knowntid2gid{$tmp[$i]}=~m/^Os/)
			{
			   push(@rapdb, $knowntid2gid{$tmp[$i]});
			}
			
		}
	}
	@msu = uniq(@msu);	@rapdb = uniq(@rapdb);
	my @msudesc;
	my @rapdbdesc;
	for(my $i = 0; $i <= $#msu; $i++)
	{
	    push(@msudesc, $gid2infor_msu{$msu[$i]});	
	}
	for(my $i = 0; $i <= $#rapdb; $i++)
	{
	    push(@rapdbdesc, $gid2infor1{$rapdb[$i]}->{"note"});	
	}

	
	my @rapdbsymbol;
	for(my $i = 0; $i <= $#rapdb; $i++)
	{
	    if($gid2infor1{$rapdb[$i]}->{"symbol"} eq "LTN1")
		{
	    push(@rapdbsymbol, "PHO2");	
		}
		else
		{
	    push(@rapdbsymbol, $gid2infor1{$rapdb[$i]}->{"symbol"});
		}
	}


	my @rapdbgenename;
	for(my $i = 0; $i <= $#rapdb; $i++)
	{
	    push(@rapdbgenename, $gid2infor1{$rapdb[$i]}->{"genename"});	
	}	
	my $msulist="NULL";
	my $rapdblist="NULL";
    if(scalar(@msu) > 0)
    {
	    $msulist=join("//", @msu);
	}
    if(scalar(@rapdb) > 0)
    {
	    $rapdblist=join("//", @rapdb);
	}

	my $msudesclist="NULL";
	my $rapdbdesclist="NULL";
    if(scalar(@msudesc) > 0)
    {
	    $msudesclist=join("//", @msudesc);
	}
    if(scalar(@rapdbdesc) > 0)
    {
	    $rapdbdesclist=join("//", @rapdbdesc);
	}	

	my $msusymbollist="NULL";
	my $rapdbsymbollist="NULL";
    if(scalar(@rapdbsymbol) > 0)
    {
	    $rapdbsymbollist=join("//", @rapdbsymbol);
	}	
	
	
	my $msugenenamelist="NULL";
	my $rapdbgenenamelist="NULL";
    if(scalar(@rapdbgenename) > 0)
    {
	    $rapdbgenenamelist=join("//", @rapdbgenename);
	}	
	print $key."\t".$gtf2exon{$key}->{locus}."\t".$rapdblist."\t".$msulist."\t".$rapdbsymbollist."\t".$rapdbgenenamelist."\t".$rapdbdesclist."\t".$msudesclist."\n";
}





