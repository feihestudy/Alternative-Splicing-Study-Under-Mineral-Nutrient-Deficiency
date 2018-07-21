use strict;
use Getopt::Long;
use vars qw($ASEdir $samplelistfile);
Getopt::Long::GetOptions(
    'ASEdir=s'    => \$ASEdir,
	'samplelist=s' =>\$samplelistfile,
);

my %RI;
my %A3SS;
my %A5SS;
my %SE;
my %MXE;

open(IN, $ASEdir."/"."AS.RI.txt");
while(<IN>)
{
    chomp;
	s/\s+$//;
	s/\"//g;
	my @arr = split(/\t/, $_);
	$RI{$arr[0]}=join("\t", @arr[1..$#arr]);
}
close(IN);

open(IN, $ASEdir."/"."AS.A3SS.txt");
while(<IN>)
{
    chomp;
	s/\s+$//;
	s/\"//g;
	my @arr = split(/\t/, $_);
	$A3SS{$arr[0]}=join("\t", @arr[1..$#arr]);
}
close(IN);


open(IN, $ASEdir."/"."AS.SE.txt");
while(<IN>)
{
    chomp;
	s/\s+$//;
	s/\"//g;
	my @arr = split(/\t/, $_);
	$SE{$arr[0]}=join("\t", @arr[1..$#arr]);
}
close(IN);

open(IN, $ASEdir."/"."AS.A5SS.txt");
while(<IN>)
{
    chomp;
	s/\s+$//;
	s/\"//g;
	my @arr = split(/\t/, $_);
	$A5SS{$arr[0]}=join("\t", @arr[1..$#arr]);
}
close(IN);


open(IN, $ASEdir."/"."AS.MXE.txt");
while(<IN>)
{
    chomp;
	s/\s+$//;
	s/\"//g;
	my @arr = split(/\t/, $_);
	$MXE{$arr[0]}=join("\t", @arr[1..$#arr]);
}
close(IN);



my %ASE;

my @samples;



open(IN, $samplelistfile);
while(<IN>)
{
    chomp;
	s/\s+$//;
	my @arrx = split(/\t/, $_);
	my $sfile = $arrx[0];
	my $sname = $arrx[1];
	#=cut;
	push(@samples,$sname);
	open(IN2, $sfile."/"."JC.RNASeq.RI.MATS.input.txt");
	my $tt = <IN2>;
	while(<IN2>)
	{
		chomp;
		s/\s+$//;
		my @arr = split(/\t/, $_);
		if(defined $RI{$arr[0]})
		{
			my @tmp = split(/\t/, $RI{$arr[0]});
			my $key = "RI"."&".$tmp[0]."&".join("&",@tmp[2..$#tmp]);
			#print $key."\t".$arr[1].",".$arr[2]."\n";
			$ASE{$key}->{$sname}=$arr[1].",".$arr[2];
		}
       	   
	}
	close(IN2);
	
	
	open(IN2, $sfile."/"."JC.RNASeq.A3SS.MATS.input.txt");
	my $tt = <IN2>;
	while(<IN2>)
	{
		chomp;
		s/\s+$//;
		my @arr = split(/\t/, $_);
		if(defined $A3SS{$arr[0]})
		{
			my @tmp = split(/\t/, $A3SS{$arr[0]});
			my $key = "A3SS"."&".$tmp[0]."&".join("&",@tmp[2..$#tmp]);
			#print $key."\t".$arr[1].",".$arr[2]."\n";
			$ASE{$key}->{$sname}=$arr[1].",".$arr[2];
		}
       	   
	}
	close(IN2);
	
	open(IN2, $sfile."/"."JC.RNASeq.A5SS.MATS.input.txt");
	my $tt = <IN2>;
	while(<IN2>)
	{
		chomp;
		s/\s+$//;
		my @arr = split(/\t/, $_);
		if(defined $A5SS{$arr[0]})
		{
			my @tmp = split(/\t/, $A5SS{$arr[0]});
			my $key = "A5SS"."&".$tmp[0]."&".join("&",@tmp[2..$#tmp]);
			#print $key."\t".$arr[1].",".$arr[2]."\n";
			$ASE{$key}->{$sname}=$arr[1].",".$arr[2];
		}
       	   
	}
	close(IN2);
	

	open(IN2, $sfile."/"."JC.RNASeq.SE.MATS.input.txt");
	my $tt = <IN2>;
	while(<IN2>)
	{
		chomp;
		s/\s+$//;
		my @arr = split(/\t/, $_);
		if(defined $SE{$arr[0]})
		{
			my @tmp = split(/\t/, $SE{$arr[0]});
			my $key = "SE"."&".$tmp[0]."&".join("&",@tmp[2..$#tmp]);
			#print $key."\t".$arr[1].",".$arr[2]."\n";
			$ASE{$key}->{$sname}=$arr[1].",".$arr[2];
		}
       	   
	}
	close(IN2);

	open(IN2, $sfile."/"."JC.RNASeq.MXE.MATS.input.txt");
	my $tt = <IN2>;
	while(<IN2>)
	{
		chomp;
		s/\s+$//;
		my @arr = split(/\t/, $_);
		if(defined $MXE{$arr[0]})
		{
			my @tmp = split(/\t/, $MXE{$arr[0]});
			my $key = "MXE"."&".$tmp[0]."&".join("&",@tmp[2..$#tmp]);
			#print $key."\t".$arr[1].",".$arr[2]."\n";
			$ASE{$key}->{$sname}=$arr[1].",".$arr[2];
		}
       	   
	}
	close(IN2);
}
close(IN);


print "ASE\t".join("\t", @samples)."\n";
foreach my $key (keys %ASE)
{
    my @arr = ();
	for(my $i = 0; $i <= $#samples; $i++)
	{
	   if(defined $ASE{$key}->{$samples[$i]})
	   {
	       push(@arr, $ASE{$key}->{$samples[$i]});
	   }
	   else
	   {
	      push(@arr, "NA");   
	   }
	
	}
	print $key."\t".join("\t", @arr)."\n";

}




