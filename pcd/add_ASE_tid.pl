use strict;
use Getopt::Long;
use List::Util qw(first max maxstr min minstr reduce shuffle sum);
use vars qw($fpkmfile);
Getopt::Long::GetOptions(
	'FPKM=s'    => \$fpkmfile,
);
my %keygid;
open(IN, $fpkmfile);
<IN>;
while(<IN>)
{
    chomp;
	s/s\+$//;
	my @arr = split(/\t/, $_);
	my $gid = shift(@arr);
	my $mean = sum(@arr)/scalar(@arr);
    $keygid{$gid}=$mean;	
}
close(IN);

open(OUT, ">RI2transcript_major.xls");
print OUT "ASEid"."\t"."Inclusive"."\t"."Exclusive"."\n";
open(IN, "RI2transcript.xls");
<IN>;
while(<IN>)
{
    chomp;
	s/\s+$//;
	s/\"//g;
	my @arr = split(/\t/, $_);
	print OUT $arr[0]."\t";
	if($arr[1] ne "---")
	{
	    my @tids;
		my @values;
	    my @tmp = split(/\,/, $arr[1]);
		for(my $i = 0; $i <= $#tmp; $i++)
		{
		    if(defined $keygid{$tmp[$i]})
			{
			    #print $arr[0]."\t".$tmp[$i]."\t"."Inclusive"."\t".$keygid{$tmp[$i]}."\n";
				push(@tids, $tmp[$i]);
				push(@values, $keygid{$tmp[$i]});
			}
		}
		if(scalar(@tids) > 0)
		{
		    #print $_."\n";
			my $maxvalue = $values[0];
			my $maxindex = 0;
			for(my $i = 0; $i <= $#values; $i++)
			{
			    if($maxvalue <= $values[$i])
				{
				    $maxindex = $i;
					$maxvalue = $values[$i];
				}
			}
			print OUT $tids[$maxindex]."\t";
		}
		else
		{
		    print OUT "---"."\t";	
		}
	
	}
	else
	{
			print OUT "---"."\t";	
	}
	if($arr[2] ne "---")
	{
	    my @tids;
		my @values;
	    my @tmp = split(/\,/, $arr[2]);
		for(my $i = 0; $i <= $#tmp; $i++)
		{
		    if(defined $keygid{$tmp[$i]})
			{
			    #print $arr[0]."\t".$tmp[$i]."\t"."Inclusive"."\t".$keygid{$tmp[$i]}."\n";
				push(@tids, $tmp[$i]);
				push(@values, $keygid{$tmp[$i]});
			}
		}
		if(scalar(@tids) > 0)
		{
		    #print $_."\n";
			my $maxvalue = $values[0];
			my $maxindex = 0;
			for(my $i = 0; $i <= $#values; $i++)
			{
			    if($maxvalue <= $values[$i])
				{
				    $maxindex = $i;
					$maxvalue = $values[$i];
				}
			}
			print OUT $tids[$maxindex]."\t";
		}
		else
		{
		    print OUT "---"."\t";	
		}
	
	}
	else
	{
			print OUT "---"."\t";	
	}

	print OUT "\n";
}
close(IN);
close(OUT);


open(OUT, ">A5SS2transcript_major.xls");
print OUT "ASEid"."\t"."Inclusive"."\t"."Exclusive"."\n";
open(IN, "A5SS2transcript.xls");
<IN>;
while(<IN>)
{
    chomp;
	s/\s+$//;
	s/\"//g;
	my @arr = split(/\t/, $_);
	print OUT $arr[0]."\t";
	if($arr[1] ne "---")
	{
	    my @tids;
		my @values;
	    my @tmp = split(/\,/, $arr[1]);
		for(my $i = 0; $i <= $#tmp; $i++)
		{
		    if(defined $keygid{$tmp[$i]})
			{
			    #print $arr[0]."\t".$tmp[$i]."\t"."Inclusive"."\t".$keygid{$tmp[$i]}."\n";
				push(@tids, $tmp[$i]);
				push(@values, $keygid{$tmp[$i]});
			}
		}
		if(scalar(@tids) > 0)
		{
		    #print $_."\n";
			my $maxvalue = $values[0];
			my $maxindex = 0;
			for(my $i = 0; $i <= $#values; $i++)
			{
			    if($maxvalue <= $values[$i])
				{
				    $maxindex = $i;
					$maxvalue = $values[$i];
				}
			}
			print OUT $tids[$maxindex]."\t";
		}
		else
		{
		    print OUT "---"."\t";	
		}
	
	}
	else
	{
			print OUT "---"."\t";	
	}
	if($arr[2] ne "---")
	{
	    my @tids;
		my @values;
	    my @tmp = split(/\,/, $arr[2]);
		for(my $i = 0; $i <= $#tmp; $i++)
		{
		    if(defined $keygid{$tmp[$i]})
			{
			    #print $arr[0]."\t".$tmp[$i]."\t"."Inclusive"."\t".$keygid{$tmp[$i]}."\n";
				push(@tids, $tmp[$i]);
				push(@values, $keygid{$tmp[$i]});
			}
		}
		if(scalar(@tids) > 0)
		{
		    #print $_."\n";
			my $maxvalue = $values[0];
			my $maxindex = 0;
			for(my $i = 0; $i <= $#values; $i++)
			{
			    if($maxvalue <= $values[$i])
				{
				    $maxindex = $i;
					$maxvalue = $values[$i];
				}
			}
			print OUT $tids[$maxindex]."\t";
		}
		else
		{
		    print OUT "---"."\t";	
		}
	
	}
	else
	{
			print OUT "---"."\t";	
	}

	print OUT "\n";
}
close(IN);
close(OUT);








open(OUT, ">A3SS2transcript_major.xls");
print OUT "ASEid"."\t"."Inclusive"."\t"."Exclusive"."\n";
open(IN, "A3SS2transcript.xls");
<IN>;
while(<IN>)
{
    chomp;
	s/\s+$//;
	s/\"//g;
	my @arr = split(/\t/, $_);
	print OUT $arr[0]."\t";
	if($arr[1] ne "---")
	{
	    my @tids;
		my @values;
	    my @tmp = split(/\,/, $arr[1]);
		for(my $i = 0; $i <= $#tmp; $i++)
		{
		    if(defined $keygid{$tmp[$i]})
			{
			    #print $arr[0]."\t".$tmp[$i]."\t"."Inclusive"."\t".$keygid{$tmp[$i]}."\n";
				push(@tids, $tmp[$i]);
				push(@values, $keygid{$tmp[$i]});
			}
		}
		if(scalar(@tids) > 0)
		{
		    #print $_."\n";
			my $maxvalue = $values[0];
			my $maxindex = 0;
			for(my $i = 0; $i <= $#values; $i++)
			{
			    if($maxvalue <= $values[$i])
				{
				    $maxindex = $i;
					$maxvalue = $values[$i];
				}
			}
			print OUT $tids[$maxindex]."\t";
		}
		else
		{
		    print OUT "---"."\t";	
		}
	
	}
	else
	{
			print OUT "---"."\t";	
	}
	if($arr[2] ne "---")
	{
	    my @tids;
		my @values;
	    my @tmp = split(/\,/, $arr[2]);
		for(my $i = 0; $i <= $#tmp; $i++)
		{
		    if(defined $keygid{$tmp[$i]})
			{
			    #print $arr[0]."\t".$tmp[$i]."\t"."Inclusive"."\t".$keygid{$tmp[$i]}."\n";
				push(@tids, $tmp[$i]);
				push(@values, $keygid{$tmp[$i]});
			}
		}
		if(scalar(@tids) > 0)
		{
		    #print $_."\n";
			my $maxvalue = $values[0];
			my $maxindex = 0;
			for(my $i = 0; $i <= $#values; $i++)
			{
			    if($maxvalue <= $values[$i])
				{
				    $maxindex = $i;
					$maxvalue = $values[$i];
				}
			}
			print OUT $tids[$maxindex]."\t";
		}
		else
		{
		    print OUT "---"."\t";	
		}
	
	}
	else
	{
			print OUT "---"."\t";	
	}

	print OUT "\n";
}
close(IN);
close(OUT);







open(OUT, ">SE2transcript_major.xls");
print OUT "ASEid"."\t"."Inclusive"."\t"."Exclusive"."\n";
open(IN, "SE2transcript.xls");
<IN>;
while(<IN>)
{
    chomp;
	s/\s+$//;
	s/\"//g;
	my @arr = split(/\t/, $_);
	print OUT $arr[0]."\t";
	if($arr[1] ne "---")
	{
	    my @tids;
		my @values;
	    my @tmp = split(/\,/, $arr[1]);
		for(my $i = 0; $i <= $#tmp; $i++)
		{
		    if(defined $keygid{$tmp[$i]})
			{
			    #print $arr[0]."\t".$tmp[$i]."\t"."Inclusive"."\t".$keygid{$tmp[$i]}."\n";
				push(@tids, $tmp[$i]);
				push(@values, $keygid{$tmp[$i]});
			}
		}
		if(scalar(@tids) > 0)
		{
		    #print $_."\n";
			my $maxvalue = $values[0];
			my $maxindex = 0;
			for(my $i = 0; $i <= $#values; $i++)
			{
			    if($maxvalue <= $values[$i])
				{
				    $maxindex = $i;
					$maxvalue = $values[$i];
				}
			}
			print OUT $tids[$maxindex]."\t";
		}
		else
		{
		    print OUT "---"."\t";	
		}
	
	}
	else
	{
			print OUT "---"."\t";	
	}
	if($arr[2] ne "---")
	{
	    my @tids;
		my @values;
	    my @tmp = split(/\,/, $arr[2]);
		for(my $i = 0; $i <= $#tmp; $i++)
		{
		    if(defined $keygid{$tmp[$i]})
			{
			    #print $arr[0]."\t".$tmp[$i]."\t"."Inclusive"."\t".$keygid{$tmp[$i]}."\n";
				push(@tids, $tmp[$i]);
				push(@values, $keygid{$tmp[$i]});
			}
		}
		if(scalar(@tids) > 0)
		{
		    #print $_."\n";
			my $maxvalue = $values[0];
			my $maxindex = 0;
			for(my $i = 0; $i <= $#values; $i++)
			{
			    if($maxvalue <= $values[$i])
				{
				    $maxindex = $i;
					$maxvalue = $values[$i];
				}
			}
			print OUT $tids[$maxindex]."\t";
		}
		else
		{
		    print OUT "---"."\t";	
		}
	
	}
	else
	{
			print OUT "---"."\t";	
	}

	print OUT "\n";
}
close(IN);
close(OUT);






open(OUT, ">MXE2transcript_major.xls");
print OUT "ASEid"."\t"."Inclusive"."\t"."Exclusive"."\n";
open(IN, "MXE2transcript.xls");
<IN>;
while(<IN>)
{
    chomp;
	s/\s+$//;
	s/\"//g;
	my @arr = split(/\t/, $_);
	print OUT $arr[0]."\t";
	if($arr[1] ne "---")
	{
	    my @tids;
		my @values;
	    my @tmp = split(/\,/, $arr[1]);
		for(my $i = 0; $i <= $#tmp; $i++)
		{
		    if(defined $keygid{$tmp[$i]})
			{
			    #print $arr[0]."\t".$tmp[$i]."\t"."Inclusive"."\t".$keygid{$tmp[$i]}."\n";
				push(@tids, $tmp[$i]);
				push(@values, $keygid{$tmp[$i]});
			}
		}
		if(scalar(@tids) > 0)
		{
		    #print $_."\n";
			my $maxvalue = $values[0];
			my $maxindex = 0;
			for(my $i = 0; $i <= $#values; $i++)
			{
			    if($maxvalue <= $values[$i])
				{
				    $maxindex = $i;
					$maxvalue = $values[$i];
				}
			}
			print OUT $tids[$maxindex]."\t";
		}
		else
		{
		    print OUT "---"."\t";	
		}
	
	}
	else
	{
			print OUT "---"."\t";	
	}
	if($arr[2] ne "---")
	{
	    my @tids;
		my @values;
	    my @tmp = split(/\,/, $arr[2]);
		for(my $i = 0; $i <= $#tmp; $i++)
		{
		    if(defined $keygid{$tmp[$i]})
			{
			    #print $arr[0]."\t".$tmp[$i]."\t"."Inclusive"."\t".$keygid{$tmp[$i]}."\n";
				push(@tids, $tmp[$i]);
				push(@values, $keygid{$tmp[$i]});
			}
		}
		if(scalar(@tids) > 0)
		{
		    #print $_."\n";
			my $maxvalue = $values[0];
			my $maxindex = 0;
			for(my $i = 0; $i <= $#values; $i++)
			{
			    if($maxvalue <= $values[$i])
				{
				    $maxindex = $i;
					$maxvalue = $values[$i];
				}
			}
			print OUT $tids[$maxindex]."\t";
		}
		else
		{
		    print OUT "---"."\t";	
		}
	
	}
	else
	{
			print OUT "---"."\t";	
	}

	print OUT "\n";
}
close(IN);
close(OUT);
