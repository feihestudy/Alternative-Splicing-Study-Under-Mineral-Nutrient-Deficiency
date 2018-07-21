use strict;
use Getopt::Long;
use vars qw($fafile);
Getopt::Long::GetOptions(
    'fa=s' => \$fafile,
);
$/=">";
open(IN, $fafile);
<IN>;
while(<IN>)
{
    my @arr = split(/\n/, $_);
	my $id = shift(@arr);
	my @tmp = split(/\s+/, $id);
	$id = $tmp[0];
	my $seq = join("", @arr);
	$seq=~s/>$//;
	next if ($id=~m/^LOC/ or $id=~m/^Os/);
	#if(defined $ids{$id})
	#{
	    print ">".$id."\n".$seq."\n";
	#}
}
close(IN);
#close(OUT);

   