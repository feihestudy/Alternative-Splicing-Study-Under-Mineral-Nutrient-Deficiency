use strict;
use List::Compare;
use List::Uniq ':all';
use Getopt::Long;
use vars qw($GOfile);
Getopt::Long::GetOptions(
    'GO=s'    =>\$GOfile,
);
open(IN, $GOfile);
my  $tt = <IN>;
print $tt;
my @lines;
while(<IN>)
{
   chomp;
   s/\s+$//;
   if(scalar(@lines) == 0)
   {
        push(@lines, $_);
        next;		
   }
   else
   {
       my @tmp=split(/\t/, $_);
	   #print $tmp[7]."\n";
       my @goids = split(/\,/, $tmp[7]);
	   my $yes = 0;
       for(my $i = 0; $i <= $#lines; $i++)
       {
	        my @tmp2=split(/\t/, $lines[$i]);
			my @goids2 = split(/\,/, $tmp2[7]);
	        my $lcsh = List::Compare->new(\@goids, \@goids2);
			my @intersection = $lcsh->get_intersection;
			if(scalar(@intersection) == scalar(@goids) and scalar(@intersection) == scalar(@goids2))
			{
			    $yes = 1;			
			}
			else
			{			
			}			
	   }
	   if($yes == 1)
	   {
	       
	   
	   }
	   else
	   {
			push(@lines, $_);	   
	   }
   }

}
close(IN);



print join("\n", @lines)."\n";

