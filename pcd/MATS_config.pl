use strict;
use Getopt::Long;
use List::Uniq ':all';
use vars qw($samplename $dir $readlen $samfile $dataType);
Getopt::Long::GetOptions(
    'dir=s' => \$dir,
	'readlen=s' => \$readlen,
	'sam=s' => \$samfile,
	'dataType=s' => \$dataType,
);
my $junctionlen = 2*($readlen - 1);
print "-readLength = $readlen
-junctionLength = $junctionlen"."\n";
print "-SE = $dir/AS.SE.txt
-MXE = $dir/AS.MXE.txt
-A5SS = $dir/AS.A5SS.txt
-A3SS = $dir/AS.A3SS.txt
-RI = $dir/AS.RI.txt
-experiment = RNASeq
-base_1 = SAMPLE_1
-base_2 = SAMPLE_2
-dataType = $dataType
-libType = fr-unstranded
-samDir = 
-outDir = ./
-input_1 = $samfile
-input_2 = $samfile
-email = yourmail\@domain.com
";
