use warnings;
use strict;

use Data::Dumper;

# Filepath, filename, file extension parsing
use File::Basename; 

use Bio::AlignIO;

my $in_msf = $ARGV[0];

my $in = Bio::AlignIO->new(-file => $in_msf ,
                        -format => 'msf');


my ($filename, $filepath, $fileext) = fileparse($in_msf, qr/\.[^.]*/);
my $out_clustalw = $filepath.$filename.".aln";

my $out = Bio::AlignIO->new(-file => ">".$out_clustalw,
                         -format => 'clustalw');

while ( my $aln = $in->next_aln ) {
    $out->write_aln($aln);
}

my $bin = "/Users/brandon/Downloads/SelectSequencesFromMSA/dist/build/SelectSequencesFromMSA/SelectSequencesFromMSA";
my $cmd = "$bin -c '$out_clustalw'  -x -o `pwd` -n 2 -i 80.0 -m 95.0";
print $cmd."\n";
