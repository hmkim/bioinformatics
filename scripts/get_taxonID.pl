use warnings;
use strict;

use Bio::DB::Taxonomy;

use Data::Dumper;

my $infile = $ARGV[0];

my $db = Bio::DB::Taxonomy->new(-source => 'entrez');
open my $fh, '<:encoding(UTF-8)', $infile or die;
while (my $row = <$fh>) {
    chomp $row;
    my $taxonid = $db->get_taxonid($row);

    if (!$taxonid){
        $taxonid = "None";
    }
    print $taxonid."\n";
}
close($fh);

