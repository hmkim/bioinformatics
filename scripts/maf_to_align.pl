use Bio::AlignIO;

$inputfilename = $ARGV[0];
$outputfilename = $ARGV[1];
$format = $ARGV[2];
$in  = Bio::AlignIO->new(-file   => $inputfilename ,
                         -format => 'maf');
$out = Bio::AlignIO->new(-file   => ">$outputfilename" ,
                         -format => "$format");

while ( my $aln = $in->next_aln() ) {
    $out->write_aln($aln);
}
