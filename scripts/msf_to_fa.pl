use Bio::AlignIO;
use Data::Dumper;

$in = Bio::AlignIO->new(-file => $ARGV[0] ,
                        -format => 'msf');

$out = Bio::AlignIO->new(-file => ">$ARGV[1]",
                         -format => 'fasta');

while ( my $aln = $in->next_aln ) {
    for my $s ($aln->each_seq){
            $s->{seq} =~ s/\-//g;
            $s->{seq} =~ s/\.//g;
            $s->{seq} = uc ($s->{seq});
    }
    $out->write_aln($aln);
}

