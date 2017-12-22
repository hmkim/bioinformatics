use warnings;
use strict;

use Bio::DB::Taxonomy;
use File::Basename; 



use Data::Dumper;

my @input_taxons = qw(
36087
);

foreach my $taxid (@input_taxons){
    # get scientific name from taxon_Id
    my $db = Bio::DB::Taxonomy->new(-source => 'entrez');
    my $taxon = $db->get_taxon(-taxonid => $taxid);
    my $scientific_name = $taxon->scientific_name;

    # Search GIs
    my $search_keyword = "(internal transcribed spacer[All Fields]) OR ITS1[All Fields] OR ITS2[All Fields] AND (".$scientific_name."[porgn])";
    #my $search_keyword = "\"".$scientific_name."\"[porgn] AND (\"ITS1\" OR \"ITS 1\" OR \"internal transcribed spacer 1\" OR \"ITS2\" OR \"ITS 2\" OR \"internal transcribed spacer 2\")";
    #my $search_keyword = "$scientific_name\[porgn] AND 500[SLEN]:1000[SLEN] AND (\"ITS1\" OR \"ITS 1\" OR \"internal transcribed spacer 1\") AND (\"ITS2\" OR \"ITS 2\" OR \"internal transcribed spacer 2\") NOT mitochondrial";
    #Ancylostoma duodenale[porgn]  AND ("ITS1" OR "ITS 1" OR "internal transcribed spacer 1") AND ("ITS2" OR "ITS 2" OR "internal transcribed spacer 2") NOT mitochondrial

    # Genbank Download From GIs
    my $gbfile = "../data/tax_$taxid.gb";
    my $cmd_downGenbankFromQuery = "python retrieve_gi.py \'$search_keyword\' $gbfile";

    if (!-f $gbfile){
        print $cmd_downGenbankFromQuery."\n";
        system($cmd_downGenbankFromQuery);
    }

    # GB to Fasta (per features)
    #my $outdir_each = "output/fa";
    #my $cmd_convertFaEach = "python parseGenBank.py $gbfile $outdir_each";
    #print $cmd_convertFaEach."\n";
    #system($cmd_convertFaEach);

    # GB to Fasta (only select ITS1, ITS2; 5.8s, 28s)  
    my $out_fa = "../data/tax_$taxid.fa";
    #my $cmd_convertFa = "python rename_convert.py $gbfile $out_fa";
    my $cmd_convertFa = "python parseGenBank_for_ITS.py $gbfile $out_fa";
    print $cmd_convertFa."\n";
    system($cmd_convertFa);

    if (!-f $out_fa){
        next;
    }
 
    my $muscle = "/Users/brandon/Downloads/bio/preprocess/tools/muscle_3.8.1551/muscle";
    my $out_muscle = "../data/tax_$taxid.msf";
    my $cmd = "$muscle -in $out_fa -out $out_muscle -msf -maxiters 1 -diags";
    print $cmd."\n";
    system($cmd);

#    my @in_fasta_files = glob("$outdir_each/tax_$taxid.*.fa");
#    foreach my $in_fa (@in_fasta_files){
#        my ($filename, $filepath, $fileext) = fileparse($in_fa, qr/\.[^.]*/);
#
#        my $out_muscle_each = "output/align/".$filename.".msf";
#        my $cmd = "$muscle -in $in_fa -out $out_muscle_each -msf -maxiters 1 -diags";
#        if (!-f $out_muscle_each){
#            #print $cmd."\n";
#            system($cmd);
#        }
#    }

}
