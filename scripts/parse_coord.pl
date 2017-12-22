#!/usr/bin/env perl

use strict;
use warnings;
use Getopt::Long;
use Time::HiRes qw(gettimeofday);

# now for each potential nucleotide accession for protein X in $tmpfile2, try to find a CDS for the 
# gene that protein X codes for

# now parse $tmpfile4
# example lines:
# CU329670.1	complement(<1..5662)	tlh1	
# CU329670.1	complement(5726..6331)	
# CU329670.1	complement(join(822332..822429,822489..822849))	smn1	

my $infile = $ARGV[0];
open(IN, "<$infile") || die "ERROR unable to open $infile for reading";

my %coords_for_prot_HA = (); # coord strings for we've found for each protein, key: protein accession, value: array of matching CDS coord strings
while(my $line = <IN>) { 
    chomp $line;
    my ($nt, $coords, $gene, $product) = split(/\s+/, $line);

    # first we have to doctor the coords string
    # examples: 
    # nt:U50746.1   coords:4..867 has to become U50746.1:4..867
    # nt:U43883.1	  join(U43876.1:608..688,U43877.1:104..175,U43878.1:118..237,U43879.1:84..284,U43880.1:69..221,U43881.1:103..198,U43882.1:53..163,209..259)
    # has to become join(U43876.1:608..688,U43877.1:104..175,U43878.1:118..237,U43879.1:84..284,U43880.1:69..221,U43881.1:103..198,U43882.1:53..163,U43883.1:209..259)

    my $head = "";
    my $tail = "";
    if($coords =~ /(^.+\()/) { $head = $1; } # save everything up to final '(', if anything
    if($coords =~ /(\).*$)/) { $tail = $1; } # save everything starting at first ')', if anything
    $coords =~ s/^.+\(//; # remove head
    $coords =~ s/\).*$//; # remove tail
    my $new_coords = $head;
    foreach my $tok (split(',', $coords)) { 
      if($tok !~ m/^\S+\:/) { # no accession.version in this token, add it
        $tok = $nt . ":" . $tok; 
      }
      $new_coords .= $tok . ",";
    }
    # we added one too many ',', remove the final one
    $new_coords =~ s/\,$//;
    $new_coords .= $tail;
    print "$nt\t$gene\t$coords\t$new_coords\n";
    $coords = $new_coords;

#    if(defined $gene) { # this indicates we found INSDQualifier_name 'gene' for this NT accession
#      foreach my $prot (@{$xref2prot_HA{$nt}}) { 
#        if(exists $prot2gene_HA{$prot}) {
#          for(my $g = 0; $g < scalar(@{$prot2gene_HA{$prot}}); $g++) { 
#            if($gene eq $prot2gene_HA{$prot}[$g]) { # match
#              my $already_printed = 0;
#              # determine CDS length from coords string
#              my $cds_len  = LengthFromCoords($coords);
#              if(! exists $prot_len_H{$prot}) { 
#                die "ERROR protein $prot somehow created in error";
#              }
#              my $prot_len = $prot_len_H{$prot};
#              if(! exists $coords_for_prot_HA{$prot}) { 
#                # initialize
#                @{$coords_for_prot_HA{$prot}} = (); 
#              }
#              else { # search for a match in the existing coords strings for this protein
#                for(my $z = 0; $z < scalar(@{$coords_for_prot_HA{$prot}}); $z++) { 
#                  if($coords_for_prot_HA{$prot}[$z] eq $coords) { 
#                    $already_printed = 1;
#                  }
#                }
#              }
#              if(! $already_printed) { 
#                # currently, output all matches to OUT5
#                print OUT5 ("$prot\t$coords\n");
#                # THIS BLOCK ALTERNATIVELY OUTPUTS TO 3 DIFFERENT FILES
#                ## check that the lengths make sense
#                #if($cds_len % 3 == 0) { # should be...
#                #  if((sprintf("%d", ($cds_len / 3)) eq $prot_len) || # matches exactly
#                #    (sprintf("%d", ($cds_len / 3)) eq ($prot_len+1))) { # CDS presumably includes stop codon
#                #    print OUT5 ("$prot\t$coords\t$prot_len\t$cds_len\n");
#                #  }
#                #  else { 
#                #    print OUT7 ("$prot\t$prot_len\t$nt\t$cds_len\t$coords\n");
#                #  }
#                #}
#                #else { # cds_len a factor of 3...
#                #  print OUT6 ("$prot\t$coords\t$prot_len\t$cds_len\n");
#                #}
#                # END OF BLOCK
#                push(@{$coords_for_prot_HA{$prot}}, $coords);
#              }
#            }
#          }
#        }
#      }
#    }

}
close(IN);
