use warnings;
use strict;

use Text::ParseWords;


## db
my $mirna_target_url = "mirna_target_url.txt";
my %db_target_url;
read_mirna_target_url($mirna_target_url,\%db_target_url);

## input

my $input_csv = "input.csv";
 
my %miRNA_list;
read_input_csv($input_csv, \%miRNA_list);

my $mirna_target_links = "mirna_target_links.txt"; # from miRBase
read_mirna_target_links($mirna_target_links,\%miRNA_list);

use Data::Dumper;
#print Dumper %miRNA_list;

# get header for target
my @db_headers;
foreach my $db (sort keys %db_target_url){
	my $db_name = $db_target_url{$db}{name};
	push @db_headers, $db_name;
}
my $db_header = join "\t", @db_headers;

my $header = "miRNA\t$db_header\tContent";
print "$header\n";

foreach my $miRNA (keys %miRNA_list){
	foreach (@{$miRNA_list{$miRNA}{content}}){
		my @db_results;
		foreach my $db (sort keys %db_target_url){
			my $db_result;

			if ($miRNA_list{$miRNA}{db_list}{$db}){
				#$db_result = $db_target_url{$db}{name};
				$db_result = $db_target_url{$db}{url};

				if ($db_result =~ /microrna\.org/){
					my $mimat_id = $miRNA_list{$miRNA}{mimat_id};
					$db_result =~ s/\<\?\>/$mimat_id/g;
				}elsif ($db_result =~ /mirdb\.org/){
					$db_result =~ s/full\=mirbase\&//g;
					$db_result =~ s/\<\?\>/$miRNA/g;
				}elsif ($db_result =~ /cm\.jefferson\.edu/){
					my $species = "MusMusculus";
					$db_result = "https://cm.jefferson.edu/rna22/Precomputed/InputController?identifier=$miRNA&minBasePairs=12&maxFoldingEnergy=-12&minSumHits=1&maxProb=.1&version=MB18E65v2&species=$species&type=mRNA";
				}elsif ($db_result =~ /DianaTools/){
					$db_result = "http://diana.imis.athena-innovation.gr/DianaTools/index.php?r=tarbase/index&mirnas=$miRNA";
				}else{
					$db_result =~ s/\<\?\>/$miRNA/g;
				}
				$db_result = "=HYPERLINK(\"".$db_result."\")";
			}else{
				$db_result = "none";
			}
			push @db_results, $db_result;
		}
		my $db_result_line = join "\t", @db_results;
		print "$miRNA\t$db_result_line\t$_\n";
	}
}

sub read_mirna_target_links{
	my ($file, $hash_ref)= @_;
	open my $fh, '<:encoding(UTF-8)', $file or die;
	while (my $row = <$fh>) {
		chomp $row;
		#my ($mature,$db,$display_name,$field1,$field2) = split /\t/, $row;
		my @l = split /\t/, $row;
		my $mature = $l[0];
		my $db = $l[1];
		my $display_name = $l[2];
		my $mimat_id = $l[3];
		if ($hash_ref->{$display_name}){
			if ($db_target_url{$db}){
				$hash_ref->{$display_name}->{db_list}->{$db}++;
				$hash_ref->{$display_name}->{mimat_id} = $mimat_id;	
			}else{ die "ERROR! no defined <$db> in <$mirna_target_url> file\n";

			}
		}
	}
	close($fh);

}

sub read_input_csv{
	my ($file, $hash_ref) = @_;
	open my $fh, '<:encoding(UTF-8)', $file or die;
	while (my $row = <$fh>) {
		chomp $row;
		if ($row =~ /^miRNA/){ next; }
		my @l = Text::ParseWords::parse_line(',', 0, $row);
		my $miRNA = $l[0];

		shift @l;
		$row = join "\t", @l;

		if ($miRNA eq "\-"){ next; }
		push @{$hash_ref->{$miRNA}->{content}}, $row;
	}
	close($fh);
}

sub read_mirna_target_url{
	my ($file, $hash_ref ) = @_;
	open my $fh, '<:encoding(UTF-8)', $file or die;
	while (my $row = <$fh>) {
		chomp $row;
		my ($db, $db_name, $db_url) =  split /\t/, $row;
		$hash_ref->{$db}->{name} = $db_name;
		$hash_ref->{$db}->{url} = $db_url;
	}
	close($fh);
}
