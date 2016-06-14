use warnings;
use strict;
use Getopt::Long;

use Text::ParseWords;

my $miR_table;

GetOptions(
	'table=s' => \$miR_table
);

if (@ARGV !=1){
	printUsage();
}

my %hash;
if ($miR_table){
	read_table($miR_table, \%hash);
}
my $input_txt = $ARGV[0];

# TargetScan
my $species_TargetScan = "Mouse"; 
my $species_code_TargetScan = "mmu_71";
my $link_TargetScan_template = "http://www.targetscan.org/cgi-bin/targetscan/$species_code_TargetScan/targetscan.cgi?species=$species_TargetScan&gid=&mir_sc=&mir_c=&mir_nc=&mir_vnc=&mirg=";
# miRBase
my $link_miRBase_template = "http://www.mirbase.org/cgi-bin/mature.pl?mature_acc=";
# MicroCosm Targets
my $genome_id_MicroCosm_Targets = 3876;
my $link_MicroCosm_Targets_template = "http://www.ebi.ac.uk/enright-srv/microcosm/cgi-bin/targets/v5/hit_list.pl?genome_id=$genome_id_MicroCosm_Targets;mirna_id=";
# microRNA.org
my $microRNA_org_organism_code = 10090;
# TarBase
my $link_TarBase_template = "http://diana.imis.athena-innovation.gr/DianaTools/index.php?r=tarbase/index&mirnas=";
# mirDB
my $link_mirDB_template = "http://mirdb.org/cgi-bin/pre_mir.cgi?name=";

open my $fh, '<:encoding(UTF-8)', $input_txt or die;
print "#Systematic\tAccession\tmiRBase\tTargetScan\tMicroCosm\tmicroRNA.org\tTarBase\tmirDB\n";
while (my $row = <$fh>) {
	chomp $row;
	if ($row =~ /^miRNA/){ next; } # pass header
	
	my @l = Text::ParseWords::parse_line(',', 0, $row);
	my $id = $l[0];
	$id =~ /(mmu)\-([a-zA-Z]+)\-([0-9a-zA-Z]+)-([0-9a-zA-Z]+)/;
	
	print STDERR $id."\n";

	my $new_id = $1."-".$2."-".$3; # mmu-let-7a
	my $new_id_2 = $1."-".$2."-".$3."-".$4; # mmu-let-7a-1
	if (!$new_id){
		die "Parse Error!! $row\n";
	}

	my $accession;

	
	my $link_TargetScan = $link_TargetScan_template.$id; # mmu-let-7a-1-3
	my $link_MicroCosm_Targets = $link_MicroCosm_Targets_template.$new_id; # mmu-let-7a
	my $link_microRNA_org = "http://www.microrna.org/microrna/getTargets.do?matureName=$new_id&organism=$microRNA_org_organism_code"; # mmu-let-7a
	my $link_TarBase = $link_TarBase_template.$id; # mmu-let-7a-1-3
	my $link_mirDB = $link_mirDB_template.$new_id_2; # mmu-let-7a-1

	if ($id !~ /p$/){
		$id .= "p"; 
	}
	if ($hash{$id}){
		$accession = $hash{$id};
	}else{ die "ERROR ! not found id in table file <$id>"; }

	my $link_miRBase = $link_miRBase_template.$accession; # only accession MIMAT0004620

	print "$id\t$accession\t$link_miRBase\t$link_TargetScan\t$link_MicroCosm_Targets\t$link_microRNA_org\t$link_TarBase\t$link_mirDB\n";
}
close($fh);

sub read_table{
	my ($input, $hash_ref) = @_;
	open my $fh_tab, '<:encoding(UTF-8)', $input or die;
	while (my $row = <$fh_tab>) {
			chomp $row;
			my @l = Text::ParseWords::parse_line(',', 0, $row);
			my $id = $l[1]; # cel-let-7

			my $mature1_acc= $l[4]; # MIMAT0000001
			my $mature1_id = $l[5]; # cel-let-7-5p
			
			my $mature2_acc = $l[7]; # MIMAT0015091
			my $mature2_id = $l[8]; # cel-let-7-3p

			if ($mature1_id eq $mature2_id){ die; }
			$hash_ref->{$mature1_id} = $mature1_acc;
			$hash_ref->{$mature2_id} = $mature2_acc;
			#print "$mature1_id => $mature1_acc";
			#print "\t$mature2_id => $mature2_acc\n";

	}
	close($fh_tab)
}

sub printUsage{
	print "Usage: perl $0 <input.txt> > <output.txt>\n";
	exit;
}
