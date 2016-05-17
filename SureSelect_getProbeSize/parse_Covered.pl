use warnings;
use strict;

my %hash;
my $file = $ARGV[0];

my @removeGene=qw(
PLEC
PCDH10
PREX2
PCLO
FLG
DMD
COL12A1
XIRP2
);

open my $fh, '<:encoding(UTF-8)', $file or die;
while (my $row = <$fh>) {
	chomp $row;
	if ($row =~ /^browser/){ next; }
	if ($row =~ /^track/){ next; }
	my ($chr, $start, $end, $gene) = split /\t/, $row;
	my $length = $end - $start;

	my $flag = 0;
	foreach (@removeGene){
		if ($gene eq $_){
			$flag = 1;
			last;
		}
	}
	if ($flag == 1){ next; }

	$hash{$gene} += $length;
}
close($fh);

#use Data::Dumper;
#print Dumper %hash;

my $sum = 0;
foreach my $gene (keys %hash){
	my $gene_len = $hash{$gene};
	$sum += $gene_len;
}
print "sum : $sum\n";
