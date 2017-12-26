use warnings;
use strict;

use Data::Dumper;


############
############
my @hiddenDirNames = ('cron', 'math', 'NAS_Backup', 'NASLog', 'NetBackup', 'Bioinformatics');
my @hiddenGroupNames = ('관리자', 'administrators', 'users', 'Bioinformatics');

if (@ARGV != 1){
	printUsage();
}

############
############

my $outputFile = $ARGV[0];

# check synonlogy tool
my $synoacltool = "/usr/syno/bin/synoacltool";

if (!-f $synoacltool){
	die "ERROR! not found <$synoacltool>\n";
}

# get list for directory in NAS
my $dirList = getDirList();

my %hash;

my %hash_group;


foreach my $path (@$dirList){
	my $groupInfoList  = getGroupsInfoForDir($path);
	my $folderName = $path;
	$folderName =~ s/^\/volume1\///;

	foreach my $groupInfo (@$groupInfoList){
		my ($group, $perm) = parse_groupInfo($groupInfo);
		#print "$folderName\t$group\t$perm\n";

		$hash{$folderName}{$group} = $perm;

		$hash_group{$group}++;
	}
}

my @groupList = sort keys %hash_group;

my $header = "Folder name\t".(join "\t", @groupList);

open my $fh, '>', $outputFile or die;
print $fh $header."\n";
foreach my $dir (sort keys %hash){

	my $res = '';
	my @arr;
	foreach my $group (@groupList){
		if ($hash{$dir}{$group}){
			push @arr, $hash{$dir}{$group};
		}else{
			push @arr, '';
		}
	}

	my $result_line = join "\t", @arr;
	print $fh "$dir\t$result_line\n";
}
close($fh);

sub parse_groupInfo{
	my $txt = shift;
	my @array = split /\:/, $txt;

	my $groupName = $array[1];
	my $groupPermissionFlag = $array[2]; # allow, deny?
	
	if ($groupPermissionFlag ne "allow"){ die "ERROR! new flag!! check the script or contact developer\n"; }
	my $groupPermission = $array[3];

	if ($groupPermission eq "rwxpdDaARWc--"){
		$groupPermission = "rw";

	}elsif ($groupPermission eq "r-x---a-R-c--"){
		$groupPermission = "r";
	}else{
		die "ERROR! new  permission information! check the script or contact developer\n"; 
	}
	return ($groupName, $groupPermission);
}

sub getGroupsInfoForDir{
	my $dir = shift;

	#my $cmd_getAllowGroupListForDir = "$synoacltool -get $dir | grep group | cut -d':' -f2,3 | awk -F':' '{ if (\$2 == \"allow\") print \$1 }'";
	my $cmd_getAllowGroupListForDir = "$synoacltool -get $dir | grep group";

	my $res = `$cmd_getAllowGroupListForDir`;
	chomp($res);

	my %hash = map { (split /\:/, $_)[1] => $_ } split(/\n/, $res) ;
	foreach my $string (@hiddenGroupNames){
		if (exists $hash{$string}){
			delete $hash{$string};
		}		
	}

	my @array = map { $_ } values %hash;
	return \@array;
}

sub getDirList{
	my $parentName = "/volume1/";

	my $cmd = "ls -al $parentName | grep \"^d-\" | awk '{ print \"$parentName\"\$NF }'";
	my $res = `$cmd` or die "ERROR! Could not get directory list\n";
	chomp($res);

	my %hash = map { $_ => undef } split(/\n/, $res);

	foreach my $string (@hiddenDirNames){
		if (exists $hash{$parentName.$string}){
			delete $hash{$parentName.$string};
		}
	}

	my @array = map { $_ } keys %hash;

	return \@array;
}

sub printUsage{
	print "Usage: perl $0 <output filename>\n\n";
	print "Hidden group: ".(join ", ", @hiddenGroupNames),"\n";
	print "Hidden directory: ".(join ", ", @hiddenDirNames),"\n\n";
	print "Example: perl $0 <result.xls>\n";

	exit;
}
