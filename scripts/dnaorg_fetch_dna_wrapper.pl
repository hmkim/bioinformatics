#!/usr/bin/env perl
# EPN, Fri Mar 13 14:17:55 2015
#
use strict;
use warnings;
use Getopt::Long;
use Time::HiRes qw(gettimeofday);

my $list_dir = "/panfs/pan1/dnaorg/share/";
my $df_ffile    = $list_dir . "all.multi-synonym.list";
my $df_sfile    = $list_dir . "all.primary-and-synonym.list";

# The definition of $usage explains the script and usage:
my $usage = "\ndnaorg_fetch_dna_wrapper.pl:\n";
$usage .= "\n"; 
$usage .= " This script fetches DNA sequences from GenBank. It can be\n";
$usage .= " run in three different modes:\n\n";

$usage .= " Mode 1: fetch CDS sequences for protein database records that link\n";
$usage .= "         to a symbol for a protein coding gene in the Gene database,\n";
$usage .= "         or if no matches exist in protein, repeat search in nuccore\n";
$usage .= "         and fetch all nucleotide records that link to the noncoding\n";
$usage .= "         symbol in the Gene database.\n";
$usage .= "\n";
$usage .= "   usage: 'perl dnaorg_fetch_dna_wrapper.pl [OPTIONS] <symbol>'\n";
$usage .= "\n";
$usage .= "\n";
$usage .= " Mode 2: fetch CDS sequences for a list of protein accessions\n";
$usage .= "\n";
$usage .= "   usage: 'perl dnaorg_fetch_dna_wrapper.pl [OPTIONS] -d <name_for_outdir> -plist <file_with_list_of_protein_accessions>'\n";
$usage .= "\n";
$usage .= "\n";
$usage .= " Mode 3: fetch nucleotide sequences for a list of nuccore accessions\n";
$usage .= "\n";
$usage .= "   usage: 'perl dnaorg_fetch_dna_wrapper.pl [OPTIONS] -d <name_for_outdir> -ntlist <file_with_list_of_nucleotide_accessions>'\n";
$usage .= "\n";
$usage .= " NOTE: NCBI Gene symbols can include whitespace (' ') but this script\n";
$usage .= "       does not allow them in the <symbol> command line argument. To\n";
$usage .= "       specify a ' ' in the <symbol> replace it with a '~' on the\n";
$usage .= "       command line. (As of 04/28/15 no '~' exist in NCBI Gene symbols\n";
$usage .= "\n";
$usage .= " BASIC OPTIONS:\n";
$usage .= "  -f      : force; if dir <symbol> exists, overwrite it\n";
$usage .= "  -v      : be verbose; output commands to stdout as they're run\n";
$usage .= "  -d <s>  : define output directory as <s>, not <symbol>\n";
$usage .= "  -nt     : search for matches to symbol in ONLY the nuccore db, not the protein db\n";
$usage .= "  -notnt  : if zero matches for symbol in the protein db, DO NOT retry using the nuccore db\n";
$usage .= "  -nosyn  : only include records for which the symbol is the primary symbol, not the synonym\n";
$usage .= "\n";
$usage .= " OPTIONS THAT ENABLE ALTERNATE MODES:\n";
$usage .= "  -plist    : <symbol> is really a list of protein accessions,    requires -d option too\n";
$usage .= "  -ntlist   : <symbol> is really a list of nucleotide accessions, requires -d option too\n";
$usage .= "  -paccn    : <symbol> is really a single protein accession,    requires -d option too\n";
$usage .= "  -ntaccn   : <symbol> is really a single nucleotide accession, requires -d option too\n";
$usage .= "  -num      : determine number of matching protein/nucleotide accessions, then exit\n";
$usage .= "  -geneinfo : output information on annotated genes in each record, then exit\n";
$usage .= "  -ftable   : output feature table for each record, then exit\n";
$usage .= "  -matpept  : output information on mat_peptide features for each record, then exit\n";
$usage .= "\n";
$usage .= " OPTIONS THAT AFFECT FAILURE/WARNING OF PRE-DETERMINED SYMBOLS:\n";
$usage .= "  -ffile <f> : fail          if a symbol listed in <f> is used as input symbol [default: $df_ffile]\n";
$usage .= "  -sfile <f> : require symbol be primary symbol for symbols listed in file <f> [default: $df_sfile]\n";
$usage .= "  -noffile   : do not fail for symbols listed in a file\n";
$usage .= "  -nosfile   : do not skip synonyms for symbols listed in a file\n";
$usage .= "\n";
$usage .= " EXPERIMENTAL/ADVANCED OPTIONS:\n";
$usage .= "  -up        : additional run experimental code for fetching non-CDS UniProt CDS via xrefs\n";
$usage .= "  -old       : use extract_fasta_multi_exon instead of esl-fetch-cds.pl\n";
$usage .= "  -ngene <n> : with -geneinfo, specify number of accessions per edirect query for gene info\n";
$usage .= "\n";
$usage .= "\n";
$usage .= " This script will create a directory called <symbol> (for modes 1 and 2)\n";
$usage .= " or called <name_for_outdir> (for modes 3 and 4) and populate it with\n";
$usage .= " several output files, including:\n";
$usage .= "\n";
$usage .= "  <symbol>.log: list of all output files and brief descriptions\n";
$usage .= "  <symbol>.cmd: list of all commands run by this script\n";
$usage .= "  <symbol>.sum: copy of all stdout from this script\n";
$usage .= "\n";
$usage .= " For modes 2 and 3, the files will start with <name_for_outdir>\n";
$usage .= " instead of <symbol>.\n";
$usage .= "\n";

my ($seconds, $microseconds) = gettimeofday();
my $start_secs = ($seconds + ($microseconds / 1000000.));
my $executable = $0;

# initialize variables that can be changed with cmdline options
my $do_force        = 0;     # set to '1' with -f, overwrite output dir if it exists
my $be_verbose      = 0;     # set to '1' with -v, output commands as they're run to stdout
my $out_dir         = undef; # set to a value <s> with -d <s>
my $do_uniprot_xref = 0;     # set to '1' if -up option used.
my $do_nt_userset   = 0;     # set to '1' if -nt option is used
my $do_nt           = 0;     # set to '1' if EITHER -nt used or we switch to searching in nuccore mid-script
my $do_not_try_nt   = 0;     # set to '1' if -notnt option is used
my $do_nosyn        = 0;     # set to '1' if -nosyn option is used
my $do_old          = 0;     # set to '1' with -old
my $do_allgene      = 0;     # set to '1' if -allgene option used.
my $ngene           = undef; # set to value if -ngene used
my $df_ngene        = 5;     # default value to use for $ngene if -ngene *not* used
my $ffile           = undef; # set to value if -ffile used
my $sfile           = undef; # set to value if -sfile used
my $no_ffile        = 0;     # set to '1' if -noffile used
my $no_sfile        = 0;     # set to '1' if -nosfile used

# variables that are indirectly changed by cmdline options
my $acclist_file    = undef; # if -plist or -ntlist is used, this is defined as $symbol
my $do_acclist_mode = 0;     # set to '1' if either -plist or -ntlist used

# different 'modes', if all are false we run in default mode
my $do_num_mode      = 0; # set to '1' if -numonly option used.
my $do_geneinfo_mode = 0; # set to '1' if -geneinfo option used.
my $do_ftable_mode   = 0; # set to '1' if -ftable option used.
my $do_matpept_mode  = 0; # set to '1' if -matpept option used.
my $do_plist_mode    = 0; # set to '1' if -plist option used.
my $do_ntlist_mode   = 0; # set to '1' if -ntlist option used.
my $do_paccn_mode    = 0; # set to '1' if -paccn option used.
my $do_ntaccn_mode   = 0; # set to '1' if -ntaccn option used.

&GetOptions( "f"        => \$do_force, 
             "v"        => \$be_verbose,
             "d=s"      => \$out_dir,
             "nt"       => \$do_nt_userset,
             "notnt"    => \$do_not_try_nt,
             "plist"    => \$do_plist_mode,
             "ntlist"   => \$do_ntlist_mode,
             "paccn"    => \$do_paccn_mode,
             "ntaccn"   => \$do_ntaccn_mode,
             "num"      => \$do_num_mode,
             "geneinfo" => \$do_geneinfo_mode,
             "ftable"   => \$do_ftable_mode,
             "matpept"  => \$do_matpept_mode,
             "ffile=s"  => \$ffile,
             "sfile=s"  => \$sfile,
             "noffile"  => \$no_ffile,
             "nosfile"  => \$no_sfile,
             "up"       => \$do_uniprot_xref,
             "nosyn"    => \$do_nosyn,
             "old"      => \$do_old, 
             "allgene"  => \$do_allgene,
             "ngene=s"  => \$ngene);

if(scalar(@ARGV) != 1) { die $usage; }
my ($symbol) = (@ARGV);
my $argv_1   = $symbol;
my $cap_symbol = Capitalize($symbol);
my $cap_no_tilde_symbol = $cap_symbol;
$cap_no_tilde_symbol =~ s/\~/ /g;

# store options used, so we can output them 
my $opts_used_short = "";
my $opts_used_long  = "";
if($do_force) { 
  $opts_used_short .= "-f ";
  $opts_used_long  .= "# option:  forcing overwrite of $symbol dir/file [-f]\n"; 
}
if($be_verbose) {
  $opts_used_short .= "-v ";
  $opts_used_long  .= "# option:  verbose mode [-v]\n"; 
}
if(defined $out_dir) { 
  $opts_used_short .= "-d $out_dir ";
  $opts_used_long  .= "# option:  output directory specified as $out_dir [-d]\n"; 
}
if($do_nt_userset) { 
  $opts_used_short .= "-nt ";
  $opts_used_long  .= "# option:  search in the nuccore database, not in the protein db [-nt]\n";
}
if($do_not_try_nt) { 
  $opts_used_short .= "-notnt ";
  $opts_used_long  .= "# option:  if no matches in the protein database, DO NOT try the nuccore database [-notnt]\n";
}
if($do_nosyn) { 
  $opts_used_short .= "-nosyn ";
  $opts_used_long  .= "# option:  require a match to the gene symbol, gene_synonyms matches (or non matches) are ignored [-nosyn]\n";
}
if($do_plist_mode) { 
  $opts_used_short .= "-plist ";
  $opts_used_long  .= "# option:  $symbol is a list of protein accessions, not a symbol [-plist]\n"; 
}
if($do_ntlist_mode) { 
  $opts_used_short .= "-ntlist ";
  $opts_used_long  .= "# option:  $symbol is a list of nucleotide accessions, not a symbol [-ntlist]\n"; 
}
if($do_paccn_mode) { 
  $opts_used_short .= "-paccn ";
  $opts_used_long  .= "# option:  $symbol is a single protein accession, not a symbol [-paccn]\n"; 
}
if($do_ntaccn_mode) { 
  $opts_used_short .= "-ntaccn ";
  $opts_used_long  .= "# option:  $symbol is a single nucleotide accession, not a symbol [-ntaccn]\n"; 
}
if($do_num_mode) { 
  $opts_used_short .= "-num ";
  $opts_used_long  .= "# option:  determining number of matching protein accessions, then exiting [-num]\n"; 
}
if($do_geneinfo_mode) { 
  $opts_used_short .= "-geneinfo ";
  $opts_used_long  .= "# option:  outputting annotated gene information, then exiting [-geneinfo]\n"; 
}
if($do_matpept_mode) { 
  $opts_used_short .= "-matpept ";
  $opts_used_long  .= "# option:  outputting mat_peptide information, then exiting [-matpept]\n"; 
}
if(defined $ffile) { 
  $opts_used_short .= "-ffile $ffile ";
  $opts_used_long  .= "# option:  list file with symbols to fail for: $ffile [-ffile]\n"; 
}
if(defined $sfile) { 
  $opts_used_short .= "-sfile $sfile ";
  $opts_used_long  .= "# option:  list file with symbols to skip synonyms for: $sfile [-sfile]\n"; 
}
if($no_ffile) { 
  $opts_used_short .= "-noffile "; 
  $opts_used_long  .= "# option:  no list file with symbols to fail for [-noffile]\n"; 
}
if($no_sfile) { 
  $opts_used_short .= "-nosfile "; 
  $opts_used_long  .= "# option:  no list file with symbols to skip synonyms for [-nosfile]\n"; 
}
if($do_uniprot_xref) { 
  $opts_used_short .= "-up ";
  $opts_used_long  .= "# option:  trying to fetch non-CDS UniProt using xrefs [-up]\n";
}
if($do_old) { 
  $opts_used_short .= "-old ";
  $opts_used_long  .= "# option:  using old strategy (extract_fasta_multi_exon) [-old]\n"; 
}
if($do_allgene) { 
  $opts_used_short .= "-allgene ";
  $opts_used_long  .= "# option:  fetching info on all linked gene entries, not just seq-specific ones [-allgene]\n"; 
}  
if(defined $ngene) { 
  $opts_used_short .= "-ngene $ngene ";
  $opts_used_long  .= "# option:  specifying number of accessions for edirect gene queries as $ngene [-ngene]\n";
}

# check for incompatible option combinations:
if($do_plist_mode || $do_ntlist_mode) { 
  if(! defined $out_dir) { 
    die "ERROR the -plist and -ntlist options require the -d option";
  }
}
if($do_plist_mode || $do_ntlist_mode || $do_paccn_mode || $do_ntaccn_mode) { 
  if($do_nt_userset) { 
    die "ERROR the -plist, -paccn, -ntlist and -ntaccn options are incompatible with -nt";
  }
  if($do_not_try_nt) { 
    die "ERROR the -plist, -paccn, -ntlist and -ntaccn options are incompatible with -notnt";
  }
  if($do_nosyn) { 
    die "ERROR the -plist, -paccn, -ntlist and -ntaccn options are incompatible with -nosyn"; 
  }
  if($do_num_mode) { 
    die "ERROR the -plist, -paccn, -ntlist and -ntaccn options are incompatible with -num";
  }
  if(defined $ffile) { 
    die "ERROR the -plist, -paccn, -ntlist and -ntaccn options are incompatible with -ffile";
  }
  if(defined $sfile) { 
    die "ERROR the -plist, -paccn, -ntlist and -ntaccn options are incompatible with -sfile";
  }
}
if(($do_num_mode + $do_geneinfo_mode + $do_ftable_mode) > 1) { 
  die "ERROR, only one of the following options can be used: -num, -geneinfo, -ftable";
}  
if(($do_num_mode + $do_geneinfo_mode + $do_matpept_mode) > 1) { 
  die "ERROR, only one of the following options can be used: -num, -geneinfo, -ftable";
}  
if($do_plist_mode && $do_ntlist_mode) { 
  die "ERROR the -plist and -ntlist options are incompatible"; 
}
if($do_nt_userset) { 
  if($do_uniprot_xref) { 
    die "ERROR the -nt option is incompatible with the -up option";
  }
  if($do_not_try_nt) { 
    die "ERROR the -nt option is incompatible with the -notnt option";
  }
}
if(defined $ffile && $no_ffile) { 
  die "ERROR the -ffile and -noffile options are incompatible"; 
}
if(defined $sfile && $no_sfile) { 
  die "ERROR the -sfile and -nosfile options are incompatible"; 
}
if(defined $ngene && (! $do_geneinfo_mode)) { 
  die "ERROR -ngene only makes sense in combination with the -gene option";
}
if($do_allgene && (! $do_geneinfo_mode)) { 
  die "ERROR -agene only makes sense in combination with the -gene option";
}

my $exec_dir  = "/panfs/pan1/dnaorg/programs";
my $idstat    = "/netopt/genbank/subtool/bin/idstat";
my $fetch_cds = "$exec_dir/esl-fetch-cds.pl";
my $idfetch   = "/netopt/ncbi_tools64/bin/idfetch";
my $query;              # an argument to a '-query' option in a edirect tool cmdline
my $cmd;                # a command to run with system() in RunCommand()
my $desc;               # description of a step
my $nlines;             # number of lines in a file
my $nlost;              # number of lost elements (often accessions) in a step
my $ncreated;           # number of created elements (often accessions) in a step, should always be 0
my $nsecs;              # number of seconds a command took
my $errmsg;             # error message
my $using_df_ffile = 0; # set to '1' if --ffile not used
my $using_df_sfile = 0; # set to '1' if --sfile not used
if(! defined $ffile) { $ffile = $df_ffile; $using_df_ffile = 1; }
if(! defined $sfile) { $sfile = $df_sfile; $using_df_sfile = 1; }
if(! defined $ngene) { $ngene = $df_ngene; }

$do_acclist_mode = ($do_plist_mode || $do_ntlist_mode || $do_paccn_mode || $do_ntaccn_mode) ? 1 : 0; # $do_acclist_mode is true if any of -plist, -ntlist, -paccn, or -ntaccn used
$do_nt           = ($do_nt_userset || ($do_ntlist_mode || $do_ntaccn_mode)) ? 1 : 0;

###############
# Preliminaries
###############
# check if our output dir $symbol exists
if(! defined $out_dir) { 
  $out_dir = $symbol . "/";
}
else { 
  if($out_dir !~ m/\/$/) { $out_dir .= "/"; } # add '/'
}
if(-d $out_dir) { 
  if($do_force) { RunCommand("rm -rf $out_dir", $be_verbose, undef); }
  else          { die "ERROR directory named $out_dir already exists. Remove it, or use -f to overwrite it."; }
}
if(-e $out_dir) { 
  if($do_force) { RunCommand("rm $out_dir", $be_verbose, undef); }
  else          { die "ERROR a file named $out_dir already exists. Remove it, or use -f to overwrite it."; }
}

if($do_acclist_mode) { 
  if((! $do_paccn_mode) && (! $do_ntaccn_mode)) { 
    $acclist_file = $symbol;
    if(! -e $acclist_file) { die "ERROR no file $acclist_file exists"; }
    if(! -s $acclist_file) { die "ERROR file $acclist_file is empty"; }
  }
  else { # -paccn or -ntaccn used
    $acclist_file = $out_dir . "/" . $symbol . ".singleacclist";
  }
  $symbol       = $out_dir; # this must be defined, we died above it if wasn't
  $symbol       =~ s/\///;
  $cap_symbol   = Capitalize($symbol);
  $no_ffile = 1;
  $no_sfile = 1;
}

# check if the symbol is in the fail file and nosyn file, if nec
if(! $no_ffile) { 
  if(checkFileForSymbol($ffile, $cap_no_tilde_symbol)) { 
    my $die_msg;
    if($using_df_ffile) { 
      $die_msg =  "FAIL: $symbol is listed in $ffile, this symbol is ambiguous because it is a synonym\n";
      $die_msg .= "for > 1 GENE primary symbols and never a primary symbol itself. Exiting.\n"; 
    }
    else { # not using the default file, unsure why we're failing:
      $die_msg = "FAIL: $symbol is listed in $ffile listed with -ffile. Exiting."; 
    }
    die $die_msg;
  }
}
my $warn_msg = undef;
if((! $no_sfile) && (! $do_nosyn)) { 
  if(checkFileForSymbol($sfile, $cap_no_tilde_symbol)) { 
    if($using_df_sfile) { 
      $warn_msg  = "# WARNING: $symbol is listed in $sfile, this symbol is ambiguous because it is a primary symbol\n";
      $warn_msg .= "# and an alias >= 1 other GENE primary symbols.\n"; 
      $warn_msg .= "# ***Accessions for which this symbol is not the primary symbol will be skipped.***\n";
      $warn_msg .= "# (To turn this behavior off, use -nosfile)\n";
    }
    else { 
      $warn_msg  = "# WARNING: $symbol is listed in $sfile supplied with -sfile.\n";
      $warn_msg .= "# ***Accessions for which this symbol is not the primary symbol will be skipped.***\n";
    }
    $do_nosyn = 1; 
  }
}

# create the dir
RunCommand("mkdir $out_dir", $be_verbose, undef);

# open the log and summary files:
my $log_file = $out_dir . $symbol . ".log";
my $sum_file = $out_dir . $symbol . ".sum";
my $cmd_file = $out_dir . $symbol . ".cmd";
my $log_FH; # file handle for log file output
my $sum_FH; # file handle for summary file output
my $cmd_FH; # file handle for command file output
open($log_FH, ">" . $log_file) || die "ERROR unable to open $log_file for writing";
open($sum_FH, ">" . $sum_file) || die "ERROR unable to open $sum_file for writing";
open($cmd_FH, ">" . $cmd_file) || die "ERROR unable to open $cmd_file for writing";
OutputFileInfo($log_file, "log file containing list and description of all output files", "", $log_FH); 
OutputFileInfo($sum_file, "summary file containing all standard output", "", $log_FH); 
OutputFileInfo($cmd_file, "command file all executed commands", "", $log_FH); 

if($do_paccn_mode || $do_ntaccn_mode) { 
  open(OUT , ">" . $acclist_file) || die "ERROR unable to open $acclist_file for writing";
  print OUT $symbol . "\n";
  close(OUT);
  $acclist_file = $out_dir . "/" . $symbol . ".singleacclist";
  OutputFileInfo($acclist_file, "list file with single accession from -paccn or -ntaccn", "", $log_FH); 
}



if(defined $warn_msg) {     
  PrintToStdoutAndFile($warn_msg, $sum_FH);
}

# output banner
my $script_name = "dnaorg_fetch_dna_wrapper.pl";
my $script_desc = "Fetch DNA sequences from GenBank";
PrintToStdoutAndFile("# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -\n", $sum_FH);
PrintToStdoutAndFile("# $script_name: $script_desc\n", $sum_FH);
PrintToStdoutAndFile("# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -\n", $sum_FH);
PrintToStdoutAndFile("# command: $executable $opts_used_short$argv_1\n", $sum_FH);
PrintToStdoutAndFile(sprintf("# date:    %s\n", scalar localtime()), $sum_FH);
if($opts_used_long ne "") { 
  PrintToStdoutAndFile("$opts_used_long", $sum_FH);
}
if(! $do_acclist_mode) { 
  PrintToStdoutAndFile("# symbol:  $symbol\n", $sum_FH);
}
PrintToStdoutAndFile("#\n", $sum_FH);

# if -ntaccn or -paccn, create an 


##############################
# Stage 1: Fetch the sequences
##############################
if($do_num_mode) { 
  PrintToStdoutAndFile(sprintf("# Special mode (-num): determining number of matching %s accessions, then exiting.\n", ($do_nt_userset) ? "nucleotide" : "protein"), $sum_FH);
}
else { 
  PrintToStdoutAndFile("# Stage 1: preparing and fetching sequences\n", $sum_FH);
}
PrintToStdoutAndFile("#\n", $sum_FH);
# create column headers
my $desc_w       = length($symbol) + 68;
my $desc_dashes = "#";
for(my $i = 0; $i < $desc_w-1; $i++) { $desc_dashes .= "-"; }
PrintToStdoutAndFile(sprintf("#%-*s  %10s  %10s  %10s  %10s  %s\n", $desc_w-1, " description", "\# output", "\# lost", "\# created", "seconds", "output-file-name"), $sum_FH);
PrintToStdoutAndFile(sprintf("%s  %10s  %10s  %10s  %10s  %s\n", $desc_dashes, "----------", "----------", "----------", "----------", "-----------------"), $sum_FH);


###################################################################################
# Step 1.1: 
# Fetch all accessions that match a query to our symbol in the Gene DB. First
# try the protein database, if 0 hits, then try the nuccore database.
# UNLESS:
# if $do_plist:  Fetch all protein    accessions listed in the accession file
# if $do_ntlist: Fetch all nucleotide accessions listed in the accession file
###################################################################################
# TODO: experiment with different queries to pull in aliases/synonyms,
# possibly with command line options. Also experiment with iterating 
# this step to pull in more accessions.
my $allacc_file     = $out_dir . $symbol . ".all.acc";
my $allacconly_file = $out_dir . $symbol . ".all.acconly";
my $keep_going      = 1; # we set this to '0' after this step UNLESS 
                         # $do_not_try_nt is FALSE and we find 0 protein accessions
                         # in this step, in that case, we set $do_nt to TRUE
                         # and retry the search in 'nuccore'.

while($keep_going) { 
  my $query_symbol = $symbol;
  if($query_symbol =~ m/\~/) { 
    $warn_msg = "# WARNING: replacing ~ characters in input symbol with spaces\n";
    PrintToStdoutAndFile($warn_msg, $sum_FH);
    $query_symbol =~ s/\~/ /g;
  }
  $query = "\"$query_symbol [GENE]\"";

  if($do_acclist_mode) { 
    $cmd  = "cat $acclist_file | sort > $allacc_file";
    $desc = "Sorted_accessions_from_file_$acclist_file";
  }
  elsif($do_nt) { 
    $cmd  = "esearch -db nuccore -query $query | efetch -format acc | sort > $allacc_file";
    $desc = "Nucleotide_accessions_fetched_from_nuccore_database";
  }
  else { # default
    $cmd  = "esearch -db protein -query $query | efetch -format acc | sort > $allacc_file";
    $desc = "Protein_accessions_fetched_from_protein_database";
  }
  $nsecs  = RunCommand($cmd, $be_verbose, $cmd_FH);
  $nlines = GetNumLinesInFile($allacc_file);
  OutputFileInfo($allacc_file, $desc, $cmd, $log_FH);
  
  $cmd    = "cat $allacc_file | sed 's/\\.[0-9]*//' > $allacconly_file";
  $nsecs += RunCommand($cmd, $be_verbose, $cmd_FH);
 
  $nlost       = 0;
  $ncreated    = $nlines;
  PrintToStdoutAndFile(sprintf("%-*s  %10d  %10d  %10d  %10.1f  %s\n", $desc_w, $desc, GetNumLinesInFile($allacc_file), $nlost, $ncreated, $nsecs, $allacc_file), $sum_FH);
  if($ncreated == 0) { # this is an okay result
   if((! $do_nt) && (! $do_not_try_nt)) { # the one case in which we retry the query in nuccore
      $do_nt      = 1;
      $keep_going = 1; 
    }
    else { 
      PrintToStdoutAndFile("#\n", $sum_FH);
      if((! $do_nt_userset) && $do_nt) { 
        PrintToStdoutAndFile("# No protein or nuccore accessions fetched.\n# Exiting.\n", $sum_FH);
      }
      elsif((! $do_nt_userset) && (! $do_acclist_mode)) { 
        PrintToStdoutAndFile("# No protein accessions fetched.\n# Maybe it's a nucleotide symbol? Use -nt or -trynt to check.\n# Exiting.\n", $sum_FH);
      }
      else { 
        PrintToStdoutAndFile("# No accessions fetched.\n# Exiting.\n", $sum_FH);
      }        
      Conclude($start_secs, -1, $sum_FH, $log_FH, $cmd_FH);
      exit 0;
    }
  }
  else { # $ncreated > 0
    $keep_going = 0;
  }
}
# End of section TODO refers to

##############################################################
# Step 1.2: Remove suppressed accessions (according to idstat) 
##############################################################
my $acc_file         = $out_dir . $symbol . ".acc";
if($do_nosyn) { 
  $acc_file = $out_dir . $symbol . ".withsym.acc";
}
my $acc_file_lost    = $acc_file . ".lost";
my $acc_file_created = $acc_file . ".created";
if($do_nt) { $desc  = "Non-suppressed_nucleotide_accessions_fetched_from_nuccore_database"; }
else       { $desc  = "Non-suppressed_protein_accessions_fetched_from_protein_database"; }
my $idstat_file = $allacconly_file . ".idstat";
$cmd = "$idstat -i PUBSEQ_OS -A $allacconly_file > $idstat_file";
$nsecs = RunCommand($cmd, $be_verbose, $cmd_FH);

# parse the idstat input to determine which accessions are not suppressed
removeGivenIdstat($allacc_file, $idstat_file, $acc_file);

# which accessions were lost because they're dead or withdrawn or suppressed?
$cmd = "awk '{ print \$1 }' $acc_file | comm -2 -3 $allacc_file - > $acc_file_lost";
RunCommand($cmd, $be_verbose, $cmd_FH);

# which accessions were created relative to the allacc file? (should be none)
$cmd = "awk '{ print \$1 }' $acc_file | comm -2 -3 - $allacc_file > $acc_file_created";
RunCommand($cmd, $be_verbose, $cmd_FH);

# enforce nothing was created
($nlost, $ncreated, $errmsg) = CheckLostAndCreated($acc_file_lost, -1, $acc_file_created, 0); 
# in above subroutine call: -1 says it's okay if we lost some sequences, '0' says it's not okay if we created some sequences

OutputFileInfo($acc_file_lost,    "Protein accessions ($nlost) that are dead|withdrawn|suppressed according to idstat.", $cmd, $log_FH);
OutputFileInfo($acc_file_created, "Protein accessions ($ncreated) added by idstat processing (in error).", $cmd, $log_FH);

# output information on this step to stdout and sum file
PrintToStdoutAndFile(sprintf("%-*s  %10d  %10d  %10d  %10.1f  %s\n", $desc_w, $desc, GetNumLinesInFile($acc_file), $nlost, $ncreated, $nsecs, $acc_file), $sum_FH);
if($errmsg ne "") { die $errmsg; }

# if -ntlist or -plist was used, create a new version of the .ntlist or .plist file with
# suppressed accessions removed (only those accessions in the $acc_file will be left)
# in the original order they appeared in the .ntlist or .plist file.
if($do_plist_mode || $do_ntlist_mode) { 
  if($do_nt) { $desc  = "Non-suppressed_nucleotide_accessions_in_original_order"; }
  else       { $desc  = "Non-suppressed_protein_accessions_in_original_order"; }
  my ($seconds, $microseconds) = gettimeofday();
  my $tmp_start_secs = ($seconds + ($microseconds / 1000000.));
  # name the new output file:
  my $acclist_file_ns = $acclist_file . ".not_suppressed";
  # but we want to make sure we put it in $out_dir, so remove dir path and prepend $out_dir
  $acclist_file_ns =~ s/^.+\///;
  $acclist_file_ns = $out_dir . $acclist_file_ns;
  
  # first create a hash of the acconly accessions in $acc_file
  my %keep_H = ();
  my $accn;
  my $nlost = 0;
  open(IN, $acc_file) || die "ERROR unable to open $acc_file for reading";
  while($accn = <IN>) { 
    if($accn =~ m/\w/) { 
      chomp $accn;
      $keep_H{$accn} = 1;
    }
  }
  close(IN);
  open(IN,  $acclist_file)          || die "ERROR unable to open $acclist_file for reading"; 
  open(OUT, ">" . $acclist_file_ns) || die "ERROR unable to open $acclist_file_ns for writing";
  while($accn = <IN>) { 
    if($accn =~ m/\w/) { 
      chomp $accn;
      my $accn_only = $accn;
      stripVersion(\$accn_only);
      if((exists $keep_H{$accn}) || (exists $keep_H{$accn_only})) { 
        print OUT $accn . "\n";
      }
      else { 
        $nlost++;
      }
    }
  }
  close(OUT);

  # output information on this step to stdout and sum file
  ($seconds, $microseconds) = gettimeofday();
  my $tmp_end_secs = ($seconds + ($microseconds / 1000000.));
  PrintToStdoutAndFile(sprintf("%-*s  %10d  %10d  %10d  %10.1f  %s\n", $desc_w, $desc, GetNumLinesInFile($acclist_file_ns), $nlost, 0, ($tmp_end_secs - $tmp_start_secs), $acclist_file_ns), $sum_FH);
  if($errmsg ne "") { die $errmsg; }

} # end of 'if($do_plist_mode || $do_ntlist_mode)'

###################################################################################
# Optional Step: If -nosyn: remove accessions which have $symbol as a synonym
# only keeping those for which $symbol is the primary symbol.
###################################################################################
if($do_nosyn) { 
  my $tmp_start_secs = ($seconds + ($microseconds / 1000000.));

  my $database = ($do_nt) ? "nuccore" : "protein";
  my $new_acc_file         = $out_dir . $symbol . ".acc";
  my $new_acc_file_lost    = $new_acc_file . ".lost";
  my $new_acc_file_created = $new_acc_file . ".created";

  my $gene_gene_file     = $out_dir . $symbol . ".gene.gene";
  $cmd = "cat $acc_file | epost -db $database -format acc | efetch -format gpc | xtract -insd gene gene | grep . | sort > $gene_gene_file";
  RunCommand($cmd, $be_verbose, $cmd_FH);
  OutputFileInfo($gene_gene_file,    "Gene:gene values for all accessions in $acc_file (this file is used to remove accessions for which symbol is a synonym).", $cmd, $log_FH);

  my $gene_gene_syn_file = $out_dir . $symbol . ".gene.gene_syn";
  $cmd = "cat $acc_file | epost -db $database -format acc | efetch -format gpc | xtract -insd gene gene_synonym | grep . | sort > $gene_gene_syn_file";
  RunCommand($cmd, $be_verbose, $cmd_FH);
  OutputFileInfo($gene_gene_syn_file,    "Gene:gene_synonym values for all accessions in $acc_file (this file is used to remove accessions for which symbol is a synonym).", $cmd, $log_FH);

  my $cds_gene_file      = $out_dir . $symbol . ".cds.gene";
  $cmd = "cat $acc_file | epost -db $database -format acc | efetch -format gpc | xtract -insd CDS gene | grep . | sort > $cds_gene_file";
  RunCommand($cmd, $be_verbose, $cmd_FH);
  OutputFileInfo($cds_gene_file,    "CDS:gene values for all accessions in $acc_file (this file is used to remove accessions for which symbol is a synonym).", $cmd, $log_FH);

  my $cds_gene_syn_file  = $out_dir . $symbol . ".cds.gene_syn";
  $cmd = "cat $acc_file | epost -db $database -format acc | efetch -format gpc | xtract -insd CDS gene_synonym | grep . | sort > $cds_gene_syn_file";
  RunCommand($cmd, $be_verbose, $cmd_FH);
  OutputFileInfo($cds_gene_syn_file,    "CDS:gene_synonym values for all accessions in $acc_file (this file is used to remove accessions for which symbol is a synonym).", $cmd, $log_FH);

  # Now, each accession should have $symbol as either a gene or a gene_synonym in >= 1 of 
  # $gene_gene_file/$gene_gene_syn_file or $cds_gene_file/$cds_gene_syn_file
  # parse them 
  my %keepme_H = ();
  foreach my $file ($gene_gene_file, $cds_gene_file) { 
    open(IN, $file) || die "ERROR unable to open $file";
    while(my $line = <IN>) { 
      chomp $line;
      my @elA = split(/\t/, $line);
      my $nel = scalar(@elA);
      if($nel > 2) { 
        die "ERROR more than one match in $file for $elA[0]"; 
      }
      elsif($nel == 2) { 
        $elA[1] =~ s/\;$//; # remove trailing semicolon if one exists;
        my $cap_el = Capitalize($elA[1]);
        if($cap_el eq $cap_no_tilde_symbol) { 
          $keepme_H{$elA[0]} = 1; 
        }
      }
      elsif($nel != 1) { 
        die "ERROR unexpected number of tokens in $file for $elA[0]";
      }
    }
    close(IN);
  }

  my %iamsyn_H = ();
  foreach my $file ($gene_gene_syn_file, $cds_gene_syn_file) { 
    open(IN, $file) || die "ERROR unable to open $file";
    while(my $line = <IN>) { 
      chomp $line;
      my @elA = split(/\s+/, $line);
      my $nel = scalar(@elA);
      for(my $i = 1; $i < $nel; $i++) { # note we start at '1', not '0'
        $elA[$i] =~ s/\;$//; # remove trailing semicolon if one exists;
        my $cap_el = Capitalize($elA[$i]);
        if($cap_el eq $cap_no_tilde_symbol) { 
          $iamsyn_H{$elA[0]} = 1; 
        }
      }
    }
    close(IN);
  }

  # create a file that is a list of accessions for which $symbol is primary
  open(IN, $acc_file) || die "ERROR unable to open $acc_file for reading";
  open(OUT, ">" . $new_acc_file) || die "ERROR unable to open $new_acc_file for writing";
  while(my $acc = <IN>) { 
    chomp $acc;
    if($keepme_H{$acc}) { 
      print OUT $acc . "\n"; 
    }
    elsif(! exists $iamsyn_H{$acc}) { 
      # this is okay, just skip it
      ;
      # alternative approach would be:
      # die "ERROR $acc does not have $symbol as a primary symbol or an alias"; 
    }
  }
  close(IN);
  close(OUT);
  ($seconds, $microseconds) = gettimeofday();
  my $tmp_end_secs = ($seconds + ($microseconds / 1000000.));

  # which accessions were lost because they have symbol as a synonym?
  $cmd = "awk '{ print \$1 }' $acc_file | comm -2 -3 $new_acc_file - > $new_acc_file_lost";
  RunCommand($cmd, $be_verbose, $cmd_FH);

  # which accessions were created here (should be none)
  $cmd = "awk '{ print \$1 }' $acc_file | comm -2 -3 - $new_acc_file > $new_acc_file_created";
  RunCommand($cmd, $be_verbose, $cmd_FH);

  # enforce nothing was created
  ($nlost, $ncreated, $errmsg) = CheckLostAndCreated($new_acc_file_lost, -1, $new_acc_file_created, 0); 
  # in above subroutine call: -1 says it's okay if we lost some sequences, '0' says it's not okay if we created some sequences

  OutputFileInfo($new_acc_file_lost,    "Accessions ($nlost) that have symbol as a synonym.", $cmd, $log_FH);
  OutputFileInfo($new_acc_file_created, "Accessions ($ncreated) added by processing synonyms.", $cmd, $log_FH);

  $desc = "Accessions_for_which_symbol_is_a_primary_symbol_and_not_a_synonymn";
  PrintToStdoutAndFile(sprintf("%-*s  %10d  %10d  %10s  %10.1f  %s\n", $desc_w, $desc, GetNumLinesInFile($new_acc_file), $nlost, $ncreated, $tmp_end_secs - $tmp_start_secs, $new_acc_file), $sum_FH);

  $acc_file = $new_acc_file; # rename $acc_file for subsequent steps
} # end of 'if($do_nosyn)'

# create the acconly file, which is the acc file minus versions
my $acconly_file = $acc_file . "only";
$cmd = "cat $acc_file | sed 's/\\.[0-9]*//' > $acconly_file";
RunCommand($cmd, $be_verbose, $cmd_FH);
OutputFileInfo($acconly_file, "Same accessions as $acc_file, but with versions removed (if any)", $cmd, $log_FH);

#########################################################################################
# Exit if we're in number-only mode, gene info mode, or we don't have any accessions left
#########################################################################################
if($do_num_mode) { 
  PrintToStdoutAndFile("#\n", $sum_FH);
  my $num_records = GetNumLinesInFile($acc_file) - $nlost;
  if($do_nt) { 
    PrintToStdoutAndFile("Number_of_nucleotide_records: $num_records\n", $sum_FH);
  }
  else { 
    PrintToStdoutAndFile("Number_of_protein_records: $num_records\n", $sum_FH);
  }    
  Conclude($start_secs, $do_nt, $sum_FH, $log_FH, $cmd_FH);
  exit 0;
}
if(GetNumLinesInFile($acc_file) == 0) { # this is an okay result
  PrintToStdoutAndFile("#\n", $sum_FH);
  PrintToStdoutAndFile("# No (non-suppressed) accessions fetched. Exiting.\n", $sum_FH);
  Conclude($start_secs, $do_nt, $sum_FH, $log_FH, $cmd_FH);
  exit 0;
}
if($do_geneinfo_mode) { 
  # we need to formulate queries for every N accessions, we can't use epost like we do elsewhere to post entire list of accessions
  my $gene_cmd;
  my $tax_cmd;
  my $gene_cmd_concat = ""; # concatenated commands for fetching gene info
  my $tax_cmd_concat  = ""; # concatenated commands for fetching taxid info
  my $naccn = 0;
  my $output_char;
  my $nsecs = 0.;

  # print first line of each file
  my $geneinfo_file = $out_dir . $symbol . ".geneinfo";
  my $taxinfo_file = $out_dir . $symbol . ".taxinfo";
  my $taxinfo_file_created = $taxinfo_file . ".created";
  my $taxinfo_file_lost    = $taxinfo_file . ".lost";
  my $geneinfo_FH;
  my $taxinfo_FH;

  open($geneinfo_FH, ">" . $geneinfo_file) || die "ERROR unable to open $geneinfo_file";
  print $geneinfo_FH "#taxid gene-id gene-name gene-aliases accession chrStart chrEnd (sep character is a tab, lines beginning with '#' are comments)\n";

  open($taxinfo_FH, ">" . $taxinfo_file) || die "ERROR unable to open $taxinfo_file";
  print $taxinfo_FH "#accession db_xref (sep character is a tab)\n";

  my $tmp_out_root   = $out_dir . $symbol . ".tmp";
  my @accn_A = ();  # array of accessions we're currently working on
  my $ncmd = 0;

  open(ACC, $acc_file) || die "ERROR unable to open $acc_file"; 
  while(my $line = <ACC>) { 
    if($line !~ m/\#/) { 
      chomp $line;
      my $accn = $line;
      $accn =~ s/^\s+//; # remove leading whitespace
      $accn =~ s/\s+$//; # remove trailing whitespace
      $naccn++;
      push(@accn_A, $accn);
      if($naccn == $ngene) { # $ngene was set above as $df_ngene, or possibly with -ngene 
        ($tax_cmd, $gene_cmd, $nsecs) = updateTaxAndGeneinfo(\@accn_A, $do_allgene, $tmp_out_root, $taxinfo_FH, $geneinfo_FH, $cmd_FH);
        $ncmd++;
        $tax_cmd_concat .= $tax_cmd . ";";
        $gene_cmd_concat .= $gene_cmd .";";
        $naccn = 0;
        @accn_A = ();
      }
    }
  } # end of 'while($line = <IN>)'
  # and don't forget to do the remaining accessions not yet processed
  if($naccn > 0) { 
    ($tax_cmd, $gene_cmd, $nsecs) = updateTaxAndGeneinfo(\@accn_A, $do_allgene, $tmp_out_root, $taxinfo_FH, $geneinfo_FH, $cmd_FH);
    $ncmd++;
    $tax_cmd_concat .= $tax_cmd . ";";
    $gene_cmd_concat .= $gene_cmd .";";
  }
  close($geneinfo_FH);
  close($taxinfo_FH);
  OutputFileInfo($taxinfo_file,   "Tabular Taxid info for all accessions in $acc_file", $tax_cmd_concat, $log_FH);
  OutputFileInfo($geneinfo_file,  "Tabular Gene info for all accessions in $acc_file", $gene_cmd_concat, $log_FH);

  # determine any lost/created tax ids 
  $cmd = "awk '{ print \$1 }' $taxinfo_file | grep -v ^\# | sed 's/\\.[0-9]*//' | sort | comm -2 -3 - $acconly_file > $taxinfo_file_created";
  RunCommand($cmd, $be_verbose, $cmd_FH);
  $cmd = "awk '{ print \$1 }' $taxinfo_file | grep -v ^\# | sed 's/\\.[0-9]*//' | sort | comm -2 -3 $acconly_file - > $taxinfo_file_lost";
  RunCommand($cmd, $be_verbose, $cmd_FH);
  
  ($nlost, $ncreated, $errmsg) = CheckLostAndCreated($taxinfo_file_lost, 0, $taxinfo_file_created, 0); # '0' are max allowed lines in each of these files
  
  $desc = "Tax_ids_for_all_accessions";
  PrintToStdoutAndFile(sprintf("%-*s  %10d  %10s  %10s  %10.1f  %s\n", $desc_w, $desc, GetNumLinesInFile($taxinfo_file)-1, $nlost, $ncreated, $nsecs, $taxinfo_file), $sum_FH);

  if($do_allgene) { 
    $desc = "Gene_database_records_which_link_to_all_accessions";
  }
  else { 
    $desc = "Gene_database_records_which_exist_in_all_accessions";
  }
  PrintToStdoutAndFile(sprintf("%-*s  %10d  %10s  %10s  %10.1f  %s\n", $desc_w, $desc, GetNumLinesInFile($geneinfo_file)-($ncmd+1), "N/A", "N/A", $nsecs, $geneinfo_file), $sum_FH);

  if($errmsg ne "") { die $errmsg; }
} # end of 'if($do_geneinfo_mode)'

my $database = ($do_nt) ? "nuccore" : "protein";
if($do_ftable_mode || $do_matpept_mode) { 
  # create a file with total lengths of each accession
  my $len_file  = $out_dir . $symbol . ".length";
  my $len_file_created = $len_file . ".created";
  my $len_file_lost    = $len_file . ".lost";
  $cmd = "cat $acconly_file | epost -db $database -format acc | efetch -format gpc | xtract -insd INSDSeq_length | grep . | sort > $len_file";
  $nsecs = RunCommand($cmd, $be_verbose, $cmd_FH);
  OutputFileInfo($len_file, "File with full lengths of sequence records with accessions in $acc_file", $cmd, $log_FH);

  # determine any lost/created accessions in the length file
  $cmd = "awk '{ print \$1 }' $len_file | sed 's/\\.[0-9]*//' | comm -2 -3 - $acconly_file > $len_file_created";
  $nsecs += RunCommand($cmd, $be_verbose, $cmd_FH);
  $cmd = "awk '{ print \$1 }' $len_file | sed 's/\\.[0-9]*//' | comm -1 -3 - $acconly_file > $len_file_lost";
  $nsecs += RunCommand($cmd, $be_verbose, $cmd_FH);
  ($nlost, $ncreated, $errmsg) = CheckLostAndCreated($len_file_lost, 0, $len_file_created, 0); # '0' are max allowed lines in each of these files
  $desc = "Lengths_of_all_accessions";
  PrintToStdoutAndFile(sprintf("%-*s  %10d  %10s  %10s  %10.1f  %s\n", $desc_w, $desc, GetNumLinesInFile($len_file), $nlost, $ncreated, $nsecs, $len_file), $sum_FH);
  if($errmsg ne "") { die $errmsg; }

  if($do_ftable_mode) { 
    # create the feature table file
    my $ft_file  = $out_dir . $symbol . ".ftable";
    my $ft_acc_file     = $ft_file . ".acc";
    my $ft_file_created = $ft_file . ".created";
    my $ft_file_lost    = $ft_file . ".lost";
    $cmd = "cat $acconly_file | epost -db $database -format acc | efetch -format ft > $ft_file";
    $nsecs = RunCommand($cmd, $be_verbose, $cmd_FH);
    OutputFileInfo($ft_file, "Feature table for all accessions in $acc_file", $cmd, $log_FH);
    
    # now create a file with a list of the accessions that are listed in the ftable,
    # so we can verify that we have an ftable for each accession
    #  $cmd = "grep ^\\> $ft_file | sed -e \'s/^>Feature\\s*[A-Za-z]*|//\' | sed -e \'s/|//g\' | sed 's/\\.[0-9]*//' | sort > $ft_acc_file";
    $cmd = "grep ^\\> $ft_file | sed -e \'s/^>Feature\\s*[A-Za-z]*|//\' | sed -e \'s/|.*//\' | sed 's/\\.[0-9]*//' | sort > $ft_acc_file";
    $nsecs += RunCommand($cmd, $be_verbose, $cmd_FH);
    OutputFileInfo($ft_acc_file,  "All accessions listed in feature table in $ft_file (used for quality checking only)", $cmd, $log_FH);

    # determine any lost/created accessions
    $cmd = "comm -2 -3 $ft_acc_file $acconly_file > $ft_file_created";
    $nsecs += RunCommand($cmd, $be_verbose, $cmd_FH);
    $cmd = "comm -1 -3 $ft_acc_file $acconly_file > $ft_file_lost";
    $nsecs += RunCommand($cmd, $be_verbose, $cmd_FH);
    
    ($nlost, $ncreated, $errmsg) = CheckLostAndCreated($ft_file_lost, 0, $ft_file_created, 0); # '0' are max allowed lines in each of these files
    
    $desc = "Feature_table_for_all_accessions";
    PrintToStdoutAndFile(sprintf("%-*s  %10d  %10s  %10s  %10.1f  %s\n", $desc_w, $desc, GetNumLinesInFile($ft_file)-1, $nlost, $ncreated, $nsecs, $ft_file), $sum_FH);
  } # end of 'if($do_ftable_mode)'

  if($do_matpept_mode) { 
    # create the mat_peptide file
    my $mp_file  = $out_dir . $symbol . ".mat_peptide";
    my $mp_acc_file     = $mp_file . ".acc";
    my $mp_file_created = $mp_file . ".created";
    my $mp_file_lost    = $mp_file . ".lost";
    # OLD AS OF 10.15.15: $cmd = "cat $acconly_file | epost -db $database -format acc | efetch -format gpc | xtract -insd mat_peptide INSDInterval_from INSDInterval_to product > $mp_file";
    $cmd = "cat $acconly_file | epost -db $database -format acc | efetch -format gpc | xtract -insd mat_peptide INSDFeature_location product > $mp_file";
    $nsecs = RunCommand($cmd, $be_verbose, $cmd_FH);
    OutputFileInfo($mp_file, "File with mat_peptide info for all accessions in $acc_file", $cmd, $log_FH);
    
    # now create a file with a list of the accessions that are listed in the mat_peptide,
    # so we can check which accessions we have info for and which we do not
    $cmd = "grep . $mp_file | awk '{ print \$1 }' | sed 's/\\.[0-9]*//' | sort | uniq > $mp_acc_file";
    $nsecs += RunCommand($cmd, $be_verbose, $cmd_FH);
    OutputFileInfo($mp_acc_file,  "All accessions for which mat_peptide info exists in $mp_file", $cmd, $log_FH);
    
    # determine any lost/created accessions
    $cmd = "comm -2 -3 $mp_acc_file $acconly_file > $mp_file_created";
    $nsecs += RunCommand($cmd, $be_verbose, $cmd_FH);
    $cmd = "comm -1 -3 $mp_acc_file $acconly_file > $mp_file_lost";
    $nsecs += RunCommand($cmd, $be_verbose, $cmd_FH);
    
    ($nlost, $ncreated, $errmsg) = CheckLostAndCreated($mp_file_lost, -1, $mp_file_created, 0); # '-1' indicates there's no restrictions on how many lines are allowed in the lost file, '0' says none can be in the created file
    
    $desc = "Mat_peptide_info_for_all_accessions";
    PrintToStdoutAndFile(sprintf("%-*s  %10d  %10s  %10s  %10.1f  %s\n", $desc_w, $desc, GetNumLinesInFile($mp_acc_file), $nlost, $ncreated, $nsecs, $mp_file), $sum_FH);
  }
  if($errmsg ne "") { die $errmsg; }
} # end of 'if($do_ftable_mode || $do_matpept_mode)'

if($do_geneinfo_mode || $do_ftable_mode || $do_matpept_mode) { 
  Conclude($start_secs, $do_nt, $sum_FH, $log_FH, $cmd_FH);
  exit 0;
}

###################################################################################
# Step 1.3: Determine the CDS sequences/coordinates that correspond to our proteins 
# from step 1.1 and fetch them.
###################################################################################
my $efa_file         = $out_dir . $symbol . ".efa";
my $efa_file_lost    = $efa_file . ".lost";
my $efa_file_created = $efa_file . ".created";
if(! $do_nt) { 
  $cmd = "cat $acc_file | epost -db protein -format acc | efetch -format gpc | xtract -insd CDS coded_by codon_start | grep . | sort ";
  if($do_old) { 
    $cmd .= " | $exec_dir/coded_by2extract_fasta_multi_exon.pl > $efa_file";
  }
  else { 
    $cmd .= " > $efa_file";
  }
  $desc = "Protein_accessions_that_have_CDS_annotation";
  $nsecs = RunCommand($cmd, $be_verbose, $cmd_FH);
  OutputFileInfo($efa_file, $desc, $cmd, $log_FH);

# which accessions were lost relative to the previous step? (we can lose some, these should have no CDS)
  if($do_old) { 
    $cmd = "awk '{ print \$NF }' $efa_file | comm -2 -3 $acc_file - > $efa_file_lost";
  }
  else {
    $cmd = "awk '{ print \$1 }' $efa_file | comm -2 -3 $acc_file - > $efa_file_lost";
  }
  RunCommand($cmd, $be_verbose, $cmd_FH);

# which accessions were created relative to the previous step? (should be none)
  if($do_old) { 
    $cmd = "awk '{ print \$NF }' $efa_file | comm -2 -3 - $acc_file > $efa_file_created";
  }
  else { 
    $cmd = "awk '{ print \$1 }' $efa_file | comm -2 -3 - $acc_file > $efa_file_created";
  }
  RunCommand($cmd, $be_verbose, $cmd_FH);

# enforce nothing was lost or created
  ($nlost, $ncreated, $errmsg) = CheckLostAndCreated($efa_file_lost, -1, $efa_file_created, 0); 
# in above subroutine call: -1 says it's okay if we lost some sequences, '0' says it's not okay if we created some sequences

  OutputFileInfo($efa_file_created, "New protein accessions ($ncreated) returned by CDS retrieval (in error).", $cmd, $log_FH);
  OutputFileInfo($efa_file_lost, "Protein accessions ($nlost) without CDS annotation.", $cmd, $log_FH);

# output information on this step to stdout and sum file
  PrintToStdoutAndFile(sprintf("%-*s  %10d  %10d  %10d  %10.1f  %s\n", $desc_w, $desc, GetNumLinesInFile($efa_file), $nlost, $ncreated, $nsecs, $efa_file), $sum_FH);
  if($errmsg ne "") { die $errmsg; }
} # end of 'if(! $do_nt)'

#############################################################
#EPN, Wed Mar 25 14:22:09 2015
# Experimental code for fetching UniProt CDS sequences, currently only run if -up
if((! $do_nt) && $do_uniprot_xref) { 
################################
# Try to fetch any uniprot sequences that we can 
  ($seconds, $microseconds) = gettimeofday();
  my $tmp_start_secs = ($seconds + ($microseconds / 1000000.));
  $nsecs = RunCommand($cmd, $be_verbose, $cmd_FH);
  my $up_root  = $out_dir . $symbol . ".up.out";
  my $tmpfile0 = $up_root . ".0";
  my $tmpfile1 = $up_root . ".1";
  my $tmpfile2 = $up_root . ".2";
  my $tmpfile3 = $up_root . ".3";
  my $tmpfile4 = $up_root . ".4";
  # get the protein lengths
  $cmd = "cat $efa_file_lost | grep -v ^WP\_ | grep -v ^YP\_ | epost -db protein -format acc | efetch -format gpc | xtract -insd INSDSeq_length | grep . | sort > $tmpfile0";
  RunCommand($cmd, $be_verbose, $cmd_FH);
  # parse the lengths
  open(IN, $tmpfile0) || die "ERROR unable to open $tmpfile0";
  my %prot_len_H = (); # lengths of each protein, value: protein accession.version, value: length in AA
  while(my $line = <IN>) { 
    chomp $line;
    my ($prot, $len) = split(/\s+/, $line);
    if($len !~ m/^\d+$/) { die "ERROR unrecognized protein length information: $line\n"; }
    $prot_len_H{$prot} = $len;
    # print "$prot $len\n";
  }
  close(IN);
  
  # get the xrefs lines
  $cmd = "cat $efa_file_lost | grep -v ^WP\_ | grep -v ^YP\_ | epost -db protein -format acc | efetch -format gpc | xtract -insd INSDSeq_source-db | grep xrefs | sort > $tmpfile1";
  RunCommand($cmd, $be_verbose, $cmd_FH);
  # parse the information we gathered on xrefs for the possible aliens in $efa_file_lost_maybe_alien_xrefs
  # Based on Alejandro's /panfs/pan1.be-md.ncbi.nlm.nih.gov/dnaorg/2015.03.19/get_xrefs_row/extract_xref_seqs.pl
  # we store the xrefs in a hash that maps the putative nucleotide accessions to their mother protein
  # accessions
  open(IN,  "<$tmpfile1") || die "ERROR unable to open $tmpfile1 for reading";
  open(OUT, ">$tmpfile2") || die "ERROR unable to open $tmpfile2 for writing";
  my $line;
  my %xref2prot_HA = (); # hash of arrays, key: xref accession, value: array of protein accessions that have this xref accession

  while(defined($line = <IN>)) { 
  # O02771.1	UniProtKB: locus SMN_CANFA, accession O02771; class: standard. created: Jul 15, 1998. sequence updated: Jul 1, 1997. annotation updated: Mar 4, 2015. xrefs: U50746.1, AAB58318.1, NP_001003226.1 xrefs (non-sequence databases): UniGene:Cfa.191, ProteinModelPortal:O02771, SMR:O02771, STRING:9615.ENSCAFP00000031438, PaxDb:O02771, PRIDE:O02771, GeneID:403896, KEGG:cfa:403896, CTD:39844, eggNOG:NOG296671, HOGENOM:HOG000232199, HOVERGEN:HBG000211, InParanoid:O02771, KO:K13129, Reactome:REACT_234653, NextBio:20817388, Proteomes:UP000002254, GO:0015030, GO:0005737, GO:0005829, GO:0097504, GO:0005654, GO:0005634, GO:0032797, GO:0034719, GO:0030018, GO:0003723, GO:0007409, GO:0007019, GO:0000387, InterPro:IPR010304, InterPro:IPR002999, Pfam:PF06003, SMART:SM00333, PROSITE:PS50304
    chomp $line;
    my $paccn = $line;
    $paccn =~ s/\s+.*$//; # remove everything after first space (giving O02771.1)
    my $xrefs = $line;
    my ($xref_info) = ($line=~m/xrefs:(.*?)xrefs/);
    my @id_A = split /\,/, $xref_info;
    foreach my $id (@id_A) { 
      my ($no_space_id) = ($id=~m/([A-Z0-9_\.]+)/);
      # we need to ensure we can map $no_space_id to exactly one protein accession
      if(! exists $xref2prot_HA{$no_space_id}) { 
        @{$xref2prot_HA{$no_space_id}} = (); # initialize
      }
      push(@{$xref2prot_HA{$no_space_id}}, $paccn);
      print OUT $paccn . " " . $no_space_id . "\n";
    }
  }
  close(IN);
  close(OUT);  
  
  # $tmpfile2 has a list of possible nucleotide accessions corresponding to the uniprot accessions
  # note this step uses protein accessions
  $cmd = "awk '{ print \$1 }' $tmpfile2 | epost -db protein -format acc | efetch -format gpc | xtract -insd Protein gene gene_synonym | grep . | sort > $tmpfile3";
  RunCommand($cmd, $be_verbose, $cmd_FH);
  
  # read $tmpfile3 to get the gene names for each protein accession:
  open(IN, $tmpfile3) || die "ERROR unable to open $tmpfile3";
  my %prot2gene_HA = (); # value is a protein accession.version, value is an array of 'gene' and 'gene_synonym's for that protein
  while(my $line = <IN>) { 
    chomp $line;
    my (@elA) = split(/\s+/, $line);
    if(scalar(@elA) >= 2) { 
      my $prot = shift @elA; # first token is protein accession
      if(exists ($prot2gene_HA{$prot})) { die "ERROR read protein $prot twice in $tmpfile3"; }
      $prot2gene_HA{$prot} = [@elA]; # remaining tokens are 'gene's or 'gene_synonym's
    }
    else { # no 'gene's or 'gene synonym's for this protein, currently we do nothing, just skip it
      ;
    }
  }
  close(IN);
  
  # now for each potential nucleotide accession for protein X in $tmpfile2, try to find a CDS for the 
  # gene that protein X codes for
  $cmd = "awk '{ print \$2 }' $tmpfile2 | epost -db nuccore -format acc | efetch -format gpc | xtract -pattern INSDSeq -ACCN INSDSeq_accession-version -group INSDFeature -match INSDFeature_key:CDS -pfx \"\\n\" -element \"\&ACCN\" -block INSDFeature -match INSDFeature_key:CDS -element INSDFeature_location -subset INSDQualifier -match INSDQualifier_name:gene -element INSDQualifier_value | grep . | cat - > $tmpfile4";
  RunCommand($cmd, $be_verbose, $cmd_FH);
  
  # now parse $tmpfile4
  # example lines:
  # CU329670.1	complement(<1..5662)	tlh1	
  # CU329670.1	complement(5726..6331)	
  # CU329670.1	complement(join(822332..822429,822489..822849))	smn1	
  my $up_efa_file = $up_root . ".uniprot-matches.efa";
  # EPN, Wed Apr  8 08:53:28 2015: decided to output all matches to $up_efa_file, we'll test them all with esl-test-cds-translate-vs-fetch.pl 
  # we may want to undo this in the future though
  #my $tmpfile6 = $up_root . ".notby3.6";
  #my $tmpfile7 = $up_root . ".lendiffer.7";
  open(IN,   "<$tmpfile4") || die "ERROR unable to open $tmpfile4 for reading";
  open(OUT5, ">$up_efa_file") || die "ERROR unable to open $up_efa_file for writing";
  #open(OUT6, ">$tmpfile6") || die "ERROR unable to open $tmpfile6 for writing";
  #open(OUT7, ">$tmpfile7") || die "ERROR unable to open $tmpfile7 for writing";
  my %coords_for_prot_HA = (); # coord strings for we've found for each protein, key: protein accession, value: array of matching CDS coord strings
  while($line = <IN>) { 
    chomp $line;
    my ($nt, $coords, $gene) = split(/\s+/, $line);
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
    $coords = $new_coords;
    if(defined $gene) { # this indicates we found INSDQualifier_name 'gene' for this NT accession
      foreach my $prot (@{$xref2prot_HA{$nt}}) { 
        if(exists $prot2gene_HA{$prot}) {
          for(my $g = 0; $g < scalar(@{$prot2gene_HA{$prot}}); $g++) { 
            if($gene eq $prot2gene_HA{$prot}[$g]) { # match
              my $already_printed = 0;
              # determine CDS length from coords string
              my $cds_len  = LengthFromCoords($coords);
              if(! exists $prot_len_H{$prot}) { 
                die "ERROR protein $prot somehow created in error";
              }
              my $prot_len = $prot_len_H{$prot};
              if(! exists $coords_for_prot_HA{$prot}) { 
                # initialize
                @{$coords_for_prot_HA{$prot}} = (); 
              }
              else { # search for a match in the existing coords strings for this protein
                for(my $z = 0; $z < scalar(@{$coords_for_prot_HA{$prot}}); $z++) { 
                  if($coords_for_prot_HA{$prot}[$z] eq $coords) { 
                    $already_printed = 1;
                  }
                }
              }
              if(! $already_printed) { 
                # currently, output all matches to OUT5
                print OUT5 ("$prot\t$coords\n");
                # THIS BLOCK ALTERNATIVELY OUTPUTS TO 3 DIFFERENT FILES
                ## check that the lengths make sense
                #if($cds_len % 3 == 0) { # should be...
                #  if((sprintf("%d", ($cds_len / 3)) eq $prot_len) || # matches exactly
                #    (sprintf("%d", ($cds_len / 3)) eq ($prot_len+1))) { # CDS presumably includes stop codon
                #    print OUT5 ("$prot\t$coords\t$prot_len\t$cds_len\n");
                #  }
                #  else { 
                #    print OUT7 ("$prot\t$prot_len\t$nt\t$cds_len\t$coords\n");
                #  }
                #}
                #else { # cds_len a factor of 3...
                #  print OUT6 ("$prot\t$coords\t$prot_len\t$cds_len\n");
                #}
                # END OF BLOCK
                push(@{$coords_for_prot_HA{$prot}}, $coords);
              }
            }
          }
        }
      }
    }
  }
  close(OUT);
  close(IN);
  ($seconds, $microseconds) = gettimeofday();
  my $tmp_end_secs = ($seconds + ($microseconds / 1000000.));
  
  # determine how many proteins we output >= 1 CDS for in $up_efa_file, and how many we didn't
  my $nxref = GetNumLinesInFile($tmpfile1);
  $nlost = $nxref - scalar(keys %coords_for_prot_HA);

  $desc = "Following_xrefs_to_find_CDS_for_accessions_with_direct_CDS_links";
  PrintToStdoutAndFile(sprintf("%-*s  %10d  %10d  %10s  %10.1f  %s\n", $desc_w, $desc, $nxref, $nlost, "?", $tmp_end_secs - $tmp_start_secs, $up_efa_file), $sum_FH);

} # end of 'if($do_uniprot_xref)'
# End of experimental code block for fetching UniProt CDS sequences, currently NOT run
################################
    
####################################################################################
# Step 1.4: Convert the extract_fasta_multi_exon input file to idfetch input format.
# We need to do this because idfetch can't read from a pipe, if it could we 
# wouldn't need to create this file.
####################################################################################
my $idfetch_file = $out_dir . $symbol . ".idfetch";
if($do_old) { 
  $cmd  = "awk '{ print \$1 }' $efa_file | sort | uniq > $idfetch_file";
  $desc = "Unique_CDS_source_sequences_(idfetch_input_file)";
  $nsecs = RunCommand($cmd, $be_verbose, $cmd_FH);
  OutputFileInfo($idfetch_file, $desc, $cmd, $log_FH); 
# don't output a summary of this step 
}

########################################################
# Step 1.5: Extract the sequences using esl-fetch-cds.pl
########################################################
my $fa_file         = $out_dir . $symbol . ".fa";
my $fa_file_names   = $fa_file . ".names";
my $fa_file_lost    = $fa_file . ".lost";
my $fa_file_created = $fa_file . ".created";
my $cds_or_nt       = ($do_nt) ? "nucleotide" : "CDS";
$desc  = "CDS_sequences_in_FASTA_format";

if($do_old) { 
  $cmd = "$idfetch -t 5 -c 1 -G $idfetch_file | $exec_dir/id_fasta.pl | $exec_dir/extract_fasta_multi_exon $efa_file > $fa_file";
}
elsif($do_nt) { 
  $cmd  = "$idfetch -t 5 -c 1 -G $acc_file | $exec_dir/id_fasta.pl > $fa_file";
  $desc = "Nucleotide_sequences_in_FASTA_format";
}
else { # normal case
  $cmd = "$fetch_cds -odir $out_dir $efa_file > $fa_file";
}
$nsecs = RunCommand($cmd, $be_verbose, $cmd_FH);
OutputFileInfo($fa_file, $desc, $cmd, $log_FH);

# get list of protein accessions we have CDS sequences (or just nucleotide sequences if $do_nt) for in fasta file
# we'll use this to make sure we fetched all of them
if($do_old) { 
  $cmd = "grep ^\\> $fa_file | sed 's/^>//' | awk '{ print \$1 }' | sed 's/^.*\://' | sort | sed 's/\\.[0-9]*//' > $fa_file_names";
}
else {
  $cmd = "grep ^\\> $fa_file | sed 's/^>//' | awk '{ print \$1 }' | sed 's/\:.*\$//' | sort | sed 's/\\.[0-9]*//' > $fa_file_names";
}
RunCommand($cmd, $be_verbose, $cmd_FH);
OutputFileInfo($fa_file_names, "Protein accessions for each $cds_or_nt sequence from FASTA file from $fa_file.", $cmd, $log_FH);

# which accessions were lost relative to the previous step? (should be none)
if($do_old) { 
  $cmd = "awk '{ print \$NF }' $efa_file | sed 's/\\.[0-9]*//' | comm -2 -3 - $fa_file_names > $fa_file_lost";
}
else { 
  if($do_nt) { $cmd = "comm -2 -3 $acconly_file $fa_file_names > $fa_file_lost"; }
  else        { $cmd = "awk '{ print \$1 }' $efa_file | sed 's/\\.[0-9]*//' | comm -2 -3 - $fa_file_names > $fa_file_lost"; }
}
RunCommand($cmd, $be_verbose, $log_FH);
OutputFileInfo($fa_file_lost, "$cds_or_nt sequences not fetched", $cmd, $log_FH);

# which accessions were created relative to the previous step? (should be none)
if($do_old) { 
  $cmd = "awk '{ print \$NF }' $efa_file | sed 's/\\.[0-9]*//' | comm -2 -3 $fa_file_names - > $fa_file_created";
}
else { 
  if($do_nt) { $cmd = "comm -2 -3 $fa_file_names $acconly_file > $fa_file_created"; }
  else        { $cmd = "awk '{ print \$1 }' $efa_file | sed 's/\\.[0-9]*//' | comm -2 -3 $fa_file_names - > $fa_file_created"; }
}
RunCommand($cmd, $be_verbose, $log_FH);
OutputFileInfo($fa_file_created, "New $cds_or_nt sequences created when fetching sequences (in error)", $cmd, $log_FH);

# enforce nothing was lost or created
($nlost, $ncreated, $errmsg) = CheckLostAndCreated($fa_file_lost, 0, $fa_file_created, 0); # '0' are max allowed lines in each of these files
if($errmsg ne "") { die $errmsg; }

# output information on this step to stdout and sum file
PrintToStdoutAndFile(sprintf("%-*s  %10d  %10d  %10d  %10.1f  %s\n", $desc_w, $desc, GetNumLinesInFile($fa_file_names), $nlost, $ncreated, $nsecs, $fa_file), $sum_FH);


###################################
# Stage 2: Validate Stage 1 results
###################################
PrintToStdoutAndFile("#\n", $sum_FH);
PrintToStdoutAndFile("# Stage 2: validating stage 1 results\n", $sum_FH);
PrintToStdoutAndFile("#\n", $sum_FH);
PrintToStdoutAndFile(sprintf("#%-*s  %10s  %10s  %10s  %10s  %s\n", $desc_w-1, " description", "\# output", "\# lost", "\# created", "seconds", "output-file-name"), $sum_FH);
PrintToStdoutAndFile(sprintf("%s  %10s  %10s  %10s  %10s  %s\n", $desc_dashes, "----------", "----------", "----------", "----------", "-----------------"), $sum_FH);

# For step 2.4 we need a de-versioned $efa_file_lost
my $efa_file_lost_acconly = $efa_file_lost . ".acconly"; # same as $efa_file_lost but with version info removed
if(! $do_nt) { 
  $cmd  = "cat $efa_file_lost | sed 's/\\.[0-9]*//' > $efa_file_lost_acconly";
  RunCommand($cmd, $be_verbose, $cmd_FH);
}

########################################################################
# Step 2.1: Validate all accessions exist in protein or nuccore database
########################################################################
my $exists_file         = $out_dir . $symbol . ".exists";
my $exists_file_lost    = $exists_file . ".lost";
my $exists_file_created = $exists_file . ".created";
my $db             = ($do_nt) ? "nuccore" : "protein";
my $protein_or_dna = ($do_nt) ? "DNA" : "protein";
$desc = sprintf("Validating_that_all_accessions_have_a_corresponding_%s_record", $protein_or_dna);
$cmd  = "cat $acc_file | epost -db $db -format acc | efetch -format acc | sed 's/\\.[0-9]*//' | sort > $exists_file";
$nsecs += RunCommand($cmd, $be_verbose, $cmd_FH);
OutputFileInfo($exists_file, $desc, $cmd, $log_FH);

# $exists_file format, multiple lines of:
# <accession>

# which accessions don't have a protein/nucleotide record? (should be none)
$cmd = "comm -2 -3 $acconly_file $exists_file > $exists_file_lost";
RunCommand($cmd, $be_verbose, $cmd_FH);
OutputFileInfo($exists_file_lost, sprintf("Accessions without %ss.", $protein_or_dna), $cmd, $log_FH);

# which accessions are created by this step? (should be none)
$cmd = "comm -2 -3 $exists_file $acconly_file > $exists_file_created";
RunCommand($cmd, $be_verbose, $cmd_FH);
OutputFileInfo($exists_file_created, "Accessions unexpectedly created by edirect command.", $cmd, $log_FH);

# enforce nothing was lost or created
($nlost, $ncreated, $errmsg) = CheckLostAndCreated($exists_file_lost, 0, $exists_file_created, 0); 
# in above subroutine call: '0's say it's not okay if we lost or created some accessions

# output information on this step to stdout and sum file
PrintToStdoutAndFile(sprintf("%-*s  %10d  %10d  %10d  %10.1f  %s\n", $desc_w, $desc, GetNumLinesInFile($exists_file), $nlost, $ncreated, $nsecs, $exists_file), $sum_FH);
if($errmsg ne "") { die $errmsg; }


################################################
# Step 2.2: Validate all accessions have a locus
################################################
my $lexists_file         = $out_dir . $symbol . ".lexists";
my $lexists_file_lost    = $lexists_file . ".lost";
my $lexists_file_created = $lexists_file . ".created";

$desc = "Validating_that_all_accessions_have_a_corresponding_locus";
$cmd  = "cat $acc_file | epost -db $db -format acc | efetch -format gpc | xtract -insd INSDSeq_locus | grep . | sort | sed 's/\\.[0-9]*//' > $lexists_file";
$nsecs = RunCommand($cmd, $be_verbose, $cmd_FH);
OutputFileInfo($lexists_file, $desc, $cmd, $log_FH);

# $lexists_file format, multiple lines of:
# <accession> <locus>

# which accessions don't have a locus? (should be none)
$cmd   = "awk '{ print \$1 }' $lexists_file | comm -2 -3 $acconly_file - > $lexists_file_lost";
RunCommand($cmd, $be_verbose, $cmd_FH);
OutputFileInfo($lexists_file_lost, "Accessions without a locus.", $cmd, $log_FH);

# which accessions are created by this step? (should be none)
$cmd   = "awk '{ print \$1 }' $lexists_file | comm -2 -3 - $acconly_file > $lexists_file_created";
RunCommand($cmd, $be_verbose, $cmd_FH);
OutputFileInfo($lexists_file_created, "Accessions unexpectedly created by edirect command.", $cmd, $log_FH);

# enforce nothing was lost or created
($nlost, $ncreated, $errmsg) = CheckLostAndCreated($lexists_file_lost, 0, $lexists_file_created, 0); 
# in above subroutine call: '0's say it's not okay if we lost or created some accessions

# output information on this step to stdout and sum file
PrintToStdoutAndFile(sprintf("%-*s  %10d  %10d  %10d  %10.1f  %s\n", $desc_w, $desc, GetNumLinesInFile($lexists_file), $nlost, $ncreated, $nsecs, $lexists_file), $sum_FH);
if($errmsg ne "") { die $errmsg; }


#############################################################################
# Step 2.3: Validate we've fetched sequences for all proteins that have a CDS
#############################################################################
my $has_cds_file         = $out_dir . $symbol . ".has_cds";
my $has_cds_file_lost    = $has_cds_file . ".lost";
my $has_cds_file_created = $has_cds_file . ".created";
if(! $do_nt) { 
  $desc = "Validating_that_we_fetched_sequences_for_all_proteins_with_a_CDS";
  $cmd  = "cat $exists_file | epost -db protein -format acc | efetch -format gpc | xtract -pattern INSDSeq -ACCN INSDSeq_accession-version -group INSDSeq_feature-table -match INSDFeature_key:CDS -element \"&ACCN\" | sort | sed 's/\\.[0-9]*//' > $has_cds_file";
  $nsecs = RunCommand($cmd, $be_verbose, $cmd_FH);
  OutputFileInfo($has_cds_file, $desc, $cmd, $log_FH);

  # $has_cds_file, multiple lines of:
  # <accession>

  # which accessions did we fetch that don't seem to have CDS? (should be none)
  $cmd = "comm -2 -3 $fa_file_names $has_cds_file > $has_cds_file_lost";
  RunCommand($cmd, $be_verbose, $cmd_FH);
  OutputFileInfo($has_cds_file_lost, "accessions we fetched which don't seem to have a CDS.", $cmd, $log_FH);

  # which accessions have CDS that we didn't fetch? (should be none)
  $cmd = "comm -2 -3 $has_cds_file $fa_file_names > $has_cds_file_created";
  RunCommand($cmd, $be_verbose, $cmd_FH);
  OutputFileInfo($has_cds_file_created, "accessions with a CDS that we failed to fetch.", $cmd, $log_FH);

  # enforce nothing was lost or created
  ($nlost, $ncreated, $errmsg) = CheckLostAndCreated($has_cds_file_lost, 0, $has_cds_file_created, 0);
  # in above subroutine call: '0's say it's not okay if we lost or created some accessions

  # output information on this step to stdout and sum file
  PrintToStdoutAndFile(sprintf("%-*s  %10d  %10d  %10d  %10.1f  %s\n", $desc_w, $desc, GetNumLinesInFile($has_cds_file), $nlost, $ncreated, $nsecs, $has_cds_file), $sum_FH);
  if($errmsg ne "") { die $errmsg; }
} # end of 'if($do_nt)'

#########################################################################################
# Step 2.4: Validate we haven't fetched any sequences for proteins that do not have a CDS
#########################################################################################
my $no_cds_file         = $out_dir . $symbol . ".no_cds";
my $no_cds_file_lost    = $no_cds_file . ".lost";
my $no_cds_file_created = $no_cds_file . ".created";
if(! $do_nt) { 
  $desc = "Validating_that_no_seqs_were_fetched_for_any_proteins_without_a_CDS";
  # same command as step 2.3 with 'avoid' replacing 'match'
  $cmd  = "cat $exists_file | epost -db protein -format acc | efetch -format gpc | xtract -pattern INSDSeq -ACCN INSDSeq_accession-version -group INSDSeq_feature-table -avoid INSDFeature_key:CDS -element \"&ACCN\" | sort | sed 's/\\.[0-9]*//' > $no_cds_file";
  $nsecs = RunCommand($cmd, $be_verbose, $cmd_FH);
  OutputFileInfo($no_cds_file, $desc, $cmd, $log_FH);

  # $no_cds_file, multiple lines of:
  # <accession>

  # which accessions did we not fetch that seem to have CDS? (should be none)
  $cmd = "comm -2 -3 $efa_file_lost_acconly $no_cds_file > $no_cds_file_lost";
  RunCommand($cmd, $be_verbose, $cmd_FH);
  OutputFileInfo($no_cds_file_lost, "accessions we did not fetch that seem to have a CDS.", $cmd, $log_FH);

  # which accessions don't have CDSs that we did fetch? (should be none)
  $cmd = "comm -2 -3 $no_cds_file $efa_file_lost_acconly > $no_cds_file_created";
  RunCommand($cmd, $be_verbose, $cmd_FH);
  OutputFileInfo($no_cds_file_created, "accessions without CDS that we did fetch.", $cmd, $log_FH);

  # enforce nothing was lost or created
  ($nlost, $ncreated, $errmsg) = CheckLostAndCreated($no_cds_file_lost, 0, $no_cds_file_created, 0);
  # in above subroutine call: '0's say it's not okay if we lost or created some accessions

  # output information on this step to stdout and sum file
  PrintToStdoutAndFile(sprintf("%-*s  %10d  %10d  %10d  %10.1f  %s\n", $desc_w, $desc, GetNumLinesInFile($no_cds_file), $nlost, $ncreated, $nsecs, $no_cds_file), $sum_FH);
  if($errmsg ne "") { die $errmsg; }
} # end of 'if(! $do_nt)'


######################################################################
# Step 2.5: Validate we can explain all proteins that don't have a CDS
######################################################################
# Look at all accessions we weren't about to fetch a CDS from, all should be either:
# (A) start with 'WP_', these are 100% identical at AA level to some other protein and have no CDS
# (B) start with 'YP_' and are dead according to idstat (return an ERROR with idstat)
# (C) are UniProt accessions (TODO: confirm we can get a CDS from these)
my $no_cds_type_a_file         = $out_dir . $symbol . ".no_cds_type_a";           # type A accessions
my $no_cds_type_b_file         = $out_dir . $symbol . ".no_cds_type_b";           # type B accessions
my $no_cds_type_c_file         = $out_dir . $symbol . ".no_cds_type_c";           # type C accessions
my $no_cds_type_a_b_c_file     = $out_dir . $symbol . ".no_cds_type_a_b_c";       # sorted concatenation of type A, B, and C accessions
my $no_cds_type_file_lost      = $out_dir . $symbol . ".no_cds_type.lost";        # accessions that are not type A, B, or C
my $no_cds_type_file_created   = $out_dir . $symbol . ".no_cds_type.created";     # accessions somehow created in this step

if(! $do_nt) { 
  ($seconds, $microseconds) = gettimeofday(); # so we can time this step
  my $step_start_secs = ($seconds + ($microseconds / 1000000.));
  
  # determine type A accessions
  $cmd   = "grep ^WP\_ $no_cds_file | cat - > $no_cds_type_a_file"; # pipe through cat to avoid grep failure (ret val of 1) if no strings found
  RunCommand($cmd, $be_verbose, $cmd_FH);
  OutputFileInfo($no_cds_type_a_file, "accessions without CDS that begin with \"WP\_\" (type A)", $cmd, $log_FH);
  
  # determine type B accessions
  $cmd   = "grep ^YP\_ $no_cds_file | cat - > $no_cds_type_b_file"; # pipe through cat to avoid grep failure (ret val of 1) if no strings found
  RunCommand($cmd, $be_verbose, $cmd_FH);
  OutputFileInfo($no_cds_type_b_file, "accessions without CDS that begin with \"YP\_\" (type B)", $cmd, $log_FH);
  
  # determine type C accessions, note that these all begin with a letter followed by a number
  # Reference: http://www.uniprot.org/help/accession_numbers
  # on that page, the regex that matches all UniProt accessions: [OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9]([A-Z][A-Z0-9]{2}[0-9]){1,2}
  # here we break it down into 3 separate grep calls:
  # Note also that these have no versions
  # ^[OPQ][0-9][A-Z0-9]{3}[0-9]$
  # ^[A-NR-Z][0-9][A-Z][A-Z0-9]{2}[0-9]$
  # ^[A-NR-Z][0-9]([A-Z][A-Z0-9]{2}[0-9]){2}$
  $cmd   = "grep '\^[OPQ][0-9][A-Z0-9]\\{3\\}[0-9]\$' $no_cds_file | cat - > $no_cds_type_c_file"; # pipe through cat to avoid grep failure (ret val of 1) if no strings found
  RunCommand($cmd, $be_verbose, $cmd_FH);
  $cmd   = "grep '\^[A-NR-Z][0-9][A-Z][A-Z0-9]\\{2\\}[0-9]\$' $no_cds_file | cat - >> $no_cds_type_c_file"; # pipe through cat to avoid grep failure (ret val of 1) if no strings found
  RunCommand($cmd, $be_verbose, $cmd_FH);
  $cmd   = "grep '\^[A-NR-Z][0-9]([A-Z][A-Z0-9]\\{2\\}[0-9])\\{2\\}\$' $no_cds_file | cat - >> $no_cds_type_c_file"; # pipe through cat to avoid grep failure (ret val of 1) if no strings found
  RunCommand($cmd, $be_verbose, $cmd_FH);
  
  # type A, type B, and type C files should comprise all the accessions, we can't explain any that are missing

  # combine the type A, type B and type C files
  $cmd = "cat $no_cds_type_a_file $no_cds_type_b_file $no_cds_type_c_file | sort > $no_cds_type_a_b_c_file";
  RunCommand($cmd, $be_verbose, $cmd_FH);
  OutputFileInfo($no_cds_type_a_b_c_file, "accessions which are missing CDS and can be classified as type A, B or C.", $cmd, $log_FH);
  
  # which accessions are not type A, type B or type C?
  $cmd = "comm -2 -3 $no_cds_file $no_cds_type_a_b_c_file > $no_cds_type_file_lost";
  RunCommand($cmd, $be_verbose, $cmd_FH);
  OutputFileInfo($no_cds_type_file_lost, "accessions which are missing CDS for reasons we don't understand.", $cmd, $log_FH);
  
  # which accessions were somehow created during this step? (should be none)
  $cmd = "comm -2 -3 $no_cds_type_a_b_c_file $no_cds_file > $no_cds_type_file_created";
  RunCommand($cmd, $be_verbose, $cmd_FH);
  OutputFileInfo($no_cds_type_file_lost, "accessions which were unexpectedly created when investigating missing CDSs.", $cmd, $log_FH);

  # enforce nothing was lost or created
  ($nlost, $ncreated, $errmsg) = CheckLostAndCreated($no_cds_type_file_lost, 0, $no_cds_type_file_created, 0);
  # in above subroutine call: '0's say it's not okay if we lost or created some accessions

  # output information on this step to stdout and sum file
  ($seconds, $microseconds) = gettimeofday(); # so we can time this step
  my $step_end_secs = ($seconds + ($microseconds / 1000000.));

  PrintToStdoutAndFile(sprintf("%-*s  %10d  %10d  %10d  %10.1f  %s\n", $desc_w, "Validating_that_proteins_without_CDS_are_explainable_(type_A,_B_or_C)", GetNumLinesInFile($no_cds_type_a_b_c_file), $nlost, $ncreated, $step_end_secs - $step_start_secs, $no_cds_type_a_b_c_file), $sum_FH);
  if($errmsg ne "") { die $errmsg; }
}

##########
# Conclude
##########
Conclude($start_secs, $do_nt, $sum_FH, $log_FH, $cmd_FH);

exit 0;

#############
# SUBROUTINES
#############
#
# Subroutine: RunCommand()
# Args:       $cmd:            command to run, with a "system" command;
#             $be_verbose:     '1' to output command to stdout before we run it, '0' not to
#             $FH:             file handle to output commands to
#
# Returns:    amount of time the command took, in seconds
# Dies:       if $cmd fails

sub RunCommand {
  my $sub_name = "RunCommand()";
  my $nargs_exp = 3;
  if(scalar(@_) != $nargs_exp) { die "ERROR $sub_name entered with wrong number of input args"; }

  my ($cmd, $be_verbose, $FH) = @_;

  if($be_verbose) { 
    print ("Running cmd: $cmd\n"); 
  }
  if(defined $FH) { 
    print $FH $cmd . "\n";
  }

  my ($seconds, $microseconds) = gettimeofday();
  my $start_time = ($seconds + ($microseconds / 1000000.));
  system($cmd);
  ($seconds, $microseconds) = gettimeofday();
  my $stop_time = ($seconds + ($microseconds / 1000000.));

  if($? != 0) { die "ERROR command failed:\n$cmd\n"; }

  return ($stop_time - $start_time);
}

# Subroutine: GetNumLinesInFile()
# Purpose:    Counts the number of non-blank lines 
#             in file $file.
# Args:       $file: file to get number of lines of

sub GetNumLinesInFile {
  my $sub_name = "GetNumLinesInFile()";
  my $nargs_exp = 1;
  if(scalar(@_) != $nargs_exp) { die "ERROR $sub_name entered with wrong number of input args"; }

  my ($file) = (@_);

  if(! -e $file) { die "ERROR file $file does not exist in CheckNumLinesInFile.\n"; }

  my $nlines = `grep -c . $file`;
  chomp $nlines;

  return $nlines;
}

# Subroutine: CheckLostAndCreated()
# Purpose:    Determine number of lines in 'lost' and 'created' files and check
#             if it exceeds our limits. 
# Args:       $lost_file:    lost file, each line represents a lost accession/id
#             $maxlost:      max number of allowed lost lines, -1 for no max
#             $created_file: created file, each line represents a newly created accession/id
#             $maxcreated:   max number of allowed created lines, -1 for no max
# 
# Returns: Three values:
#          $nlost:    number of lines in $lostfile
#          $ncreated: number of lines in $createdfile
#          $errmsg:   "" if no error, else filled with reason caller should die
#                     will be != "" if we find > maximum number of lines in 
#                     $lost_file and $created_file
#     
# Dies:    if $lostfile or $createdfile do not exist
sub CheckLostAndCreated { 
  my $sub_name = "OutputLostAndCreated()";
  my $nargs_exp = 4;
  if(scalar(@_) != $nargs_exp) { die "ERROR $sub_name entered with wrong number of input args"; }

  my ($lost_file, $maxlost, $created_file, $maxcreated) = (@_);

  if(! -e $lost_file)    { die "ERROR in $sub_name, lost file $lost_file does not exist"; }
  if(! -e $created_file) { die "ERROR in $sub_name, created file $created_file does not exist"; }

  $nlost    = GetNumLinesInFile($lost_file);
  $ncreated = GetNumLinesInFile($created_file);
  
  my $errmsg = "";
  if(($maxlost    != -1) && ($nlost    > $maxlost))    { $errmsg = "FAILED: #lost ($nlost) > $maxlost in file $lost_file"; }
  if(($maxcreated != -1) && ($ncreated > $maxcreated)) { $errmsg = "FAILED: #created ($ncreated) > $maxcreated in file $created_file"; }

  return ($nlost, $ncreated, $errmsg);
}

# Subroutine: PrintToStdoutAndFile()
# Purpose:    Output a string to stdout and to a file handle.
# Args:       $string:    string to print to stdout and $FH
#             $FH:        file handle to print to, unless undef
#
# Returns: void
sub PrintToStdoutAndFile { 

  my $sub_name = "PrintToStdoutAndFile()";
  my $nargs_exp = 2;
  if(scalar(@_) != $nargs_exp) { die "ERROR $sub_name entered with wrong number of input args"; }

  my ($string, $FH) = (@_);

  print $string;
  if(defined $FH) { print $FH $string; }

  return;
}

# Subroutine: OutputFileInfo()
# Purpose:    Add output file name to an array and description
#             to a hash, and output it to stdout and FH (if defined)
# Args:       $name:    name of output file
#             $desc:    description of file
#             $cmd:     command used to generate file $name
#             $FH:      file handle to output to
# Returns: void
sub OutputFileInfo {
  my $sub_name  = "OutputFileInfo()";
  my $nargs_exp = 4;
  if(scalar(@_) != $nargs_exp) { die "ERROR $sub_name entered with wrong number of input args"; }

  my ($name, $desc, $cmd, $FH) = (@_);

  if(! defined $FH) { die "ERROR file handle not defined in $sub_name"; }

  print $FH ("--------------------------------------------------\n");
  print $FH ("filename:        $name\n");
  print $FH ("description:     $desc\n");
  if($cmd ne "") { 
    print  $FH ("command:         $cmd\n");
    printf $FH ("number-of-lines: %d\n", GetNumLinesInFile($name));
  }
  else { 
    print $FH ("command:         N/A\n");
  }
  print $FH ("--------------------------------------------------\n");

  return;
}

# Subroutine: LengthFromCoords()
# Purpose:    Determine the length of a CDS given a coordinate 
#             string $coords.
# Args:       $coords:  the coords string
# Returns:    the length (in nt) of the CDS described by $coords
sub LengthFromCoords {
  my $sub_name  = "LengthFromCoords()"; 
  my $nargs_exp = 1;
  if(scalar(@_) != $nargs_exp) { die "ERROR $sub_name entered with wrong number of input args"; }

  my ($coords) = (@_);
  #print("IN LengthFromCoords(): $coords\n");

  my $orig_coords = $coords; # only nec so we can output it if there's an error
  # example $coords values:
  # U77714.1:1..867
  # complement(U77714.1:1..867)
  # join(U43876.1:608..688,U43877.1:104..175)
  # complement(join(U43876.1:608..688,U43877.1:104..175))
  # 
  # the 'join(' and 'complement' are useless
  $coords =~ s/^complement\(//;
  $coords =~ s/^join\(//;
  $coords =~ s/\)$//;
  $coords =~ s/\)$//;
  my @exon_A = split(",", $coords);

  if(scalar(@exon_A) == 0) { # single exon CDS
    push(@exon_A, $coords);
  }
  my $len = 0;
  foreach my $exon (@exon_A) { 
    if($exon =~ /^\S+\.?\d*\:\<?(\d+)\.\.\>?(\d+)$/) { 
      my ($start, $end) = ($1, $2);
      #print("\t$start $end\n");
      if($end < $start) { 
        die "ERROR found end coordinate ($end) less than start coord ($start) in coords string $orig_coords";
      }
      $len += ($end-$start+1);
    }
    else { 
      die "ERROR unable to parse exon $exon in coords string $orig_coords when trying to get length in LengthFromCoords()";
    }
  }
  # printf("LengthFromCoords(): $orig_coords $len\n");
  return $len;
}

# Subroutine: removeGivenIdstat()

# Purpose: Given an output file from 'idstat' ($idstat_file) and a
#             file that lists all accessions that 'idstat' was run on
#             (*in order*) ($allacc_file) output a new file with
#             some accessions removed ($out_file). Remove all accessions
#             that meet any of the following criteria:
#             - have 'Dead' in the 'State' column.
#             - have any value except 'No' in the 'Wdrw' (withdrawn) column
#             - have 'Yes' or 'Perm' in the 'Supp' (suppressed) column
#
# Args:       $idstat_file: the idstat output file
#             $allacc_file: the accession file containing all accession+versions 
#                           that idstat was run on the accession-onlys of to produce
#                           idstat_file
#             $out_file:    file to output non-suppressed accession+versions to
# Returns:    number of seconds required to do this parsing
# Dies:       if we have a problem parsing $idstat_file
sub removeGivenIdstat {
  my $sub_name  = "removeSuppressedGivenIdstat()";
  my $nargs_exp = 3;
  if(scalar(@_) != $nargs_exp) { die "ERROR $sub_name entered with wrong number of input args"; }
  
  my ($allacc_file, $idstat_file, $out_file) = (@_);
  my ($seconds, $microseconds) = gettimeofday();
  my $start_secs = ($seconds + ($microseconds / 1000000.));

  my @accn_A  = (); # array of accessions, in order read from idstate fule
  my @state_A = (); # array of 'State' column values same order as accn_A
  my @wdrw_A  = (); # array of 'Wdrw' column values same order as accn_A
  my @supp_A  = (); # array of 'Supp' column values same order as accn_A

  # first store 'suppressed' value for all accessions
  open(IDSTAT, $idstat_file) || die "ERROR unable to open $idstat_file"; 
  while(my $line = <IDSTAT>) { 
    #example:
    # Accession        GI         SatKey    SatName   State Owner                 Wdrw Supp Date Loaded           Username
    # ---------------- ---------- --------- --------- ----- --------------------- ---- ---- --------------------- --------
    # AAH45158         28175331   122024391 NCBI      Alive NCBI-MGC              No   No   9/23/2014 3:10:49PM   kans

    # chew up everything until we hit a line that starts with '---'
    if($line =~ m/^\-\-\-/) { # next line is special, its the most recent entry for an accession
      my $spec_line = <IDSTAT>;
      # AAH45158         28175331   122024391 NCBI      Alive NCBI-MGC              No   No   9/23/2014 3:10:49PM   kans
      my @elA = split(/\s+/, $spec_line);
      my ($accn, $state, $withdrawn, $suppressed) = ($elA[0], $elA[4], $elA[6], $elA[7]);
      push(@accn_A,  $accn);
      push(@state_A, $state);
      push(@wdrw_A,  $withdrawn);
      push(@supp_A,  $suppressed);
    }
  }
  close(IDSTAT);

  open(ACC,       $allacc_file) || die "ERROR unable to open $allacc_file for reading";
  open(OUT, ">" . $out_file)    || die "ERROR unable to open $out_file for writing";
  my $ctr = 0;
  while(my $line = <ACC>) { 
    if($line =~ m/\w/) { # skip blank lines
      if((! defined $accn_A[$ctr]) || (! defined $supp_A[$ctr])) { 
        die "ERROR ran out of parsed idstat accessions prematurely!"; 
      }
      chomp $line;
      my $accn_version = $line;
      my $accn2match = $accn_A[$ctr];
      if($accn_version !~ m/$accn2match\.?\d*/) { 
        die "ERROR, problem parsing $idstat_file and $allacc_file in $sub_name, trying to match $accn2match to $accn_version"; 
      }

      my $do_remove = 0; # set to '1' below if removed
      # parse 'state' column
      if($state_A[$ctr]   eq "Dead") { $do_remove = 1; } # dead, remove it

      # parse 'wdrw' column
      if($wdrw_A[$ctr]    ne "No")   { $do_remove = 1; }

      # parse 'supp' column
      if($supp_A[$ctr]    eq "Yes")  { $do_remove = 1; } # suppressed, remove it
      elsif($supp_A[$ctr] eq "Perm") { $do_remove = 1; } # permanently suppressed, remove it
      elsif($supp_A[$ctr] eq "No")   { ; } # not suppressed
      elsif($supp_A[$ctr] eq "Eblk") { ; } # 'edit-blocked', not suppressed as far as I understand, so we do not remove it
      else                           { die "ERROR unable to parse Supp column for $accn2match, read $supp_A[$ctr]"; }

      if(! $do_remove) { 
        print OUT $accn_version . "\n";
      }

      $ctr++;
    }
  }
  close ACC;
  close OUT;

  ($seconds, $microseconds) = gettimeofday(); # so we can time this step
  my $end_secs = ($seconds + ($microseconds / 1000000.));

  return ($end_secs - $start_secs);
}
  
# Subroutine: Capitalize()
# Purpose:    Return the capitalized version of a string.
# Args:       $str: string to capitalize
# Returns:    $cap_str: $str capitalize
sub Capitalize {
  my $sub_name  = "Capitalize()";
  my $nargs_exp = 1;
  if(scalar(@_) != $nargs_exp) { die "ERROR $sub_name entered with wrong number of input args"; }

  my ($str) = (@_);

  my $ret_str = $str;
  $ret_str =~ tr/a-z/A-Z/;

  return $ret_str;
}

# Subroutine: checkFileForSymbol()
# Purpose:    Given a file that is one line per symbol, return '1'
#             if $symbol exists as one of the symbols in the file,
#             case insensitive!
# Args:       $in_file: the file to check in for $in_symbol
#             $in_symbol: string to look for
# Returns:    '1' if $in_symbol is found in $in_file, else '0'
sub checkFileForSymbol {
  my $sub_name  = "checkFileForSymbol()";
  my $nargs_exp = 2;
  if(scalar(@_) != $nargs_exp) { die "ERROR $sub_name entered with wrong number of input args"; }

  my ($in_file, $in_symbol) = (@_);

  my $cap_in_symbol = Capitalize($in_symbol);

  open(IN, $in_file) || die "ERROR unable to open file $in_file in subroutine $sub_name"; 
  while(my $line = <IN>) { 
    chomp $line;
    my $symbol = $line;
    $symbol =~ s/^\s+//; # remove leading whitespace
    $symbol =~ s/\s+$//; # remove trailing whitespace
    $symbol =~ s/\t.*$//; # remove anything after first tab
    if(Capitalize($symbol) eq $cap_in_symbol) { 
      close(IN);
      return 1; # found a match
    }
  }
  close(IN);
  
  return 0; # didn't find a match
}

# Subroutine: stripVersion()
# Purpose:    Given a ref to an accession.version string, remove the version.
# Args:       $accver_R: ref to accession version string
# Returns:    Nothing, $$accver_R has version removed
sub stripVersion {
  my $sub_name  = "stripVersion()";
  my $nargs_exp = 1;
  if(scalar(@_) != $nargs_exp) { die "ERROR $sub_name entered with wrong number of input args"; }
  
  my ($accver_R) = (@_);

  $$accver_R =~ s/\.[0-9]*$//; # strip version

  return;
}

# Subroutine: updateTaxAndGeneinfo()
# Purpose:    Given a ref to an array of accessions, use edirect to output Gene db information
#             linked to those accessions, as well as taxid information.
# Args:       $accn_AR:       ref to array of accessions
#             $do_allgene:    '1' if -allgene flag used, we'll output all linked Gene info
#                             '0' to only output info for Genes annotated in accessions in $accnAR
#             $tmp_out_root: name for temporary output files we'll delete before leaving
#             $taxinfo_FH:   file handle to output taxid info to
#             $geneinfo_FH:  file handle to output gene info to 
#             $cmd_FH:       file handle to output command info to
# Returns:    2 values:
#             1) geneinfo command performed in this subroutine
#             3) number of seconds spent in this subroutine including time required for 
#                edirect commands
sub updateTaxAndGeneinfo {
  my $sub_name  = "updateTaxAndGeneinfo()";
  my $nargs_exp = 6;
  if(scalar(@_) != $nargs_exp) { die "ERROR $sub_name entered with wrong number of input args"; }
  
  my ($accn_AR, $do_allgene, $tmp_out_root, $taxinfo_FH, $geneinfo_FH, $cmd_FH) = (@_);
  my ($seconds, $microseconds) = gettimeofday();
  my $start_secs = ($seconds + ($microseconds / 1000000.));

  my %accn_H = (); # key accession, used to lookup accessions
  my %tax_H  = (); # key accession, value taxid for that accession
  my $query  = "\""; # all accessions concatenated together
  my $tmp_tax_info_file = $tmp_out_root . ".taxinfo";
  my $tmp_gene_info_file = $tmp_out_root . ".geneinfo";

  my $naccn = scalar(@{$accn_AR});
  for(my $i = 0; $i < $naccn; $i++) { 
    my $accn = $accn_AR->[$i];
    stripVersion(\$accn);
    $tax_H{$accn} = "";
    $accn_H{$accn} = 1;
    if($i > 0) { $query .= " OR "; }
    $query .= "$accn [accession]";
  }
  $query .= "\"";

  # first, determine the tax ids of each of the accessions
  my $tax_cmd = "esearch -query $query -db nuccore | efetch -format gbc | xtract -insd source db_xref | grep . > $tmp_tax_info_file";
  RunCommand($tax_cmd, $be_verbose, $cmd_FH);
  # we need to parse this file to figure out what the allowable tax ids are for the geneinfo file
  open(IN,  $tmp_tax_info_file)    || die "ERROR unable to open $tmp_tax_info_file for reading"; 
  while(my $taxline = <IN>) { 
    chomp $taxline;
    #printf("taxinfo line: $taxline");
    # A07116.1	taxon:11676
    # OR
    # NC_002501.1	ATCC:VR-896	taxon:114102
    my @elA = split(/\t/, $taxline);
    my $found_taxon = 0;
    my $taxid = "";
    my $tax_accn = $elA[0];
    stripVersion(\$tax_accn);
    for(my $i = 1; $i < scalar(@elA); $i++) { 
      if($elA[$i] =~ /^taxon\:(\d+)/) { 
          $taxid = $1;
      }
    }
    if($taxid eq "") { 
      die "ERROR unable to parse line when trying to find taxid info: $taxline\n";
    }
    if(! exists $tax_H{$tax_accn}) { die "ERROR found unexpected accession $tax_accn returned by edirect command: $tax_cmd"; }
    if($tax_H{$tax_accn} eq "") { 
      # we only want to print this line if it's not a duplicate; some records have taxid's stored twice 
      # and the esearch | efetch | xtract command above will return 2 identical lines.
      print $taxinfo_FH $taxline . "\n"; # output this line to our master file of all taxinfo
    }
    $tax_H{$tax_accn} = $taxid;
    } # end of 'while($taxline = <IN>)'
  close(IN);
  
  # ensure we have tax ids for all accessions
  foreach my $tax_accn (keys %tax_H) { 
    if($tax_H{$tax_accn} eq "") { 
      die "ERROR did not find taxid information for accession $tax_accn using edirect command: $tax_cmd"; 
    }
  }

  # now get the gene information
  my $gene_cmd = "esearch -query $query -db gene | efetch -format docsum | xtract -pattern DocumentSummary -group DocumentSummary -element Status -group Organism -element TaxID -group DocumentSummary -element Id -element Name -block GenomicInfoType -element ChrAccVer -element ChrStart -element ChrStop | grep -P \"\^0\\t\" | sort -k 2 > $tmp_gene_info_file"; 
  # the 'grep -P "^0\\t"' command toward the end specifies that only those records with status==0 (live gene records) will be kept
  # the final 'sort' makes it so the command doesn't return an error status if nothing survives past the grep

  $nsecs += RunCommand($gene_cmd, $be_verbose, $cmd_FH);
  print $geneinfo_FH "# Running cmd: $gene_cmd\n";
  # we need to parse this file to keep only allowable tax ids
  open(IN,  $tmp_gene_info_file)    || die "ERROR unable to open $tmp_gene_info_file for reading"; 
  while(my $geneline = <IN>) { 
    if($geneline =~ s/^0\t//) { 
      chomp $geneline;
      my $output_this_gene = 0;
      if($do_allgene) { 
        $output_this_gene = 1;
      }
      else { # ! $do_allgene
        # determine the accession for this gene
        my @elA = split(/\t/, $geneline);
        ##taxid gene-id gene-name gene-aliases accession chrStart chrEnd (sep character is a tab, lines beginning with '#' are comments)
        #10506	10971099	A329bR	NC_000852.5	166039	166224
        if(scalar(@elA) < 6) { die "ERROR parsing geneinfo line $geneline, at least 6 tab-delimited tokens expected: $geneline"; }
        my $cur_accn = $elA[3];
        stripVersion(\$cur_accn);
        if(exists $accn_H{$cur_accn}) { 
          $output_this_gene = 1;
        }
      }
      if($output_this_gene) { 
        print $geneinfo_FH $geneline . "\n";
      }
    } 
    else { 
      die "ERROR geneinfo line did not begin with '0': $geneline\n";
    }
  }

  ($seconds, $microseconds) = gettimeofday();
  my $end_secs = ($seconds + ($microseconds / 1000000.));

  if(-e $tmp_tax_info_file)  { unlink $tmp_tax_info_file; }
  if(-e $tmp_gene_info_file) { unlink $tmp_gene_info_file; }

  return($tax_cmd, $gene_cmd, ($end_secs - $start_secs));
}

# Subroutine: Conclude()
# Purpose:    Print out conclusion text and close file handles in preparation for exit.
# Args:       $start_secs: number of seconds since epoch and start of this script running
#             $do_nt:      '1' if we're in nucleotide mode, '0' if in protein mode, -1 if no records fetched
#             $sum_FH:     file handle for summary file
#             $log_FH:     file handle for log file
#             $cmd_FH:     file handle for command file
# Returns:    void
sub Conclude {
  my $sub_name  = "Conclude()";
  my $nargs_exp = 5;
  if(scalar(@_) != $nargs_exp) { die "ERROR $sub_name entered with wrong number of input args"; }
  
  my ($start_seconds, $do_nt, $sum_FH, $log_FH) = (@_);
  ($seconds, $microseconds) = gettimeofday();
  my $end_secs = ($seconds + ($microseconds / 1000000.));

  my $protein_or_nt;
  if   ($do_nt == -1) { $protein_or_nt = "N/A"; }
  elsif($do_nt ==  1) { $protein_or_nt = "nucleotide"; }
  elsif($do_nt ==  0) { $protein_or_nt = "protein"; }

  PrintToStdoutAndFile("#\n", $sum_FH);
  PrintToStdoutAndFile("# Output files created by this script with brief descriptions listed in log file:  $log_file\n",    $sum_FH);
  PrintToStdoutAndFile("# This output printed to stdout written to summary file:                           $sum_file\n",    $sum_FH);
  PrintToStdoutAndFile("# All commands executed by this script listed in cmd file:                         $cmd_file\n",    $sum_FH);
  PrintToStdoutAndFile("# All output files created in directory:                                           \.\/$out_dir\n", $sum_FH);
  PrintToStdoutAndFile(sprintf("# Total seconds elapsed:                                                           %.1f\n", $end_secs-$start_secs), $sum_FH);
  PrintToStdoutAndFile("#\n",     $sum_FH);
  PrintToStdoutAndFile("# Type of records: $protein_or_nt\n", $sum_FH);
  PrintToStdoutAndFile("#[ok]\n", $sum_FH);
  
  close $log_FH;
  close $sum_FH;
  close $cmd_FH;
  return;
}


