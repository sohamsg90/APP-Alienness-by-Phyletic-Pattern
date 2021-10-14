#!/usr/bin/env perl
#FILE: APP.pl
#AUTH: Soham Sengupta (sohamsengupta@my.unt.edu;sohamsengupta90@gmail.com)
#DATE: August 15, 2021
#VERS: 1.0
##Link to manual: https://github.com/sohamsg90

use strict;
use Data::Dumper;
use Getopt::Long;
use List::MoreUtils qw(uniq);
### Get my parameters for program ###
my %param = (
    query                 => undef,
    taxonomyFile          => undef,
    fileType              => undef,
    output_file           => undef,
    threads               => 0,
    markerGeneEnrichment  => 0,
    generateGeneMap      => 0,
    verbose               => 1,
    verboseDetailed       => 0,
    expert                => 0,
    help                  => undef
);

Getopt::Long::Configure('bundling');#https://perldoc.perl.org/Getopt::Long
GetOptions (
    'q|query=s'                => \$param{query},
    't|taxonomyFile=s'         => \$param{taxonomyFile},
    'f|fileType=s'             => \$param{fileType},
    'o|output_file=s'          => \$param{output_file},
    'n|threads=i'              => \$param{threads},
    'm|markerGeneEnrichment=i' => \$param{markerGeneEnrichment},
    'g|generateGeneMap=i'      => \$param{generateGeneMap},
    'v|verbose=i'              => \$param{verbose},
    'd|verboseDetailed=i'      => \$param{verboseDetailed},
    'e|expert=i'               => \$param{expert},
    'h|help'                   => \$param{help}
);
### Validate input parameters ###
## print user manual ##
if (defined($param{help}))
  {
      print_usage();
      exit(0);
  }
## turn on detailed program output ##
my $verbose = 1; #Basic verbose always
my $verboseDetailed = 0;
if ($param{verboseDetailed} == 1)
  {
    $verboseDetailed = 1;
  }
## check required input parameters ##
if (!(defined($param{query}))
    || !(defined($param{fileType}))
    || !defined($param{taxonomyFile}))
  {
      print_usage();
      exit(1);
  }
## Input file processing or set to download ##
if ($param{fileType} eq "multifasta")
  {
    print "\nInput file type: MultiFASTA.\n\n" if ($verbose == 1);
    if (-e "$param{query}.modified")
      {
        print "\nMultiFASTA located\n";
      }
    else 
      {
        process_fasta($param{query});    
      }
  }
elsif ($param{fileType} eq "accession")
  {
    print "\nUsing NCBI e-utilities FAA file will be downloaded.\n\n" if ($verbose == 1);
    obtain_faa_eutilities($param{query});
  }
else
  {
    print "\nError in declaring input file type. See instructions below\n\n" if ($verbose == 1);
    print_usage();
    exit(2);
  }

## obtain total number of CPU cores ##
my $num_threads = $param{threads};
unless ($param{threads})
  {
    # $num_threads = $param{threads};
    if (-e "/proc/cpuinfo") {
        $num_threads = `grep -c ^processor /proc/cpuinfo`;
    } elsif ($^O eq "darwin") {
        $num_threads = `sysctl -n hw.ncpu`;
    } elsif ($^O eq "MSWin32") {
        $num_threads = `echo %NUMBER_OF_PROCESSORS%`;
    }
    $num_threads =~ s/\s+$//;
    $num_threads = 0 unless $num_threads =~ /^\d+$/;
    unless ($num_threads) {
        print "Cannot determine the number of CPUs. Do single threading.\n";
        $num_threads = 1;
    }
  }
print "\nTotal number of CPU cores detected: $num_threads\n" if ($verbose == 1);

## check if marker gene enrichment is enabled
my $FLAG_markerGeneEnrichment = 0;#Marker gene enrichment is by user choice. By default set to 0
if ($param{markerGeneEnrichment} == 1)
  {
    $FLAG_markerGeneEnrichment = 1;
  }
## check if gene map is enabled
my $FLAG_generateGeneMap = 0;
if ($param{generateGeneMaps} == 1)
  {
    $FLAG_generateGeneMap = 1;
  }
##### Main program #####

### Default input files ###
my $f1 = "1.REFSEQ_all_kingdoms_ftplinks_lineage.txt";#taxidlist_bacteria_complete_full_lineage.txt
open IN, $f1 or die;
print "\nTaxonomic database imported....\n" if ($verbose == 1);
my @refseq_dataset = <IN>;

my $f2 = "$param{query}.modified";#genes test_gene.fasta. Multi-fasta sequences
open IN2, $f2 or die;
print "\nQuery sequence file imported....\n" if ($verbose == 1);
my @genes = <IN2>;
my @f = split(/\./,$f2);
my $f4;#Filename for all output files
if (defined($param{output_file}))
  {
    $f4 = $param{output_file};
  }
else 
  {
    $f4 = $f[0];
  }
chomp $f4;
print "\nOutput filename prefix: $f4\n" if ($verbose == 1);


my $f3 = $param{taxonomyFile};#all genomes species ID and taxon ID extacted from NCBI taxdump using taxonkit
open IN3, $f3 or die;
print "\nQuery taxonomic lineage imported....\n" if ($verbose == 1);
my @gID_speciesID_taxID = <IN3>;

my $taxonomic_level_to_blast = "all";
chomp $taxonomic_level_to_blast;

###Following are log-files for tracking various informations
open ERROR1, ">Error1.txt" or die;
open ERROR2, ">Error2.txt" or die;
open ERROR3, ">Error3.txt" or die;
open ERROR4, ">Error4.txt" or die;
open ERROR5, ">Error5_FTP_links_not_included.txt";

###Count how many genomes used at each level
open ERROR6, ">$f4\_Mega_Database_genome_counter.txt";
print ERROR6 "Gene\tstrain\tgenus\tfamily\torder\n";
###Determine which level to blast to
my $FLAG = 0;
if ($taxonomic_level_to_blast eq "all")
	{
		my $s = "species";
		my $g = "genus";
		my $f = "family";
    my $type = 1; #Header files for type 1 inclusion
    if ($type == 1)
      {
        my $accession_file1 = "$s\_$type\_Accession_no_species.txt";
        my $accession_file2 = "$g\_$type\_Accession_no_species.txt";
        my $accession_file3 = "$f\_$type\_Accession_no_species.txt";
        if (-e $accession_file1) {print "$s\_$type\_Accession_no_species.txt exists\n";$FLAG = 0;} else {open OUTPUT5, ">$s\_$type\_Accession_no_species.txt" or die;	$FLAG = 1;print "\nSpecies tracking file created....\n" if ($verboseDetailed == 1);}
        if (-e $accession_file2) {print "$g\_$type\_Accession_no_species.txt exists\n";$FLAG = 0;} else {open OUTPUT6, ">$g\_$type\_Accession_no_species.txt" or die;	$FLAG = 1;print "\nGenus tracking file created....\n" if ($verboseDetailed == 1);}
        if (-e $accession_file3) {print "$f\_$type\_Accession_no_species.txt exists\n";$FLAG = 0;}	else {open OUTPUT7, ">$f\_$type\_Accession_no_species.txt" or die;	$FLAG = 1;print "\nFamily tracking file created....\n" if ($verboseDetailed == 1);}
      }
    $type = 2; #Header files for type 2 exclusion
    if ($type == 2)
      {
        my $accession_file1 = "$s\_$type\_Accession_no_species.txt";
        my $accession_file2 = "$g\_$type\_Accession_no_species.txt";
        my $accession_file3 = "$f\_$type\_Accession_no_species.txt";
        if (-e $accession_file1) {print "$s\_$type\_Accession_no_species.txt exists\n";$FLAG = 0;} else {open OUTPUT8, ">$s\_$type\_Accession_no_species.txt" or die;	$FLAG = 1;print "\nSpecies tracking file created....\n" if ($verboseDetailed == 1);}
        if (-e $accession_file2) {print "$g\_$type\_Accession_no_species.txt exists\n";$FLAG = 0;} else {open OUTPUT9, ">$g\_$type\_Accession_no_species.txt" or die;	$FLAG = 1;print "\nGenus tracking file created....\n" if ($verboseDetailed == 1);}
        if (-e $accession_file3) {print "$f\_$type\_Accession_no_species.txt exists\n";$FLAG = 0;}	else {open OUTPUT10, ">$f\_$type\_Accession_no_species.txt" or die;	$FLAG = 1;print "\nFamily tracking file created....\n" if ($verboseDetailed == 1);}
      }
		
	}
elsif ($taxonomic_level_to_blast eq "species")
	{
		my $s = "species";
		my $accession_file1 = "$s\_Accession_no_species.txt";
		if (-e $accession_file1) {print "$s\_Accession_no_species.txt exists\n";$FLAG = 0;}	else {open OUTPUT5, ">$s\_Accession_no_species.txt" or die;	$FLAG = 1;}
	}
elsif ($taxonomic_level_to_blast eq "genus")
	{
		my $g = "genus";
		my $accession_file2 = "$g\_Accession_no_species.txt";
		if (-e $accession_file2) {print "$g\_Accession_no_species.txt exists\n";$FLAG = 0;}	else {open OUTPUT6, ">$g\_Accession_no_species.txt" or die;	$FLAG = 1;}
	}
elsif ($taxonomic_level_to_blast eq "family")
	{
		my $f = "family";
		my $accession_file3 = "$f\_Accession_no_species.txt";
		if (-e $accession_file3) {print "$f\_Accession_no_species.txt exists\n";$FLAG = 0;}	else {open OUTPUT7, ">$f\_Accession_no_species.txt" or die;	$FLAG = 1;}
	}

###Generate a Hash table of all species/tax IDs and their corresponding ftp-links
my %speciesID_taxID_ftplink;#All from NCBI
foreach my $line2 (@refseq_dataset)
	{
		chomp $line2;
		my @arr = split("\t", $line2);
		my $sID = $arr[12];
		my $tID = $arr[13];
		my $g_name = $arr[14];
		my $ftp_link = $arr[15];
		chomp($sID, $tID, $g_name, $ftp_link);
		if (defined($speciesID_taxID_ftplink{$sID}{$tID}{$g_name}))
			{
				print ERROR5 $sID,"\t",$tID,"\t",$g_name,"\t",$ftp_link,"\n";
			}
		else
			{
				$speciesID_taxID_ftplink{$sID}{$tID}{$g_name} = $ftp_link;
			}
	}

###Generate a Hash table of all species ID and taxonomic ID from file 3
my %genomeID_speciesID_taxID;
foreach my $line3 (@gID_speciesID_taxID)
	{
		chomp $line3;
		my @arr1 = split("\t", $line3);
		# print $arr1[1];
		my $a = $arr1[0];
		my $b = $arr1[1];
		my $c = $arr1[2];
		chomp ($a, $b, $c);
		my $str1 = $b.":".$c;
		$genomeID_speciesID_taxID{$a} = $str1;
	}

###Database of full NCBI complete genomes taxonomy to ftp-link. Here only taxonomic IDs are used.
my %global_taxonomyID_database;
my %global_species_to_full_lineage_database;#Databse of full NCBI complete genomes taxonomy with speciesID as Keys. Each record will have the full NCBI taxonomic classification
foreach my $line (@refseq_dataset)
	{
		chomp $line;
		my @arr = split("\t", $line);
		my $kingdom_rank = $arr[0];		my $phylum_rank = $arr[1];
		my $class_rank = $arr[2];		my $order_rank = $arr[3];
		my $family_rank = $arr[4];		my $genus_rank = $arr[5];
		my $kingdom_id = $arr[6];		my $phylum_id = $arr[7];
		my $class_id = $arr[8];		my $order_id = $arr[9];
		my $family_id = $arr[10];		my $genus_id = $arr[11];
		my $SPECIES_ID = $arr[12];		my $TAX_ID = $arr[13];
		my $Species_name = $arr[14];		my $FTP_link = $arr[15];
		my @ftp_links = $speciesID_taxID_ftplink{$SPECIES_ID};
		###For a given species all associated ftp-links are stored as values of array.
		$global_taxonomyID_database{$kingdom_id}{$phylum_id}{$class_id}{$order_id}{$family_id}{$genus_id}{$SPECIES_ID} = [@ftp_links];
    ###Traceback hash table generated for taxonomy extraction using just species ID
		$global_species_to_full_lineage_database{$SPECIES_ID}{kingdom} = $kingdom_id;
		$global_species_to_full_lineage_database{$SPECIES_ID}{phylum} = $phylum_id;
		$global_species_to_full_lineage_database{$SPECIES_ID}{class} = $class_id;
		$global_species_to_full_lineage_database{$SPECIES_ID}{order} = $order_id;
		$global_species_to_full_lineage_database{$SPECIES_ID}{family} = $family_id;
		$global_species_to_full_lineage_database{$SPECIES_ID}{genus} = $genus_id;
		$global_species_to_full_lineage_database{$SPECIES_ID}{species} = $Species_name;

	}
print "\nTaxonomic database logged....\n" if ($verboseDetailed == 1);
print ERROR1 Dumper \%global_taxonomyID_database if ($verbose==1);

##Extract speciesID and taxID for the genome
my ($species_ID, $tax_ID);
my $s_t_IDs = $genomeID_speciesID_taxID{$f4};
my @arr2 = split(":", $s_t_IDs);
$species_ID = $arr2[0];
$tax_ID = $arr2[1];
print "\nExtracted speciesID and taxID for the genome....\n" if ($verboseDetailed == 1);

#Subroutine selection for taxonomic level of choice
my $FLAG_pdt = 0;
if ($taxonomic_level_to_blast eq "species")	{my $level = "species"; my @ftp_links = species_level_blast($species_ID, $tax_ID, $level);}
elsif ($taxonomic_level_to_blast eq "genus"){my $level = "genus"; my @ftp_links = genus_level_blast($species_ID, $tax_ID, $level);}
elsif ($taxonomic_level_to_blast eq "family"){my $level = "family"; my @ftp_links = family_level_blast($species_ID, $tax_ID, $level);}
elsif ($taxonomic_level_to_blast eq "all"){my $level = "all"; my @ftp_links = all_level_blast($species_ID, $tax_ID, $level);	}

###Blasthits output processing post BLASTP analysis
# $FLAG_pdt = 1;
if ($taxonomic_level_to_blast eq "all")
	{
		my $type = 1;
    if ($type == 1)
      {
        my $input_blasthit_species = "$f4\_species_$type\_output_blasthits";
        my $input_blasthit_genus = "$f4\_genus_$type\_output_blasthits";
        my $input_blasthit_family = "$f4\_family_$type\_output_blasthits";
        # my $input_blasthit_order = "$f4\_order_output_blasthits";
        ###Check whether file exists and assign thresholds
        my $threshold_pident_species = 60; 
        my $threshold_pident_genus = 50; 
        my $threshold_pident_family = 25; 
        # my $threshold_pident_order = 25; my $cluster_separation_threshold_order = 15;#Values not checked
        if (-e $input_blasthit_species) {print "Species level hits exists\n";} else {print "Type 1 (inclusion) Species level hits doesn't exists or file not found\n";}
        if (-e $input_blasthit_genus) {print "Genus level hits exists\n";} else {print "Type 1 (inclusion) Genus level hits doesn't existsor file not found\n";}
        if (-e $input_blasthit_family) {print "Family level hits exists\n";} else {print "Type 1 (inclusion) Family level hits doesn't existsor file not found\n";}
        # if (-e $input_blasthit_order) {print "Order level hits exists\n";} else {print "Order level hits doesn't exists\n";}
        my $f6 = "species_$type\_Accession_no_species.txt";
        my $f7 = "genus_$type\_Accession_no_species.txt";
        my $f8 = "family_$type\_Accession_no_species.txt";
        # my $f9 = "order_Accession_no_species.txt";
        
        if ($FLAG_pdt == 1)
          {
            print "\nProcessing BLAST hits of species level of Type 1 (inclusion). Please wait....\n" if ($verbose == 1);
            output_blasthit_processor($type, $input_blasthit_species, $f6, $threshold_pident_species,  $f2, $f4);
          }
        print "\nProcessing BLAST hits of genus level of Type 1 (inclusion). Please wait....\n" if ($verbose == 1);
        output_blasthit_processor($type, $input_blasthit_genus, $f7, $threshold_pident_genus, $f2, $f4);
        print "\nProcessing BLAST hits of family level of Type 1 (inclusion). Please wait....\n" if ($verbose == 1);
        output_blasthit_processor($type, $input_blasthit_family, $f8, $threshold_pident_family, $f2, $f4);
        # output_blasthit_processor($input_blasthit_order, $f9, $threshold_pident_order, $cluster_separation_threshold_order, $f2);
      }
    $type = 2;
    if ($type == 2)
      {
        my $input_blasthit_species = "$f4\_species_$type\_output_blasthits";
        my $input_blasthit_genus = "$f4\_genus_$type\_output_blasthits";
        my $input_blasthit_family = "$f4\_family_$type\_output_blasthits";
        # my $input_blasthit_order = "$f4\_order_output_blasthits";
        ###Check whether file exists and assign thresholds
        my $threshold_pident_species = 60; 
        my $threshold_pident_genus = 50; 
        my $threshold_pident_family = 25; 
        # my $threshold_pident_order = 25; my $cluster_separation_threshold_order = 15;#Values not checked
        if (-e $input_blasthit_species) {print "Species level hits exists\n";} else {print "Type 2 (exclusion) Species level hits doesn't exists or file not found\n";}
        if (-e $input_blasthit_genus) {print "Genus level hits exists\n";} else {print "Type 2 (exclusion) Genus level hits doesn't existsor file not found\n";}
        if (-e $input_blasthit_family) {print "Family level hits exists\n";} else {print "Type 2 (exclusion) Family level hits doesn't existsor file not found\n";}
        # if (-e $input_blasthit_order) {print "Order level hits exists\n";} else {print "Order level hits doesn't exists\n";}
        my $f6 = "species_$type\_Accession_no_species.txt";
        my $f7 = "genus_$type\_Accession_no_species.txt";
        my $f8 = "family_$type\_Accession_no_species.txt";
        # my $f9 = "order_Accession_no_species.txt";
        if ($FLAG_pdt == 0)
          {
            print "\nProcessing BLAST hits of species level of Type 2 (exclusion). Please wait....\n" if ($verbose == 1);
            output_blasthit_processor($type, $input_blasthit_species, $f6, $threshold_pident_species, $f2, $f4);
          }
        print "\nProcessing BLAST hits of genus level of Type 2 (exclusion). Please wait....\n" if ($verbose == 1);
        output_blasthit_processor($type, $input_blasthit_genus, $f7, $threshold_pident_genus, $f2, $f4);
        print "\nProcessing BLAST hits of family level of Type 2 (exclusion). Please wait....\n" if ($verbose == 1);
        output_blasthit_processor($type, $input_blasthit_family, $f8, $threshold_pident_family, $f2, $f4);
        # output_blasthit_processor($input_blasthit_order, $f9, $threshold_pident_order, $cluster_separation_threshold_order, $f2);
      }

	}
elsif ($taxonomic_level_to_blast eq "species")
	{
		my $input_blasthit_species = "$f4\_species_output_blasthits";
		my $threshold_pident_species = 60; my $cluster_separation_threshold_species = 80;
		if (-e $input_blasthit_species) {print "Species level hits exists\n";} else {print "Species level hits doesn't exists or file not found\n";}
		my $f6 = "species_Accession_no_species.txt";
		output_blasthit_processor($input_blasthit_species, $f6, $threshold_pident_species, $cluster_separation_threshold_species, $f2, $f4);
	}
elsif ($taxonomic_level_to_blast eq "genus")
	{
		my $input_blasthit_genus = "$f4\_genus_output_blasthits";
		my $threshold_pident_genus = 50; my $cluster_separation_threshold_genus = 70;
		if (-e $input_blasthit_genus) {print "Genus level hits exists\n";} else {print "Genus level hits doesn't existsor file not found\n";}
		my $f7 = "genus_Accession_no_species.txt";
		output_blasthit_processor($input_blasthit_genus, $f7, $threshold_pident_genus, $cluster_separation_threshold_genus, $f2, $f4);
	}
elsif ($taxonomic_level_to_blast eq "family")
	{
		my $input_blasthit_family = "$f4\_family_output_blasthits";
		# my $input_blasthit_order = "$f4\_order_output_blasthits";
		my $threshold_pident_family = 25; my $cluster_separation_threshold_family = 40;
		# my $threshold_pident_order = 25; my $cluster_separation_threshold_order = 15;#Values not checked
		if (-e $input_blasthit_family) {print "Family level hits exists\n";} else {print "Family level hits doesn't existsor file not found\n";}
		# if (-e $input_blasthit_order) {print "Order level hits exists\n";} else {print "Order level hits doesn't exists\n";}
		my $f8 = "family_Accession_no_species.txt";
		# my $f9 = "order_Accession_no_species.txt";
		output_blasthit_processor($input_blasthit_family, $f8, $threshold_pident_family, $cluster_separation_threshold_family, $f2, $f4);
		# output_blasthit_processor($input_blasthit_order, $f9, $threshold_pident_order, $cluster_separation_threshold_order, $f2);
	}

###Alieness by phyletic pattern; Detect alien genes
##Import type 1 files with phyletic distribution values
my $f10 = "$f4\_60_calculations_species\_1.txt";
my $f11 = "$f4\_50_calculations_genus_1.txt";
my $f12 = "$f4\_25_calculations_family_1.txt";

##Import type 2 files with phyletic distribution values
my $f14 = "$f4\_60_calculations_species_2.txt";
my $f15 = "$f4\_50_calculations_genus_2.txt";
my $f16 = "$f4\_25_calculations_family_2.txt";


print "\nFinding alien genes. Please wait....\n" if ($verbose == 1);
alien_genes_finder();

###Perform maker gene enrichment if user wants
print "\nFinding marker genes. Please wait....\n" if ($FLAG_markerGeneEnrichment == 1);
if ($FLAG_markerGeneEnrichment == 1)
  {
    marker_gene_enrichment();
  }

###Generate gene maps
print "\nGenerating gene map using CGViewer. Please wait ...\n" if ($FLAG_generateGeneMap == 1);
if ($FLAG_generateGeneMap == 1)
  {
    get_genbank();
  }
##Program End; close FILEHANDLES ##
close(IN); close(IN2); close(IN3);
close(ERROR1);close(ERROR2);close(ERROR3);
close(ERROR4);close(ERROR5);close(ERROR6);

# Check expert option; Prevent from deleting temporary files ##
if ($param{expert} == 0)
  {
    system("rm Error*.txt");
    system("rm *.pdb *.phr *.pin *.pog *.pos *.pot *.psq *.ptf *.pto *.psi *.psd");
    system("rm FOR_INPUT*");
    system("rm *_raw_database.fasta");
    system("rm *_Accession_no_*");
    # system("rm ");
    # system("rm ");
    # system("rm ");
  }

#### Relevant subroutines ####
###All three levels of comparison performed. If no Family level database found, goes to Order level.
sub all_level_blast
	{
		my ($w, $x, $y) = @_;
		my $species_ID = $w;
		my $tax_ID = $x;
		my $level = $y;
		chomp ($species_ID, $tax_ID, $level);
		my @ftp_links;
    my $type; ##Determines inclusion and exclusion
 		my $kingdom = $global_species_to_full_lineage_database{$species_ID}{kingdom};
		my $phylum = $global_species_to_full_lineage_database{$species_ID}{phylum};
		my $class = $global_species_to_full_lineage_database{$species_ID}{class};
		my $order = $global_species_to_full_lineage_database{$species_ID}{order};
		my $family = $global_species_to_full_lineage_database{$species_ID}{family};
		my $genus = $global_species_to_full_lineage_database{$species_ID}{genus};
		my $order_under_consideration = $global_species_to_full_lineage_database{$species_ID}{order};
		my $family_under_consideration = $global_species_to_full_lineage_database{$species_ID}{family};
		my $genus_under_consideration = $global_species_to_full_lineage_database{$species_ID}{genus};
		my $species_under_consideration = $species_ID;
		my $taxa_under_consideration = $tax_ID;
    print "\nGenome download started....\n" if ($verbose == 1);
    if ($verbose==1)
      {
        print ERROR2 "kingdom under consideration: $kingdom\n"; print ERROR2 "phylum under consideration: $phylum\n";
        print ERROR2 "clase under consideration: $class\n"; print ERROR2 "order under consideration: $order_under_consideration\n";
        print ERROR2 "family under consideration: $family_under_consideration\n"; print ERROR2 "genus under consideration: $genus_under_consideration\n";
        print ERROR2 "species under consideration: $species_under_consideration\n"; print ERROR2 "taxID under consideration: $taxa_under_consideration\n";
      }
    ####Type1-inclusion####
    $type = 1;
		#####Species#####
    if ($type == 1)
      {
          print "\nAt species level (inclusion)\n" if ($verbose == 1);
          my %temp_array1 =  %{$global_taxonomyID_database{$kingdom}{$phylum}{$class}{$order}{$family}{$genus}{$species_under_consideration}[0]};
          my @keys1 = sort keys %temp_array1;
          my $MAIN_taxID_of_interest = join('',@keys1);
          print ERROR2 "MAIN_taxID_of_interest ",$MAIN_taxID_of_interest,"\n";
          my %database_of_genomes_to_blast_against1;
          foreach my $taxid (sort keys %{$global_taxonomyID_database{$kingdom}{$phylum}{$class}{$order}{$family}{$genus}{$species_under_consideration}[0]})
            {
                my %temp_array5 = %{$global_taxonomyID_database{$kingdom}{$phylum}{$class}{$order}{$family}{$genus}{$species_under_consideration}[0]{$taxid}};
                foreach my $genome (sort keys %temp_array5)
                  {
                    $database_of_genomes_to_blast_against1{$genome} = $temp_array5{$genome};
                  }
            }
          print ERROR3 Dumper \%database_of_genomes_to_blast_against1;
          my $size_of_strain = keys %database_of_genomes_to_blast_against1;
          print ERROR6 "$f4\t$size_of_strain";
          my $taxa_level_for_send1 = $species_under_consideration."_".$taxa_under_consideration;
          my $level2 = "species";
          if ($size_of_strain >= 1)
            {
              download_links_blastp_all($type, $taxa_level_for_send1, $level2, $size_of_strain, %database_of_genomes_to_blast_against1);
            }
          else {print "\nBlast not performed at species level (inclusion)\n" if ($verbose == 1);}
          

          #####Genus#######
          print "\nAt genus level (inclusion)\n" if ($verbose == 1);
          my %temp_array2 =  %{$global_taxonomyID_database{$kingdom}{$phylum}{$class}{$order}{$family}{$genus_under_consideration}};
          my %database_of_genomes_to_blast_against2;
          foreach my $species (sort keys %temp_array2)
            {
              chomp $species;
              my %temp_array3 = %{$temp_array2{$species}[0]};
              foreach my $taxID (sort keys %temp_array3)
                {
                  foreach my $bacteria (sort keys %{$temp_array3{$taxID}})
                    {
                      $database_of_genomes_to_blast_against2{$bacteria} = $temp_array3{$taxID}{$bacteria};
                    }
                }
            }
          print ERROR3 Dumper \%database_of_genomes_to_blast_against2;
          my $size_of_genus = keys %database_of_genomes_to_blast_against2;
          print ERROR6 "\t$size_of_genus";
          my $taxa_level_for_send2 = $genus_under_consideration."_".$species_under_consideration."_".$taxa_under_consideration;
          my $level3 = "genus";
          if ($size_of_genus >= 1)
            {
              download_links_blastp_all($type, $taxa_level_for_send2, $level3, $size_of_genus, %database_of_genomes_to_blast_against2);
            }
          else {print "\nBlast not performed at genus level (inclusion)\n" if ($verbose == 1);}

          #####Family######
          print "\nAt family level (inclusion)\n" if ($verbose == 1);
          my %temp_array3 =  %{$global_taxonomyID_database{$kingdom}{$phylum}{$class}{$order}{$family_under_consideration}};
          my %database_of_genomes_to_blast_against3;
          foreach my $genus (sort keys %temp_array3)
            {
              chomp $genus;
              foreach my $species (sort keys %{$temp_array3{$genus}})
                {
                  chomp $species;
                  my %temp_array4 = %{$temp_array3{$genus}{$species}[0]};
                  foreach my $taxID (sort keys %temp_array4)
                    {
                      foreach my $bacteria (sort keys %{$temp_array4{$taxID}})
                        {
                          $database_of_genomes_to_blast_against3{$bacteria} = $temp_array4{$taxID}{$bacteria};
                        }
                    }
                }
            }
          # print ERROR3 "Databse family\n";
          print ERROR3 Dumper \%database_of_genomes_to_blast_against3;
          my $size_of_family = keys %database_of_genomes_to_blast_against3;
          print ERROR6 "\t$size_of_family";
          if ($size_of_family > 0)
            {
              my $taxa_level_for_send3 = $family_under_consideration."_".$genus_under_consideration."_".$species_under_consideration."_".$taxa_under_consideration;
              my $level4 = "family";
              print ERROR6 "\t0";
              download_links_blastp_all($type, $taxa_level_for_send3, $level4, $size_of_family, %database_of_genomes_to_blast_against3);
            }
          else
            {
              ###Go to Order level only when no Family level database available
              print ERROR3 "No family database found;hence going to order level\n";
              print ERROR3 "\nAt order level\n";
              my %temp_array4 =  %{$global_taxonomyID_database{$kingdom}{$phylum}{$class}{$order_under_consideration}};
              print ERROR4 "\nAt order level\n";
              my %database_of_genomes_to_blast_against4;
              foreach my $family (sort keys %temp_array4)
                {
                  print ERROR4 $family,"\n";
                  foreach my $genus (sort keys %{$temp_array4{$family}})
                    {
                      print ERROR4 "\t$genus\n";
                      foreach my $species (sort keys %{$temp_array4{$family}{$genus}})
                        {
                          print ERROR4 "\t\t$species\n";
                          foreach my $taxid (sort keys %{$temp_array4{$family}{$genus}{$species}[0]})
                            {
                              my %temp_array6 = %{$temp_array4{$family}{$genus}{$species}[0]{$taxid}};
                              foreach my $genome (sort keys %temp_array6)
                                {
                                  print ERROR4 "\t\t\t$genome\n";
                                  $database_of_genomes_to_blast_against4{$genome} = $temp_array6{$genome};
                                }
                            }
                        }
                    }
                }
              print ERROR3 Dumper \%database_of_genomes_to_blast_against4;
              my $size_of_order = keys %database_of_genomes_to_blast_against4;
              print ERROR6 "\t$size_of_order";
              my $taxa_level_for_send4 = $order_under_consideration."_".$family_under_consideration."_".$genus_under_consideration."_".$species_under_consideration."_".$taxa_under_consideration;
              my $level5 = "order";
              print ERROR6 "\t0";
              download_links_blastp_all($type, $taxa_level_for_send4, $level5, $size_of_order, %database_of_genomes_to_blast_against4);
            }
            print ERROR6 "\n";
      }
    ####Type2-Exclusion####
    $type = 2; 
    if ($type == 2)
      {
        print "\nAt species level (exclusion)\n" if ($verbose == 1);
        my %temp_array1 =  %{$global_taxonomyID_database{$kingdom}{$phylum}{$class}{$order}{$family}{$genus}{$species_under_consideration}[0]};
        my @keys1 = sort keys %temp_array1;
        my $MAIN_taxID_of_interest = join('',@keys1);
        print ERROR2 "MAIN_taxID_of_interest ",$MAIN_taxID_of_interest,"\n";
        my %database_of_genomes_to_blast_against1;
        foreach my $taxid (sort keys %{$global_taxonomyID_database{$kingdom}{$phylum}{$class}{$order}{$family}{$genus}{$species_under_consideration}[0]})
          {
              next if $taxid == $taxa_under_consideration;
              my %temp_array5 = %{$global_taxonomyID_database{$kingdom}{$phylum}{$class}{$order}{$family}{$genus}{$species_under_consideration}[0]{$taxid}};
              foreach my $genome (sort keys %temp_array5)
                {
                  $database_of_genomes_to_blast_against1{$genome} = $temp_array5{$genome};
                }
          }
        print ERROR3 Dumper \%database_of_genomes_to_blast_against1;
        my $size_of_strain = keys %database_of_genomes_to_blast_against1;
        print ERROR6 "$f4\t$size_of_strain";
        if ($size_of_strain == 0) {$FLAG_pdt = 1;}
        my $taxa_level_for_send1 = $species_under_consideration."_".$taxa_under_consideration;
        my $level2 = "species";
        if ($size_of_strain >= 1)
          {
            download_links_blastp_all($type, $taxa_level_for_send1, $level2, $size_of_strain, %database_of_genomes_to_blast_against1);
          }
        else {print "\nBlast not performed at species level (exclusion)\n" if ($verbose == 1);}

        #####Genus#######
        print "\nAt genus level (exclusion)\n" if ($verbose == 1);
        my %temp_array2 =  %{$global_taxonomyID_database{$kingdom}{$phylum}{$class}{$order}{$family}{$genus_under_consideration}};
        my %database_of_genomes_to_blast_against2;
        foreach my $species (sort keys %temp_array2)
          {
            chomp $species;
            next if $species == $species_under_consideration;
            my %temp_array3 = %{$temp_array2{$species}[0]};
            foreach my $taxID (sort keys %temp_array3)
              {
                next if $taxID == $taxa_under_consideration;
                foreach my $bacteria (sort keys %{$temp_array3{$taxID}})
                  {
                    $database_of_genomes_to_blast_against2{$bacteria} = $temp_array3{$taxID}{$bacteria};
                  }
              }
          }
        print ERROR3 Dumper \%database_of_genomes_to_blast_against2;
        my $size_of_genus = keys %database_of_genomes_to_blast_against2;
        print ERROR6 "\t$size_of_genus";
        my $taxa_level_for_send2 = $genus_under_consideration."_".$species_under_consideration."_".$taxa_under_consideration;
        my $level3 = "genus";
        if ($size_of_genus >= 1)
          {
            download_links_blastp_all($type, $taxa_level_for_send2, $level3, $size_of_genus, %database_of_genomes_to_blast_against2);
          }
        else {print "\nBlast not performed at genus level (exclusion)\n" if ($verbose == 1);}

        #####Family######
        print "\nAt family level (exclusion)\n" if ($verbose == 1);
        my %temp_array3 =  %{$global_taxonomyID_database{$kingdom}{$phylum}{$class}{$order}{$family_under_consideration}};
        my %database_of_genomes_to_blast_against3;
        foreach my $genus (sort keys %temp_array3)
          {
            chomp $genus;
            next if $genus == $genus_under_consideration ;
            foreach my $species (sort keys %{$temp_array3{$genus}})
              {
                chomp $species;
                next if $species == $species_under_consideration;
                my %temp_array4 = %{$temp_array3{$genus}{$species}[0]};
                foreach my $taxID (sort keys %temp_array4)
                  {
                    next if $taxID == $taxa_under_consideration;
                    foreach my $bacteria (sort keys %{$temp_array4{$taxID}})
                      {
                        $database_of_genomes_to_blast_against3{$bacteria} = $temp_array4{$taxID}{$bacteria};
                      }
                  }
              }
          }
        # print ERROR3 "Databse family\n";
        print ERROR3 Dumper \%database_of_genomes_to_blast_against3;
        my $size_of_family = keys %database_of_genomes_to_blast_against3;
        print ERROR6 "\t$size_of_family";
        if ($size_of_family > 0)
          {
            my $taxa_level_for_send3 = $family_under_consideration."_".$genus_under_consideration."_".$species_under_consideration."_".$taxa_under_consideration;
            my $level4 = "family";
            print ERROR6 "\t0";
            download_links_blastp_all($type, $taxa_level_for_send3, $level4, $size_of_family, %database_of_genomes_to_blast_against3);
          }
        else
          {
            ###Go to Order level only when no Family level database available
            print "No family database found;hence going to order level\n";
            print "\nAt order level\n";
            my %temp_array4 =  %{$global_taxonomyID_database{$kingdom}{$phylum}{$class}{$order_under_consideration}};
            print ERROR4 "\nAt order level\n";
            my %database_of_genomes_to_blast_against4;
            foreach my $family (sort keys %temp_array4)
              {
                print ERROR4 $family,"\n";
                next if $family == $family_under_consideration;
                foreach my $genus (sort keys %{$temp_array4{$family}})
                  {
                    next if $genus == $genus_under_consideration;
                    print ERROR4 "\t$genus\n";
                    foreach my $species (sort keys %{$temp_array4{$family}{$genus}})
                      {
                        next if $species == $species_under_consideration;
                        print ERROR4 "\t\t$species\n";
                        foreach my $taxid (sort keys %{$temp_array4{$family}{$genus}{$species}[0]})
                          {
                            next if $taxid == $taxa_under_consideration;
                            my %temp_array6 = %{$temp_array4{$family}{$genus}{$species}[0]{$taxid}};
                            foreach my $genome (sort keys %temp_array6)
                              {
                                print ERROR4 "\t\t\t$genome\n";
                                $database_of_genomes_to_blast_against4{$genome} = $temp_array6{$genome};
                              }
                          }
                      }
                  }
              }
            print ERROR3 Dumper \%database_of_genomes_to_blast_against4;
            my $size_of_order = keys %database_of_genomes_to_blast_against4;
            print ERROR6 "\t$size_of_order";
            my $taxa_level_for_send4 = $order_under_consideration."_".$family_under_consideration."_".$genus_under_consideration."_".$species_under_consideration."_".$taxa_under_consideration;
            my $level5 = "order";
            print ERROR6 "\t0";
            download_links_blastp_all($type, $taxa_level_for_send4, $level5, $size_of_order, %database_of_genomes_to_blast_against4);
          }
          print ERROR6 "\n";
    }
	}

###Downloads NCBI genomes from FTP website based on previous selection. Common sub-routine irrespective of taxonomic level of interest. This subroutine downloads and renames files for convenience and generates a sequence header file for keeping track of accession numbers. This is primarily done to take into account identical proteins scenario.
sub download_links_blastp_all
  {
      my ($t, $x, $y, $w, %z) = @_;
      my $type = $t;
      my $taxa_under_consideration = $x;#Full lineage
      my $level = $y;#specific taxonomic level for other namings
      my $size = $w;
      my %database_links = %z;

      # print $size,"\n";
      print ERROR3 "inside download: $level";
      my $count = 1;
      my $target_database_file = "$level\_$type\_raw_database.fasta";
      # print $target_database_file,"\n";

      if (-e $target_database_file)
          {
              # system(`makeblastdb -in $level\_raw_database.fasta -parse_seqids -dbtype prot`);
              print "Database exists\n";
              system(`blastp -num_threads $num_threads -query $f2 -db $level\_$type\_raw_database.fasta -out $f4\_$level\_$type\_output_blasthits  -outfmt "6 qseqid sseqid pident qcovs length mismatch gapopen qstart qend sstart send evalue bitscore" -max_target_seqs $size`);

          }
      else
          {
              my $random_counter = 1;
              print "\nDownloading $level level genomes. Please wait....\n" if ($verbose == 1);
              print "\nTotal no. of genomes to download: $size \n" if ($verbose == 1);
              # my $count_d = 1;
              foreach my $genome (sort keys %database_links)
                  {
                      chomp $genome;
                      print "$random_counter of $size genomes downloaded\n" if ($verboseDetailed == 1);
                      my $link = $database_links{$genome};
                      chomp $link;
                      $genome =~ s/=| /_/g;#for blast database filename purposes
                      $genome =~ s/\.|\(|\)|\\|\/|;|,|\'|\"//g;#for blast database filename purposes
                      my $filename = $genome;
                      print ERROR3 $filename,"\n";
                      my @arr2 = split("/", $link);
                      my $file = $arr2[scalar(@arr2)-1]."_protein.faa.gz";
                      my $wget_link = $link."/".$file;
                      system (`wget -c -q $wget_link`);
                      if (-e $file)
                          {
                              system (`mv $file temp.faa.gz`);
                              system (`gzip -d temp.faa.gz`);
                              ###Following done to accomodate identical proteins in makeblastdb -parse_seqids scenario
                              if ($FLAG == 1)
                                  {
                                      open IN6, "temp.faa" or die;
                                      open TEMPOUT, ">$filename.faa" or die;
                                      while (<IN6>)
                                          {
                                              if (/^>/)
                                                  {
                                                      my @arr1 = split(">", $_);
                                                      my @arr2 = split(" ", $arr1[1]);
                                                      my $acc_no = $arr2[0];
                                                      my $new_header = $acc_no.":g".$random_counter;#for blast database filename purposes
                                                      print TEMPOUT ">$new_header\n";
                                                      ###For Accession_no_species files
                                                      if ($type == 1)
                                                        {
                                                          if ($level eq "species")	{print OUTPUT5 $new_header,"\t",$filename,"\n";}
                                                          elsif ($level eq "genus") {print OUTPUT6 $new_header,"\t",$filename,"\n";}
                                                          elsif ($level eq "family") {print OUTPUT7 $new_header,"\t",$filename,"\n";}
                                                        }
                                                      if ($type == 2)
                                                        {
                                                          if ($level eq "species")	{print OUTPUT8 $new_header,"\t",$filename,"\n";}
                                                          elsif ($level eq "genus") {print OUTPUT9 $new_header,"\t",$filename,"\n";}
                                                          elsif ($level eq "family") {print OUTPUT10 $new_header,"\t",$filename,"\n";}
                                                        }
                                                                                                        }
                                                  else
                                                  {
                                                      print TEMPOUT $_;
                                                  }
                                          }
                                          close (IN6);
                                          close (TEMPOUT);
                                          system (`rm temp.faa`);
                                  }
                          }
                      else
                          {
                              print ERROR3 "$file doesn't exists\n";
                          }
                    $random_counter++;
                  }
              system(`cat *.faa > $level\_$type\_raw_database.fasta`);
              # system(`ls *.faa | wc -l`);
              system(`rm *.faa`);
              print "\nConstructing local blast database. Please wait....\n" if ($verbose == 1);
              system(`makeblastdb -in $level\_$type\_raw_database.fasta -parse_seqids -dbtype prot`);
              print "\nLocal blast database generated. Performing BLAST. Please wait....\n" if ($verbose == 1);
              system(`blastp -num_threads $num_threads -query $f2 -db $level\_$type\_raw_database.fasta -out $f4\_$level\_$type\_output_blasthits -outfmt "6 qseqid sseqid pident qcovs length mismatch gapopen qstart qend sstart send evalue bitscore" -max_target_seqs $size`);
          }
      # system (`mkdir $f4\_main_output_blasthits_folder\_$level`);
      # system (`mv $f4\_$level\_output_blasthits $f4\_main_output_blasthits_folder\_$level`);


  }
###Process blasthits, based on given thresholds, count absence or presence of genes, generate phyletic pattern
sub output_blasthit_processor
	{
		my ($t, $a, $b, $c, $d, $e) = @_;
    my $type = $t;
		my $blasthits = $a;
		my $acc_no_list = $b;
		my $threshold_pident = $c;
		my $query_file = $d;#Input multifasta file
		my $f4 = $e;#output filename

		open IN7, $blasthits or die;
		my @main_blasthits_file = <IN7>;

		open IN8, $acc_no_list or die;
		my @accno_species = <IN8>;

		open IN9, $query_file or die;
		my @sequences = <IN9>;

		my $level = taxonomic_level($threshold_pident);

		open ERROR7, ">Error.txt" or die;
		open OUTPUT1, ">$f4\_Main_blast_output_counts_pident_$threshold_pident\_$level\_$type.txt" or die;
		open OUTPUT2, ">$f4\_$threshold_pident\_calculations\_$level\_$type.txt" or die;
		# open OUTPUT4, ">$f4\_$threshold_pident\_$cluster_separation_threshold\_gene_counts\_$level\_$type.txt" or die;
		print OUTPUT2 "Gene\tno. of genomes\tcount in which present\tfraction of genomes with the gene\n";
		open OUTPUT3, ">FOR_INPUT\_$level\_$type.txt" or die;

		#####Post processing files of blast _output_blasthits
		###Generate database of accession nos and their species with strain information
		my @names_all_genomes;
		my %accno_genome_name_database;
		foreach my $line1 (@accno_species)
			{
				chomp $line1;
				my @arr1 = split("\t", $line1);
				my $accno = $arr1[0];
				my $genome_name = $arr1[1];
				chomp ($accno, $genome_name);
				push (@names_all_genomes, $genome_name);
				$accno_genome_name_database{$accno} = $genome_name;
			}
		my @uniq_all_genomes = uniq @names_all_genomes;
		print ERROR7 $level,"\n";
		# print ERROR7 Dumper \%accno_genome_name_database;
		# print ERROR7 Dumper \@names_all_genomes;
		###Extract all JSCB/Query accession numbers
		my @JSCB_accno_all;
		foreach my $line2 (@sequences)
			{
				chomp $line2;
				if ($line2 =~ /^>/)
					{
						my @arr2 = split(">", $line2);
						push (@JSCB_accno_all, $arr2[1]);
					}
			}

		##Assign 0 counts to all gene vs genomes combos
		my %blasthits_counter;
		foreach my $gene (@JSCB_accno_all)
			{
				chomp $gene;
				foreach my $genome (@uniq_all_genomes)
					{
						chomp $genome;
						$blasthits_counter{$gene}{$genome} = 0;
					}
			}
		# print ERROR7 Dumper \%blasthits_counter;
		##Fill out counts as 1 for actual hits
		foreach my $hit (@main_blasthits_file)
			{
				chomp $hit;
				# $hit =~ s/gb|emb|dbj|ref//g; #remove leading gb|emb|dj
				# $hit =~ s/\|//g;#remove | symbol
				my @arr = split("\t", $hit);
				my $query = $arr[0];
				my $subject = $arr[1];
				my $pident = $arr[2];
				my $qcoverage = $arr[3];
				# my @arr2 = split(":", $file);
				my $target_species = $accno_genome_name_database{$subject};
				# if (defined! $target_species){print ERROR7 $hit,"\n";}
				if ($pident >= $threshold_pident)
					{
						if ($qcoverage >= 70)
							{
										$blasthits_counter{$query}{$target_species} = 1;
							}
					}
				# print ERROR7 $hit,"\n";
			}
		# print ERROR7 Dumper \%blasthits_counter;

		####All of next to convert multidimensional hash to table.txt
		my @final_blasthits_output;
		my @keys;
		foreach my $query (sort keys %blasthits_counter)
			{
			@keys = sort keys %{$blasthits_counter{$query}};
			}
		# print Dumper \@keys;
		push (@final_blasthits_output, "genes","\t");
		foreach my $out (@keys)
		{
			chomp $out;
				my $species = $out;
			push (@final_blasthits_output, $species,"\t");
		}
		push(@final_blasthits_output, "\n");
		# print ERROR7 Dumper \@final_blasthits_output;

		foreach my $query (sort keys %blasthits_counter)
			{
			push (@final_blasthits_output, $query,"\t");
			foreach my $filename (sort keys %{$blasthits_counter{$query}})
					{
						my $count = $blasthits_counter{$query}{$filename};
						# my @arr = split("_", $filename);
				# my $arr_size = scalar(@arr);
				# my $subject_taxID = $arr[$arr_size-3];
				push(@final_blasthits_output,$count,"\t");
					}
			push(@final_blasthits_output,"\n");
			}
		foreach my $l (@final_blasthits_output)
			{
				print OUTPUT1 $l;
				print OUTPUT3 $l;
			}
		close(OUTPUT1);
		close(OUTPUT3);
		my $f9 = "FOR_INPUT\_$level\_$type.txt";
		open IN10, $f9 or die;
		my @final_blasthits_output2 = <IN10>;
		# print OUTPUT2 Dumper \@final_blasthits_output2;
		#####Calculations of NATIVE or ALIEN genes###########
		my @calculations;
		# my $type_of_gene = 0;
		foreach my $l (@final_blasthits_output2)
			{
				chomp $l;
				next if $l =~ /^genes/;
				# print  $l,"\n";
				my @arr = split("\t", $l);
				my $gene = $arr[0];
				chomp ($gene);
				my $gene_counter = 0;
				for (my $i = 0; $i <= scalar(@arr)-2; $i++)
					{
						if ($arr[$i+1] == 1)
							{
								$gene_counter++;
							}
					}
				my $size = scalar(@arr)-1;
				my $fraction_of_genomes = ($gene_counter/$size)*100;
				# if ($fraction_of_genomes > $cluster_separation_threshold)
				# 	{
				# 		$type_of_gene = 1;
				# 	}
				# else
				# 	{
				# 		$type_of_gene = 0;
				# 	}
				my $str = $gene."\t".$size."\t".$gene_counter."\t".$fraction_of_genomes."\n";
				push (@calculations,$str);
			}
		# my $core_gene_counter = 0;
		# my $total_gene_counter = 0;
		# foreach my $l (@calculations)
		# 	{
		# 		chomp $l;
		# 		my @arr = split("\t", $l);
		# 		my $type_of_gene = $arr[4];
		# 		$core_gene_counter = $core_gene_counter + $type_of_gene;
		# 		$total_gene_counter++;
		# 	}
		# my $foreign_genes = $total_gene_counter - $core_gene_counter;
		# print OUTPUT4 "$total_gene_counter\t";
		# print OUTPUT4 "$core_gene_counter\t";
		# print OUTPUT4 "$foreign_genes\n";
		# print "No. of total genes = $total_gene_counter\n";
		# print "No. of native genes = $core_gene_counter\n";
		# print "No. of foreign genes = $foreign_genes\n";
		# print "If the cluster is native or of foreign origin: "
		foreach my $l (@calculations)
			{
				# chomp $l;
        print OUTPUT2 $l;
			}
    close(IN7); close(IN8); close(IN9); close(OUTPUT2);
	}

 sub taxonomic_level
	{
		my ($a) = @_;
		my $threshold_pident = $a;
		my $level;

		if ($threshold_pident == 60)		{$level = "species";}
		elsif ($threshold_pident == 50)	{$level = "genus";}
		elsif ($threshold_pident == 25)	{$level = "family";}
		
		return($level);
	}


sub alien_genes_finder
  {
    # my ($a, $b, $c, $d, $e, $f) = @_;
    # print Dumper \@_;

    #Calculation phyletic distribution values: Species, Type1
    my %species_level_type1;
    if ($FLAG_pdt == 1)
      {
        my $f1 = $f10; 
        open IN11, $f1 or die;
        my @species_type1 =<IN11>;
        foreach my $line2 (@species_type1)
        {
            chomp $line2;
            next if $line2 =~ /^Gene/;
            my @arr2 = split("\t", $line2);
            my $acc = $arr2[0];
            my $phyletic_distribution_threshold = $arr2[3];
            my $type = $arr2[4];
            chomp ($acc, $phyletic_distribution_threshold, $type);
            $species_level_type1{$acc}{"phyletic_distribution_threshold"} = $phyletic_distribution_threshold;
            # $species_level_type1{$acc}{"type"} = $type;

        }
      }
    
    my $f2 = $f11; #Calculation phyletic distribution values: Genus, Type1
    open IN12, $f2 or die;
    my @genus_type1 = <IN12>;
    my $f3 = $f12; #Calculation phyletic distribution values: Family, Type1
    open IN14, $f3 or die;
    my @family_type1 = <IN14>;

    #Calculation phyletic distribution values: Species, Type2
    my %species_level_type2;
    if ($FLAG_pdt == 0)
      {
        my $f4 = $f14;
        open IN15, $f4 or die;
        my @species_type2 =<IN15>;
        
        foreach my $line2 (@species_type2)
            {
                chomp $line2;
                next if $line2 =~ /^Gene/;
                my @arr2 = split("\t", $line2);
                my $acc = $arr2[0];
                my $phyletic_distribution_threshold = $arr2[3];
                my $type = $arr2[4];
                chomp ($acc, $phyletic_distribution_threshold, $type);
                $species_level_type2{$acc}{"phyletic_distribution_threshold"} = $phyletic_distribution_threshold;
                # $species_level_type2{$acc}{"type"} = $type;

            }
      }
    my $f5 = $f15; #Calculation phyletic distribution values: Genus, Type2
    open IN16, $f5 or die;
    my @genus_type2 = <IN16>;
    my $f6 = $f16; #Calculation phyletic distribution values: Family, Type2
    open IN17, $f6 or die;
    my @family_type2 = <IN17>;

    ###Extract all accession numbers in sequential order from faa file. Order checked and matched with gbk
    my %accession_nos;
    my $random_counter = 1;
    my $f7 = "$param{query}.modified";
    open IN18, $f7 or die;
    while(<IN18>)
      {
        if (/^>/)
          {
            my @arr1 = split(">", $_);
            my $g = $arr1[1];
            chomp $g;
            $accession_nos{$random_counter} = $g;
            $random_counter++;
          }
      }
    # open OUTPUT11, ">$f4\_phyletic_pattern_analysis_native_aliene_genes_$phyletic_distribution_threshold_native\_$phyletic_distribution_threshold_alien.txt";
    open OUTPUT11, ">$f4\_phyletic_pattern_analysis_native_alien_genes_full_record.txt";
    open OUTPUT12, ">$f4\_APP_Alien_genes.txt";
    
    ####Type1 data
    # print OUTPUT11 Dumper \%species_level;
    my %genus_level_type1;
    foreach my $line2 (@genus_type1)
        {
            chomp $line2;
            next if $line2 =~ /^Gene/;
            my @arr2 = split("\t", $line2);
            my $acc = $arr2[0];
            my $phyletic_distribution_threshold = $arr2[3];
            my $type = $arr2[4];
            chomp ($acc, $phyletic_distribution_threshold, $type);
            $genus_level_type1{$acc}{"phyletic_distribution_threshold"} = $phyletic_distribution_threshold;
            # $genus_level_type1{$acc}{"type"} = $type;

        }

    my %family_level_type1;
    foreach my $line2 (@family_type1)
        {
            chomp $line2;
            next if $line2 =~ /^Gene/;
            my @arr2 = split("\t", $line2);
            my $acc = $arr2[0];
            my $phyletic_distribution_threshold = $arr2[3];
            my $type = $arr2[4];
            chomp ($acc, $phyletic_distribution_threshold, $type);
            $family_level_type1{$acc}{"phyletic_distribution_threshold"} = $phyletic_distribution_threshold;
            # $family_level_type1{$acc}{"type"} = $type;

        }
    #Type2 data  
    
    # print OUTPUT11 Dumper \%species_level_type2;
    my %genus_level_type2;
    foreach my $line2 (@genus_type2)
        {
            chomp $line2;
            next if $line2 =~ /^Gene/;
            my @arr2 = split("\t", $line2);
            my $acc = $arr2[0];
            my $phyletic_distribution_threshold = $arr2[3];
            my $type = $arr2[4];
            chomp ($acc, $phyletic_distribution_threshold, $type);
            $genus_level_type2{$acc}{"phyletic_distribution_threshold"} = $phyletic_distribution_threshold;
            # $genus_level_type2{$acc}{"type"} = $type;

        }

    my %family_level_type2;
    foreach my $line2 (@family_type2)
        {
            chomp $line2;
            next if $line2 =~ /^Gene/;
            my @arr2 = split("\t", $line2);
            my $acc = $arr2[0];
            my $phyletic_distribution_threshold = $arr2[3];
            my $type = $arr2[4];
            chomp ($acc, $phyletic_distribution_threshold, $type);
            $family_level_type2{$acc}{"phyletic_distribution_threshold"} = $phyletic_distribution_threshold;
            # $family_level_type2{$acc}{"type"} = $type;

        }                    
    my %final_native_alien_genes;  
    my %species_level_data;
    if ($FLAG_pdt == 0) {%species_level_data = %species_level_type2;}
    elsif($FLAG_pdt == 1) {%species_level_data = %species_level_type1;}
    foreach my $acc (sort keys %species_level_data)
        {
            chomp $acc;

            #Type1 Inclusion data
            my $phyletic_distribution_threshold_species_type1 = $species_level_type1{$acc}{"phyletic_distribution_threshold"};
            my $phyletic_distribution_threshold_genus_type1 = $genus_level_type1{$acc}{"phyletic_distribution_threshold"};
            my $phyletic_distribution_threshold_family_type1 = $family_level_type1{$acc}{"phyletic_distribution_threshold"};

            #Type2 Exclusion data
            my $phyletic_distribution_threshold_species_type2 = $species_level_type2{$acc}{"phyletic_distribution_threshold"};
            my $phyletic_distribution_threshold_genus_type2 = $genus_level_type2{$acc}{"phyletic_distribution_threshold"};
            my $phyletic_distribution_threshold_family_type2 = $family_level_type2{$acc}{"phyletic_distribution_threshold"};
            # print "$phyletic_distribution_threshold_species_type1\t$phyletic_distribution_threshold_genus_type1\t$phyletic_distribution_threshold_family_type1\t
            # $phyletic_distribution_threshold_species_type2\t$phyletic_distribution_threshold_genus_type2\t$phyletic_distribution_threshold_family_type2\n";

            ##Decision rules; 
            #"type" => Native = 1; Alien = 0
            #"mode" => Ancient = 1; Recent = 0; Native = 2; Native values assigned to avoid error warning
            #"PDT_species_T1_or_T2" => When exlusion data is not available, type 1 data is used at species level
            my $PDT_species_T1_or_T2;
            if ($FLAG_pdt == 1)
              {
                print "\nNo exclusion data, using inclusion data for species level\n" if ($verbose == 1);
                $PDT_species_T1_or_T2 = $phyletic_distribution_threshold_species_type1;
              }
            elsif ($FLAG_pdt == 0)
              {
                print "\nUsing exclusion data for species level\n" if ($verbose == 1);
                $PDT_species_T1_or_T2 = $phyletic_distribution_threshold_species_type2;
              }
            if ($PDT_species_T1_or_T2 <= 30)
                {
                  $final_native_alien_genes{$acc}{"type"} = 0;
                  $final_native_alien_genes{$acc}{"mode"} = 0;#Recent
                }
            elsif ($PDT_species_T1_or_T2 >= 80)
                {
                    if ($phyletic_distribution_threshold_genus_type2 >= 70)
                        {
                            if ($phyletic_distribution_threshold_family_type2 >= 40)
                                {
                                    $final_native_alien_genes{$acc}{"type"} = 1; 
                                    $final_native_alien_genes{$acc}{"mode"} = 2;#Native hence ignore
                                }
                            elsif ($phyletic_distribution_threshold_family_type2 <= 30)
                                {
                                    $final_native_alien_genes{$acc}{"type"} = 0; 
                                    $final_native_alien_genes{$acc}{"mode"} = 1;#Ancient
                                }
                            elsif ($phyletic_distribution_threshold_family_type2 > 30 && $phyletic_distribution_threshold_family_type2 < 40)
                                {
                                    $final_native_alien_genes{$acc}{"type"} = "Ambiguous1";
                                    $final_native_alien_genes{$acc}{"mode"} = 2;#Native hence ignore
                                }
                            else 
                                {
                                    $final_native_alien_genes{$acc}{"type"} = "Error1";
                                }
                        }
                    elsif ($phyletic_distribution_threshold_genus_type2 <= 70)
                        {
                            if ($phyletic_distribution_threshold_genus_type2 >= 30)
                                {
                                    if ($phyletic_distribution_threshold_family_type1 <= 30)
                                        {
                                            $final_native_alien_genes{$acc}{"type"} = 0;
                                            $final_native_alien_genes{$acc}{"mode"} = 1;#Ancient
                                        }
                                    elsif ($phyletic_distribution_threshold_family_type1 >= 40 )
                                        {
                                            $final_native_alien_genes{$acc}{"type"} = 1;
                                            $final_native_alien_genes{$acc}{"mode"} = 2;#Native hence ignore
                                        }
                                    elsif ($phyletic_distribution_threshold_family_type1 > 30 && $phyletic_distribution_threshold_family_type1 < 40)
                                        {
                                            $final_native_alien_genes{$acc}{"type"} = "Ambiguous2";
                                            $final_native_alien_genes{$acc}{"mode"} = 2;#Native hence ignore
                                        } 
                                    else 
                                        {
                                            $final_native_alien_genes{$acc}{"type"} = "Error2";
                                        }
                                }
                            elsif ($phyletic_distribution_threshold_genus_type2 < 30)
                                {
                                    $final_native_alien_genes{$acc}{"type"} = 0;
                                    $final_native_alien_genes{$acc}{"mode"} = 1;#Ancient
                                }
                            else 
                                {
                                    $final_native_alien_genes{$acc}{"type"} = "Error3";
                                }
                        }
                }
            elsif ($PDT_species_T1_or_T2 < 80 && $PDT_species_T1_or_T2 > 30)
                {
                    if ($phyletic_distribution_threshold_genus_type1 <= 30)
                        {
                            $final_native_alien_genes{$acc}{"type"} = 0;
                            $final_native_alien_genes{$acc}{"mode"} = 1;#Ancient
                        }
                    elsif ($phyletic_distribution_threshold_genus_type1 >= 70)
                        {
                            $final_native_alien_genes{$acc}{"type"} = 1;
                            $final_native_alien_genes{$acc}{"mode"} = 2;#Native hence ignore
                        }
                    elsif ($phyletic_distribution_threshold_genus_type1 > 30 && $phyletic_distribution_threshold_genus_type1 < 70)
                        {
                            if ($phyletic_distribution_threshold_family_type1 <= 30)
                                {
                                    $final_native_alien_genes{$acc}{"type"} = 0;
                                    $final_native_alien_genes{$acc}{"mode"} = 1;#Ancient
                                }
                            elsif ($phyletic_distribution_threshold_family_type1 >= 40 )
                                {
                                    $final_native_alien_genes{$acc}{"type"} = 1;
                                    $final_native_alien_genes{$acc}{"mode"} = 2;#Native hence ignore
                                }
                            elsif ($phyletic_distribution_threshold_family_type1 > 30 && $phyletic_distribution_threshold_family_type1 < 40)
                                {
                                    $final_native_alien_genes{$acc}{"type"} = "Ambiguous3";
                                    $final_native_alien_genes{$acc}{"mode"} = 2;#Native hence ignore
                                } 
                            else 
                                {
                                    $final_native_alien_genes{$acc}{"type"} = "Error4";
                                }
                        }
                }
            else
                {
                            $final_native_alien_genes{$acc}{"type"} = "Error";
                }
        }
    # print OUTPUT11 Dumper \%final_native_alien_genes;
    foreach my $a (sort keys %final_native_alien_genes)
        {
            my $b = $species_level_type1{$a}{"phyletic_distribution_threshold"};
            my $c = $genus_level_type1{$a}{"phyletic_distribution_threshold"};
            my $d = $family_level_type1{$a}{"phyletic_distribution_threshold"};
            my $e = $species_level_type2{$a}{"phyletic_distribution_threshold"};
            my $f = $genus_level_type2{$a}{"phyletic_distribution_threshold"};
            my $g = $family_level_type2{$a}{"phyletic_distribution_threshold"};
            my $h = $final_native_alien_genes{$a}{"type"};
            my $t = $final_native_alien_genes{$a}{"mode"};
            my $string = "$a\t$b\t$c\t$d\t$e\t$f\t$g\t$h\t$t\n";
            print OUTPUT11 "$string";
            if ($h == 0)
              {
                if ($t == 1)
                  {
                    my $mode = "Ancient";
                    my $string = "$a\t$t\t$mode\n";
                    print OUTPUT12 "$string";
                  }
                elsif ($t == 0)
                  {
                    my $mode = "Recent";
                    my $string = "$a\t$t\t$mode\n";
                    print OUTPUT12 "$string";
                  }
              }
        }
    ## Selecting only genes and listing as recent or ancient.

    close(IN11);
    close(IN12);
    close(IN14);
    close(IN15);
    close(IN16);
    close(IN17);
    close(IN18);
    close(OUTPUT11);
    close(OUTPUT12);
  }

sub marker_gene_enrichment 
  {
      my $file1 = "$f4\_APP_Alien_genes.txt";
      open IN19, $file1 or die;
      my @alien_genes = <IN19>;
      # print Dumper \@alien_genes;

      my $file2 = glob("*.modified");
      open IN20, $file2 or die;
      my @all_genes = <IN20>;
      # print Dumper \@all_genes;
      # open OUTPUT14, ">FOR_INPUT_all_gene_sequence.fa" or die;

      my $total_genes_across_whole_genome = 0;
      my @all_gene_list_enrichment;
      for(my $i = 0; $i <= scalar(@all_genes)-2; $i = $i + 2)
          {
              my $l1 = $all_genes[$i];
              my @arr1 = split(">", $l1);
              my $header = $arr1[1];
              chomp $header;
              push (@all_gene_list_enrichment, $header);
              # my $seq = $all_genes[$i+1];
              # chomp $seq;
              # $gene_database{$header} = $seq;
              $total_genes_across_whole_genome++;
          }
      # print Dumper \%gene_database;
      
      my $total_alien_genes = 0;
      my @alien_gene_list_enrichment;
      foreach my $g (@alien_genes)    
          {
              chomp $g;
              my @arr = split("\t", $g);
              my $gene = $arr[0];
              # print $gene,"\n";
              # print OUTPUT14 ">$gene\n";
              push(@alien_gene_list_enrichment, $gene);
              # my $seq = $gene_database{$gene};
              # print $seq,"\n";
              # print OUTPUT14 "$seq\n";
              $total_alien_genes++;
          }
      print $total_alien_genes,"\n";
      print $total_genes_across_whole_genome,"\n";    
      print "\nPerforming HMMSCAN. Please wait....\n" if ($FLAG_markerGeneEnrichment == 1);    
      # system (`hmmscan -o $f4\:hmm.txt --tblout $f4\:table.txt reserve/cafe_database_table $file2`);
      my $file3 = glob("*table.txt");
      open IN21, $file3 or die;
      my @hmmscan_data = <IN21>;
      open OUTPUT15, ">$f4\_marker_genes.txt";
      open OUTPUT16, ">$f4\_marker_enrichment_statistics.txt";
      print "\nExtracting marker genes. Please wait....\n" if ($FLAG_markerGeneEnrichment == 1);   
      my $total_marker_gene_count = 0;
      my @marker_gene_list_enrichment;
      for (my $i = 0; $i < scalar(@hmmscan_data); $i++)
          {
              chomp $hmmscan_data[$i];
              next if $hmmscan_data[$i] =~ /^#/;
              my @arr2 = split('\s+', $hmmscan_data[$i]);
              my $gene_accession = $arr2[2];
              my $evalue = $arr2[4];
              my $gene_description = join(' ', @arr2[18..$#arr2]);
              my @arr3 = split ('\s+', $hmmscan_data[$i-1]);
              my $gene_accession1 = $arr3[2];
              next if ($gene_accession eq $gene_accession1 && defined($gene_accession1));
              if ($evalue < 0.05)
                  {
                      $total_marker_gene_count++;
                      print OUTPUT15 $gene_accession, "\n";
                      push (@marker_gene_list_enrichment, $gene_accession);
                  }
              # else {print $evalue,"\n";}

          }
      print "\nPerforming enrichment. Please wait ....\n" if ($FLAG_markerGeneEnrichment == 1);   
      my @union; my @isect; 
      my %union; my %isect;
      foreach my $a_gene (@alien_gene_list_enrichment)
          {
              $union{$a_gene} = 1;
          }
      foreach my $m_gene (@marker_gene_list_enrichment)
          {
              if ($union{$m_gene})
                  {
                      $isect{$m_gene} = 1;
                  }
              $union{$m_gene} = 1;
          }
      my $count_marker_genes_also_marked_alien_by_APP = keys %isect;
      print OUTPUT16 "Out of $total_genes_across_whole_genome genes in the genome, we found $total_marker_gene_count marker genes. APP identified $total_alien_genes genes as of alien origin. Out of $total_alien_genes alien genes, $count_marker_genes_also_marked_alien_by_APP were found to be marker genes.\n";
      my $overlap = $count_marker_genes_also_marked_alien_by_APP;#overlap of set 1 & 2
      my $total = $total_genes_across_whole_genome;#Total no of genes in genome
      my $m = $total_alien_genes;#set 1; no. of alien genes
      my $n = $total - $m;#Total genes - set 1
      my $k = $total_marker_gene_count;#set 2; no. of marker genes across entiregenome
      my $expected_no_of_genes = ($m)*($k/$total);
      my $enrichment = $overlap/$expected_no_of_genes;
      my $x;#overlap calculation for phyper formulae 
      if ($overlap < $expected_no_of_genes)
          {
              $x = $overlap;
          }
      elsif ($overlap > $expected_no_of_genes)
          {
              $x = $overlap-1;
          }
      my $R = Statistics::R->new();
      $R->run( qq`x = phyper($x,$m,$n,$k, lower.tail=FALSE)` );
      my $squares = $R->get('x');
      print OUTPUT16 "Enrichment = $enrichment\n";
      if ($squares > 0.5)
          {
              my $sig = 1 - $squares;
              print OUTPUT16 "p-value = $sig\n";
              
          }
      else 
          {
              print OUTPUT16 "p-value = ", $squares,"\n";
          }
      close(IN19); close(IN20); close(IN21);
      close(OUTPUT15); close(OUTPUT16);
  }

sub get_genbank
  {
      print ("\nDownloading GenBank record from NCBI. Please wait....\n");
      system(`efetch -db nuccore -format gbwithparts -id $f4 > $f4.gb`);
      my $op_g = "$f4.gb";
      if (-e $op_g)
          {
              print "$op_g has been found and downloaded. Post-processing .....\n";
              genbank_to_ptt($op_g);
          }
      else
          {
              print "\n$op_g not found\n";
              print "\nRecheck accession number. Make sure NCBI e-utilities are installed and in path\n";
              exit(0);
          }
      
  }
sub genbank_to_ptt
  {
      my ($a) = @_;
      my $gbk = $a;
      chomp $gbk;
      if (-e $gbk)
          {
              open IN, $gbk or die;
          }
      else {print "\n$gbk doesn't exists\n\n"; exit(0);}
      open TEMPOUT, ">$f4.tbl" or die;
      my $random_counter = 0;
      my $gene_counter = 1;
      while (<IN>)
          {
              if (/^     CDS             /)
                  {
                      $_ =~ s/     CDS             //g;
                      $_ =~ s/<|>//g;
                      chomp $_;
                      if ($_ =~ /complement/)
                          {
                              if ($_ =~ /join/)
                                  {
                                      $_ =~ s/\(|\)|join|complement//g;
                                      my @a = split(/\./, $_);
                                      print TEMPOUT "g$gene_counter\t$a[0]\t$a[scalar(@a)-1]\t-\n"; 
                                  }
                              else 
                                  {
                                      $_ =~ s/\(|\)|complement//g;
                                      my @a = split(/\./, $_);
                                      print TEMPOUT "g$gene_counter\t$a[0]\t$a[2]\t-\n"; 
                                  }
                          }
                      else
                          {
                              if ($_ =~ /join/)
                                  {
                                      $_ =~ s/\(|\)|join//g;
                                      my @a = split(/\./, $_);
                                      print TEMPOUT "g$gene_counter\t$a[0]\t$a[scalar(@a)-1]\t+\n";
                                      my $c1 = $a[0];

                                  }
                              else 
                                  {
                                      my @a = split(/\./, $_);
                                      print TEMPOUT "g$gene_counter\t$a[0]\t$a[2]\t+\n";
                                  }
                          }
                  
                      $gene_counter++;
                  }
          }
      print ("\nFeature table computed for gene map. Please wait....\n");
      close(TEMPOUT);
      gene_map_builder();

  }
sub gene_map_builder
  {
      my $file1 = "$f4\_APP_Alien_genes.txt";
      open IN21, $file1 or die;
      my @alien_genes = <IN21>;

      my $file2 = "$f4.tbl";
      open IN22, $file2 or die;
      my @ptt = <IN22>;
      # print scalar(@ptt),"\n";
      my %ptt_db;
      foreach my $g (@ptt)
          {
              chomp $g;
              my @a = split(/\t/, $g);
              # print $a[0],"\n";
              $ptt_db{$a[0]}{"start"} = $a[1];
              $ptt_db{$a[0]}{"end"} = $a[2];
              $ptt_db{$a[0]}{"strand"} = $a[3];
          } 
      # print Dumper \%ptt_db;
      print "\nGenerating plotting files for CGViewer. Please wait...\n";
      open OUTPUT17, ">cgview_alien_genes.txt";
      print OUTPUT17"seqname	source	feature	start	end	score	strand	frame\n";
      # open OUTPUT18, ">cgview_Recent_alien_genes.txt";
      # print OUTPUT18"seqname	source	feature	start	end	score	strand	frame\n";
      foreach my $l (@alien_genes)
          {
              chomp $l;
              next if $l !~ /:g/;
              my @a = split("\t", $l);
              my @b = split(":", $a[0]);
              my $gene = $b[scalar(@b)-1];
              
              my $gene_start = $ptt_db{$gene}{"start"};
              my $gene_end = $ptt_db{$gene}{"end"};
              my $gene_strand = $ptt_db{$gene}{"strand"};
              my $transfer = $a[2];
              if ($gene_strand eq '+')
                  {
                      if ($transfer eq "Ancient")
                          {
                              print OUTPUT17 "$gene\t\tF\t$gene_start\t$gene_end\t\t+\n";
                          }
                      elsif ($transfer eq "Recent")
                          {
                              print OUTPUT17 "$gene\t\tA\t$gene_start\t$gene_end\t\t+\n";
                          }
                  }
              elsif ($gene_strand eq '-')
                  {
                      if ($transfer eq "Ancient")
                          {
                              print OUTPUT17 "$gene\t\tF\t$gene_start\t$gene_end\t\t-\n";
                          }
                      elsif ($transfer eq "Recent")
                          {
                              print OUTPUT17 "$gene\t\tA\t$gene_start\t$gene_end\t\t-\n";
                          }
                  }
              
          }
      close(OUTPUT17);
      print "\nRunning CGViewer. Please wait...\n";
      if ($verboseDetailed == 1)
          {
              system (`perl cgview_xml_builder.pl -sequence $f4.gb -out $f4.xml -genes cgview_alien_genes.txt  -details F -show_sequence_feature F -gene_labels F -custom moveInnerLabelsToOuter=true tickLength=20 featureOpacity=1 featureOpacityOther=1 backboneRadius=1000 -legend T `);
          }
      else 
          {
              system (`perl cgview_xml_builder.pl -verbose F -sequence $f4.gb -out $f4.xml -genes cgview_alien_genes.txt  -details F -show_sequence_feature F -gene_labels F -custom moveInnerLabelsToOuter=true tickLength=20 featureOpacity=1 featureOpacityOther=1 backboneRadius=1000 -legend T `);
          }
      
      my $file3 = "$f4.xml";
      open IN23, $file3 or die;
      my @xml = <IN23>;
      open OUTPUT18, ">$f4\_modified.xml";
      foreach my $l (@xml)
          {
              # chomp $l;
              if ($l =~ /A COG/) {$l =~ s/A COG/Recent/g;}
              if (($l =~ /F COG/)) {$l =~ s/F COG/Ancient/g;}
              next if $l =~ /B COG|J COG|K COG|L COG|D COG|O COG|M COG|N COG|P COG|T COG|U COG|V COG|W COG|Y COG|Z COG|C COG|G COG|E COG|H COG|CDS|I COG|Q COG|R COG|S COG|Unknown COG|legendItem text="tRNA"|legendItem text="rRNA"|legendItem text="Other"/;
              print OUTPUT18 $l;
          }
      close(OUTPUT18);
      system(`java -jar -Xmx2000m cgview.jar -i $f4\_modified.xml -o $f4\_HGT_map.png -f png -D 48`);
      close(IN21);close(IN22);close(IN23);

  }



sub process_fasta
  {
    # print "process fasta\n";
    my ($a) = @_;
    my $fna = $a;
    chomp $fna;
    if (-e $fna)
      {
        open IN, $fna or die;
      }
    else {print "\n$fna doesn't exists\n\n"; exit(0);}
    my $random_counter = 0;
    open TEMPOUT, ">$fna\.modified" or die;
    while (<IN>)
        {
            if (/^>/)
                {
                    if ($random_counter > 0) {print TEMPOUT "\n"};
                    $random_counter++;
                    $_ =~ s/(-->)|(->)//g;
                    my @arr1 = split(">", $_);
                    my @arr2 = split(/(\]\ \[)/, $arr1[1]);
                    my $FLAG = 0;
                    my $real_header = $_; #for default headers
                    foreach my $key (@arr2)
                        {
                            chomp $key;
                            if ($key =~ /(protein_id=)/)
                                {
                                    my $acc_no = $key;#Protein ID
                                    $acc_no =~ s/(protein_id=)//g;
                                    if ($acc_no)
                                        {
                                           my $new_header = $acc_no.":g".$random_counter;#for blast database filename purposes
                                           print TEMPOUT ">$new_header\n";
                                        }
                                    $FLAG = 1;
                                }
                        }
                    if ($FLAG == 0)
                                {
                                    chomp $real_header;
                                    my @arr3 = split(">", $real_header);
                                    my $fname = $arr3[1];
                                    my $new_header = $fname.":g".$random_counter;#for blast database filename purposes
                                    print TEMPOUT ">$new_header\n";
                                }
                }
                else
                  {
                      chomp $_;
                      print TEMPOUT $_;
                  }
        }
        close (TEMPOUT);
        close (IN);
      print "Modified $fna has been generated with revised fasta headers.\n\n";
  }
sub obtain_faa_eutilities
  {
    my ($a) = @_;
    my $acc_no = $a;
    system(`efetch -db nuccore -format fasta_cds_aa -id $acc_no > $acc_no.faa`);
    my $op_f = "$acc_no\.faa";
    if (-e $op_f)
        {
            print "$op_f has been found and downloaded. Post-processing .....\n";
            process_fasta($op_f);
        }
    else
        {
           print "\n$op_f not found\n";
           print "\nRecheck accession number. Make sure NCBI e-utilities are installed and in path\n";
           exit(0);

        }
  }

sub print_usage {
    print <<BLOCK;

APP.pl - performs phyletic pattern analysis to predict horizontally acquired genes.
This is version 1.0, created by Soham Sengupta.

DISPLAY HELP AND EXIT:

usage:

  perl APP.pl --help

PERFORM ALIEN GENE DETECTION (PHYLETIC PATTERN ANALYSIS)

usage:

  perl APP.pl -q <query fileName> -t <query taxonomy file> -o <Output fileName> -f <fileType> [Options]

required arguments:

-q - multifasta amino acid file or genome acession number (NCBI only).

-t - query specific taxonomic id file.

-o - Output file name to create.

-f - File type as 'multifasta' or 'accession'.

optional arguments:

-n - No. of CPU cores to use for performing blast. By default, uses all available cores.

-m - Marker gene enrichment. Default set to 0. Use 1 only when whole genome is analyzed.

-g - Generate circular gene map with alien genes. Default set to 0. Use 1 only when whole genome is analyzed.

-e - turn on Expert option (1; keep temporary and intermediate files). Default set to 0. 

-v - Provide basic progress messages. Default set to 1.

-d - Provide detailed progress messages. Default set to 0.

example usage:
1. Input multifasta file
  perl APP.pl -q example/NC_004088.faa -t 3.query_speciesID_taxID.txt -o NC_004088 -f multifasta
2. Accession number of genomes
  perl APP.pl -q NC_004088 -t 3.query_speciesID_taxID.txt -o NC_004088 -f accession
3. Input multifasta file with marker gene enrichment
  perl APP.pl -q example/NC_004088.faa -t 3.query_speciesID_taxID.txt -o NC_004088 -f multifasta -m 1
BLOCK
}