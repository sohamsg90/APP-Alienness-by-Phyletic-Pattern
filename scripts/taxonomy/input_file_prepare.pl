use strict;
use warnings;


my $f1 = $ARGV[0];
if ($f1)
	{
		system(`sh scripts/taxonomy/extract_taxIDs.sh $f1`);
	}
else {print "\nUsage: perl 2.input_file_prep.pl acc_no.txt\n";}
open IN, $f1 or die "\nNo input file provided\n";
my @acc_nos = <IN>;

my $f2 = "tax_ID_of_interest_lineage_ID_format.txt";
open IN2, $f2 or die;
my @ncbi_taxids = <IN2>;

open OUT, ">3.query_speciesID_taxID.txt";
print OUT "Acc\tspecies\ttaxid\n";
my $i = 0;
foreach my $l (@ncbi_taxids)
	{
		chomp $l;
		# print $l,"\n";
		my $FLAG = 0;
		my ($taxid, $speciesid);
		if ($l !~ /^taxid|\-/)
			{
				# print $l,"\n";
				my @arr = split("\t", $l);
				$taxid = $arr[0];
				$speciesid = $arr[7];
				$FLAG = 1;
			}
		if ($FLAG == 1) {chomp $acc_nos[$i];print OUT "$acc_nos[$i]\t$speciesid\t$taxid\n";$i++;}
	}
system(`rm tax_ID*`);






