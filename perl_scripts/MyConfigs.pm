package MyConfigs;
use Carp ;
use POSIX ;
require Exporter;
use MyUtils;
no warnings 'redefine';
my $EPSILON = 0.01;

local $SIG{__WARN__} = sub {};

my $RAY = 1 ; 

@ISA = qw(Exporter);
@EXPORT = qw( 
	Config_Helix
	Config_GetUncharStrings
	Config_getCodonTable
	Config_tableforRevComl
	Config_getCountStr
	Config_getTissueList
	Config_NTFastaUnknown
	    );

use strict ;
use FileHandle ;
use Getopt::Long;

sub Config_Helix{
   my $tableATOM = {};
   $tableATOM->{LYS} = "NZ";
   $tableATOM->{ARG} = "NH1";
   $tableATOM->{HIS} = "ND1";
   $tableATOM->{SER} = "OG";
   $tableATOM->{THR} = "OG1";
   $tableATOM->{CYS} = "SG";
   $tableATOM->{ASN} = "OD1";
   $tableATOM->{GLN} = "NE2";
   $tableATOM->{TYR} = "CB";
   $tableATOM->{PHE} = "CB";
   $tableATOM->{TRP} = "CB";
   $tableATOM->{ASP} = "OD1";
   $tableATOM->{GLU} = "OE1";
   $tableATOM->{GLY} = "CA";
   $tableATOM->{VAL} = "CB";
   $tableATOM->{ALA} = "CB";
   $tableATOM->{LEU} = "CB";
   $tableATOM->{MET} = "CB";
   $tableATOM->{ILE} = "CB";
   $tableATOM->{PRO} = "CB";

    my $HYDROVAL = {};
    ## from "Computer programs to identify and classify amphipathic a helical domains"
    $HYDROVAL->{MET} = 0.975;
    $HYDROVAL->{ILE} = 0.913;
    $HYDROVAL->{LEU} = 0.852;
    $HYDROVAL->{VAL} = 0.811;
    $HYDROVAL->{CYS} = 0.689;
    $HYDROVAL->{ALA} = 0.607;
    $HYDROVAL->{THR} = 0.525;
    $HYDROVAL->{GLY} = 0.484;
    $HYDROVAL->{SER} = 0.402;
    $HYDROVAL->{HIS} = 0.333;
    $HYDROVAL->{PRO} = 0.239;
    $HYDROVAL->{PHE} = 1.036;
    $HYDROVAL->{TRP} = 0.668;
    $HYDROVAL->{TYR} = 0.137;
    $HYDROVAL->{GLN} = -0.558;
    $HYDROVAL->{ASN} = -0.701;
    $HYDROVAL->{GLU} = -1.396;
    $HYDROVAL->{LYS} = -1.518;
    $HYDROVAL->{ASP} = -1.600;
    $HYDROVAL->{ARG} = -2.233;

    my $colortable = {};
    my $value = {};
    my $chargedtable = {};
    my @NP = qw(G P A V L I M );
    foreach my $k (@NP){
    	$colortable->{$k} = "red";
    	$value->{$k} = 100 ;
    }
    
    my @pos = qw (H K R);
    foreach my $k (@pos){
    	$chargedtable->{$k} = 1 ;
    	$colortable->{$k} = "blue";
    	$value->{$k} = 100 ;
    }
    
    my @neg = qw (E D); 
    foreach my $k (@neg){
    	$chargedtable->{$k} = -1 ;
    	$colortable->{$k} = "blue";
    	$value->{$k} = 60 ;
    }
    my @amide = qw(Q N );
    foreach my $k (@amide){
    	$colortable->{$k} = "blue";
    	$value->{$k} = 10 ;
    }
    
    my @misc = qw(F Y W C S T);
    foreach my $k (@misc){
    	$colortable->{$k} = "red";
    	$value->{$k} = 75 ;
    }

	return ($tableATOM,$HYDROVAL,$colortable,$value,$chargedtable);
}


sub Config_GetUncharStrings{
my $UNCHARSTRINGS = "unknown|chromosome|hypothe|unnamed|uncharacterized|Predicted protein";
return $UNCHARSTRINGS ;
}

sub Config_getCodonTable{
my %DNA_code = (
'GCT' => 'A', 'GCC' => 'A', 'GCA' => 'A', 'GCG' => 'A', 'TTA' => 'L',
'TTG' => 'L', 'CTT' => 'L', 'CTC' => 'L', 'CTA' => 'L', 'CTG' => 'L',
'CGT' => 'R', 'CGC' => 'R', 'CGA' => 'R', 'CGG' => 'R', 'AGA' => 'R',
'AGG' => 'R', 'AAA' => 'K', 'AAG' => 'K', 'AAT' => 'N', 'AAC' => 'N',
'ATG' => 'M', 'GAT' => 'D', 'GAC' => 'D', 'TTT' => 'F', 'TTC' => 'F',
'TGT' => 'C', 'TGC' => 'C', 'CCT' => 'P', 'CCC' => 'P', 'CCA' => 'P',
'CCG' => 'P', 'CAA' => 'Q', 'CAG' => 'Q', 'TCT' => 'S', 'TCC' => 'S',
'TCA' => 'S', 'TCG' => 'S', 'AGT' => 'S', 'AGC' => 'S', 'GAA' => 'E',
'GAG' => 'E', 'ACT' => 'T', 'ACC' => 'T', 'ACA' => 'T', 'ACG' => 'T',
'GGT' => 'G', 'GGC' => 'G', 'GGA' => 'G', 'GGG' => 'G', 'TGG' => 'W',
'CAT' => 'H', 'CAC' => 'H', 'TAT' => 'Y', 'TAC' => 'Y', 'ATT' => 'I',
'ATC' => 'I', 'ATA' => 'I', 'GTT' => 'V', 'GTC' => 'V', 'GTA' => 'V',
'GTG' => 'V',);

return \%DNA_code ;

}

sub Config_getCountStr{
return " TRS   CE   CI   CK   EM   FL   HC   HL   HP   HU   IF   LE   LM   LY   PK   PL   PT   RT   SE   TZ   VB";
}

sub Config_getTissueList{
my @l = qw (CE   CI   CK   EM   FL   HC   HL   HP   HU   IF   LE   LM   LY   PK   PL   PT   RT   SE   TZ   VB);
return @l ;
}


sub Config_NTFastaUnknown{
return "R|Y|K|M|S|W|B|D|H|V|N";
}


sub Config_tableforRevComl{
    my $tableforRevComl = {};
    $tableforRevComl->{"A"} = "T";
    $tableforRevComl->{"T"} = "A";
    $tableforRevComl->{"G"} = "C";
    $tableforRevComl->{"C"} = "G";
    return $tableforRevComl ;
}



