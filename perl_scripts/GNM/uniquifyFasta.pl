#!/usr/bin/perl -w 
use strict ;
use FileHandle ;
use Getopt::Long;
use Cwd ;
use MyUtils;
use POSIX qw(floor);
my $commandline = util_get_cmdline("",\@ARGV) ;
my ($all,$infile,$outfile,$or,$silent,$groupinfo);
my ($tag,$size,$DIR,$listfile,$ignorefile,$mergedir);
my $howmany = 600000 ; 
my $cutofflength = 0 ; 
my @types = (); 
my $protein ;
my @motifs = (); 
my $createnewname = 0 ;
my $writedata = 0 ;
my $reverse = 0 ;
my $printbande = 0 ;
use Memory::Usage;
my $mu = Memory::Usage->new();
# Record amount of memory used by current process
$mu->record('starting work');

GetOptions(
            "all"=>\$all ,
            "groupinfo"=>\$groupinfo ,
            "silent"=>\$silent ,
            "infile=s"=>\$infile ,
            "mergedir=s"=>\$mergedir ,
            "dir=s"=>\$DIR ,
            "listfile=s"=>\$listfile ,
            "ignorefile=s"=>\$ignorefile ,
            "createnewname=s"=>\$createnewname ,
            "writedata=i"=>\$writedata ,
            "printbande=i"=>\$printbande ,
            "reverse=i"=>\$reverse ,
            "size=i"=>\$size ,
            "or=i"=>\$or ,
            "cutofflength=i"=>\$cutofflength ,
            "type=s"=>\@types,
            "protein=s"=>\$protein,
            "motif=s"=>\@motifs,
            "outfile=s"=>\$outfile 
           );
die "Dont recognize command line arg @ARGV " if(@ARGV);
usage( "Need to give a input file name => option -infile ") if(!defined $infile);
usage( "Need to give a output file name => option -outfile ") if(!defined $outfile);
my $ofh = util_write("$outfile");
my ($RESULTDIR,$PDBDIR,$FASTADIR,$APBSDIR,$FPOCKET,$SRC,$MATCH3D,$ANNDIR, $UNIPROT) = util_SetEnvVars();
print "INFO: Parameters are writedata=$writedata, createnewname=$createnewname and reverse=$reverse,FASTADIR =$FASTADIR\n";
print STDERR "Info: parsing file $infile - might take some time\n";
my ($info,$infoSeq2PDB,$mapChainedName2Name) = util_parsePDBSEQRESNEW($infile,0,$writedata,$reverse);

my $ofhmappingname = util_write("$outfile.mappingname");
my $ofhmappingsame = util_write("$outfile.mappingsame");
my $ofhlength = util_write("$outfile.length");
my $ofhlist = util_write("$outfile.list");

my $ofhb100 = util_write("$outfile.begin.100");
my $ofhe100 = util_write("$outfile.end.100");
my $ofhb30 = util_write("$outfile.begin.30");
my $ofhe30 = util_write("$outfile.end.30");

## map seq and names to numbers
my $TRSID = 0 ;
my $MAPCNT2NAME = {};
my $MAPCNT2SEQ = {};
my $MAPSEQ2CNT = {};

foreach my $seq (keys %{$infoSeq2PDB}){
	$TRSID++;
	my $len = length($seq);
	my @l = sort @{$infoSeq2PDB->{$seq}};
	my $name = shift @l ;

	## create new name? else mappingname remains the same
	my $ID = $createnewname ? $createnewname. "_" . $TRSID: $name ;

	print $ofhlength "$ID $len \n";
	print $ofhlist "$ID\n";
	print $ofh ">$ID\n";
	print $ofh "$seq\n";
	print $ofhmappingname "$ID $name \n";

	## mutiple
	if(@l){
	   print $ofhmappingsame "$ID ";
	   foreach my $l (@l){
	       print $ofhmappingsame "$l\n";
	   }
	}

	if($printbande){
	    PrintBandE($ID,$seq,100,$ofhb100,$ofhe100);
	    PrintBandE($ID,$seq,30,$ofhb30,$ofhe30);
	}

}

# Record amount in use afterwards
$mu->record('after something_memory_intensive()');

# Spit out a report
$mu->dump();


sub PrintBandE{
	my ($ID,$seq,$N,$ofhb,$ofhe) = @_ ;
	my ($begin,$end) = util_findLength($seq,$N);
	print $ofhb ">$ID\n";
	print $ofhb "$begin\n";
	print $ofhe ">$ID\n";
	print $ofhe "$end\n";
}

print "Uniquifying $infile to $outfile, wrote $TRSID\n";

system("sortOnLast.pl -rev -in $outfile.length");
system("wc -l $outfile.*");


sub usage{
    my ($msg) = @_ ;
    print $msg , "\n" ; 
print << "ENDOFUSAGE" ; 
ENDOFUSAGE
    die ;
}
