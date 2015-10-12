#!/usr/bin/perl -w 
use strict ;
use FileHandle ;
use Getopt::Long;
use Cwd ;
use MyUtils;
use MyGeom;
use PDB;
use ConfigPDB;
use MyGNM;
use Math::Geometry ;
use Math::Geometry::Planar;
 no warnings "all";


use Time::HiRes qw( usleep ualarm gettimeofday tv_interval clock_gettime clock_getres  clock);
use POSIX qw(floor);
my $commandline = util_get_cmdline("",\@ARGV) ;
my ($blastdir,$force,$justchecking,$checkforsame,$infile,$p1,$p2,$cutoff,$which_tech,$listfile,$protein);
my ($promoters,$scaffoldfasta,$ignorefile);
my $howmany = 100000 ;
my $verbose = 1 ;

my $percentlength = 10;
my $percentmatched = 70;
my $percentidentity = 30;


GetOptions(
            "force"=>\$force ,
            "protein=s"=>\$protein ,
            "infile=s"=>\$infile ,
            "blastdir=s"=>\$blastdir ,
            "scaffoldfasta=s"=>\$scaffoldfasta ,
            "checkforsame"=>\$checkforsame ,
            "p1=s"=>\$p1 ,
            "p2=s"=>\$p2 ,
            "listfile=s"=>\$listfile ,
            "ignorefile=s"=>\$ignorefile ,
            "promoters"=>\$promoters ,
            "howmany=i"=>\$howmany ,
            "verbose=i"=>\$verbose ,
            "percentlength=i"=>\$percentlength ,
            "percentmatched=i"=>\$percentmatched ,
            "percentidentity=i"=>\$percentidentity ,
            "cutoff=f"=>\$cutoff ,
           );

my ($RESULTDIR,$PDBDIR,$FASTADIR,$APBSDIR,$FPOCKET,$SRC,$MATCH3D,$ANNDIR, $UNIPROT,$PREMONITION,$HELIXDIR,$DSSP,$CONFIGGRP,$BLASTOUT,$BLASTDB) = util_SetEnvVars();
my $HARD_DIFFFORMATCH = 5;
my $HARD_LONGDISTANCEFACTOR = 10;


## hardcoded 
$blastdir = "BLASTOUT_SCAFFOLD/";
my $SCAFFOLDDIR  = "SCAFFOLDDIR/";
my $blasterrs = "blast.err";:
my $infodir = "INFO/";
#usage( "Need to give a input file name => option -blastdir ") if(!defined $blastdir);
#usage( "Need to give a input file name => option -infile ") if(!defined $infile);
#usage( "Need to give a input file name => option -scaffoldfasta ") if(!defined $scaffoldfasta);



die "blast.err" if( ! -e "blast.err");
my ($errtable,$N) = util_maketablefromfile("blast.err");



#my $ifh = util_read($infile);
my $ERRSTATE = 0 ;
print "Checking for BLAST OUTS \n";
if(1){

     #next if(/^\s*$/);
	 #next if(/^\s*#/);
     my ($sname, @trss) = @ARGV;

   foreach my $trs (@trss){
       my $infile = "$blastdir/$trs.$sname.blast";
       if(! -e $infile){
	       warn " $infile does not exist ";
	       $ERRSTATE = 1 ;
	       next ;
       }
  }
}
die "Errors : see log " if($ERRSTATE);
print "Checking for BLAST OUTS done .......\n";
  
my $ofhlog = util_open_or_append("logfile");
if(1){

     my ($sname, @trss) = @ARGV;
     $scaffoldfasta = "$SCAFFOLDDIR/$sname.ALL.1.fasta" ;
     my ($infoScaff,$infoSeq2PDB,$mapChainedName2Name) = util_parsePDBSEQRESNEW("$scaffoldfasta",0,0);
     my ($strScaffold) = $infoScaff->{$sname};
     my $scafflen = length($strScaffold);


        my $outfile = "$infodir/$sname.info";
		if(-e $outfile && !defined $force){
               print "Already done $outfile\n" ;
		       exit ;
        }

     print "Processing $sname $scafflen \n";
     my $strScaffold = $infoScaff->{$sname};
     foreach my $trs (@trss){

        my $infile = "$blastdir/$trs.$sname.blast";
		if(exists $errtable->{$infile}){
			print "Doing next for $infile since it has errors\n";
			next ;
		}

        my $ofh = util_open_or_append($outfile);
        my $str = GNM_MapOneTRS2Scaffold($trs,$sname,$infile,$strScaffold,$scafflen,$ofh,$ofhlog,$HARD_DIFFFORMATCH,$HARD_LONGDISTANCEFACTOR);
		close($ofh);
  }

}



sub usage{
    my ($msg) = @_ ;
    print $msg , "\n" ; 
print << "ENDOFUSAGE" ; 
ENDOFUSAGE
    die ;
}

