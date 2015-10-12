#!/usr/bin/perl -w 
use strict ;
use FileHandle ;
use Getopt::Long;
use Cwd ;
use MyUtils;
use MyGeom;
use PDB;
use ConfigPDB;
use Math::Geometry ;
use Math::Geometry::Planar;


use Time::HiRes qw( usleep ualarm gettimeofday tv_interval clock_gettime clock_getres  clock);
use POSIX qw(floor);
my $commandline = util_get_cmdline("",\@ARGV) ;
my ($nostderr,$infile,$p1,$p2,$outfile,$cutoff,$sleeptime,$listfile,$interval);
my (@expressions);
my $howmany = 100000 ;
my $verbose = 1 ;
GetOptions(
            "sleeptime=s"=>\$sleeptime ,
            "nostderr"=>\$nostderr ,
            "interval=s"=>\$interval ,
            "infile=s"=>\$infile ,
            "p1=s"=>\$p1 ,
            "p2=s"=>\$p2 ,
            "listfile=s"=>\$listfile ,
            "outfile=s"=>\$outfile ,
            "expr=s"=>\@expressions,
            "howmany=i"=>\$howmany ,
            "verbose=i"=>\$verbose ,
            "cutoff=f"=>\$cutoff ,
           );
die "Dont recognize command line arg @ARGV " if(@ARGV);
$outfile = "Sch.". $infile  if(!defined $outfile);
my $ofh = util_write($outfile);

usage( "Need to give a input file name => option -infile ") if(!defined $infile);
my $ifh = util_read($infile);
usage( "Need to give a interval pdb id -option -interval  ") if(!defined $interval);
#my ($RESULTDIR,$PDBDIR,$FASTADIR,$APBSDIR,$FPOCKET,$SRC) = util_SetEnvVars();
 my ($RESULTDIR,$PDBDIR,$FASTADIR,$APBSDIR,$FPOCKET,$SRC,$MATCH3D,$ANNDIR, $UNIPROT,$PREMONITION,$HELIXDIR,$DSSP,$CONFIGGRP,$BLASTOUT,$BLASTDB) = util_SetEnvVars();


my $CNT = 0 ;

my $strtoadd = defined $nostderr ? " & " : "> & ! /dev/null &" ;
while(<$ifh>){
     chomp ;
     next if(/^\s*$/);
     if(!/.dev.null/){
        $_ = $_ . $strtoadd ;
     }
     	print $ofh "$_ \n" ;
     	$CNT++;
	 if($CNT > $interval){
	 	print $ofh "sleep $sleeptime\n";
		$CNT = 0 ;
	 }
}

print "Wrote to outfile $outfile\n";



print "Wrote to outfile $outfile\n";

sub usage{
    my ($msg) = @_ ;
    print $msg , "\n" ; 
print << "ENDOFUSAGE" ; 
ENDOFUSAGE
    die ;
}
