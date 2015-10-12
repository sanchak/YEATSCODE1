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
my ($fastadir,$infile,$p1,$p2,$orfdir,$outfile,$cutoff,$listfile,$orf,$trs);
my ($fastafile,@expressions);
my $howmany = 100000 ;
my $verbose = 0 ;
GetOptions(
            "listfile=s"=>\$listfile ,
            "trs=s"=>\$trs ,
            "infile=s"=>\$infile ,
            "p1=s"=>\$p1 ,
            "p2=s"=>\$p2 ,
            "orf=s"=>\$orf ,
            "fastadir=s"=>\$fastadir ,
            "outfile=s"=>\$outfile ,
            "orfdir=s"=>\$orfdir ,
            "expr=s"=>\@expressions,
            "howmany=i"=>\$howmany ,
            "verbose=i"=>\$verbose ,
            "cutoff=f"=>\$cutoff ,
           );
die "Dont recognize command line arg @ARGV " if(@ARGV);
usage( "Need to give a output file name => option -outfile ") if(!defined $outfile);
usage( "Need to give a output file name => option -orfdir ") if(!defined $orfdir);
my $ofh = util_write($outfile);
usage( "Need to give a listfile -option -listfile  ") if(!defined $listfile);
usage( "Need to give a fastadir -option -fastadir  ") if(!defined $fastadir);

my $info = {};

my $ignoredduetounknown = 0 ;
my $ignoredduetosize = 0 ;
my $actuallprocessed = 0 ;
my $IFH = util_read($listfile);
my $exactcodingseqs  = {};
while(<$IFH>){
	next if(/^\s*$/);
	next if(/^\s*#/);
	 my ($trs,$orf) = split ; 
	 my $orfile = "$orfdir/$trs.orf" ;
	 my $fastafile = "$fastadir/$trs.ALL.1.fasta";
	 my ($NT,$AA,$exactcodingseq,$X,$Y,$Z) = util_ProcessOneORF($trs,$orfile,$fastafile,$orf,$info,$cutoff);
	 $exactcodingseqs->{$exactcodingseq} = 1;
	 $actuallprocessed = $actuallprocessed + $X;
	 $ignoredduetounknown = $ignoredduetounknown + $Y;
	 $ignoredduetosize = $ignoredduetosize + $Z;
}

my $N = (keys %{$exactcodingseqs});
print "There are $N exactcodingseqs \n";

util_EmitCodonBiasfromInfo($ofh,$info);



print "actuallprocessed ignoredduetounknown ignoredduetosize \n";
print "$actuallprocessed $ignoredduetounknown $ignoredduetosize \n";

sub usage{
    my ($msg) = @_ ;
    print $msg , "\n" ; 
print << "ENDOFUSAGE" ; 
ENDOFUSAGE
    die ;
}
