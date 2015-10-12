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
my ($idx,$infile,$p1,$p2,$outfile,$cutoff,$which_tech,$listfile,$protein);
my ($datadir,$ignorefile,@expressions);
my $howmany = 100000 ;
my $verbose = 1 ;
GetOptions(
            "which_tech=s"=>\$which_tech ,
            "datadir=s"=>\$datadir ,
            "protein=s"=>\$protein ,
            "infile=s"=>\$infile ,
            "p1=s"=>\$p1 ,
            "p2=s"=>\$p2 ,
            "listfile=s"=>\$listfile ,
            "ignorefile=s"=>\$ignorefile ,
            "outfile=s"=>\$outfile ,
            "expr=s"=>\@expressions,
            "howmany=i"=>\$howmany ,
            "verbose=i"=>\$verbose ,
            "idx=i"=>\$idx ,
            "cutoff=f"=>\$cutoff ,
           );
die "Dont recognize command line arg @ARGV " if(@ARGV);
usage( "Need to give a input file name => option -infile ") if(!defined $infile);
usage( "Need to give a input file name => option -outfile ") if(!defined $outfile);
usage( "Need to give a input file name => option -datadir ") if(!defined $datadir);
my $ifh = util_read($infile);
my $ofh = util_write($outfile);

my $CNT = 0 ; 
#my ($RESULTDIR,$PDBDIR,$FASTADIR,$APBSDIR,$FPOCKET,$SRC) = util_SetEnvVars();
 my ($RESULTDIR,$PDBDIR,$FASTADIR,$APBSDIR,$FPOCKET,$SRC,$MATCH3D,$ANNDIR, $UNIPROT,$PREMONITION,$HELIXDIR,$DSSP,$CONFIGGRP,$BLASTOUT,$BLASTDB) = util_SetEnvVars();


my $info = {};
my $doneID = {};
my $table =  {};
while(<$ifh>){
     next if(/^\s*$/);
	 next if(/^\s*#/);

	 last if($CNT > $howmany) ;
	 my ($nm,$id,$sub,$lineage) = split ; 
	 my @l ;
	 my $NM ;
	 if(defined $idx){
	     $lineage =~ s/;/   /g;
	     @l = split " " , $lineage ;
	     $NM = $l[$idx];
	     $NM =~ s/\//./g;
	 }
	 else{
	     $NM = $sub;
	 }
	 $table->{$NM} = 1 ;

	 if(exists $doneID->{$id}){
	 	 die if($doneID->{$id} ne $NM);
		 next ;
	 }
	 $doneID->{$id} = $NM ;
	 if(! defined $info->{$NM}){
	 	$info->{$NM} = util_write("$datadir/$NM.$outfile");
	 }
	 my $fh = $info->{$NM} ;
	 print $fh $_ ;
}

foreach my $k (keys %{$table}){
	print $ofh "$k\n";
}

#system ("wc -l $datadir/*.$outfile");

sub usage{
    my ($msg) = @_ ;
    print $msg , "\n" ; 
print << "ENDOFUSAGE" ; 
ENDOFUSAGE
    die ;
}
