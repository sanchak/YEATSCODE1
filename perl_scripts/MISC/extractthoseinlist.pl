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
my ($duplicate,$blastout,$infile,$p1,$p2,$outfile,$trs,$cutoff,$which_tech,$listfile,$protein);
my ($chooseone,$ignorefile,$inorderoflist,$tag,$donone,@expressions);
my $howmany = 100000 ;
my $verbose = 0 ;
GetOptions(
            "which_tech=s"=>\$which_tech ,
            "trs=s"=>\$trs ,
            "tag=s"=>\$tag ,
            "blastout=s"=>\$blastout ,
            "infile=s"=>\$infile ,
            "p1=s"=>\$p1 ,
            "p2=s"=>\$p2 ,
            "listfile=s"=>\$listfile ,
            "ignorefile=s"=>\$ignorefile ,
            "outfile=s"=>\$outfile ,
            "inorderoflist"=>\$inorderoflist ,
            "donone"=>\$donone ,
            "duplicate"=>\$duplicate ,
            "expr=s"=>\@expressions,
            "howmany=i"=>\$howmany ,
            "verbose=i"=>\$verbose ,
            "cutoff=f"=>\$cutoff ,
           );
die "Dont recognize command line arg @ARGV " if(@ARGV);

usage( "Need to give a listfile -option -listfile  ") if(!defined $listfile);
usage( "Need to give a infile -option -infile  ") if(!defined $infile);
$tag = "ttt" if(!defined $tag);
$outfile = "$tag.$infile" if(!defined $outfile);
$outfile =~ s/\///g;
my $ofh = util_write($outfile);
my $ifh = util_read($infile);

my $table = {};
my @list= util_read_list_words($listfile);
map {  my @ll = split ;  $table->{$ll[0]} = 1 ; } @list ;

my $cnt = 0 ;
my $missed = 0 ;
my $printed = {};
while(<$ifh>){

   my $orog = $_;

   # just hardcoded, essentially we take the first
   s/ln -s//;
   s/.ORF.*//;
   my ($trs) = split ;
   if(exists $table->{$trs}){
   	   $cnt++;
   	   print $ofh "$orog" if(!defined $inorderoflist);
	   $printed->{$trs} = $orog;
   }
}

if(defined $inorderoflist){
	foreach my $i (@list){
	   if(exists $printed->{$i}){
	      my $orog = $printed->{$i};
   	      print $ofh "$orog" ;
	   }
	}
}

foreach my $k (keys %{$table}){
	if(! exists $printed->{$k}){
		$missed++;
	}
}

print "Wrote $cnt to $outfile, missed $missed\n";
   
sub usage{
    my ($msg) = @_ ;
    print $msg , "\n" ; 
print << "ENDOFUSAGE" ; 
ENDOFUSAGE
    die ;
}
