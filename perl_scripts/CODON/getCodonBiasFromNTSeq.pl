#!/usr/bin/perl -w 
use strict ;
use FileHandle ;
use Getopt::Long;
use Cwd ;
use MyUtils;
use POSIX qw(floor);
my $commandline = util_get_cmdline("",\@ARGV) ;
my ($dbfile,$all,$infile,$outfile,$or,$silent,$groupinfo);
my ($size,$DIR,$listfile,$ignorefile);
my $howmany = 600000 ; 
my $cutofflength = 0 ; 
my @types = (); 
my $protein ;
my @motifs = (); 
my $verbose = 0 ;
GetOptions(
            "all"=>\$all ,
            "groupinfo"=>\$groupinfo ,
            "silent"=>\$silent ,
            "infile=s"=>\$infile ,
            "dbfile=s"=>\$dbfile ,
            "dir=s"=>\$DIR ,
            "listfile=s"=>\$listfile ,
            "ignorefile=s"=>\$ignorefile ,
            "howmany=i"=>\$howmany ,
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

my ($info,$infoSeq2PDB,$mapChainedName2Name) = util_parsePDBSEQRESNEW($infile,0);
die if(keys %{$infoSeq2PDB} ne 1);

## get the protein name and sequence...
my ($protseq) = (keys %{$infoSeq2PDB});
my @QQ = @{$infoSeq2PDB->{$protseq}};
my @RR = split " ", $QQ[0];
my $protname = $RR[0];


my $codoninfo = {};
my $AASEQ = util_GetCodonBiasFromNucleotide($codoninfo,$protseq);


my $OFH = util_write("codonbias");
util_EmitCodonBiasfromInfo($OFH,$codoninfo);
system("processCodonBias.pl -outf codonbias.1 -inf codonbias -getAA");
system("processCodonBias.pl -outf codonbias.2 -inf codonbias.1 ");


sub usage{
    my ($msg) = @_ ;
    print $msg , "\n" ; 
print << "ENDOFUSAGE" ; 
ENDOFUSAGE
    die ;
}
