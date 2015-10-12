#!/usr/bin/perl -w 
use strict ;
use FileHandle ;
use Getopt::Long;
use MyUtils;
use MyConfigs;
use Memory::Usage;
use MyGNM;
my $mu = Memory::Usage->new();
$mu->record('');


my $debugone = 0 ;
my $savetable = 0 ; ## applies only to seqeunce, and not genome

my ($infile,$p1,$p2,$outfile,$cutoff,$which_tech,$listfile,$sequencefile);
my ($ignorefile,@expressions);
my $ksize = 30 ;
my $verbose =0  ;
GetOptions(
            "sequencefile=s"=>\$sequencefile ,
            "infile=s"=>\$infile ,
            "p1=s"=>\$p1 ,
            "p2=s"=>\$p2 ,
            "listfile=s"=>\$listfile ,
            "ignorefile=s"=>\$ignorefile ,
            "which_tech=s"=>\$which_tech ,
            "outfile=s"=>\$outfile ,
            "expr=s"=>\@expressions,
            "ksize=i"=>\$ksize ,
            "savetable=i"=>\$savetable ,
            "verbose=i"=>\$verbose ,
            "debugone=i"=>\$debugone ,
            "cutoff=f"=>\$cutoff ,
           );
die "Dont recognize command line arg @ARGV " if(@ARGV);

my ($RESULTDIR,$PDBDIR,$FASTADIR,$APBSDIR,$FPOCKET,$SRC,$MATCH3D,$ANNDIR, $UNIPROT,$PREMONITION,$HELIXDIR,$DSSP,$CONFIGGRP,$BLASTOUT,$BLASTDB) = util_SetEnvVars();
usage( "Need to give a input file name => option -infile ") if(!defined $infile);
usage( "Need to give a input file name => option -sequencefile ") if(!defined $sequencefile);
my $info = {};
system ("mkdir -p $FASTADIR");

print STDERR "Info: Loggin in tmp2monitor\n";
my $ofhtmp = util_write("tmp2monitor");

my $STR = "";
my $NAME;

my ($GENOMETABLE,$NAMESgenome) = GNM_GenomeBreakHoppingOrSliding($savetable,$debugone,$ofhtmp,$infile,$ksize,1);
my $N2 = (keys %{$GENOMETABLE}) ;
print "There are $N2 in genome \n";
#my @l = (keys %{$GENOMETABLE});
#$, = " KKK \n" ;
#print @l ;
#die ;



my ($seqtable,$NAMESseq) = GNM_GenomeBreakHoppingOrSliding($savetable,$debugone,$ofhtmp,$sequencefile,$ksize,0,$GENOMETABLE);
my $N1 = (keys %{$seqtable}) ;
print "There are $N1 in seq \n";


my $done = {};
foreach my $s (keys %{$seqtable}){
	if(exists $GENOMETABLE->{$s}){
		my @l = @{$seqtable->{$s}};
		foreach my $l (@l){
			$done->{$l} = 1;
		}
	}
}
my ($common,$inAbutnotinB,$inBbutnotinA) = util_table_diff($NAMESseq,$done);
print "There are $common,$inAbutnotinB,$inBbutnotinA  not done \n";
foreach my $x (keys %{$NAMESseq}){
	if(! exists $done->{$x}){
		print $ofhtmp "$x notdone\n";
	}
}

## testing memory
#foreach my $k (keys %{$GENOMETABLE}){
   #delete $GENOMETABLE->{$k};
#}

$mu->record('after something_memory_intensive()');
$mu->dump();


####################



sub usage{
    my ($msg) = @_ ;
    print $msg , "\n" ; 
print << "ENDOFUSAGE" ; 
ENDOFUSAGE
    die ;
}
