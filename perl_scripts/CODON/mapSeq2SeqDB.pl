#!/usr/bin/perl -w 
use strict ;
use FileHandle ;
use Getopt::Long;
use Cwd ;
use MyUtils;
use POSIX qw(floor);
my $commandline = util_get_cmdline("",\@ARGV) ;
my ($dbfile,$all,$infile,$outfile,$randomize,$or,$silent,$groupinfo);
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
            "randomize"=>\$randomize ,
            "listfile=s"=>\$listfile ,
            "ignorefile=s"=>\$ignorefile ,
            "howmany=i"=>\$howmany ,
            "size=i"=>\$size ,
            "verbose=i"=>\$verbose ,
            "or=i"=>\$or ,
            "cutofflength=i"=>\$cutofflength ,
            "type=s"=>\@types,
            "protein=s"=>\$protein,
            "motif=s"=>\@motifs,
            "outfile=s"=>\$outfile 
           );
die "Dont recognize command line arg @ARGV " if(@ARGV);
usage( "Need to give a input file name => option -infile ") if(!defined $infile);
usage( "Need to give a input file name => option -dbfile ") if(!defined $dbfile);
usage( "Need to give a output file name => option -outfile ") if(!defined $outfile);
usage( "Need to give a protein -option -cutofflength  ") if(!defined $cutofflength);

my $ofh = util_write($outfile);
my $ofhlog = util_write("log");

print STDERR "Info: parsing file $infile - might take some time\n";
my ($info,$infoSeq2PDB,$mapChainedName2Name) = util_parsePDBSEQRESNEW($infile,0);
my ($info1,$infoSeq2PDB1,$mapChainedName2Name1) = util_parsePDBSEQRESNEW($dbfile,0);

die if(keys %{$infoSeq2PDB} ne 1);

## get the protein name and sequence...
my ($protseq) = (keys %{$infoSeq2PDB});
my @QQ = @{$infoSeq2PDB->{$protseq}};
my @RR = split " ", $QQ[0];
my $protname = $RR[0];


## make the full DB into one string - warning...
my $fulldb = "";

my @keysofPDB1 ;
if(defined $randomize){
    @keysofPDB1 = keys %{$infoSeq2PDB1} ;
}
else{
    @keysofPDB1 = sort keys %{$infoSeq2PDB1} ;
}

foreach my $db (@keysofPDB1){
	$fulldb = $fulldb . "B" . $db;
}


### start from the whole protein...
if(!defined $size){
    $size = length($protseq);
}

## track which indices are used...
my $doneindices = {};

print "Iterating downwards from size $size, log in log\n";
my $CNTFRAG = 0 ;
while($size){
	print "Processing for $size\n" if($verbose);

    my @strings = util_splitstring($protseq,1,$size);

	my $CNT = 0 ;
    foreach my $string (@strings){
		$CNT++;

		next if(exists $doneindices->{$CNT});

		## do individual only if there is match in the whole
	    $_ = $fulldb ;
	    if(/$string/){
           foreach my $db (@keysofPDB1){
		       if($db =~ /$string/){
		           next if(exists $doneindices->{$CNT});


				   ## find the length before the match
				   $_ = $db ;
				   my ($strprev) = (/^(.*)$string/);
				   my @PP = split "", $strprev ;
				   my $PREVLEN = @PP ;




				   # get the name - and split it into TRS and ORF
				   my $W = $infoSeq2PDB1->{$db};
				   my @l = @{$W};
				   $_ = $l[0];

				   my ($trs,$orf) = /(.*)\.(ORF.*)/;
				   my $orfile = "ORF/$trs.orf" ;
				   my $fastafile = "FASTADIR/$trs.ALL.1.fasta";
				   my $info ;
				   my ($NT,$AA) = util_ProcessOneORF($trs,$orfile,$fastafile,$orf,$info);
				   my @NT = @{$NT} ;
				   my @AA = @{$AA} ;

			       print $ofhlog "found $string for $size and $CNT and ID @l PREVLEN = $PREVLEN $trs $orf \n";
				   $CNTFRAG++;
				   ## extract slice of size from the ORF 
				   my $TTT = "";
				   my @codons ;
				   {
				      my $END = $PREVLEN+$size -1;

				      foreach my $idx ($PREVLEN..$END){
					  	$TTT = $TTT. $AA[$idx];
						push @codons, $NT[$idx];
				   	  
				      }
				   }
				   die "$TTT and $string ne"  if($TTT ne $string);




				   my $END = $CNT+$size -1;

				   ## map the indices in the protein to the corresponding codon
				   foreach my $idx ($CNT..$END){
				   		my $codon = shift @codons ;
				   	    $doneindices->{$idx} = $codon;
				   }
				   die if(@codons);
		       }
	       }
        }
	}
	## try one lenght lesser
	$size--;
}

my $LEN  = length($protseq);
print "Mapped $LEN into $CNTFRAG fragments\n";

my $TRANSLATED = "";;
foreach my $x (1..$LEN){
	if(! exists $doneindices->{$x}){
		print "Did not find $x...\n";
	}
	else{
		$TRANSLATED = $TRANSLATED . $doneindices->{$x} ;
		#print $ofh "$doneindices->{$x}\n";
	}
}
print $ofh ">${protname}TRANS\n";
print $ofh "$TRANSLATED \n";


my $codoninfo = {};
my $SEQ = util_GetCodonBiasFromNucleotide($codoninfo,$TRANSLATED);


my $OFH = util_write("codonbias");
util_EmitCodonBiasfromInfo($OFH,$codoninfo);

sub usage{
    my ($msg) = @_ ;
    print $msg , "\n" ; 
print << "ENDOFUSAGE" ; 
ENDOFUSAGE
    die ;
}
