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
my ($duplicate,$blastout,$infile,$p1,$p2,$outfile,$trs,$cutoff,$orfdir,$listfile,$protein);
my ($chooseone,$ignorefile,$donone,@expressions);
my $howmany = 100000 ;
my $verbose = 0 ;
GetOptions(
            "orfdir=s"=>\$orfdir ,
            "trs=s"=>\$trs ,
            "blastout=s"=>\$blastout ,
            "infile=s"=>\$infile ,
            "p1=s"=>\$p1 ,
            "p2=s"=>\$p2 ,
            "listfile=s"=>\$listfile ,
            "ignorefile=s"=>\$ignorefile ,
            "outfile=s"=>\$outfile ,
            "donone"=>\$donone ,
            "duplicate"=>\$duplicate ,
            "expr=s"=>\@expressions,
            "howmany=i"=>\$howmany ,
            "verbose=i"=>\$verbose ,
            "cutoff=f"=>\$cutoff ,
           );
die "Dont recognize command line arg @ARGV " if(@ARGV);
usage( "Need to give a output file name => option -outfile ") if(!defined $outfile);
my $ofh = util_write($outfile);

usage( "Need to give a trs -option -trs  ") if(!defined $trs);
usage( "Need to give a orfdir -option -orfdir  ") if(!defined $orfdir);
usage( "Need to give a blastout -option -blastout  ") if(!defined $blastout);
usage( "Need to give a cutoff -option -cutoff  ") if(!defined $cutoff);


my $EVALFORCOMMON = 1E-15 ;
my $BLASTEXT = ".blast.nt";

system("mkdir -p ANN/ERR");
system("mkdir -p ANN/WARN");
system("mkdir -p ANN/ANNOTATE");
system("mkdir -p ANN/COMMANDS");
system("mkdir -p FASTADIR_ORFBEST");


my $ofhwarn = util_write("ANN/WARN/$trs.w");
my $ofhannotate = util_write("ANN/ANNOTATE/$trs.a");
my $ofhcommands = util_write("ANN/COMMANDS/$trs.c");


my $trsfasta = "$trs.ALL.1.fasta";

my @l = <$blastout/$trs.*>;


$, = "\n";
print "Files are @l \n" if($verbose);

my $warnstr = "";
if(@l > 3){
	$warnstr = "repeats?";
	print $ofhwarn "$trs repeats?\n";
	print "$trs repeats?\n" if($verbose);
	shift @l ;
}

if(@l < 3){
        my $ofhERR = util_write("ANN/ERR/$trs.e");
	print $ofhERR "Not all found for $trs\n";
	exit ;
}
foreach my $infile (@l){
	if(-z $infile){
            my $ofhERR = util_write("ANN/ERR/$trs.e");
	    print $ofhERR "Not all found $trs\n";
	    exit ;
	}
}


my $sort = {};
my $sortEVALS2SCORES = {};
my $sortEVALS2ALLSTRS = {};
my $sortEVALS2INFILE = {};
my $MAXSCORE = 0 ;
my $ALLSTRS = "";

my $sortTRSLen = {};
foreach my $infile (@l){
   my ($fulllength,$TRSnm,$STRS,$LLLs,$SCORES,$IDENTITIES,$EVALUES) = util_ParseWebBlast($infile);

	my $orfile = $infile;
	$orfile =~ s/$BLASTEXT/.ALL.1.fasta/;
	$orfile =~ s/$blastout.//;
	my $OOO = "$orfdir/$orfile";
	if( ! -e $OOO){
	  # ignore missing TRS file if there are more than 3 - repeats
	  next if($warnstr ne "");

          my $ofhERR = util_write("ANN/ERR/$trs.e");
   	  print $ofhERR "Did not find $OOO. \n";
	  exit ;

	}
	my ($str,$firstline) = util_readfasta($OOO);
	my $len = length($str);
	$sortTRSLen->{$orfile} = $len;

   if(!defined $fulllength){
          my $ofhERR = util_write("ANN/ERR/$trs.e");
   	  print $ofhERR "File not complete, rerun $trs\n";
   	  exit ;
   }

   if(!defined $TRSnm){
      next ;
   }

   my @TRSnm = @{$TRSnm};
   my @STRS = @{$STRS};
   my @LLLs = @{$LLLs};
   my @SCORES = @{$SCORES};
   my @EVALUES = @{$EVALUES};
   my @IDENTITIES = @{$IDENTITIES};

   $infile =~ s/^$blastout.//;
   $ALLSTRS = $ALLSTRS . "\n" . "\t$infile ". $STRS[0];



   my $EEE = $EVALUES[0] ;
   while(exists $sort->{$EEE}){
	  print "incrementing $EEE $trs \n"  if($verbose);
   	  $EEE  = $EEE + $EEE/10000 + 1E-100 ;
   }

   $sort->{$EEE} = $STRS ;
   $sortEVALS2SCORES->{$EEE} = $SCORES[0] ;
   $sortEVALS2ALLSTRS->{$EEE} = $STRS ;
   $sortEVALS2INFILE->{$EEE} = $infile ;
   if($MAXSCORE < $SCORES[0]){
   	  $MAXSCORE = $SCORES[0];
   }
   
}

print "Max score = $MAXSCORE   \n" if($verbose);

my $UNCHARSTRINGS = "chromosome|hypothe|unnamed|uncharacterized|Predicted protein";

my $cnt = 0 ;
my $trsrealname  ;
my $firsteval = {};

foreach my $k (sort { $a <=> $b} keys %{$sort}){

	my $SCORE = $sortEVALS2SCORES->{$k};
	if(!defined $SCORE){
        my $ofhERR = util_write("ANN/ERR/$trs.e");
	    print $ofhERR "$k $trs has not SCORE " ;
		die ;
	}
	my $infile = $sortEVALS2INFILE->{$k};
	$infile =~ s/$BLASTEXT/.ALL.1.fasta/;
	$infile =~ s/\.blast/.ALL.1.fasta/;
	$infile =~ s/\._/.ORF_/;
	


	my $len = $sortTRSLen->{$infile};
	print "$infile blastscore=$SCORE len=$len eval=$k\n" if($verbose);


	if(($k < $cutoff && $SCORE > $MAXSCORE/2)){
		## only for the first 
		if(!defined $trsrealname){
	        $trsrealname = $infile ;
			$firsteval = $k ;
		}

		#if you have one, be more stringent for the next match
		if($cnt){
			if($k > $EVALFORCOMMON){
				print "$infile Ignoring $k since we aleady have one match with $firsteval\n";
				next ;
			}
		}


		$cnt++;
	    my @l = @{$sort->{$k}};

		## try to get an charaterization - only if the evalue does not fall off
		foreach my $i (@l){
			my @l = split " ", $i ;
			my $EEE = getEvalue($i);
			next if($EEE > $cutoff);

			if(!($i =~ /($UNCHARSTRINGS)/i)){
				last ;


             }
		}

	}
}

print "cnt is $cnt \n" if($verbose);
## not annotated - even at the low level
if($cnt eq 0){
		## choose the longest ....
		my @sorted = sort {$sortTRSLen->{$b} <=> $sortTRSLen->{$a}} (keys %{$sortTRSLen});

		my $trsrealname = $sorted[0];
		my $len = $sortTRSLen->{$trsrealname};
	    print $ofhcommands "\\cp -f  $orfdir/$trsrealname FASTADIR_ORFBEST/$trsfasta # none \n";
	    print $ofhannotate "$trs #none none \n";
	    print "None - choosing $trsrealname with len = $len\n" if($verbose);
            exit ;
}

my @sorted = (sort {$a <=> $b} keys %{$sortEVALS2INFILE});
if(!@sorted){
        my $ofhERR = util_write("ANN/ERR/$trs.e");
	print $ofhERR "Expect a value here \n";
	die ;
}

## Evalue of best match
my $X = shift @sorted;
if($cnt eq 1){		

	my $SSS = "uniqhigh";
	if($X > 1E-10){
		$SSS = "uniqlow";
	}
	my @SCORESA = @{$sortEVALS2ALLSTRS->{$X}};
	my $STRA  = $SCORESA[0];

	print $ofhannotate "$trs #$SSS $STRA \n";
	print $ofhcommands "\\cp -f  $orfdir/$trsrealname FASTADIR_ORFBEST/$trsfasta # uniq \n";
	print "$SSS $trsrealname\n" if($verbose);

}
elsif($cnt > 1 ){

		my $Y = shift @sorted;

		my $Z ;

		my @SCORESA = @{$sortEVALS2ALLSTRS->{$X}};
		my @SCORESB = @{$sortEVALS2ALLSTRS->{$Y}};
		my @SCORESC ;

		my $STRA  = $SCORESA[0];
		my $STRB  = $SCORESB[0];
		my $STRC  ;


		my $Atable = getTableofUniprotIds(\@SCORESA,$EVALFORCOMMON);
		my $Btable = getTableofUniprotIds(\@SCORESB,$EVALFORCOMMON);
		my $Ctable = {};
		if($cnt eq 3){
		    $Z = shift @sorted;
		    @SCORESC = @{$sortEVALS2ALLSTRS->{$Z}};
		    $STRC = $SCORESC[0];
		    $Ctable = getTableofUniprotIds(\@SCORESC,$EVALFORCOMMON);
		}

		my $foundcommon = 0 ;
		foreach my $uniprotID (keys %{$Atable}){
			if(exists $Btable->{$uniprotID} || exists $Ctable->{$uniprotID} ){
				$foundcommon = 1 ;
				print "Common (error) $l[1]  $uniprotID \n" if($verbose);
	            print $ofhannotate "$trs #uniqfix $STRA  \n";
	            print $ofhcommands "\\cp -f  $orfdir/$trsrealname FASTADIR_ORFBEST/$trsfasta #uniqfix \n";
				last ;
			}
		}


	    my $trsorig = $trs.  ".ALL.1.fasta";
		if(!$foundcommon){
	        print "Duplicate $X $Y \n" if($verbose);
		    my $a = $sortEVALS2INFILE->{$X};
		    my $b = $sortEVALS2INFILE->{$Y};
	        my $trsaonlyname = $trs. "_A" ;
	        my $trsbonlyname = $trs. "_B" ;
	        my $trsa = $trsaonlyname.  ".ALL.1.fasta";
	        my $trsb = $trsbonlyname.  ".ALL.1.fasta";
            #$a =~ s/\./.ORF/;
            #$b =~ s/\./.ORF/;

            $a =~ s/$BLASTEXT/\.ALL.1.fasta/;
            $b =~ s/$BLASTEXT/\.ALL.1.fasta/;

		    print $ofhcommands "\\cp -f  $orfdir/$a FASTADIR_ORFBEST/$trsa  \n";
		    print $ofhcommands "\\cp -f  $orfdir/$b FASTADIR_ORFBEST/$trsb  \n";
		    print $ofhcommands "\\cp -f  FASTADIR_NT/$trsorig FASTADIR_NT/$trsa # duplicate 2 $a \n";
		    print $ofhcommands "\\cp -f  FASTADIR_NT/$trsorig FASTADIR_NT/$trsb # duplicate 2 $b \n";


	        print $ofhannotate "$trsaonlyname #dupli $STRA  \n";
	        print $ofhannotate "$trsbonlyname #dupli $STRB  \n";

		    if($cnt eq 2){
		    }
		    else{
			    die if($cnt ne 3);
		        my $c = $sortEVALS2INFILE->{$Z};
                $c =~ s/$BLASTEXT/\.ALL.1.fasta/;
	            my $trsconlyname = $trs . "_C" ;
	            my $trsc = $trsconlyname . ".ALL.1.fasta";

		        print $ofhcommands "\\cp -f  FASTADIR_NT/$trsorig FASTADIR_NT/$trsc # duplicate 3 $c\n";
		        print $ofhcommands "\\cp -f  $orfdir/$c FASTADIR_ORFBEST/$trsc # duplicate 3 $c\n";
	            print $ofhannotate "$trsconlyname #dupli $STRC  \n";
		    }
		    print $ofhcommands "\n";
		}
	}

sub getUniprotID{
	my ($str) = @_ ;
	#$str =~ s/\|/ /g;
	my @l = split " ", $str ;
	return $l[0];
}
sub getEvalue{
	my ($str) = @_ ;
	$str =~ s/\|/ /g;
	my @l = split " ", $str ;
	my $N = @l - 1;
	my $EEE = $l[$N];
	return $EEE ;
}
   
sub getTableofUniprotIds{
		    ## make table ofA loAwer uniprot IDs 
		    my ($l,$CUTOFF) = @_ ;
		    my @SCORESB = @{$l};
		    my $Atable = {};
		    foreach my $str (@SCORESB){
			    #$str =~ s/\|/ /g;
			    #my @l = split " ", $str ;
			    my $uniprotID = getUniprotID($str);
			    my $EEE = getEvalue($str);
                next if($EEE > $CUTOFF);


			    $Atable->{$uniprotID} = $EEE ;
		    }
			return $Atable ;
}

sub usage{
    my ($msg) = @_ ;
    print $msg , "\n" ; 
print << "ENDOFUSAGE" ; 
ENDOFUSAGE
    die ;
}
