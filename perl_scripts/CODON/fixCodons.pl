#!/usr/bin/perl -w 
use strict ;
use FileHandle ;
use Getopt::Long;
use Cwd ;
use MyUtils;
use MyGeom;
use PDB;
use ConfigPDB;

use Time::HiRes qw( usleep ualarm gettimeofday tv_interval clock_gettime clock_getres  clock);
use POSIX qw(floor);
my $commandline = util_get_cmdline("",\@ARGV) ;
my ($infile,$outfile,$which_tech,$listfile,$protein);
my (@expressions,$config);
my $aminoacid = 1;
my $verbose = 1 ;
GetOptions(
            "which_tech=s"=>\$which_tech ,
            "protein=s"=>\$protein ,
            "config=s"=>\$config ,
            "infile=s"=>\$infile ,
            "listfile=s"=>\$listfile ,
            "outfile=s"=>\$outfile ,
            "expr=s"=>\@expressions,
            "aminoacid=i"=>\$aminoacid ,
           );
die "Dont recognize command line arg @ARGV " if(@ARGV);
usage( "Need to give a output file name => option -outfile ") if(!defined $outfile);
my $ofh = util_write($outfile);
usage( "Need to give a input file name => option -infile ") if(!defined $infile);

my ($RESULTDIR,$PDBDIR,$FASTADIR,$APBSDIR,$FPOCKET,$SRC) = util_SetEnvVars();


my ($string,$firstline) = util_readfasta($infile);


my $STARTIDX=5;
my $ENDIDX=10;
my $CNT = CountRepeatsOfBases($string,$STARTIDX,$ENDIDX);
while($CNT){
    $string = FixCodonOnePass($string);
    $CNT = CountRepeatsOfBases($string,$STARTIDX,$ENDIDX);
    print "There were $CNT base repeats\n";
}

print $ofh "$firstline";
print $ofh "$string \n";

sub FixCodonOnePass{
	my ($STR) = @_ ;
		 my @NT = ($STR =~ /(...)/g);
		 my $codontable = util_getCodonTable();
		 my $NNN = @NT -1; 
		 my $SSS = "";
		 my $FIXED = "";

		 foreach my $idx (0..$NNN){
		 	    my $prev = "";
		 	    my $this = "";
		 	    my $next = "";

		 	    $prev = $NT[$idx-1] if($idx ne 0);
		 	    $this = $NT[$idx];
		 	    $next = $NT[$idx+1] if($idx ne $NNN);

				my $AAA = $prev. $this. $next ;

				my $s1 = $prev . $this ;
                $CNT = CountRepeatsOfBases($s1,$STARTIDX,$ENDIDX);
				if($CNT || $s1 =~ /CGCG/){
				    print " ALL = $AAA \n";
					print "fddd $s1 \n";
					my $replace = ReplaceCodonWithAnyOther($this);
					print "Relace prev $this with $replace at index $idx\n";
					$this =$replace ;
					$NT[$idx] = $replace ;
				}
				

				my $s2 = $this . $next ;
                $CNT = CountRepeatsOfBases($s2,$STARTIDX,$ENDIDX);
				if($CNT || $s2 =~ /CGCG/){
				    print " ALL = $AAA \n";
					print "fddd $s2 \n";
					my $replace = ReplaceCodonWithAnyOther($this);
					print "Relace next $this with $replace at index $idx\n";
					$this =$replace ;
					$NT[$idx] = $replace ;
				}

				$FIXED = $FIXED . $this ;

		 	    #my $a = $codontable->{$i} ;
			    #$SSS = $SSS . $a ;
		 }
		 return $FIXED ;

}

sub ReplaceCodonWithAnyOther{
	my ($codon)=@_ ;
		 my $codontable = util_getCodonTable();
	my $a = $codontable->{$codon} ;

	my @codonsforaa = util_getCodonForAA($a);
	my $CC ;
	while(1){
	   $CC = util_pick_random_from_list(\@codonsforaa);
	   if($CC ne $codon){
	   	   last ;
	   }
	}
	return $CC ;
}


sub util_getCodonForAA{
	my ($aa) = @_ ;
	my @ret ;
		 my $codontable = util_getCodonTable();
	foreach my $i (keys %{$codontable}){
		my $AA = $codontable->{$i} ;
		if($AA eq $aa){
			push @ret, $i ;
		}
	}
	die if(!@ret);
	return @ret ;

}

sub CountRepeatsOfBases{
  my ($STR,$start,$end) = @_ ;
  my @bases = qw(A T G C);
  my $cnt = 0 ;
  foreach my $i ($start..$end){
	 my $str = "";
	 foreach my $base (@bases){
		 ## create query string
	     foreach my $x (1..$i){
		     $str = $str . $base ;
	     }
	     if($STR =~ /$str/){
	         #print "len=$i str=$str matched in CountRepeatsOfBases \n";
			 $cnt++;
	     }
	 }
  }
  
  if($STR =~ /CGCG/){
  	$cnt++;
  }
  return $cnt ;
}



sub usage{
    my ($msg) = @_ ;
	    print $msg , "\n" ;
print << "ENDOFUSAGE" ; 
ENDOFUSAGE
		    die ;
			}
