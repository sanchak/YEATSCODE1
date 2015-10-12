#!/bin/csh -f

if($#argv != 2  ) then 
  echo "Usage : ./generaterep.csh  <impl dir> <file_having_list_of_designs> <tech - eg altera> <mode> <dirfortech - eg stratixii> "
  echo "You need to set ,  BENCH_HOME , BIN_HOME &  MGC_HOME "
  exit 
endif 

set TRS=$1
set SCAFFOLD=$2

setenv FASTADIR /home/sandeepc/DATA/rafael/walnut/2/NEW/FASTADIR_NT
#setenv BLASTDB /data2/export/tsbutter/NEW/NEW/DB
#set    GENOMEDB=//home/sandeepc/DATA/rafael/walnut/2/NEW/DB/wgs.5d.scafSeq200+.trimmed
set  SCAFFOLDDIR=/home/sandeepc/DATA/rafael/walnut/2/NEW/SCAFFOLDDIR

set BLASTDIR=BLASTOUT_SCAFFOLD

if( -z "$BLASTDIR/$TRS.$SCAFFOLD.blast") then
  \rm $BLASTDIR/$TRS.$SCAFFOLD.blast 
endif
if( -e "$BLASTDIR/$TRS.$SCAFFOLD.blast") then
 echo "done $BLASTDIR/$TRS.$SCAFFOLD.blast "
 exit 
endif 





set SCAFFOLDFASTA=$SCAFFOLDDIR/$SCAFFOLD.ALL.1.fasta
#unsetenv FORCERUN 
#if( ! -e $SCAFFOLDDIR/$SCAFFOLD.ALL.1.fasta) then 
    #extractsinglefastafromfile.pl -inf $BLASTDB/$GENOMEDB -wh $SCAFFOLD
#else
	#echo "$TRS and $SCAFFOLD seems to be processed. "
#endif



#\rm -f $SCAFFOLD.ALL.1.fasta
#\cp -f $SCAFFOLDFASTA $SCAFFOLD.ALL.1.fasta

mkdir -p $BLASTDIR
mkdir -p INFO

set REVFILE=$TRS.$SCAFFOLD.isrev
\rm -rf $REVFILE
myblastcompare2Fastafiles.csh $FASTADIR/$TRS.ALL.1.fasta $SCAFFOLDFASTA $BLASTDIR/$TRS.$SCAFFOLD.blast N



#parseBlastLatestpairwiseOneScaffold.pl -outf INFO/$TRS.$SCAFFOLD.tmp -inf $BLASTDIR/$TRS.$SCAFFOLD.blast -trs $TRS  -justchecking 1 -scaff $SCAFFOLDFASTA -sname $SCAFFOLD -blastdir $BLASTDIR
checkBLASTtrs2scaffolddirection.pl -inf $BLASTDIR/$TRS.$SCAFFOLD.blast -outf $REVFILE 




if( -e $REVFILE) then 
	echo "======== Reverese ==========="
	if(-e "$SCAFFOLDFASTA.comp.fasta") then 
	else
	      echo "doing comp"
              writeCompFasta.pl -in $SCAFFOLDFASTA 
	endif

        myblastcompare2Fastafiles.csh $FASTADIR/$TRS.ALL.1.fasta $SCAFFOLDFASTA.comp.fasta  $BLASTDIR/$TRS.$SCAFFOLD.blast N
	#unlink $SCAFFOLD.ALL.1.fasta
else
	echo "======== SAME ==========="
	#\rm -f TMPDIR/$TRS.$SCAFFOLD.fasta
	#ln -s $SCAFFOLDFASTA TMPDIR/$TRS.$SCAFFOLD.fasta
endif



\rm -rf $REVFILE
