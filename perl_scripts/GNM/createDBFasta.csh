#!/bin/csh 

if($#argv != 3  ) then 
  echo "Usage : ./generaterep.csh  <impl dir> <file_having_list_of_designs> <tech - eg altera> <mode> <dirfortech - eg stratixii> "
  echo "You need to set ,  BENCH_HOME , BIN_HOME &  MGC_HOME "
  exit 
endif 
#echo "PERLLIB = $PERLLIB ,  $BENCH_HOME = BENCH_HOME , BIN_HOME = $BIN_HOME , MGC_HOME = $MGC_HOME "


set list = $1 
set FASTADIR = $2 
set out = $3 

 
 if(! -e $out) then 
newfile.csh $out 
foreach i (`cat $list`)
	cat $FASTADIR/$i.ALL.1.fasta >> $out
end 
else
	echo "Info: File $out exists"
endif

