#!/bin/csh -f

if($#argv != 1  ) then 
  echo "Usage : ./generaterep.csh  <impl dir> <file_having_list_of_designs> <tech - eg altera> <mode> <dirfortech - eg stratixii> "
  echo "You need to set ,  BENCH_HOME , BIN_HOME &  MGC_HOME "
  exit 
endif 
#echo "PERLLIB = $PERLLIB ,  $BENCH_HOME = BENCH_HOME , BIN_HOME = $BIN_HOME , MGC_HOME = $MGC_HOME "

set list = $1

mkdir -p ORF
mkdir -p ORFTMP

newfile.csh ORFCommands.csh 
foreach i (`cat $list`)
     echo orfprocessSingle.csh $i >>  ORFCommands.csh 
end 

scheduleprocessInsertingsleep.pl -int 100 -sleep 2 -inf ORFCommands.csh 
runInBack.csh Sch.ORFCommands.csh 
