#!/bin/csh -f

if($#argv != 3  ) then 
  echo "Usage : ./generaterep.csh  <impl dir> <file_having_list_of_designs> <tech - eg altera> <mode> <dirfortech - eg stratixii> "
  echo "You need to set ,  BENCH_HOME , BIN_HOME &  MGC_HOME "
  exit 
endif 
#echo "PERLLIB = $PERLLIB ,  $BENCH_HOME = BENCH_HOME , BIN_HOME = $BIN_HOME , MGC_HOME = $MGC_HOME "

set PWD = ` pwd`
set run = $1
set base = $2
set resolution = $3

cd $PWD/$run
ls *.png | nawk -F'-' '{print $NF, $0}' | sort -n | cut -d ' ' -f2- > ! list
#echo Sharpening
#$SRC/SHELLSCRIPTS/sharpen.csh list

cd $PWD/$base
ls *.png > ! list 


set resultsdir = $PWD/RESULTS.$run
mkdir -p $resultsdir

set testfile = "out.peri"

cd $resultsdir
foreach i (`cat $PWD/$base/list`)
	mkdir -p $i 
	cd $i 

	\cp -f $PWD/$base/$i base.png
	ln -s $PWD/$run/* .  > & ! /dev/null
	mkdir -p found

	#if(! -e $testfile || -z $testfile) then 
	if(! -e $testfile ) then 
		echo Processing $i 
	    time $SRC/IMAGE/3DprocessAmba.pl -in base.png -lis list -spec red  -color white -dist 1.5  -resolution $resolution -justide red > & ! log 
	else
		echo $i has already been processed
	endif

	cd ../
end 

cd $PWD
