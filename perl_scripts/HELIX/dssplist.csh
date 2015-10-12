#!/bin/csh -f

if($#argv != 3  ) then 
  echo "Usage : ./generaterep.csh  <impl dir> <file_having_list_of_designs> <tech - eg altera> <mode> <dirfortech - eg stratixii> "
  echo "You need to set ,  BENCH_HOME , BIN_HOME &  MGC_HOME "
  exit 
endif 
#echo "PERLLIB = $PERLLIB ,  $BENCH_HOME = BENCH_HOME , BIN_HOME = $BIN_HOME , MGC_HOME = $MGC_HOME "

set list = $1 
set outfile = $2 
set writeind = $3 
touch $outfile

echo DSSP=$DSSP and HELIXDIR=$HELIXDIR
touch HTH
foreach i (`cat $list`)
	if(! -e "$DSSP/$i.dssp") then
	    echo mkdssp -i $PDBDIR/$i.pdb -o $DSSP/$i.dssp
	    mkdssp -i $PDBDIR/$i.pdb -o $DSSP/$i.dssp
	else
		echo Info: DSSP $DSSP/$i.dssp file exists
	endif
	grep -w $i $outfile > & ! /dev/null 
	if  ($status != 0 ) then 
	   helixparseDSSPoutput.pl -outf $outfile -p $i -dssp $DSSP/$i.dssp -what BETA -writeind $writeind
	   helixparseDSSPoutput.pl -outf $outfile -p $i -dssp $DSSP/$i.dssp -what HELIX -writeind $writeind
	else
		echo $i already processed in output file $outfile
	endif 

end 

 

