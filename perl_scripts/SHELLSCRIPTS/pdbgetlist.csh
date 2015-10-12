#!/bin/csh -f

if(! -e $1) then
	set list = /tmp/tttttt
	echo $1 > ! $list
else
    set list = $PWD/$1 
endif 

### First move into PDBDIR - so set it properly
cd $PDBDIR 

newfile.csh $list.new
foreach i (`cat $list`)
   if(! -e ${i}A.pdb) then 
   	  if(! -e $i.pdb) then
         wget    http://www.rcsb.org/pdb/files/$i.pdb 
	  endif 
      splitpdbIntChains.pl -p $i  -outf $list.new
   endif 
   #echo ls ${i}?.pdb
   ls ${i}?.pdb >> $list.new 
end 


$SRC/SHELLSCRIPTS/removePDBfromname.csh $list.new
UNIQ $list.new

getPDBModel1ChainAlist.csh $list.new

\mv -f $list.new $PWD
if(! -e $1) then
	\rm $list
endif 


#foreach i (`cat list.new`)
   #\mv -f $i.pdb $PDBDIR
#end 
#apbs.csh list.new
 

