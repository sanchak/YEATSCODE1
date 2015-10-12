#!/bin/csh  -f

if($#argv != 4  ) then 
    echo "Wrong args, required 1 " ; exit 
endif 

set list = $1 
set DBfile = $2 
set FASTADIR = $3 

if(! -e $DBfile) then 
    newfile.csh $DBfile  
    foreach i (`cat $list`)
        cat $FASTADIR/$i.ALL.1.fasta >> $DBfile
    end 
    \rm -f DBFILE
    ln -s $DBfile DBFILE
	makeblastdbP.csh DBFILE
endif 


\rm -f FASTADIRSUBSET
ln -s $FASTADIR FASTADIRSUBSET 
mkdir -p BLASTOUTSUBSET

setupCommandsFromListAndSchedule.pl -lis $list -out blast.exec.csh -sleep 2 -inter 100 -string 'blastp -db DBFILE -query FASTADIRSUBSET/$i.ALL.1.fasta -out BLASTOUTSUBSET//$i.blast.nt '
 

setupCommandsFromListAndSchedule.pl -lis $list -out group.exec.csh -sleep 2 -inter 50 -string 'parseBlastLatest.pl -outf group.nr -inf BLASTOUTSUBSET/$i.blast.nt -tr $i'
