#!/bin/csh -f
mkdir -p ORF
mkdir -p ORFTMP
set i=$1 

if(! -e FASTADIR_ORF/$i.list) then 
   getorf $FASTADIR/$i.ALL.1.fasta ORFTMP/$i.orf 
   fixORFnames.pl -inf ORFTMP/$i.orf -out ORF/$i.orf
   findrepeatedorfs.pl -trs $i -orfdir ORF/ -write 4 
endif 

