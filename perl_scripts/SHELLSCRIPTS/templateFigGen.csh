#!/bin/csh  -f

if($#argv != 2  ) then 
    echo "<name> <pdf>"
endif 

set nm=$1 
set file=$2

cp fig.sample.tex fig.$nm.tex 
replacestring.pl -in fig.$nm.tex -wh figlabel -with $nm -same -outf kkkkk
replacestring.pl -in fig.$nm.tex -wh microbio -with $nm -same -outf kkkkk

echo "\input{fig.$nm.tex}" >> fig.tex


\cp -f $2 $nm.pdf 

 

