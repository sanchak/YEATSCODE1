package MyGNM;
use Carp ;
use POSIX ;
require Exporter;
use Algorithm::Combinatorics qw(combinations) ;
use Math::NumberCruncher;
use Math::MatrixReal;  # Not required for pure vector math as above
use Math::Geometry ; 
use Math::VectorReal qw(:all);  # Include O X Y Z axis constant vectors
#use Math::Trig;
#use Math::Trig ':radial';
no warnings 'redefine';
my $EPSILON = 0.01;
use MyPymol;
use MyConfigs;
use MyUtils;

my $verbose = 0 ;
no warnings 'recursion'; 

  local $SIG{__WARN__} = sub {};

@ISA = qw(Exporter);
@EXPORT = qw( 
GNM_PARSEBLAST
GNM_DirectionBlast
GNM_MakeGroupOfNR_OrAnnotate
GNM_MapOneTRS2Scaffold
GNM_parseCountSummed
GNM_parseANNFile
GNM_KMERIZE_sliding
GNM_KMERIZE_hopping
GNM_GenomeBreakHoppingOrSliding
        );

use strict ;
use FileHandle ;
use Getopt::Long;



## this functions writes in a group file format - (:$trs $name $blastscore}  (forWGS = 0)
## or write the annotation 

sub GNM_MakeGroupOfNR_OrAnnotate{
    my ($info,$querylength,$Subjectname,$isNT,$findcharstring,$forWGS,$trs,$infile,$ofh,$ofhanno,$blastcutoff,$MAPSscafftoTRS,$checkforsame) = @_ ;
    #my ($info,$querylength,$Subjectname) = GNM_PARSEBLAST($infile);
    #die "$trs - querylength not defined in $infile" if(!defined $querylength);
    
    my $sorted = {};
    my $DONEBSCORE = {};


	#print STDERR "Sorting based on BLAST SCORE\n";
    foreach my $k (@{$info}){
        my $org = $k ;
        my ($name,$description,$blastscore,$subjectlength,$iden_percent,$subjectmatched,$querystart,$queryend,$subjectstart,$subjectend,$expect) = split " ",$k ;
        while(exists $DONEBSCORE->{$blastscore}){
            $blastscore = $blastscore - 0.1;
        }
        $DONEBSCORE->{$blastscore} = 1 ;
        $sorted->{$k} = $blastscore ;
    }
     my @sorted = sort { $sorted->{$b} <=> $sorted->{$a} } (keys %{$sorted});
     my $NNN = @sorted ;
    
    my $done = {};
    my @others ;
    my $found100percent = 0;
    foreach my $k (@sorted){
        my $org = $k ;
        my ($name,$description,$blastscore,$subjectlength,$iden_percent,$subjectmatched,$querystart,$queryend,$subjectstart,$subjectend,$expect) = split " ",$k ;
        print "$name,$description,$blastscore,$subjectlength,$iden_percent,$subjectmatched,$querystart,$queryend,$subjectstart,$subjectend,$expect\n" if($verbose);
        $name =~ s/;//g;
        if($name eq $trs && $iden_percent eq 100){
            $found100percent  = 1 ;
            $querylength = $subjectlength;
            $done->{$name} = 1 ;
        }
        else{
            push @others, $k ;
        }

        #while(exists $DONEBSCORE->{$blastscore}){
            #$blastscore = $blastscore - 0.1;
        #}
        #$DONEBSCORE->{$blastscore} = 1 ;
        #$sorted->{$k} = $blastscore ;
    }
    die if (!$found100percent && defined $checkforsame);
    die if (!defined $querylength && defined $checkforsame);

    

     if(!@sorted){
          print $ofhanno "$trs Warning:Did not find any characterized\n";
     }



	### Annotate ###
    my $UNCHARSTRINGS = Config_GetUncharStrings();

    my $first ; 
    my $percent ; 
    foreach my $k (@sorted){
        my $org = $k ;
        my ($name,$description,$blastscore,$subjectlength,$iden_percent,$subjectmatched,$querystart,$queryend,$subjectstart,$subjectend,$expect) = split " ",$k ;
    
        ## choose the first percent
        $description =~ s/ZZZ/ /g;
        my $str2print ;
        if(!defined $first){
    	    if($blastscore < $blastcutoff){
               print $ofhanno "$trs Warning:Did not find any characterized\n";
			   last ;
	        }

			## note when matching to scaffolds this percent might be very low
            $percent = int(100 * ($querylength/$subjectlength));
            $str2print = "$trs\t$name\t$description $percent $blastscore $expect\n";
            $first = $str2print ;
        }
    
        $str2print = "$trs\t$name\t$description isNT=$isNT $percent $blastscore $expect\n";
    
        next if($findcharstring && $description =~ /$UNCHARSTRINGS/);

		if($blastscore < $blastcutoff){
           print $ofhanno "$first";
        }
		else{
             print $ofhanno "$str2print";
		}
        last ;
    }
        
	########## Make groups ############
    my $grouping = {};
    my $foundone = 0 ;
    my $cntfound = 0 ;
    if($forWGS){
       print $ofh "$trs ";    
     }
	my $ofhCommands = util_open_or_append("pairwiseTRS2scaffold.csh");
	my $last ; 

    my $DONENAMES = {};
    foreach my $k (@sorted){
        my $org = $k ;
        my ($name,$description,$blastscore,$subjectlength,$iden_percent,$subjectmatched,$querystart,$queryend,$subjectstart,$subjectend,$expect) = split " ",$k ;
    
        $name =~ s/;//g;
        $name =~ s/.ORF.*//;
        next if($trs eq $name);
    
        print STDERR "This one has $name $blastscore $subjectlength $subjectmatched $iden_percent.  querylength = $querylength\n" if($verbose);
        if($forWGS){
          
          if(exists $DONENAMES->{$name}){
		  		## this is happening as we are doing collison correction
                #warn "$trs has $blastscore > $DONENAMES->{$name}" if($blastscore > $DONENAMES->{$name});
          }
          if($blastscore > $blastcutoff){
              print $ofh " $name $blastscore " ;

			  ## HARDCODED
			  if(! -e "BLASTOUT_SCAFFOLD/$trs.$name.blast"){
                 print $ofhCommands " parseBlastLatestpairwiseStep1.csh $trs $name\n";
			  }
			  $MAPSscafftoTRS->{$name} = {} if(!defined $MAPSscafftoTRS->{$name});
			  $MAPSscafftoTRS->{$name}->{$trs} = 1 ;
              $DONENAMES->{$name} = $blastscore ;
          }
        }
        elsif($blastscore > $blastcutoff){
            $foundone = 1 ;
            if(! exists $done->{$name}){
                print $ofh "$trs $name $blastscore\n" ;
                $done->{$name} = 1;
                $cntfound++;
             }
    
        }
    }
    if($forWGS){
       print $ofh "\n";    
     }
     else{
         print "Info: found $cntfound\n";
         if(!$foundone){
             print $ofh "$trs $trs 0 # only one \n";
         }
    }
}



sub GNM_MapOneTRS2Scaffold{
    my ($trs,$sname,$infile,$strScaffold,$scafflen,$ofh,$ofhlog,$HARD_DIFFFORMATCH,$HARD_LONGDISTANCEFACTOR) = @_ ;
    my ($info,$querylength,$Subjectname) = GNM_PARSEBLAST($infile);
    
    die "$trs - querylength not defined in $infile" if(!defined $querylength);
    
    my $done = {};
    my @others ;
    my $found100percent = 0;
    my @allQ ;
    my @allS ;
    
    
    ## Sort based on blastscore
                                   
    my $sort = {};
    my $DONEBSCORE = {};
    foreach my $k (@{$info}){
        my $org = $k ;
        my ($name,$description,$blastscore,$subjectlength,$iden_percent,$subjectmatched,$querystart,$queryend,$subjectstart,$subjectend,$expect) = split " ",$k ;
        while(exists $DONEBSCORE->{$blastscore}){
            $blastscore = $blastscore - 0.1;
        }
        $DONEBSCORE->{$blastscore} = 1 ;
        $sort->{$k} = $blastscore ;
    }
    my @sorted = sort {$sort->{$b} <=> $sort->{$a}} (keys %{$sort});

    
    
    ## keep a count of where we start
    my $START ;
    my $ignoredfar = 0 ;
    my $ignoredrev = 0 ;
    my $TOTALMATCHEDINSUB = 0 ;
    
    my $MATCHEDScale = {};
    my $bestblastscore ;
    my $MAPSTARTOFBOTH = {};
    
    my $ERRSTATE = 0 ;
    foreach my $k (@sorted){
        my $org = $k ;
        my ($name,$description,$blastscore,$subjectlength,$iden_percent,$subjectmatched,$querystart,$queryend,$subjectstart,$subjectend,$expect) = split " ",$k ;
        print "$name,$description,$blastscore,$subjectlength,$iden_percent,l=$subjectmatched,$querystart,$queryend,$subjectstart,$subjectend,$expect\n" if($verbose);
        $bestblastscore = $blastscore if(!defined $bestblastscore);
    
        if (!defined $START){
           if($subjectstart > $subjectend){
               # this was a rev seq problem
                 print $ofhlog "Error: kk Should not have got reverse here for $trs and $infile.";
				 return;
           }
           $START = $subjectstart ;
        }
        else{
           if($subjectstart > $subjectend){
               print "Info: ignoring since it is reversed, and not the first entry\n" if($verbose);
               $ignoredrev++;
               next ;
           }
        }
    
        ### Ignore something which is really far
        if(abs($START - $subjectstart) > $HARD_LONGDISTANCEFACTOR*$querylength){
            print "Info: Ignored since far away \n" if($verbose);
            $ignoredfar++;
            next ;
        }
    
    
        push @allQ, $querystart;
        push @allQ, $queryend;
        my $diffQ = $queryend - $querystart + 1 ;
    
        push @allS, $subjectstart;
        push @allS, $subjectend;
        $MAPSTARTOFBOTH->{$querystart} = $subjectstart ;
    
        # Annotate a scale
        foreach my $i ($subjectstart..$subjectend){
            $MATCHEDScale->{$i} = 1 ;
        }
    
    
        $TOTALMATCHEDINSUB = $TOTALMATCHEDINSUB + $subjectmatched ;
        ## if we find a match which is long, why bother further
        if(abs($diffQ -$querylength) < $HARD_DIFFFORMATCH){
            print "Info: $diffQ -$querylength, equal or almost equal match found \n" if($verbose);
            last ;
        }
    }
    if($ERRSTATE){
        die "Killing due to previous errors";
    }
    
    
    print "$querylength is querylength and TOTALMATCHEDINSUB is $TOTALMATCHEDINSUB \n" if($verbose);
    
    ## this was just meant for writing $trs.$sname.isrev
    ## please note to check this
    
    my @sortQ = sort {$a <=> $b} @allQ ;
    my @sortS = sort {$a <=> $b} @allS ;
    
    my $N = @sortQ ;
    die if ($N ne @sortS);
    
    my $qS = $sortQ[0];
    my $qE = $sortQ[$N-1];
    
    my $sS = $sortS[0];
    my $sE = $sortS[$N-1];
    
    
    ## ensure the starts are the same 
    die if(!exists $MAPSTARTOFBOTH->{$qS});
    my $MAPPEDEND = $MAPSTARTOFBOTH->{$qS};
    my $orderscrewed = 0;
    if($MAPPEDEND ne $sS){
        print "Order screwed up \n" if($verbose);
        $orderscrewed = 1;
    }
    
    my $foundNNNinProximity= 0 ;
                                                            
    {
       my @llll = split "",$strScaffold;
       my $NNN = @llll  ;
       my $postend = $sE+500 ;
       $postend = $NNN  if($postend >= $NNN);
       my $poststr = util_extractSliceFromFastaString($strScaffold,$sE,$postend);
       $foundNNNinProximity= 1 if($poststr =~ /NNNNNNNN/);
    
       my $prestart = $sS-500 ;
       $prestart = 0 if($prestart < 0);
       $poststr = util_extractSliceFromFastaString($strScaffold,$prestart,$sS);
       $foundNNNinProximity= 1 if($poststr =~ /NNNNNNNNN/);
    }
    
    my $diffQ = $qE - $qS + 1 ;
    my $diffS = $sE - $sS + 1  ;
    
    my $diffMatchedQ = $querylength - $diffQ ;

    ## this can be misleading, esp when orderscrewed is true, but is necesssary when there are X's...
    my $percentmatched = int(100 *($diffQ/$querylength));
    ## this is the real thing matched
    my $percentmatchedreal = int(100 *($TOTALMATCHEDINSUB/$querylength));
    if($diffMatchedQ < $HARD_DIFFFORMATCH){
        print "Exact match (diff=$diffMatchedQ) for $trs and $Subjectname\n";
    }
    
    my $diffMatchedQ2S = $diffS - $querylength ;
    
    
    my $introncount = 0 ;
    my $exoncount = 1 ;

    my $maxintronlength = 0 ;
    my $foundNNNinIntron = 0 ;
    if(1)
    {
       unlink ("INTRONS/$trs.intron");
       die if(!exists $MATCHEDScale->{$sS});
       die if(!exists $MATCHEDScale->{$sE});
       my $STATE = 0 ;
       my $intronstart = 0;
       foreach my $i ($sS..$sE){
           if($STATE eq 0 && !exists $MATCHEDScale->{$i}){
                  # switch to intron state 
                  $STATE  = 1;
              $exoncount++;
               $intronstart = $i ;
           }
           else{
                if($STATE eq 1 && exists $MATCHEDScale->{$i}){
                     # switch to exon state 
                    $STATE = 0 ;
                    my $intronend = $i - 1;
                    my $diff = abs($intronstart - $intronend);
                    $maxintronlength = $diff if($diff > $maxintronlength);
                    my $intronstr = util_extractSliceFromFastaString($strScaffold,$intronstart,$intronend);
                    $foundNNNinIntron = 1 if($intronstr =~ /NNNNNNNNN/);
                    my $ofhintron = util_open_or_append("INTRONS/$trs.intron");
                    print $ofhintron ">$sname.$intronstart.$intronend\n";
                    print $ofhintron ">$intronstr\n";
                    $introncount++;
                }
           }
       }
    }
    my $str = "$trs $scafflen $bestblastscore $exoncount $introncount $diffMatchedQ $percentmatched% $querylength $qS $qE $sS $sE $diffQ $diffS $diffMatchedQ $diffMatchedQ2S $ignoredfar $ignoredrev $maxintronlength $orderscrewed $foundNNNinProximity $foundNNNinIntron $Subjectname";
    print $ofh "$str\n";
    return $str ;
}

## Note that this parses not the first lines (after "sequences") - but the complete file.
sub GNM_PARSEBLAST{
    my ($infile,$verbose) =@_ ;
	$verbose = 0 if(!defined $verbose);
	print "GNM_PARSEBLAST for $infile \n" if($verbose);
    
    my $ifh = util_read($infile);
    my $junk ;
    my $Query ;
    my $Querylength ;
    my $Subjectlength ;
    my @SAVED ;
    my $found = 0 ;
	my $Subjectname ; 
    while(<$ifh>){
        chomp ;
        next if(/^BLAST/);
        next if(/^\s*$/);
        push @SAVED, $_ ;
        if(/Subject=/){
            chomp;
            s/Subject=//;
            $Subjectname = $_ ;
            next ;
        }
        if(/Query= /){
            ($junk,$Query) = split ;
			$Query =~ s/;.*//;
            next ;
        }
        if(/Length=/){
            s/Length=//;
            if(!defined $Querylength){
                $Querylength =$_ ;
            }
            else{
                ## for pw only
                $Subjectlength =$_ ;
            }
            next ;
        }
        if(/Sequences producing significant alignments/){
            $found = 1 ;
            last ;
        }

    }   

    my @alllines ;
    if(!$found){
        print "Processing pairwise\n" if($verbose);
        @alllines = @SAVED ;
        my $CNT = 0 ;
        my @newalllines ;
        my $subject ;
        my $subjectlen ;
        while(@alllines){
            $_ = shift @alllines ;
            if(/^\s*Subject/){
               s/\s*//g;
               ($subject) = (/Subject=(.*)/);
               $_ = shift @alllines ;
               s/\s*//g;
               ($subjectlen) = (/Length=(.*)/);
            }
            if(/Score/){
                push @newalllines, ">$CNT";
                push @newalllines, "Length=$subjectlen";
                $CNT++;
            }
            push @newalllines, $_ ;
        }

        @alllines = @newalllines ;
    }
    else{
    
        while(<$ifh>){
           chomp ;
           push @alllines, $_ ;
        }
    }

    
    my $info = [];
    while(@alllines){
        $_ = shift @alllines ;
        if(/^\s*>/){
            ParseSingleMatch($ifh,$_,\@alllines,$info);
            $found = 1 ;
            last ;
        }
    }


    return ($info,$Querylength,$Subjectname);
}


sub ParseSingleMatch{
    my ($IFH,$line,$alllines,$ALLINFO) = @_ ;
    my $N = @{$alllines};
    #print "$alllines $N kkkkkkkkk\n";
    $line =~ s/>//;


	if(0){
	if($line !~ /\]\s*$/){
        my $NN  = shift @{$alllines};
		chomp $NN ;
		$line = $line . $NN ;
	}
	}
	#if($line !~ /\]\s*$/){
		#print "BAAAD $line \n";
	#}



    my ($name,@l) = split " " ,$line ;
    my $description = join "ZZZ", @l;


    my ($junk,$subjectlength,$blastscore,$expect,$iden_percent,$ratio,$subjectmatched) ;
    my ($querystart,$queryend) ;
    my ($subjectstart,$subjectend) ;

    while(@{$alllines}){
        ## recurse
        $_ = shift @{$alllines};

        if(/^\s*>/){
            ParseSingleMatch($IFH,$_,$alllines,$ALLINFO);
        }
        else{
           if(/Length=/){
               s/Length=//;
               $subjectlength =$_ ;
           }
           elsif(/Score/){
                   if(defined $querystart){
                     $description = "UUU" if($description eq "");
                   my $str = "$name $description $blastscore $subjectlength $iden_percent $subjectmatched  $querystart $queryend $subjectstart $subjectend $expect";
                   #print "XXXXXXX $str\n";
                   push @{$ALLINFO}, $str ;
                }
                   undef $querystart;
                   undef $queryend;
                   undef $subjectstart;
                   undef $subjectend;
                   ($junk,$junk,$blastscore,$junk,$junk,$junk,$junk,$expect)= split ;
                $expect =~ s/,//;
            }
           elsif(/Identities/){
                     ($junk,$junk,$ratio,$iden_percent) = split ;
                   ($junk,$subjectmatched) = split "/", $ratio;
                   $iden_percent =~ s/\%//;
                   $iden_percent =~ s/,//;
                   $iden_percent =~ s/\(//;
                   $iden_percent =~ s/\)//;
           }
           elsif(/Query/){
                 my (@l) = split ;
              my $N= @l -1;
              $querystart = $l[1] if(!defined $querystart);;
              $queryend = $l[$N] ;
           }
           elsif(/Sbjct/){
                 my (@l) = split ;
              my $N= @l -1;
              $subjectstart = $l[1] if(!defined $subjectstart);;
              $subjectend = $l[$N] ;
              #print "$subjectstart $subjectend \n";
           }


        }
    }

    $description = "UUU" if($description eq "");
    my $str = "$name $description $blastscore $subjectlength $iden_percent $subjectmatched  $querystart $queryend $subjectstart $subjectend $expect";
    #print "XXXXXXX $str\n";
    push @{$ALLINFO}, $str ;

}


sub GNM_DirectionBlast{
	my ($infile) = @_ ;
    my $ifh = util_read($infile);
    while(<$ifh>){
		if(/^\s*Strand/){
			if(!/Strand=Plus\/Plus/){
				return 0 ;
			}
			else{
		       return 1 ;
			}
		}
    }
	## no hits
	return -1 ;
}

sub GNM_parseCountSummed{
	my ($annofile,$cutoff) = @_ ;
	my $ifh = util_read($annofile);
	my $table = {};
	my $allinfocounts = {};
    while(<$ifh>){
         next if(/^\s*$/);
	     next if(/^\s*#/);
		 my @l = split ;
		 my $N = @l -1 ; 
		 my $nm = $l[0] ;
		 my $mean = $l[3] ;
		 if($mean > $cutoff){
		     $table->{$nm} = $mean ;
		     $allinfocounts->{$nm} =  {};
		     $allinfocounts->{$nm}->{TOTAL} = $l[1];
		     $allinfocounts->{$nm}->{NUMBER} = $l[2];
		     $allinfocounts->{$nm}->{MEAN} = $l[3];
		     $allinfocounts->{$nm}->{SD} = $l[4];
		     $allinfocounts->{$nm}->{RATION} = $l[5];
		 }
		 
	}
	return ($table,$allinfocounts) ;
    
}
sub GNM_parseANNFile{
	my ($annofile,$tableRRNA) = @_ ;
	my $ifh = util_read($annofile);
	my $table = {};
    while(<$ifh>){
         next if(/^\s*$/);
	     next if(/^\s*#/);
		 my @l = split ;
		 my $N = @l -1 ; 
		 my $nm = $l[0];
		 next if(exists $tableRRNA->{$nm});
		 my $eval = $l[$N ] ;
		 my $blastscore = $l[$N - 1] ;
		 my $percent = $l[$N - 2] ;
		 $table->{$nm} = {};
		 $table->{$nm}->{EVAL} = $eval ;
		 $table->{$nm}->{BLASTSCORE} = $blastscore ;
		 $table->{$nm}->{PERCENT} = $percent ;
		 $table->{$nm}->{FULLLINE} = $_ ;
		 if($nm =~ /_A/){
		 	$nm =~ s/_A//;
		    $table->{$nm} = {};
		    $table->{$nm}->{EVAL} = $eval ;
		    $table->{$nm}->{BLASTSCORE} = $blastscore ;
		    $table->{$nm}->{PERCENT} = $percent ;
		    $table->{$nm}->{FULLLINE} = $_ ;
		 }
	}
	return $table ;
    
}

sub GNM_KMERIZE_sliding{
	my ($savetable,$debugone,$ofhtmp,$str,$table,$LEN,$name,$del,$genome) = @_ ;
	my $len = length($str);
	print "Processing $name \n" if($verbose);

	return if($len < $LEN);
	

	my $start = 0 ;
	my $end = $LEN ;
	my $mapped = 0 ;
	while(1){
	    my $s1 = substr ($str , $start, $LEN)   ;
		if(!defined $s1){
			die "$name substr ( , $start, $LEN)  \n";
		}
		my $lll = length($s1);

		if(!defined $genome || $savetable){ ## save if genome is not defined
		    if(! exists $table->{$s1}){
			    $table->{$s1} = [];
		    }
		    push @{$table->{$s1}}, $name ;
		}
		else{
			if(exists $genome->{$s1}){
				my $genomemap = $genome->{$s1};
				print $ofhtmp "$name $genomemap $s1 QQQ\n"; 
				$mapped = 1;
			}
		}

		if($debugone){
		    my $lll = length($s1);
		    print "$lll $start $end \n";
		    print "$s1\n";
		}
		last if($mapped);

		last if($end eq $len );
		$start = $start + $del ;
		$end = $start + $LEN ;
	}
	die "debugone = $debugone" if($debugone);
	return $mapped ;

}

sub GNM_KMERIZE_hopping{
	my ($savetable,$debugone,$ofhtmp,$str,$table,$LEN,$name) = @_ ;
	my $len = length($str);
	return if($len < $LEN);
	print "Processing $name \n" if($verbose);

	my $start = 0 ;
	my $end = $LEN ;
	while(1){

	    my  $s1 = $end > $len ? substr ($str , $start) : substr ($str , $start, $LEN)   ;

		if(0) {
		    my $AA = util_ConvertNucleotide2Amino($s1,1);
		    $s1 = $AA ;
		}

		### Add data
		if(1){
		    if(! exists $table->{$s1}){
			    #$table->{$s1} = [];
			    #$table->{$s1} = $name;
			    $table->{$s1} = 0 ;

			    #$table->{$s1} = {};
	           #$table->{$s1}->{$name} = 1 ;
		    }
		    else{
		       #push @{$table->{$s1}}, $name ;
	           #$table->{$s1} = $table->{$s1} . ";" . $name ;
	           $table->{$s1} +=1 ;
	           #$table->{$s1}->{$name} = 1 ;
			    #$table->{$s1} = $name;
		    }
		}



		if($debugone){
		    my $lll = length($s1);
		    print "$lll $start $end \n";
		    print "$s1\n";
		}

		last if($end > $len );

		$start = $start + $LEN ;
		$end = $start + $LEN ;
	}
	die "debugone = $debugone" if($debugone);
}

## genometable might be undefined
sub GNM_GenomeBreakHoppingOrSliding{
	my ($savetable,$debugone,$ofhtmp,$INFILE,$ksize,$hopping,$genometable) = @_; 
    my $kmertable = {};
    my $NAMES = {};
    my $ifh = util_read($INFILE);
	my $NAME ;
	my $STR ;
    while(<$ifh>){
	     next if(/^\s*$/);
         if(/^\s*>/){
	 	    ## process one fasta
	 	    if(defined $NAME){
	           _ProcessSingleFasta($kmertable,$NAMES,$NAME,$STR,$savetable,$debugone,$ofhtmp,$INFILE,$ksize,$hopping,$genometable) ;

                $STR = "";
		    }
    
		    my @l = split ;
		    $NAME =  shift @l  ;
		    $NAME =~ s/>//;
	     }
	     else{
	        die if(/>/);
		    ## remove spaces 
		    s/ //g;
		    chomp;
		    $STR = $STR . $_ ;
	     }
    }
	_ProcessSingleFasta($kmertable,$NAMES,$NAME,$STR,$savetable,$debugone,$ofhtmp,$INFILE,$ksize,$hopping,$genometable) ;
    close($ifh);
    return ($kmertable,$NAMES);
}

sub _ProcessSingleFasta{
	my ($kmertable,$NAMES,$NAME,$STR,$savetable,$debugone,$ofhtmp,$INFILE,$ksize,$hopping,$genometable) = @_ ;
	$NAMES->{$NAME} = 1;
	my $len = length($STR);
	if($hopping){
       GNM_KMERIZE_hopping($savetable,$debugone,$ofhtmp,$STR,$kmertable,$ksize,$NAME);
       #my $rev = util_getComplimentaryString($STR) ;
       #GNM_KMERIZE_hopping($savetable,$debugone,$ofhtmp,$rev,$kmertable,$ksize,$NAME);
	}
	else{
        my $mapped = GNM_KMERIZE_sliding($savetable,$debugone,$ofhtmp,$STR,$kmertable,$ksize,$NAME,1,$genometable);
		if(!$mapped){
            my $rev = util_getComplimentaryString($STR) ;
            GNM_KMERIZE_sliding($savetable,$debugone,$ofhtmp,$rev,$kmertable,$ksize,$NAME,1,$genometable);
			}
	}
	print $ofhtmp "$NAME $len kkk\n";
}
