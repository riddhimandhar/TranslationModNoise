#!/usr/bin/perl -w 

## Simulate how bursty translation impacts gene expression noise 
## Single mRNA tracking and identify all events 
## Protien birth and degradation 
## Collect ribo availability with translation initiation rate 
### Parallelization 
## FN16 riboav calculation 

use warnings;

use Statistics::R;

$timcut=200; 

$tottime=1000; ## Min
$uot=1;  ## Min

$numevt=10000; 
$bas=0; 

$mrnarate=4;  ## Transcription initiation rate per min 
$mrnaddtime=2; ## in min (for mRNA synthesis + transport)
$mrnahl=10; ## mRNA half-life in min 
#$degm=0.693/$mrnahl; ## per mRNA molecule per unit time per min 

$onrate=0.02;  ## per min  ## 0.5
$offrate=0.2; ##  per min 
$tlon=0.5;
$tloff=0.5;

$tlinrate=5; ## Translation initiation rate per min  ## default 0.5 
$khl=30; ## Km quivalent
$aatrav=100;    ## per mRNA molecule per unit time (codons per min) ## 100  
$orflen=300; ## No. of codons/AA

$trvtime=($orflen-1)/$aatrav+log($orflen)*$khl/$aatrav; ## Hill fn.  
$t10=(10-1)/$aatrav+log(10)*$khl/$aatrav; 

$trvtime=sprintf("%0.3f",$trvtime); 
$t10=sprintf("%0.3f",$t10); 

$prothl=100; ## Protein halflife in min 
#$degp=0.693/$prothl; ## per protein moleculer per unit time (per min) 

$inimrna=0; 
$iniprot=0; 
$prevmr=0; 
$prevrb=0; 


## -------------------------------------------------------------------
##  Noise calculation in protein levels from 10000 cells 
## -------------------------------------------------------------------

## Setting Parameter ranges 
## Setting Parameter ranges
### $mrnarate 0.1 to 20
## $onrate 0.01 to 10
### $offrate 0.01 to 10
### $tlon 0.1 to 10
### $tloff 0.1 to 10

$ont[0]=8;  # $onrate 0.01 to 10
$ont[1][0]=0.01; $ont[1][1]=0.05; $ont[1][2]=0.1; $ont[1][3]=0.5; $ont[1][4]=1; $ont[1][5]=2; $ont[1][6]=5; $ont[1][7]=10; 

$oft[0]=8;   ## $offrate 0.01 to 10
$oft[1][0]=0.01; $oft[1][1]=0.05; $oft[1][2]=0.1; $oft[1][3]=0.5; $oft[1][4]=1; $oft[1][5]=2; $oft[1][6]=5; $oft[1][7]=10;  

$ton[0]=6;  ## $tlon 0.01 to 10  
$ton[1][0]=0.01; $ton[1][1]=0.05; $ton[1][2]=0.1; $ton[1][3]=0.5; $ton[1][4]=1; $ton[1][5]=5;  

$tof[0]=6;  ## $tloff 0.01 to 10 
$tof[1][0]=0.01; $tof[1][1]=0.05; $tof[1][2]=0.1; $tof[1][3]=0.5; $tof[1][4]=1; $tof[1][5]=5;  

$mr[0]=8;  ## $mrnarate 0.1 to 20 
$mr[1][0]=0.1; $mr[1][1]=0.5; $mr[1][2]=1; $mr[1][3]=2; $mr[1][4]=5; $mr[1][5]=10; $mr[1][6]=15; $mr[1][7]=20; 

$tlr[0]=8; ## $tlinrate 0.1 to 10
$tlr[1][0]=0.1; $tlr[1][1]=0.5; $tlr[1][2]=1; $tlr[1][3]=2; $tlr[1][4]=4; $tlr[1][5]=6; $tlr[1][6]=8;  $tlr[1][7]=10;      

$aatrav=100;    ## per mRNA molecule per unit time (codons per min) ## 100
$orflen=300; ## No. of codons/AA


$runfl=0;
$tcl=0;
$modrun=1;

$mp=0;

while($runfl==0) 
{
   $mrnarate=4; ## $mr 
   $mrnahl=10;  ## $mhl
   $onrate=0.02;  ## $ont
   $offrate=0.2; ##  $oft
   $tlon=0.5;  ##  $ton
   $tloff=0.5; ##  $tof 
   $tlinrate=5; ## $tlr
   $khl=30; ##  $kk 
   $aatrav=100;    ## per mRNA molecule per unit time (codons per min) ## 100
   $orflen=300; ## No. of codons/AA
   $trvtime=($orflen-1)/$aatrav+log($orflen)*$khl/$aatrav; ## Hill fn.
   $t10=(10-1)/$aatrav+log(10)*$khl/$aatrav;
   $prothl=100;  ## $phl

   $trvtime=sprintf("%0.3f",$trvtime); 
   $t10=sprintf("%0.3f",$t10); 

   if($tcl==0) 
   {
     for($jj=0;$jj<$ont[0];$jj++)
     {
        $onrate=$ont[1][$jj]; 

        for($cc=0;$cc<$tlr[0];$cc++)
        {
          $tlinrate=$tlr[1][$cc]; 
  	  $khl=20+20*5/(5+$tlinrate);
          $trvtime=($orflen-1)/$aatrav+log($orflen)*$khl/$aatrav; ## Hill fn.  
          $t10=(10-1)/$aatrav+log(10)*$khl/$aatrav; 

	  $khl=sprintf("%0.3f",$khl); 
          $trvtime=sprintf("%0.3f",$trvtime); 
          $t10=sprintf("%0.3f",$t10); 

          $modlist[$mp][0]=$modrun; 
          $modlist[$mp][1]=$mrnarate; 
          $modlist[$mp][2]=$mrnahl; 
          $modlist[$mp][3]=$onrate; 
          $modlist[$mp][4]=$offrate; 
          $modlist[$mp][5]=$tlon; 
          $modlist[$mp][6]=$tloff; 
          $modlist[$mp][7]=$tlinrate; 
          $modlist[$mp][8]=$khl; 
          $modlist[$mp][9]=$trvtime; 
          $modlist[$mp][10]=$t10; 
          $modlist[$mp][11]=$prothl; 

	  $mp++; 
	  $modrun++; 
       }
     }
     $tcl++; 
     next; 
   }
   if($tcl==1) 
   {
     for($jj=0;$jj<$oft[0];$jj++)
     {
        $offrate=$oft[1][$jj]; 
        for($cc=0;$cc<$tlr[0];$cc++)
        {
          $tlinrate=$tlr[1][$cc]; 
  	  $khl=20+20*5/(5+$tlinrate);
          $trvtime=($orflen-1)/$aatrav+log($orflen)*$khl/$aatrav; ## Hill fn.  
          $t10=(10-1)/$aatrav+log(10)*$khl/$aatrav; 

	  $khl=sprintf("%0.3f",$khl); 
          $trvtime=sprintf("%0.3f",$trvtime); 
          $t10=sprintf("%0.3f",$t10); 

          $modlist[$mp][0]=$modrun; 
          $modlist[$mp][1]=$mrnarate; 
          $modlist[$mp][2]=$mrnahl; 
          $modlist[$mp][3]=$onrate; 
          $modlist[$mp][4]=$offrate; 
          $modlist[$mp][5]=$tlon; 
          $modlist[$mp][6]=$tloff; 
          $modlist[$mp][7]=$tlinrate; 
          $modlist[$mp][8]=$khl; 
          $modlist[$mp][9]=$trvtime; 
          $modlist[$mp][10]=$t10; 
          $modlist[$mp][11]=$prothl; 

	  $mp++; 
	  $modrun++; 
       }
     }
     $tcl++; 
     next; 
   }
   if($tcl==2) 
   {
     for($jj=0;$jj<$ton[0];$jj++)
     {
        $tlon=$ton[1][$jj]; 

        for($cc=0;$cc<$tlr[0];$cc++)
        {
          $tlinrate=$tlr[1][$cc]; 
  	  $khl=20+20*5/(5+$tlinrate);
          $trvtime=($orflen-1)/$aatrav+log($orflen)*$khl/$aatrav; ## Hill fn.  
          $t10=(10-1)/$aatrav+log(10)*$khl/$aatrav; 

	  $khl=sprintf("%0.3f",$khl); 
          $trvtime=sprintf("%0.3f",$trvtime); 
          $t10=sprintf("%0.3f",$t10); 

          $modlist[$mp][0]=$modrun; 
          $modlist[$mp][1]=$mrnarate; 
          $modlist[$mp][2]=$mrnahl; 
          $modlist[$mp][3]=$onrate; 
          $modlist[$mp][4]=$offrate; 
          $modlist[$mp][5]=$tlon; 
          $modlist[$mp][6]=$tloff; 
          $modlist[$mp][7]=$tlinrate; 
          $modlist[$mp][8]=$khl; 
          $modlist[$mp][9]=$trvtime; 
          $modlist[$mp][10]=$t10; 
          $modlist[$mp][11]=$prothl; 

	  $mp++; 
	  $modrun++; 
       }
     }
     $tcl++; 
     next; 
   }
   if($tcl==3) 
   {
     for($jj=0;$jj<$tof[0];$jj++)
     {
        $tloff=$tof[1][$jj]; 

        for($cc=0;$cc<$tlr[0];$cc++)
        {
          $tlinrate=$tlr[1][$cc]; 
  	  $khl=20+20*5/(5+$tlinrate);
          $trvtime=($orflen-1)/$aatrav+log($orflen)*$khl/$aatrav; ## Hill fn.  
          $t10=(10-1)/$aatrav+log(10)*$khl/$aatrav; 

	  $khl=sprintf("%0.3f",$khl); 
          $trvtime=sprintf("%0.3f",$trvtime); 
          $t10=sprintf("%0.3f",$t10); 

          $modlist[$mp][0]=$modrun; 
          $modlist[$mp][1]=$mrnarate; 
          $modlist[$mp][2]=$mrnahl; 
          $modlist[$mp][3]=$onrate; 
          $modlist[$mp][4]=$offrate; 
          $modlist[$mp][5]=$tlon; 
          $modlist[$mp][6]=$tloff; 
          $modlist[$mp][7]=$tlinrate; 
          $modlist[$mp][8]=$khl; 
          $modlist[$mp][9]=$trvtime; 
          $modlist[$mp][10]=$t10; 
          $modlist[$mp][11]=$prothl; 

	  $mp++; 
	  $modrun++; 
       }
     }
     $tcl++; 
     next; 
   }
   if($tcl==4) 
   {
     for($jj=0;$jj<$mr[0];$jj++)
     {
        $mrnarate=$mr[1][$jj]; 

        for($cc=0;$cc<$tlr[0];$cc++)
        {
          $tlinrate=$tlr[1][$cc]; 
  	  $khl=20+20*5/(5+$tlinrate);
          $trvtime=($orflen-1)/$aatrav+log($orflen)*$khl/$aatrav; ## Hill fn.  
          $t10=(10-1)/$aatrav+log(10)*$khl/$aatrav; 

	  $khl=sprintf("%0.3f",$khl); 
          $trvtime=sprintf("%0.3f",$trvtime); 
          $t10=sprintf("%0.3f",$t10); 

          $modlist[$mp][0]=$modrun; 
          $modlist[$mp][1]=$mrnarate; 
          $modlist[$mp][2]=$mrnahl; 
          $modlist[$mp][3]=$onrate; 
          $modlist[$mp][4]=$offrate; 
          $modlist[$mp][5]=$tlon; 
          $modlist[$mp][6]=$tloff; 
          $modlist[$mp][7]=$tlinrate; 
          $modlist[$mp][8]=$khl; 
          $modlist[$mp][9]=$trvtime; 
          $modlist[$mp][10]=$t10; 
          $modlist[$mp][11]=$prothl; 

	  $mp++; 
	  $modrun++; 
       }
     }
     $tcl++; 
     last; 
   }
}

print "Total no. of models: $mp\n"; 

## 

use LWP::Simple;
use Parallel::ForkManager;

$numcell=1000; 

$mc=0;
for($mc=169;$mc<$mp;$mc++) ##
{
   $modrun=$modlist[$mc][0];
   $mrnarate=$modlist[$mc][1];
   $mrnahl=$modlist[$mc][2];
   $onrate=$modlist[$mc][3];
   $offrate=$modlist[$mc][4];
   $tlon=$modlist[$mc][5];
   $tloff=$modlist[$mc][6];
   $tlinrate=$modlist[$mc][7];
   $khl=$modlist[$mc][8];
   $trvtime=$modlist[$mc][9];
   $t10=$modlist[$mc][10];
   $prothl=$modlist[$mc][11];

   print "$modrun mRNArate $mrnarate mRNAhl $mrnahl Onrate $onrate Offrate $offrate TLon $tlon TLoff $tloff TlInrate $tlinrate Ktr $khl TrvTi $trvtime Trv10 $t10 ProtHL $prothl\n";

   $pm = Parallel::ForkManager->new(100);
   $pm->run_on_finish (
   sub {
     local($pid, $exit_code, $ident, $exit_signal, $core_dump, $dataref) = @_;
     local($lcl,$lcell,@darr,$li);

     @darr=@{$dataref}; 

     $lcl=@{$darr[1]}; 
     $lcell=$darr[4];
     
     #print "$lcell $lcl $darr[1]\n"; 
     for($li=0;$li<$lcl;$li++) 
     {
	$timelist[$li]=$li; 
	$trlist[$lcell][$li]=$darr[3][$li];
	$mrnalist[$lcell][$li]=$darr[0][$li];
	$protlist[$lcell][$li]=$darr[1][$li];
	#print "  $li $darr[1][$li] $darr[0][$li] $darr[2][$li]\n"; 
     }
     undef(@darr); 
     undef $dataref; 
     undef(@trackmrn); 
   }); 

   LINKS:
   for($cn=0;$cn<$numcell;$cn++)
   {
	   ##print "$cn\n"; 
     $pm->start and next LINKS; # do the fork

     $ref=simulate_cell($cn);

     $pm->finish(0,$ref); 
   }
 
   $pm->wait_all_children; 
   undef $pm; 


   $cl=@timelist; 

   open(WR,">../results/parameter_expl_waittime3/$modrun"."_NoiseData.txt") or die; 
   print WR "$modrun mRNArate $mrnarate mRNAhl $mrnahl Onrate $onrate Offrate $offrate TLon $tlon TLoff $tloff TlInrate $tlinrate Ktr $khl TrvTi $trvtime Trv10 $t10 ProtHL $prothl\n";
   print WR "Timepoint\tAvg_TR\tSd_TR\tCV_TR\tAvg_mRNA\tSd_mRNA\tCV_mRNA\tAvg_Prot\tSd_Prot\tCV_Prot\n"; 

  $avmep=0;  
  $avcvp=0; $cnc=0; 

  for($cc=$timcut;$cc<$cl;$cc++)
  {
    $avgtr=0; $sdtr=0; 
    $avgm=0; $sdm=0; 
    $avgp=0; $sdp=0; 

    for($cn=0;$cn<$numcell;$cn++)
    {
      $avgtr+=$trlist[$cn][$cc];
      $sdtr+=($trlist[$cn][$cc])**2;
      $avgm+=$mrnalist[$cn][$cc];
      $sdm+=($mrnalist[$cn][$cc])**2;
      $avgp+=$protlist[$cn][$cc];
      $sdp+=($protlist[$cn][$cc])**2;
    }

    $sdm=(($sdm/$numcell)-(($avgm/$numcell)**2));
    if($sdm<0) { $sdm=0; }
    $sdm=sprintf("%0.4f",sqrt($sdm));
    $avgm=sprintf("%0.4f",$avgm/$numcell);
    if($avgm!=0) 
    { 
       $cvm=sprintf("%0.4f",$sdm/$avgm);
    }
    else 
    {   
       $cvm="NA"; 
    } 

    $sdp=(($sdp/$numcell)-(($avgp/$numcell)**2));
    if($sdp<0) { $sdp=0; }
    $sdp=sprintf("%0.4f",sqrt($sdp));
    $avgp=sprintf("%0.4f",$avgp/$numcell);
    if($avgp!=0) 
    { 
       $cvp=sprintf("%0.4f",$sdp/$avgp);
       $avmep+=$avgp; 
       $avcvp+=$cvp; 
       $cnc++; 
    }
    else
    {
       $cvp="NA"; 
    }

    $sdtr=(($sdtr/$numcell)-(($avgtr/$numcell)**2));
    if($sdtr<0) { $sdtr=0; }
    $sdtr=sprintf("%0.4f",sqrt($sdtr));
    $avgtr=sprintf("%0.4f",$avgtr/$numcell);
    if($avgtr!=0) 
    { 
       $cvtr=sprintf("%0.4f",$sdtr/$avgtr);
    }
    else
    {
       $cvtr="NA"; 
    }

    if($cc>$timcut) { print WR "$cc\t$avgtr\t$sdtr\t$cvtr\t$avgm\t$sdm\t$cvm\t$avgp\t$sdp\t$cvp\n"; }

    undef $avgtr;
    undef $sdtr;
    undef $cvtr; 
    undef $avgm;
    undef $sdm;
    undef $cvm;
    undef $avgp;
    undef $sdp;
    undef $cvp;

  }
  close(WR); 


  if($cnc!=0) 
  {
     $avcvp=sprintf("%0.4f",$avcvp/$cnc); 
     $avmep=sprintf("%0.4f",$avmep/$cnc); 
  }
  else { $avcvp="NA"; $avmep="NA";  } 

  print "AvgProtExp: $avmep  AvgProtCV: $avcvp\n"; 


  open(WR,">../results/parameter_expl_waittime3/$modrun"."_Burst_SizeFreq.txt") or die; 
  print WR "$modrun mRNArate $mrnarate mRNAhl $mrnahl Onrate $onrate Offrate $offrate TLon $tlon TLoff $tloff TlInrate $tlinrate Ktr $khl TrvTi $trvtime Trv10 $t10 ProtHL $prothl\n";
  print WR "CellID\tBurstSize\tBurstFreq\n"; 

  for($cn=0;$cn<$numcell;$cn++)
  {
    $bsize=0; $bfreq=0; 
    for($cc=0;$cc<$cl;$cc++)
    {
      if($trlist[$cn][$cc]>$bas) 
      {
        $bsize+=$trlist[$cn][$cc];
        $bfreq++; 
      }
    }
    if($bfreq!=0) 
    {
      $bsize=sprintf("%0.4f",$bsize/$bfreq); 
    }
    else { $bsize="NA"; } 
    $bfreq=sprintf("%0.4f",$bfreq/$cl); 
    print WR "$cn\t$bsize\t$bfreq\n"; 
    undef $bsize; 
    undef $bfreq; 
  }
  close(WR); 

  undef $cl;
  undef(@trlist);
  undef(@timelist); 
  undef(@mrnalist); 
  undef(@protlist); 
  undef(@trackmrn); 
}

undef(@modlist); 


sub simulate_cell
{
  local($lcn)=@_; 

  $R=Statistics::R->new();
  $R->start_sharedR;

  ##print "  Cell no: $lcn\n"; 

  $R->run(qq'di=rexp($numevt,rate=$onrate)');
  $ref1=$R->get('di');
  $R->run(q'rm(di)');
  @arr1=@{$ref1};

  $R->run(qq'di=rexp($numevt,rate=$offrate)');
  $ref2=$R->get('di');
  $R->run(q'rm(di)');
  @arr2=@{$ref2};

  $no=@arr1;

  $time=0;

  for($i=0;$i<$no && $time<=$tottime;$i++)
  {
    $rn=0;
    $time+=$arr1[$i];
    $ntime=sprintf("%0.0f",$time+$mrnaddtime);
    $ttlist{$ntime}=$mrnarate;
    $time+=$arr2[$i];
    $ntime=sprintf("%0.0f",$time+$mrnaddtime);
    $ttlist{$ntime}=$bas;
    undef $rn;
  }
  
  undef $ref1;
  undef $ref2;
  undef $no;
  undef(@arr1);
  undef(@arr2);

  $prev=$bas;
  $mrna=$inimrna;
  $prevmr=0; 

  $mcc=0; 
  $incr=0; 

  for($i=0;$i<$tottime;$i+=$uot) 
  {
    $trackmrn[0][$i]=$inimrna; 
    $trackmrn[1][$i]=$iniprot;
    $trackmrn[2][$i]=0;
    $trackmrn[3][$i]=0;
  }
  $trackmrn[4]=$cn; 

  for($i=0;$i<$tottime;$i+=$uot)
  {
    $i=sprintf("%0.0f",$i);
    ## mRNA numbers (active) 
  
    if(exists $ttlist{$i})
    {
      $trackmrn[3][$i]=$ttlist{$i}; 
      $mrna+=($ttlist{$i}); 
      $incr=$ttlist{$i}; 

      $prev=$ttlist{$i};
    }
    else
    {
      $trackmrn[3][$i]=$prev; 
      $mrna+=($prev); 
      $incr=$prev; 
    }

    for($k=$mcc;$k<$mcc+$incr;$k++) 
    {
        ## $trackmrn[$k][0]=$i-$mrnaddtime;  ## mRNA birth time 
        ## $trackmrn[$k][1]=$i-$mrnaddtime+$mrnahl*2; ## mRNA death time; 
    
        $lb=$mrnahl*0.1; 
        $off=int(rand($lb)); ## Offset  
        $sig=rand(1); 
        if($sig>0.5) { $off*=-1; }  ## Sign of the offset 
        undef $lb; 
        undef $sig; 

        $fin=$i-$mrnaddtime+$mrnahl+$off-$trvtime; 
        if($fin>$tottime) { $fin=$tottime; } 

        $tln=sprintf("%0.0f",$mrnahl*1000);
        $R->run(qq'di=rexp($tln,rate=$tlon)');
        $ref1=$R->get('di');
        $R->run(q'rm(di)');
        @arr1=@{$ref1};

        $R->run(qq'di=rexp($tln,rate=$tloff)');
        $ref2=$R->get('di');
        $R->run(q'rm(di)');
        @arr2=@{$ref2};

        $tlt=0;
        for($m=0;$m<$tln && $tlt<=$tottime;$m++)
        {
           $tlt+=$arr1[$m];
           $tlt=sprintf("%0.0f",$tlt);
           $tldat{$tlt}=0;
           #print "$tlt $tldat{$tlt} $arr1[$m]\n";
           $tlt+=$arr2[$m];
           $tlt=sprintf("%0.0f",$tlt);
           $tldat{$tlt}=2;
           #print "$tlt $tldat{$tlt} $arr2[$m]\n";
        }
        ##print "$tln $tlt\n";
        undef $tln;
        undef $ref1;
        undef $ref2;
        undef(@arr1); 
        undef(@arr2); 
        undef $tlt; 

        $tlprev=2; 

        for($m=$i;$m<$fin;$m+=$uot) 
        {
	  $m=sprintf("%0.0f",$m);
          $trackmrn[0][$m]+=1; 

	  if(exists $rnafl{$m}) 
	  {
             $cnt=$rnafl{$m}[0]; 
	     if(exists $tldat{$m})
	     {
	        $rnafl{$m}[1][$cnt]=$tldat{$m};
	        $tlprev=$tldat{$m};
	     }
	     else
	     {
	        $rnafl{$m}[1][$cnt]=$tlprev;
	     }
	     $rnafl{$m}[2][$cnt]=$k; 
	     $rnafl{$m}[0]++; 
	     undef $cnt; 
	  }
	  else
	  {
	     if(exists $tldat{$m})
	     {
	        $rnafl{$m}[1][0]=$tldat{$m};
	        $tlprev=$tldat{$m};
	     }
	     else
	     {
	        $rnafl{$m}[1][0]=$tlprev;
	     }
	     $rnafl{$m}[2][0]=$k; 
	     $rnafl{$m}[0]=1; 
	  }
        }

        for($m=$fin;$m<$fin+$trvtime;$m+=$uot)
        {
	   $m=sprintf("%0.0f",$m);
	   $trackmrn[0][$m]+=1;
	}
	undef $off; 
	undef $fin; 
	undef $tlprev; 
	undef %tldat; 
    }

    $riboav=(1+$prevrb)/(1+$trackmrn[2][$i])*(1+$prevmr)/(1+$trackmrn[0][$i]);
    $lkc=(1+5/(1+$prevrb+$prevmr));
    $templin=$tlinrate*($lkc**10)/(($lkc**10)+(1/$riboav)**10);
    undef $lkc; 
    #    $riboav=(1+$prevmr)/(1+$trackmrn[0][$i]);
    #   $kc=(1+6/(1+$prevmr));
    #$templin=$tlinrate*($kc**10)/(($kc**10)+(1/$riboav)**10);


    if(exists $rnafl{$i}) 
    {
      for($k=0;$k<$rnafl{$i}[0];$k++)
      {
        if($rnafl{$i}[1][$k]==0)
        {
	  $rnafl{$i}[1][$k]=0;
	  $tlref=1/$templin-0.001;

     	  for($j=$i;$j<$i+$uot;$j+=$tlref)
	  {
	    $st=sprintf("%0.0f",$j); 
            if($st>=$tottime) { next; } 

            $flt=$st+$t10; 
            if($flt>=$tottime) { $flt=$tottime; } 

            for($u=$st;$u<$flt;$u+=$uot) 
            {
              $u=sprintf("%0.0f",$u);

	      if(exists $rnafl{$u}) 
	      {
	         $lfl=0; 
	         for($m=0;$m<$rnafl{$u}[0];$m++)
	         {
	           if($rnafl{$u}[2][$m] == $rnafl{$i}[2][$k])
	           {
		     $lfl=1; 
	             $rnafl{$u}[1][$m]=1; 
	           }
	         }
	      }
            }

	    $occ=$st+$trvtime; 
	    if($occ>$tottime) { $occ=$tottime; } 
	    for($u=$st;$u<$occ;$u+=$uot) 
    	    {
		$u=sprintf("%0.0f",$u);
                $trackmrn[2][$u]+=1;            
	    }

            $prin=sprintf("%0.0f",$st+$trvtime); 

            $lb=$prothl*0.1; 
            $off=int(rand($lb)); ## Offset
            $sig=rand(1); 
            if($sig>0.5) { $off*=-1; } ## Sign of the offset 
            undef $lb; 
            undef $sig; 

            $pren=$st+$trvtime+$prothl+$off; 
            undef $off; 

            if($prin>$tottime) { $prin=$tottime; } 
            if($pren>$tottime) { $pren=$tottime; } 

            for($u=$prin;$u<$pren;$u+=$uot) 
            {
              $u=sprintf("%0.0f",$u);
              $trackmrn[1][$u]+=1; 
            }
	    undef $st; 
	    undef $occ; 
            undef $flt; 
            undef $prin; 
            undef $pren; 
            }
            undef $tlref; 
           }
         }
       }
       undef $templin; 

       $mcc+=$incr; 
       $prevmr=$trackmrn[0][$i]; 
       $prevrb=$trackmrn[2][$i]; 
  }

  undef %rnafl;
  undef %ttlist; 

  undef $prev;
  undef $ref1;
  undef $ref2;
  undef $no;
  undef(@arr1);
  undef(@arr2);

  $R->stopR; 

  $lref=\@trackmrn;
  return $lref; 
}

