#!/usr/bin/perl -w 

## Transcriptional burst but no translational burst 
## Explore the impact of mRNA copy number and protein synthesis rates on gene expression noise 
## Also quantify effect of mRNA and protein degradation rates  


use Statistics::R;

print "Calculating and plotting time-variant dynamics ....\n"; 

$numevt=10000; ## Num of events 
$bas=0; ## Basal transcription rate  
$tottime=1000; ## Total time of simulation 
$uot=1;  ## Time interval 

$inimrna=0; 
$iniprot=0; 

## Variables - synthesis rates, onrate, offrate, degradation rates 
## --------------------------
 
$mrnarate=10; ## mRNA production rate  

$onrate=0.02; ## rate parameter exponential distribution for onrate 
$offrate=0.2; ## rate parameter exponential distribution for offrate 

$protrate=100; ## protein synthesis rate per mRNA molecule per unit time 

$degm=0.07; ## mRNA degradation rate per mRNA molecule per unit time 
$degp=0.007; ## protein degradation rate per protein molecule per unit time 


## Parameter exploration - one at a time 
#

$mrn[0]=9; # $mrnarate 
$mrn[1][0]=1; $mrn[1][1]=5; $mrn[1][2]=10; $mrn[1][3]=20; $mrn[1][4]=50; $mrn[1][5]=100; $mrn[1][6]=200; $mrn[1][7]=500; $mrn[1][8]=1000;  
 
$ont[0]=9;  # $onrate 0.001 to 10
$ont[1][0]=0.001; $ont[1][1]=0.005; $ont[1][2]=0.01; $ont[1][3]=0.05; $ont[1][4]=0.1; $ont[1][5]=0.5; $ont[1][6]=1; $ont[1][7]=5; $ont[1][8]=10;

$oft[0]=9;   ## $offrate 0.001 to 10
$oft[1][0]=0.001; $oft[1][1]=0.005; $oft[1][2]=0.01; $oft[1][3]=0.05; $oft[1][4]=0.1; $oft[1][5]=0.5; $oft[1][6]=1; $oft[1][7]=5; $oft[1][8]=10;

$prt[0]=9; ## $protrate 
$prt[1][0]=1; $prt[1][1]=5; $prt[1][2]=10; $prt[1][3]=20; $prt[1][4]=50; $prt[1][5]=100; $prt[1][6]=200; $prt[1][7]=500; $prt[1][8]=1000;  

$dm[0]=8; ## $degm
$dm[1][0]=0.01; $dm[1][1]=0.02; $dm[1][2]=0.05; $dm[1][3]=0.1; $dm[1][4]=0.2; $dm[1][5]=0.3; $dm[1][6]=0.4; $dm[1][7]=0.5;  ## 0.01 to 0.5;

$dp[0]=9; ## $degp; 
$dp[1][0]=0.001; $dp[1][1]=0.005; $dp[1][2]=0.01; $dp[1][3]=0.05; $dp[1][4]=0.1; $dp[1][5]=0.2; $dp[1][6]=0.3; $dp[1][7]=0.4; $dp[1][8]=0.5; ## 0.001 to 0.5; 


$runfl=0; 
$tcl=0;
$modrun=1; 

$mp=0; 

while($runfl==0)
{
  ## Parameter value reset 
  
  $mrnarate=10; ## mRNA production rate  
  $onrate=0.02; ## rate parameter exponential distribution for onrate 
  $offrate=0.2; ## rate parameter exponential distribution for offrate 
  $protrate=100; ## protein synthesis rate per mRNA molecule per unit time 
  $degm=0.07; ## mRNA degradation rate per mRNA molecule per unit time 
  $degp=0.007; ## protein degradation rate per protein moleculer per unit time 

   if($tcl==0)
   {
      for($li=0;$li<$mrn[0]; $li++)
      {
	 $mrnarate=$mrn[1][$li]; 

	 $modlist[$mp][0]=$modrun; 
	 $modlist[$mp][1]=$mrnarate; 
	 $modlist[$mp][2]=$onrate; 
	 $modlist[$mp][3]=$offrate; 
	 $modlist[$mp][4]=$protrate; 
	 $modlist[$mp][5]=$degm; 
	 $modlist[$mp][6]=$degp; 
	 $mp++;

	 $modrun++;
      }
      $tcl++; 
      next; 
   } 
   if($tcl==1) 
   {
      for($li=0;$li<$ont[0]; $li++)
      {
         $onrate=$ont[1][$li];

	 $modlist[$mp][0]=$modrun; 
	 $modlist[$mp][1]=$mrnarate; 
	 $modlist[$mp][2]=$onrate; 
	 $modlist[$mp][3]=$offrate; 
	 $modlist[$mp][4]=$protrate; 
	 $modlist[$mp][5]=$degm; 
	 $modlist[$mp][6]=$degp; 
	 $mp++;

         $modrun++;
      }
      $tcl++; 
      next; 
   } 
   if($tcl==2) 
   {
      for($li=0;$li<$oft[0]; $li++)
      {
         $offrate=$oft[1][$li];

	 $modlist[$mp][0]=$modrun; 
	 $modlist[$mp][1]=$mrnarate; 
	 $modlist[$mp][2]=$onrate; 
	 $modlist[$mp][3]=$offrate; 
	 $modlist[$mp][4]=$protrate; 
	 $modlist[$mp][5]=$degm; 
	 $modlist[$mp][6]=$degp; 
	 $mp++;

         $modrun++;
      }
      $tcl++; 
      next; 
   } 
   if($tcl==3) 
   {
      for($li=0;$li<$prt[0]; $li++)
      {
         $protrate=$prt[1][$li];

	 $modlist[$mp][0]=$modrun; 
	 $modlist[$mp][1]=$mrnarate; 
	 $modlist[$mp][2]=$onrate; 
	 $modlist[$mp][3]=$offrate; 
	 $modlist[$mp][4]=$protrate; 
	 $modlist[$mp][5]=$degm; 
	 $modlist[$mp][6]=$degp; 
	 $mp++;

         $modrun++;
      }
      $tcl++; 
      next; 
   } 
   if($tcl==4) 
   {
      for($li=0;$li<$dm[0]; $li++)
      {
         $degm=$dm[1][$li];

	 $modlist[$mp][0]=$modrun; 
	 $modlist[$mp][1]=$mrnarate; 
	 $modlist[$mp][2]=$onrate; 
	 $modlist[$mp][3]=$offrate; 
	 $modlist[$mp][4]=$protrate; 
	 $modlist[$mp][5]=$degm; 
	 $modlist[$mp][6]=$degp; 
	 $mp++;

         $modrun++;
      }
      $tcl++; 
      next; 
   } 
   if($tcl==5) 
   {
      for($li=0;$li<$dp[0]; $li++)
      {
         $degp=$dp[1][$li];

	 $modlist[$mp][0]=$modrun; 
	 $modlist[$mp][1]=$mrnarate; 
	 $modlist[$mp][2]=$onrate; 
	 $modlist[$mp][3]=$offrate; 
	 $modlist[$mp][4]=$protrate; 
	 $modlist[$mp][5]=$degm; 
	 $modlist[$mp][6]=$degp; 
	 $mp++;

         $modrun++;
      }
      $tcl++; 
      next; 
   } 
   if($tcl==6) 
   {
      for($li=0;$li<$ont[0]; $li++)
      {
         $onrate=$ont[1][$li];
    
	 for($lj=0;$lj<$prt[0];$lj++)
	 {
           $protrate=$prt[1][$lj]; 

	   $modlist[$mp][0]=$modrun; 
	   $modlist[$mp][1]=$mrnarate; 
	   $modlist[$mp][2]=$onrate; 
	   $modlist[$mp][3]=$offrate; 
	   $modlist[$mp][4]=$protrate; 
	   $modlist[$mp][5]=$degm; 
	   $modlist[$mp][6]=$degp; 
	   $mp++;

           $modrun++;
	 }
      }
      $tcl++; 
      next; 
   } 
   if($tcl==7) 
   {
      for($li=0;$li<$oft[0]; $li++)
      {
         $offrate=$oft[1][$li];
    
	 for($lj=0;$lj<$prt[0];$lj++)
	 {
           $protrate=$prt[1][$lj]; 

	   $modlist[$mp][0]=$modrun; 
	   $modlist[$mp][1]=$mrnarate; 
	   $modlist[$mp][2]=$onrate; 
	   $modlist[$mp][3]=$offrate; 
	   $modlist[$mp][4]=$protrate; 
	   $modlist[$mp][5]=$degm; 
	   $modlist[$mp][6]=$degp; 
	   $mp++;

           $modrun++;
	 }
      }
      $tcl++; 
      next; 
   } 
   if($tcl==8) 
   {
      for($li=0;$li<$mrn[0]; $li++)
      {
	 $mrnarate=$mrn[1][$li]; 
    
	 for($lj=0;$lj<$prt[0];$lj++)
	 {
           $protrate=$prt[1][$lj]; 

	   $modlist[$mp][0]=$modrun; 
	   $modlist[$mp][1]=$mrnarate; 
	   $modlist[$mp][2]=$onrate; 
	   $modlist[$mp][3]=$offrate; 
	   $modlist[$mp][4]=$protrate; 
	   $modlist[$mp][5]=$degm; 
	   $modlist[$mp][6]=$degp; 
	   $mp++;

           $modrun++;
	 }
      }
      $tcl++; 
      last; 
   } 
}

use LWP::Simple;
use Parallel::ForkManager;

$pm = Parallel::ForkManager->new(20);

LINKS:

for($k=215;$k<217;$k++)
{
   $pm->start and next LINKS; # do the fork

   $modrun=$modlist[$k][0]; 
   $mrnarate=$modlist[$k][1];
   $onrate=$modlist[$k][2];
   $offrate=$modlist[$k][3];
   $protrate=$modlist[$k][4];
   $degm=$modlist[$k][5];
   $degp=$modlist[$k][6];

   print "$modrun mRNArate $mrnarate Onrate $onrate Offrate $offrate Protrate $protrate DegmRNA $degm DegProt $degp\n";

   simulate($modrun); 

   $pm->finish; # do the exit in the child process
}

$pm->wait_all_children;
print "Done...\n"; 

undef(@modlist); 


sub simulate
{

$R=Statistics::R->new();
$R->startR;

local($modno)=@_; 

local $li; 
## -------------------------------------------------------------------
##  Noise calculation in protein levels from 1000 cells 
## -------------------------------------------------------------------


print "Noise calculation across cells ...\n"; 

$numcell=1000; 

for($cn=0;$cn<$numcell;$cn++)
{
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

  for($li=0;$li<$no && $time<=$tottime;$li++)
  {
    $rn=0;
    $time+=$arr1[$li];
    $time=sprintf("%0.0f",$time); ##
    $ttlist{$time}=$mrnarate;
    $time+=$arr2[$li];
    $time=sprintf("%0.0f",$time); ##
    $ttlist{$time}=$bas;
    undef $rn;
  }
  
  $prev=$bas;
  $mrna=$inimrna;
  $prot=$iniprot; 

  $cl=0; 
  for($li=0;$li<$tottime;$li+=$uot)
  {
    $li=sprintf("%0.0f",$li); ##
    ##print "$li\n"; 
    $timelist[$cl]=$li;
    if(exists $ttlist{$li})
    {
      $trlist[$cn][$cl]=$ttlist{$li}; 

      if($prev==$ttlist{$li})
      {
        $prot+=($protrate*$mrna-$degp*$prot);
        $mrna+=($mrnarate-$degm*$mrna); 
      }
      else
      {
        $prot+=($protrate*$mrna-$degp*$prot);
        $mrna+=($bas-$degm*$mrna); 
      }
      if($mrna<0) { $mrna=0; }
      if($prot<0) { $prot=0; }

      $mrnalist[$cn][$cl]=$mrna;
      $protlist[$cn][$cl]=$prot;

      $prev=$ttlist{$li};
    }
    else
    {
      if($prev>$bas)
      {
        $prot+=($protrate*$mrna-$degp*$prot);
        $mrna+=($mrnarate-$degm*$mrna); 
      }
      else
      {
        $prot+=($protrate*$mrna-$degp*$prot);
        $mrna+=($bas-$degm*$mrna); 
      }

      if($mrna<0) { $mrna=0; }
      if($prot<0) { $prot=0; }

      $trlist[$cn][$cl]=$prev; 
      $mrnalist[$cn][$cl]=$mrna;
      $protlist[$cn][$cl]=$prot;
    }
    $cl++;
  }
  undef $prev;
  undef $ref1;
  undef $ref2;
  undef $no;
  undef(@arr1);
  undef(@arr2);
  
  undef %ttlist; 
}


$pp=0; 
for($i=0;$i<$cl;$i++)
{
  $avgtr=0; $sdtr=0; 
  $avgm=0; $sdm=0; 
  $avgp=0; $sdp=0; 

  for($cn=0;$cn<$numcell;$cn++)
  {
    $avgtr+=$trlist[$cn][$i];
    $sdtr+=($trlist[$cn][$i])**2;
    $avgm+=$mrnalist[$cn][$i];
    $sdm+=($mrnalist[$cn][$i])**2;
    $avgp+=$protlist[$cn][$i];
    $sdp+=($protlist[$cn][$i])**2;
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

  if($i>200) 
  {
      ## print WR "$i\t$avgtr\t$sdtr\t$cvtr\t$avgm\t$sdm\t$cvm\t$avgp\t$sdp\t$cvp\n"; }
      $silist[$pp][0]=$i; 
      $silist[$pp][1]=$avgtr; 
      $silist[$pp][2]=$sdtr; 
      $silist[$pp][3]=$cvtr; 
      $silist[$pp][4]=$avgm; 
      $silist[$pp][5]=$sdm; 
      $silist[$pp][6]=$cvm; 
      $silist[$pp][7]=$avgp; 
      $silist[$pp][8]=$sdp; 
      $silist[$pp][9]=$cvp; 
      $pp++;
  }

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


for($cn=0;$cn<$numcell;$cn++)
{
  $bsize=0; $bfreq=0; 
  for($i=0;$i<$cl;$i++)
  {
    if($trlist[$cn][$i]>$bas) 
    {
      $bsize+=$trlist[$cn][$i];
      $bfreq++; 
    }
  }
  if($bfreq!=0)
  {
     $bsize=sprintf("%0.4f",$bsize/$bfreq); 
  }
  else
  {
     $bsize="NA"; 
  }
  $bfreq=sprintf("%0.4f",$bfreq/$cl); 
  $silist2[$cn][0]=$bsize; 
  $silist2[$cn][1]=$bfreq; 
  undef $bsize; 
  undef $bfreq; 
}

undef $cl;
undef(@trlist);
undef(@timelist); 
undef(@mrnalist); 
undef(@protlist); 


## Competitive TFs
$ff="../results/parameter_expl_nonburstytl/$modno"."_NoiseData.txt";
open(WR1,">$ff") or die; 

print WR1 "$modrun mRNArate $mrnarate Onrate $onrate Offrate $offrate Protrate $protrate DegmRNA $degm DegProt $degp\n";

print WR1 "Timepoint\tAvg_TR\tSd_TR\tCV_TR\tAvg_mRNA\tSd_mRNA\tCV_mRNA\tAvg_Prot\tSd_Prot\tCV_Prot\n"; 

for($i=0;$i<$pp;$i++)
{
  for($lo=0;$lo<10;$lo++)
  {
	  print WR1 "$silist[$i][$lo]\t"; 
  }
  print WR1 "\n"; 
}
close(WR1);


$ff="../results/parameter_expl_nonburstytl/$modno"."_Burst_SizeFreq.txt";
open(WR2,">$ff") or die; 

print WR2 "$modrun mRNArate $mrnarate Onrate $onrate Offrate $offrate Protrate $protrate DegmRNA $degm DegProt $degp\n";

print WR2 "CellID\tBurstSize\tBurstFreq\n"; 

for($i=0;$i<$numcell;$i++)
{
  print WR2 "$i\t$silist2[$i][0]\t$silist2[$i][1]\n"; 
}
close(WR2); 

undef(@silist);
undef(@silist2);

undef $cl;
undef(@trlist);
undef(@timelist);
undef(@mrnalist);
undef(@protlist);

$R->stopR;

}
