#!/usr/bin/perl -w

#system "rm -r Output/";
#system "mkdir Output/";

system "ls ../Input_FCS/*.fcs > FCSLIST";

open(FCS,"FCSLIST") or die;
while($fcs=<FCS>)
{
  chomp($fcs);
  print "$fcs\n";

  @wr=split(/\//,$fcs);
  $lno=@wr;
  $infile=$wr[$lno-1];
  undef(@wr);
  undef $lno;
  @re=split(/\./,$infile);
  @arr=split(/\_/,$re[0]);
  $no=@arr;


  if($no==3) 
  { 
    $dr=$arr[1];
    $outfile="../Output/$dr/DATA/$dr".".txt"; 
    $imgout="../Output/$dr/PLOTS/$arr[1]"; 
    $txtout="../Output/$dr/DATA/$arr[1]"; 
  }  
  elsif($no==4) 
  { 
    $dr=$arr[1]."_".$arr[2];
    $outfile="../Output/$dr/DATA/$dr".".txt"; 
    $imgout="../Output/$dr/PLOTS/$arr[1]"; 
    $txtout="../Output/$dr/DATA/$arr[1]"; 
  }  
  elsif($no==5) 
  { 
    $dr=$arr[1]."_".$arr[2]."_".$arr[3];
    $outfile="../Output/$dr/DATA/$dr".".txt"; 
    $imgout="../Output/$dr/PLOTS/$arr[1]"; 
    $txtout="../Output/$dr/DATA/$arr[1]"; 
  }  
  elsif($no==6) 
  { 
    $dr=$arr[1]."_".$arr[2]."_".$arr[3]."_".$arr[4];
    $outfile="../Output/$dr/DATA/$dr".".txt"; 
    $imgout="../Output/$dr/PLOTS/$arr[1]"; 
    $txtout="../Output/$dr/DATA/$arr[1]"; 
  }  
  ##print "$no $dr $infile $outfile $imgout\n";
  ##next; 
  if(-d "../Output/$dr") { system "rm -r ../Output/$dr"; }
  system "mkdir ../Output/$dr"; 
  system "mkdir ../Output/$dr/DATA"; 
  system "mkdir ../Output/$dr/PLOTS"; 
  undef $no;
  

  $dir="../Input_FCS/";
  system "R --slave --no-save --no-restore --no-environ --silent --args $infile $outfile $dir $imgout $txtout< FCS_R_code.R";
  
  undef $infile;
  undef $outfile;
  undef $dr;
  undef $dir;
  undef $imgout;
  undef(@arr);
  undef(@re);
}
close(FCS);

system "rm FCSLIST";

