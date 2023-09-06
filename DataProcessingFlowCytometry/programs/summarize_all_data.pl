#!/usr/bin/perl -w 

system "ls ../Input_FCS/*.fcs > FCSLIST";

open(FCS,"FCSLIST") or die;
open(WR,">../Output/All_data_summary.txt") or die; 
while($fcs=<FCS>)
{
  chomp($fcs);

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
  }
  elsif($no==4)
  {
    $dr=$arr[1]."_".$arr[2];
    $outfile="../Output/$dr/DATA/$dr".".txt";
  }
  elsif($no==5)
  {
    $dr=$arr[1]."_".$arr[2]."_".$arr[3];
    $outfile="../Output/$dr/DATA/$dr".".txt";
  }
  elsif($no==6)
  {
    $dr=$arr[1]."_".$arr[2]."_".$arr[3]."_".$arr[4];
    $outfile="../Output/$dr/DATA/$dr".".txt";
  }
  print "$dr.txt\n";
  print WR "$dr.txt\n";
  ##print "$no $dr $infile $outfile $imgout\n";
  ##next;
  $cn=0;

  $flag=0;
  open(FP,$outfile) or die;
  while($fp=<FP>)
  {
    chomp($fp);
    @pw=split(/\s+/,$fp);
    if($fp=~/\[1\]/ && $flag==0) 
    {
      $par1=$pw[1];
      $flag=1; 
      next;
    }
    if($flag==1 && $fp=~/myEllipsoidGate/)
    {
      $event=$pw[2]; 
      $totevent=$pw[4]; 
      $flag=2;  
      next; 
    }
    if($flag==2 && $fp=~/\[1\]/)
    {
      $par2=$pw[1]; 
      $flag=3;  
      next; 
    }
    if($flag==3 && $fp=~/\[1\]/)
    {
      $mean=$pw[1]; 
      $flag=4;  
      next; 
    }
    if($flag==4 && $fp=~/\[1\]/)
    {
      $sd=$pw[1]; 
      $list[$cn][0]=$par1;
      $list[$cn][1]=$par2;
      $list[$cn][2]=$event;
      $list[$cn][3]=$totevent;
      $list[$cn][4]=$mean;
      $list[$cn][5]=$sd;
      $list[$cn][6]=sprintf("%0.4f",$sd/$mean);
      $cn++;

      undef $par2;
      undef $event;
      undef $totevent;
      undef $mean;
      undef $sd;
      $flag=1;  
      next; 
    }
    if($flag==1 && $fp=~/\-+/)
    {
       undef $par1; 
       $flag=0;
    }
   
    undef(@pw); 
  }
  close(FP);
 
  for($k=0;$k<7;$k++)
  {
    for($l=0;$l<$cn;$l++)
    {
      print WR "$list[$l][$k]\t";
    }
    print WR "\n";
  }
  print WR "\n";
  undef(@list);
  
  undef $no;
  undef $lno; 
  undef $infile;
  undef $outfile;
  undef $dr;
  undef(@arr);
  undef(@re);
}
close(WR);
close(FCS); 

system "rm FCSLIST";

