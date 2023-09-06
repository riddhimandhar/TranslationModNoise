1.  master_R.pl 

    - reads *.fcs files in the "Input_FCS" directory 
    - calls FCS_R_code.R for processing 


2. FCS_R_code.R 
  
   - processes FCS files 
   - marks ellipses of different sizes and shapes 
   - generates noise data 



3. summarize_all_data.pl 

   - summarizes data generated in the "Output" folder 