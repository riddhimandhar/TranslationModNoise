Maximum Liklihood approach for estimation of burst parameters from single-cell RNA-seq data


burst_kinetics_ML.py - estimates the parameters of bursty mode of trancription i.e switching 'on' rate - Kon, switching 'off' rate- Koff, systhesis rate of mRNA - βm,  transcriptional efficieny- burst size (BS)


Installation
Following python libraries are required for running burst_kinetics_ML.py in ipython console

a) pandas
b) numpy
c) scipy
d) joblib

Usage

usage: burst_kinetics_ML.py  'file_name.csv'

Input 
A file containing UMI counts in 'csv' format

Output
A csv file containing all the burst estimates for each gene and a boolean entry ‘predict’ indicating whether a gene qualified a basic filtering step based on the quality of the inference procedure. 
