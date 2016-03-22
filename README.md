# ISC_final_project
This is the GitHub repository of my ISC final project, which includes the R codes as well as a brief README file.

## R packages required
ggplot2

## Running the codes
This GitHub repository contains mainly two sets of files -- data files (seq_data_\*.txt) and R codes (\*.R). And in order to test the codes, you need to pull all files into your local laptop (Mac or Linux system only).

After downloading all files (excluding files in parallel_version_alpha folder), you can simply see the result by running run.R file, which will load data from 5 files (seq_data_1.txt to seq_data_5.txt), import functions from functions_codes.R file, and calculate the likelihood of different selection coefficients. The log-likelihoods of different selection coefficients are saved in a data.frame called result, which is then used to plot the data using ggplot2.

If everything works well (which it should, as I have tested it on both Mac and CentOS systems), you shall be able to see the following plot, indicating that the most probably selection coefficient is s = 0.8.

