# ISC_final_project
This is the GitHub repository of Wenhan Chang's ISC final project, which includes the R codes as well as a brief README file.

## R packages required
ggplot2

## Running the codes
This GitHub repository contains mainly two sets of files -- data files (seq_data_\*.txt) and R codes (\*.R). And in order to test the codes, you need to pull all files into your local laptop (Mac or Linux system only).

After downloading all files (excluding files in parallel_version_alpha folder), you can simply see the result by running run.R file, which will load data from 5 files (seq_data_1.txt to seq_data_5.txt), import functions from functions_codes.R file, and calculate the likelihood of different selection coefficients. The log-likelihoods of different selection coefficients are saved in a data.frame called result, which is then used to plot the data using ggplot2.

If everything works well (which it should, as I have tested it on both Mac and CentOS systems), you shall be able to see the following plot, indicating that the most probably selection coefficient is s = 0.8.

![Correct Results](https://github.com/wenhan-uchicago/ISC_final_project/blob/master/result.png)

## Introcution
This is a Maximum-Likelihood estimator for deducing the most probable selection coefficient from Evolution-Resequencing experiment (ER) data [1]. As a simple ER experiment example, suppose we are selecting 100 fruitflies for their body size for 100 generations. And as well all know, body size is a classic example of Quantitative Trait that is contributed by hundreds of Quantitattive Trait Loci (QTL). Therefore, unlike classic Mendelian genetics, the specific genes that are responsible for body size could not be easily discovered. One way to do so is by QTL-mapping, which is tedious and less-efficient; another way is to conduct an Evolution-Resequencing (ER) experiment, since those genes responsible for body size would increase their allelic frequencies after ~100 generations of selection. 

However, since the population size is normally small in ER experiments [2], the effect of genetic drift, gene frequency change caused by randomness, cannot be omitted. Past method relied on using Chi-squared test to discover genomic regions that have undergone dramatic differentiation after selection [2], which was proven to be useless by simulation data. Therefore, the availability of cheap Next-Generation Sequencing (NGS) data urges for an sensitive way to detect loci under selection, which I present here by using a model based maximum likelihood estimator.

## Materials and Methods

### Input files
The five data files (seq_data_1.txt to seq_data_5.txt) were generated by another ad hoc R script that follows the classic Wright-Fisher model, with a user-defined selection coefficient s = 0.8. And the first column means the sequence depth from NGS,  second column means absolute frequency of allele A, and the third column means the number of generations. For example, the first line in seq_data_1.txt is "97      52      1", which means that in Genration 1, the sequence depth is 97 and the frequency of allele A is 52. There are 5 different files, because in ER experiment it is required to have multiple replicates. Therefore, these five data files represent five independent replications under exactly the same conditions.

### R Functions
In functions_codes.R file, there are several funcitons I defined. 

First, there are two functions for reading the data and processing them. The output is a data.frame that have all input data.

Second, there are ~5 functions that calculate the likelihood of selection coefficient given the data, i.e. P(Data|Selection Coefficient). The output of these functions is the log-likelihood of given selection coefficient.

Third, there are two "hidden" functions that use Metropolis–Hastings algorithm [3] to sample the most probable selection coefficient. These two functions are not used since the advantage of using MH algorithm is not large enough and is usually more time-consuming.

In run.R file, I firstly import all functions in functions_codes.R file and then read all input data. Then, the output of log-likehoods are saved into a data.frame called result, which is then used for plotting the data with ggplot2.

## Results
As a simple demonstration of the power of this estimator, it is obvious from the figure above that it correctly estimated the user-defined selection coefficient, 0.8. Given the limited time, more sophisticated benchmarking results are not shown, but this model-based estimator clearly surpass Chi-square based test proposed in [2].

## Discussion
As a painful process, I started writing this as a parallel version which utilized the "snow" package of R. The legacy of which could still be found under parallel_version_alpha folder, and it works great in RCC-midway cluster. However, when trying to test it in local laptop, espically Mac, the required MPI (Message-Passing Interface) is very hard to deal with. Therefore, I have to give up writing a easy-to-use parallel version, but instead writing a simple prototype.

However, the time and efforts I spent on writing the parallel version definitely paid back, as it not only trained my ability to write workable codes but also showed me the amazing power of computer clusters. And as a quick "showoff", the following picture shows the results from a 1MB genome region (with ~25,000 markers) that was analyzed by using 48 CPU nodes and finished within ~6 hours (it would take days to finish analyzing the same amount of data without using the power of parallel). Thank you all for the fantastic teaching and I really enjoyed the learning process that taught me how to take R to the limit! Thank you!


![Results](https://github.com/wenhan-uchicago/ISC_final_project/blob/master/parallel_version_alpha/parallel_results.png)
(The above figure is a Forward-Simulation conducted by SLiM. And the only SNP that is under selection is the center-most 500 region. Therefore, it is obvious that though it did not precisely found the SNP under selection, the most probable one estimated is very close to the correct one (~20kb away). And hundreds of benchmarks have been done, which are not shown here, all support the idea that this estimator is a better option than Chi-square based test.)



## Reference
1. Kofler R, Schlötterer C. A guide for the design of evolve and resequencing studies. Molecular biology and evolution. 2014 Feb 1;31(2):474-83.
 
2.	Baldwin-Brown JG, Long AD, Thornton KR. The power to detect quantitative trait loci using resequenced, experimentally evolved populations of diploid, sexual organisms. Molecular biology and evolution. 2014;31(4):1040-55.

3. Metropolis N, Rosenbluth AW, Rosenbluth MN, Teller AH, Teller E. Equation of state calculations by fast computing machines. The journal of chemical physics. 1953 Jun 1;21(6):1087-92.
