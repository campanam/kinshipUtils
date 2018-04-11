# kinshipUtils
Michael G. Campana, 2018  
Smithsonian Conservation Biology Institute  
Contact: <campanam@si.edu>  

R code to perform bootstrapping analysis and generate triangular kinship matrices from [SNPRelate](https://github.com/zhengxwen/SNPRelate) ([Zheng et al. 2012. Bioinformatics. 28: 3326-3328.](https://academic.oup.com/bioinformatics/article/28/24/3326/245844)) Identity-by-Descent results.  

## License  
The software is made available under the Smithsonian Institution [terms of use](https://www.si.edu/termsofuse).  

## Installation  
In the terminal:  
`git clone https://github.com/campanam/kinshipUtils`  
In your R environment:  
`source("/path/to/kinshipUtils.R")`  

## Usage  
### bootstrap.kinship  
bootstrap.kinship is a wrapper function around snpgdsLDpruning, snpgdsIBDMoM and snpgdsIBDMLE. Unless specified below, all parameters correspond to the parameters in those SNPRelate functions.  
`ibdmethod`: "MoM" specifies to use the snpgdsIBDMoM function. "MLE" specifies to use the snpgdsIBDMLE function.  
`mlemethod`: Method to perform MLE if snpgdsIBDMLE specified. Passes the value to the method parameter of snpgdsIBDMLE.  
`resample`: Number of bootstrap replicates. Default is 100.  
`ldmethod`: Method to perform LD pruning. Passes the value to the method parameter of snpgdsLDpruning.  

bootstrap.kinship generates a table with the following columns:  
Sample1&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;Sample2&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;KinshipRep1&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;...&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;KinshipRepN&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;MeanKinship&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;StandardError&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;LowCI&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;HighCI

MeanKinship and StandardError list the mean kinship estimate and the corresponding standard errors for each sample pair. LowCI and HighCI give the lower and upper bounds of the 95% Confidence interval respectively.  

### write.kinship.matrix  
write.kinship.matrix generates triangular matrices for each sample pair from the bootstrap.kinship results table. Results are written as text files in the current working directory.  
`x`: bootstrap.kinship results file to process.  
`meanfile`: Name of output file to include only mean kinship values.  
`cifile`: Name of output file to include mean and confidence interval values.  
`sep`: Separator for values. Default is a comma (for CSV).  
`digits`: Number of digits for values rounding. Default is 3.  
