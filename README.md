# koeken
Linear Discriminant Analysis (LEfSe) on A Longitudinal Microbial Dataset.

### Introduction
Koeken is a tool developed to speed up the process of analyzing microbial community data with LEfSe, when starting with a QIIME .biom OTU table. A typical workflow of processing a .biom file into a lefse-ready format is as follows.

1. Run summarize_taxa.py on the biom table with your mapping file to create a relative abundance table with all sample metadata.
2. Open excel and remove the metadata columns you do not need and leave only metadata for class and subclass
3. Convert ; to | within the taxonomy names
4. Upload to Galaxy-Lefse and analyze (format, run, graph)

Now this workflow may seem simple, but when you have longitudinal dataset with 12 different timepoint and 4 different groups, its becomes a very tedious task. You will need to repeat this process over and over manually subsetting the groups. This is where koeken comes into play. koeken is just a simple fancy for loop which iterates over all the different groups found in you Time variable and runs all the processing steps found above.

##### Features
Koeken has some additional features built in that allow you to subset you data in various ways. Not only can you iterate over different timepoints with the ```--split``` parameter, but you can choose to compare whichever groups you want to with the ```--compare``` command. If you have 4 different groups, but are only interested in running LEfSe with 2 or 3 groups, you can list the names of the groups, and koeken will subset the tables and iterate over each day with only those specified groups.  
  
#### Outputs
The outputs from koeken are many different text files. These files are normal output from a LEfSe analysis, but there are one for each timepoint. Within the one output folder there are 2 addional folders. One for the output from running summarize_taxa.py and the other is from running LEfSe. They are named summarize_taxa_L# and lefse_output. The summarize_taxa folder is made up of relative abundance plots for each timepoint and includes metadata for the chosen class/sublcass. The lefse_output has the files for the two stages of LEfSe analysis; formatting and running. These files can be used on the galaxy server. Specifically, the files wihin run_lefse can be uploaded to the galaxy browser and used to plot the cladogram and bar charts. You will just need to choose ```lefse_res``` as the type of file when uploading.
  
Output file structure:   
output_folder/
  ---summarize_taxa
  ---lefse_output
    ---format_lefse
    ---run_lefse


#### Credits
Big thank you to the QIIME group:  
"QIIME allows analysis of high-throughput community sequencing data"  
    J Gregory Caporaso, Justin Kuczynski, Jesse Stombaugh, Kyle Bittinger, Frederic D Bushman, Elizabeth K Costello, Noah Fierer, Antonio Gonzalez Pena, Julia K Goodrich, Jeffrey I Gordon, Gavin A Huttley, Scott T Kelley, Dan Knights, Jeremy E Koenig, Ruth E Ley, Catherine A Lozupone, Daniel McDonald, Brian D Muegge, Meg Pirrung, Jens Reeder, Joel R Sevinsky, Peter J Turnbaugh, William A Walters, Jeremy Widmann, Tanya Yatsunenko, Jesse Zaneveld and Rob Knight; Nature Methods, 2010;   doi:10.1038/nmeth.f.303  

Big thank you to the LEfSe group: 
"Metagenomic biomarker discovery and explanation"  
Nicola Segata, Jacques Izard, Levi Waldron, Dirk Gevers, Larisa Miropolsky, Wendy S Garrett, and Curtis Huttenhower  
Genome Biology, 12:R60, 2011  

### Installing Dependencies
Koeken is heavily dependent upon the QIIME and LEfSe tools. Because koeken is basicaly a wrapper for the connection between these two tools, it requires the depencies found in both. If you have MacQIIME already installed you have everything you need for the QIIME side. LEfSe on the other hand requires a few R packages. If you dont already have R, you will need to go install it first before running the commands below. (). After you have R installed you can run the command below, which will install the necessary dependencies required by LEfSe. First you will need to open a new Terminal window. Then copy and paste the two commands below. You may see some errors regarding the 'splines' and 'stats4' packages. You can ignore these warnings. Once the command is run and everything is properly installed, exit out of the Terminal window and open a new one to install Koeken.

```shell
R --quiet
install.packages(c('splines','stats4','survival','mvtnorm','modeltools','coin','MASS'), repos = "http://cran.us.r-project.org")
```  
  
  
### Installing Koeken
Because Koeken needs scripts found within the QIIME package, it is easiest to use when you are in a MacQIIME session. If you have MacQIIME installed, you must first initialize it before installing Koeken, so that it is placed in the correct location. If you do not have macqiime installed, you can still run koeken. You will just need to have the scripts available in your path, or if you dont have qiime installed and you also dont have MacQIIME, installing koeken will install qiime for you, but be aware that it can take quite a bit of time to install all of qiime and it's dependencies using 'pip'.  
```shell
macqiime
```
  
Once you have initialized a new macqiime session, you can now go ahead and install koeken! Run the command below while in a macqiime session and it should install relatively quickly. To verify the installation was completed, you can run the second command found below which runs the scripts' help command. If you see options/parameters, you're good to go! 
```shell
pip install https://github.com/twbattaglia/koeken/zipball/master
koeken.py --help
```
  
### Simple Example
In this example we have a typical wokflow using koeken. We have an otu table in .biom format, a mapping text file, and an output folder to place the results. The ```--class``` parameter corresponds to the column name which describes the different groups in your mapping file. The ```--split``` parameter corresponds to the column name which describes the different timepoints in your data.  
```shell
koeken.py --input otu_table.biom --output output_folder/ --map mapping_file.txt --class Treatment --split Day 
```


### Usage
```shell
usage: koeken.py [-h] [-v] -i INPUT -o OUTPUT -m MAP [-l {2,3,4,5,6,7}] -cl
                 CLASSID [-sc SUBCLASSID] [-su SUBJECTID] [-p P_CUTOFF]
                 [-e LDA_CUTOFF] [-str {0,1}] [-c COMPARE [COMPARE ...]] -sp
                 SPLIT

Performs Linear Discriminant Analysis (LEfSe) on A Longitudinal Dataset.

optional arguments:
  -h, --help            show this help message and exit
  -v, --version         show program's version number and exit
  -i INPUT, --input INPUT
                        Location of the OTU Table for main analysis. (Must be
                        .biom format)
  -o OUTPUT, --output OUTPUT
                        Location of the folder to place all resulting files.
                        If folder does not exist, the program will create it.
  -m MAP, --map MAP     Location of the mapping file associated with OTU Table.
  -l {2,3,4,5,6,7}, --level {2,3,4,5,6,7}
                        Level for which to use for summarize_taxa.py. [default
                        = 6]
  -cl CLASSID, --class CLASSID
                        Location of the OTU Table for main analysis. (Must be
                        .biom format)
  -sc SUBCLASSID, --subclass SUBCLASSID
                        Directory to place all the files.
  -su SUBJECTID, --subject SUBJECTID
                        Only change if your Sample-ID column names differs.
                        [default] = #SampleID.
  -p P_CUTOFF, --pval P_CUTOFF
                        Change alpha value for the Anova test (default 0.05)
  -e LDA_CUTOFF, --effect LDA_CUTOFF
                        Change the cutoff for logarithmic LDA score (default
                        2.0).
  -str {0,1}, --strict {0,1}
                        Change the strictness of the comparisons. Can be
                        changed to less strict(1). [default = 0](more-strict).
  -c COMPARE [COMPARE ...], --compare COMPARE [COMPARE ...]
                        Which groups should be kept to be compared against one
                        another. [default = all factors]
  -sp SPLIT, --split SPLIT
                        The name of the timepoint variable in you mapping
                        file. This variable will be used to split the table
                        for each value in this variable.
```
