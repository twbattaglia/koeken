# koeken
A Linear Discriminant Analysis (LEfSe) wrapper.  


### Installing Koeken
Because Koeken needs scripts found within the QIIME package, it is easiest to use when you are in a MacQIIME session. If you have MacQIIME installed, you ***must*** first initialize it before installing Koeken. If you do not have macqiime installed, you can still run koeken as long as you have the scripts available in your path. See http://qiime.org/install/install.html for more information.  
```shell
macqiime
```

Run the command below while in a ***MacQIIME*** session. The command will install new R packages. If you do not have R installed, please see https://cran.r-project.org/, and install the package first before running.  
```shell
pip install https://github.com/twbattaglia/koeken/zipball/master
koeken.py --help
```

### Simple Example
In this example we have a typical workflow using koeken. We have an otu-table in .biom format, a mapping text file, and an output folder to place the results. The ```--class``` parameter corresponds to the column name which describes the different groups in your mapping file. The ```--split``` parameter corresponds to the column name which describes the different timepoints in your data. ```--clade``` will produce the cladograms popular with LEfSe.

```bash
koeken.py \
--input otu_table.biom \
--output output_folder/ \
--map mapping_file.txt \
--class Treatment \
--split Day \
--clade 
```

### Introduction
Koeken is a wrapper around QIIME and LEfSe to create a more convienient and easier workflow from a QIIME OTU table to biomarker microbial species. It was specifically developed for the experiments which have more than one time point and many groups. Koeken allows the user to specify the variable which correspond to the 'Time' of a users dataset to run LEfSe on each 'Time' factor. Additionally it has the capabilities to specify particular groups to compare to one another, without any additionally subsetting. 

### Features
* Run LEfSe over many time points
* Specify groups (2+) to compare to eachother
* Automatically generate cladograms for each comparison.
* Clean GreenGenes naming conventions (e.g k__,g__)
* Find significant biomarker PICRUST Level 3 features.
* Combine the data from many timepoints into a 'heatmap-like' table. (See 'pretty_lefse.py -h' for more information)
* Modify thresholds for testing parameters (e.g effect size, p-value, strictness)
* Collapse OTU data on different levels [Default = L6(Genus)]

### Outputs
Koeken generates many intermediate files as it iterates over each time point. For each timepoint, koeken runs a QIIME command (summarize_taxa.py), LEfSe formatting and LEfSe significance tests, but each of the steps are separated into their respective folders to minimize confusion. A breakdown of the folders are shown below.
  
```
output_folder/  
|--summarize_taxa/ (Output from running QIIME commands)  
|--lefse_output/ (Root folder for storing LEfSe outputs)
|----format_lefse/ (Output from running LEfSe formatting command)
|----run_lefse/ (Output from running main LEfSe command) 

```
### Using Output Files on Galaxy-LEfSe Application
If you want to use the galaxy version of LEfSe to generate figures, you can use the provided output from koeken. During file upload you must select ```lefse_res``` for the files in the 'run_lefse/' folder and must select ```lefse_for``` for the files in the 'format_lefse/'.


### Credits
Big thank you to the QIIME group:  
"QIIME allows analysis of high-throughput community sequencing data"  
    J Gregory Caporaso, Justin Kuczynski, Jesse Stombaugh, Kyle Bittinger, Frederic D Bushman, Elizabeth K Costello, Noah Fierer, Antonio Gonzalez Pena, Julia K Goodrich, Jeffrey I Gordon, Gavin A Huttley, Scott T Kelley, Dan Knights, Jeremy E Koenig, Ruth E Ley, Catherine A Lozupone, Daniel McDonald, Brian D Muegge, Meg Pirrung, Jens Reeder, Joel R Sevinsky, Peter J Turnbaugh, William A Walters, Jeremy Widmann, Tanya Yatsunenko, Jesse Zaneveld and Rob Knight; Nature Methods, 2010;   doi:10.1038/nmeth.f.303  

Big thank you to the LEfSe group:
"Metagenomic biomarker discovery and explanation"  
Nicola Segata, Jacques Izard, Levi Waldron, Dirk Gevers, Larisa Miropolsky, Wendy S Garrett, and Curtis Huttenhower  
Genome Biology, 12:R60, 2011  


### Parameters
```bash
usage: koeken.py [-h] [-v] [-d] -i INPUT_BIOM -o OUTPUTDIR -m MAP_FP
                 [-l {2,3,4,5,6,7}] -cl CLASSID [-sc SUBCLASSID]
                 [-su SUBJECTID] [-p P_CUTOFF] [-e LDA_CUTOFF] [-str {0,1}]
                 [-c COMPARE [COMPARE ...]] -sp SPLIT [-pc]
                 [-it {png,pdf,svg}] [-dp DPI] [-pi]

Performs Linear Discriminant Analysis (LEfSe) on A Longitudinal Dataset.

optional arguments:
  -h, --help            show this help message and exit
  -v, --version         show program's version number and exit
  -d, --debug           Enable debugging
  -i INPUT_BIOM, --input INPUT_BIOM
                        Location of the OTU Table for main analysis. (Must be
                        .biom format)
  -o OUTPUTDIR, --output OUTPUTDIR
                        Location of the folder to place all resulting files.
                        If folder does not exist, the program will create it.
  -m MAP_FP, --map MAP_FP
                        Location of the mapping file associated with OTU
                        Table.
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
                        changed to less strict (1). [default = 0](more-
                        strict).
  -c COMPARE [COMPARE ...], --compare COMPARE [COMPARE ...]
                        Which groups should be kept to be compared against one
                        another. [default = all factors]
  -sp SPLIT, --split SPLIT
                        The name of the timepoint variable in you mapping
                        file. This variable will be used to split the table
                        for each value in this variable.
  -pc, --clade          Plot Lefse Cladogram for each output time point.
                        Outputs are placed in a new folder created in the
                        lefse results location.
  -it {png,pdf,svg}, --image {png,pdf,svg}
                        Set the file type for the image create when using
                        cladogram setting
  -dp DPI, --dpi DPI    Set DPI resolution for cladogram
  -pi, --picrust        Run analysis with PICRUSt biom file. Must use the
                        cateogirze by function level 3. Next updates will
                        reflect the difference levels.
```


### Features for the future
- [ ] install LEfSe files separately and much better
- [ ] add results table for a combined all time points analysis
- [ ] add more package dependency verification
- [ ] export2graphlan.py support
- [ ] add more run_lefse.py parameters
- [ ] add more plot_cladogram.py parameters
