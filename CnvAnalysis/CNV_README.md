Clustering of samples by CNV
============================

We use the software package [FACETS](https://github.com/mskcc/facets) to
perform copy number variation (CNV) analysis on the whole genome sequencing
data of our samples. It consists of several steps, starting from the aligned
BAM files. Although we cannot include the BAM files in this repository because
they are from human subjects and are under controlled access, we will describe
all the analysis steps, and provide results when there is no longer a privacy
concern, and the data size has been reduced sufficiently to be included in
GitHub repositories.

Briefly, our CNV analysis and sample clustering pipeline involves the following
steps

  1. extraction of sequening depth at known polymorphism sites
  2. CNV analysis using FACETS
      * we provide our analysis in script `facets.r`
  3. FACETS results post processing (described in manuscript)
      * we provide the result of this step in file `binnedLogR.txt`
  4. Clustering and results visualization
      * we provide script `plotylyCnvHeatmap.py` to reproduce the CNV heatmap
      * we provide script `CNVcluster.R` to reproduce the clustering results

Extraction of sequencing depth
------------------------------

We used the command line utility `snp-pileup` as part of FACETS to perform this
step. We refer you to the [FACETS](https://github.com/mskcc/facets) repository
for more information.

CNV analysis using FACETS
-------------------------

Although FACETS is a self contained R package, it does not by default save the
output we need for downstream analysis. We developed a wrapper script
`facets.r` to perform the analysis and save the results. It can be run in R
batch mode

`$ R CMD BATCH facets.r`

Outputs from this step include the FACETS CNV figures, the `Rout` file that
contains the ploidy and tumor purity estimations, as well as a txt file that
includes _logRatio_ values which we use to visualize the copy number events
using heatmap

FACETS results post-processing
------------------------------

CNV events identified from different tumor samples need to be integrated for
further analysis. Briefly, we binned the segments produced by FACETS into 1Mb
windows, and compute the _logR_ value for each bin. This is an interactive,
manual analysis step that we describe in the methods section. We provide the
result of this step in file `binnedLogR.txt` that is the input to the following
steps

CNV events visualization via heatmap
------------------------------------

We provide the script `plotlyCnvHeatmap.py` to reproduce the CNV heatmap showed
in **Figure 2A**. To run the script, first install the python package
`chart_studio` via your python package management system, e.g. `$ pip install
chart_studio`. Then invoke the script in the following manner

`$ cat binnedLogR.txt | python3 plotlyCnvHeatmap.py 1000000`

The parameter `1000000` is the window size for the bins. 

Sample clustering based on CNV events
-------------------------------------

To reproduce the result on sample clustering via CNV events, we provide the R
script `CNVcluster.R`. The output is a cluster dendrogram as shown in **Fig
S3B**. 
