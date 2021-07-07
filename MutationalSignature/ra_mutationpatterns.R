#if (!requireNamespace("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")

#BiocManager::install("MutationalPatterns")
library(MutationalPatterns)
library(rtracklayer)
library(BSgenome)
require(ggplot2)
require(reshape);
ref_genome <- "BSgenome.Hsapiens.UCSC.hg19"
library(ref_genome, character.only = TRUE)

vcf_files = list.files(pattern = ".vcf")
vcf_files
#samples_names need to be the same order as vcf names
sample_names = vcf_files

vcfs = read_vcfs_as_granges(vcf_files,sample_names,ref_genome)

########### 96 mutational profile ##########
mut_mat <- mut_matrix(vcf_list = vcfs, ref_genome = ref_genome)
plot_96_profile(mut_mat,condensed = TRUE)

########## Mutational signatures ##########
apply(mut_mat,2,sum)
mut_mat <- mut_mat + 0.0001
library("NMF")
estimate <- nmf(mut_mat, rank=2:5, method="brunet", nrun=10, seed=123456)
sp_url = paste("https://cancer.sanger.ac.uk/cancergenome/assets/","signatures_probabilities.txt",sep="")
cancer_signatures = read.table(sp_url,sep="\t",header=TRUE)
new_order = match(row.names(mut_mat), cancer_signatures$Somatic.Mutation.Type)
cancer_signatures = cancer_signatures[as.vector(new_order),]
row.names(cancer_signatures) = cancer_signatures$Somatic.Mutation.Type
cancer_signatures = as.matrix(cancer_signatures[,4:33])

cos_sim_samples_signatures = cos_sim_matrix(mut_mat, cancer_signatures)
fit_res <- fit_to_signatures(mut_mat, cancer_signatures)
select <- which(rowSums(fit_res$contribution) > 300)
plot_contribution(fit_res$contribution[select,],cancer_signatures[,select],coord_flip = FALSE,mode = "absolute")

plot_contribution_pie = function(contribution,
                                 samplename,
                                 signatures,
                                 index=c(),
                                 palette=c())
{
  # optional subsetting if index parameter is provided
  if(length(index > 0)){contribution = contribution[,index]}
  
  # These variables will be available at run-time, but not at compile-time.
  # To avoid compiling trouble, we initialize them to NULL.
  Contribution = NULL
  Signature = NULL
  
  m_contribution = melt(contribution)
  colnames(m_contribution) = c("Contribution")
  m_contribution$Signature = c("1)Signature 1","2)Signature 2","3)Signature 3",
                               "4)Signature 5","5)Signature 8","6)signature 9",
                               "7)Signature 13","8)Signature 18")
  rownames(m_contribution) = c(1,2,3,4,5,6,7,8)
  print(m_contribution)
  
  plot = ggplot(m_contribution,
                aes(x = samplename,
                    y = Contribution,
                    fill = Signature)) +
    geom_bar(width=1,size=1,stat="identity", colour="white")  +
    # ylabel
    labs(x = "", y = "Relative contribution") +
    # white background
    theme_bw() +
    # no gridlines
    theme(panel.grid.minor.x=element_blank(),
          panel.grid.major.x=element_blank()) +
    theme(panel.grid.minor.y=element_blank(),
          panel.grid.major.y=element_blank()) +
    theme(axis.line = element_blank(),
          axis.text = element_blank(),
          axis.ticks = element_blank()) +
    coord_polar("y")
  
  # Allow custom color palettes.
  if (length(palette) > 0)
    plot = plot + scale_fill_manual(name="Signature", values=palette)
  else
    plot = plot + scale_fill_discrete(name="Signature")
  
  return(plot)
}
palette = c("#66c2a5",
            "#f46d43",
            "#3288bd",
            "#abdda4",
            "#f1b6da",
            "#4d9221",
            "#fdae61",
            "#762a83")


plot_contribution_pie(fit_res$contribution[select,][,12],
                      samplename="Truncal",
                      cancer_signatures[,select],
                      palette=palette)
plot_contribution_pie(fit_res$contribution[select,][,1],
                      samplename="G1 private",
                      cancer_signatures[,select],
                      palette=palette)
plot_contribution_pie(fit_res$contribution[select,][,2],
                      samplename="G1 shared",
                      cancer_signatures[,select],
                      palette=palette)
plot_contribution_pie(fit_res$contribution[select,][,3],
                      samplename="G2 private",
                      cancer_signatures[,select],
                      palette=palette)
plot_contribution_pie(fit_res$contribution[select,][,4],
                      samplename="G2 shared",
                      cancer_signatures[,select],
                      palette=palette)
plot_contribution_pie(fit_res$contribution[select,][,5],
                      samplename="G3 private",
                      cancer_signatures[,select],
                      palette=palette)
plot_contribution_pie(fit_res$contribution[select,][,6],
                      samplename="G3 shared",
                      cancer_signatures[,select],
                      palette=palette)
plot_contribution_pie(fit_res$contribution[select,][,7],
                      samplename="G3 shared by >=2 samples",
                      cancer_signatures[,select],
                      palette=palette)
plot_contribution_pie(fit_res$contribution[select,][,8],
                      samplename="G3 shared by >=5 samples",
                      cancer_signatures[,select],
                      palette=palette)
plot_contribution_pie(fit_res$contribution[select,][,9],
                      samplename="G3 shared by >=7 samples",
                      cancer_signatures[,select],
                      palette=palette)
plot_contribution_pie(fit_res$contribution[select,][,10],
                      samplename="G4 private",
                      cancer_signatures[,select],
                      palette=palette)
plot_contribution_pie(fit_res$contribution[select,][,11],
                      samplename="G4 shared",
                      cancer_signatures[,select],
                      palette=palette)
