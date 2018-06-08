

# A treemap R script produced by the REVIGO server at http://revigo.irb.hr/
# If you found REVIGO useful in your work, please cite the following reference:
# Supek F et al. "REVIGO summarizes and visualizes long lists of Gene Ontology
# terms" PLoS ONE 2011. doi:10.1371/journal.pone.0021800

# author: Anton Kratz <anton.kratz@gmail.com>, RIKEN Omics Science Center, Functional Genomics Technology Team, Japan
# created: Fri, Nov 02, 2012  7:25:52 PM
# last change: Fri, Nov 09, 2012  3:20:01 PM

# -----------------------------------------------------------------------------
# If you don't have the treemap package installed, uncomment the following line:
# install.packages( "treemap" );
library(treemap) 								# treemap package by Martijn Tennekes

# Set the working directory if necessary
# setwd("C:/Users/username/workingdir");

# --------------------------------------------------------------------------
# Here is your data from REVIGO. Scroll down for plot configuration options.

revigo.names <- c("term_ID","description","freqInDbPercent","abslog10pvalue","uniqueness","dispensability","representative");
revigo.data <- rbind(c("GO:0001227","transcriptional repressor activity, RNA polymerase II transcription regulatory region sequence-specific binding",0.039,1.7020,0.973,0.000,"transcriptional repressor activity, RNA polymerase II transcription regulatory region sequence-specific binding"),
c("GO:0004871","signal transducer activity",2.778,5.3565,0.987,0.000,"signal transducer activity"),
c("GO:0004872","receptor activity",2.694,7.8239,0.695,0.000,"receptor activity"),
c("GO:0005085","guanyl-nucleotide exchange factor activity",0.244,2.8601,0.961,0.000,"guanyl-nucleotide exchange factor activity"),
c("GO:0005198","structural molecule activity",3.268,4.9208,0.987,0.000,"structural molecule activity"),
c("GO:0005200","structural constituent of cytoskeleton",0.079,4.2596,0.972,0.000,"structural constituent of cytoskeleton"),
c("GO:0097493","structural molecule activity conferring elasticity",0.002,2.9208,0.974,0.368,"structural constituent of cytoskeleton"),
c("GO:0008307","structural constituent of muscle",0.006,2.1643,0.973,0.392,"structural constituent of cytoskeleton"),
c("GO:0005215","transporter activity",8.494,1.9115,0.988,0.000,"transporter activity"),
c("GO:0005261","cation channel activity",0.299,11.2441,0.628,0.000,"cation channel activity"),
c("GO:0005509","calcium ion binding",0.967,7.6990,0.937,0.000,"calcium ion binding"),
c("GO:0005543","phospholipid binding",0.270,1.4516,0.940,0.138,"calcium ion binding"),
c("GO:0042166","acetylcholine binding",0.003,2.0453,0.954,0.207,"calcium ion binding"),
c("GO:0043169","cation binding",15.599,4.0506,0.930,0.224,"calcium ion binding"),
c("GO:0070405","ammonium ion binding",0.018,2.1630,0.950,0.236,"calcium ion binding"),
c("GO:0030898","actin-dependent ATPase activity",0.002,4.9586,0.975,0.000,"actin-dependent ATPase activity"),
c("GO:0008296","3'-5'-exodeoxyribonuclease activity",0.006,1.9307,0.974,0.113,"actin-dependent ATPase activity"),
c("GO:0000146","microfilament motor activity",0.004,4.1487,0.974,0.238,"actin-dependent ATPase activity"),
c("GO:0003774","motor activity",0.399,3.6990,0.972,0.321,"actin-dependent ATPase activity"),
c("GO:0060089","molecular transducer activity",2.707,7.8239,0.987,0.000,"molecular transducer activity"),
c("GO:0008146","sulfotransferase activity",0.093,4.5086,0.976,0.016,"sulfotransferase activity"),
c("GO:0016782","transferase activity, transferring sulfur-containing groups",0.427,3.8861,0.975,0.198,"sulfotransferase activity"),
c("GO:0004842","ubiquitin-protein transferase activity",0.352,2.6716,0.975,0.222,"sulfotransferase activity"),
c("GO:0019787","ubiquitin-like protein transferase activity",0.378,2.3288,0.975,0.223,"sulfotransferase activity"),
c("GO:0042165","neurotransmitter binding",0.003,3.2218,0.964,0.034,"neurotransmitter binding"),
c("GO:1901338","catecholamine binding",0.004,2.6840,0.959,0.035,"catecholamine binding"),
c("GO:0071813","lipoprotein particle binding",0.004,2.2118,0.938,0.035,"lipoprotein particle binding"),
c("GO:0001540","beta-amyloid binding",0.007,2.9626,0.954,0.037,"beta-amyloid binding"),
c("GO:0060090","binding, bridging",0.052,1.7433,0.959,0.043,"binding, bridging"),
c("GO:0005516","calmodulin binding",0.072,6.2076,0.825,0.044,"calmodulin binding"),
c("GO:0034185","apolipoprotein binding",0.002,2.1445,0.854,0.333,"calmodulin binding"),
c("GO:0031690","adrenergic receptor binding",0.002,2.9208,0.837,0.333,"calmodulin binding"),
c("GO:0050431","transforming growth factor beta binding",0.003,2.0453,0.845,0.340,"calmodulin binding"),
c("GO:0030506","ankyrin binding",0.003,5.9586,0.843,0.341,"calmodulin binding"),
c("GO:0043548","phosphatidylinositol 3-kinase binding",0.005,1.6803,0.848,0.349,"calmodulin binding"),
c("GO:0043236","laminin binding",0.006,1.8557,0.847,0.350,"calmodulin binding"),
c("GO:0097110","scaffold protein binding",0.007,1.4484,0.845,0.355,"calmodulin binding"),
c("GO:0017147","Wnt-protein binding",0.009,3.2291,0.844,0.360,"calmodulin binding"),
c("GO:0043394","proteoglycan binding",0.011,1.5262,0.837,0.365,"calmodulin binding"),
c("GO:0019838","growth factor binding",0.034,1.9570,0.832,0.392,"calmodulin binding"),
c("GO:0030674","protein binding, bridging",0.045,1.4806,0.830,0.399,"calmodulin binding"),
c("GO:0005539","glycosaminoglycan binding",0.107,3.0362,0.952,0.045,"glycosaminoglycan binding"),
c("GO:0005540","hyaluronic acid binding",0.011,2.4776,0.949,0.182,"glycosaminoglycan binding"),
c("GO:0033218","amide binding",0.413,1.3684,0.954,0.051,"amide binding"),
c("GO:0044877","macromolecular complex binding",0.740,2.6073,0.952,0.054,"macromolecular complex binding"),
c("GO:0043565","sequence-specific DNA binding",2.222,2.4908,0.921,0.074,"sequence-specific DNA binding"),
c("GO:0043167","ion binding",33.492,3.5229,0.937,0.106,"sequence-specific DNA binding"),
c("GO:0008143","poly(A) binding",0.008,1.7001,0.950,0.184,"sequence-specific DNA binding"),
c("GO:0003676","nucleic acid binding",21.226,1.6799,0.930,0.197,"sequence-specific DNA binding"),
c("GO:0097159","organic cyclic compound binding",41.137,1.5290,0.936,0.262,"sequence-specific DNA binding"),
c("GO:1990837","sequence-specific double-stranded DNA binding",0.279,1.8642,0.918,0.375,"sequence-specific DNA binding"),
c("GO:0003690","double-stranded DNA binding",0.508,1.4740,0.929,0.399,"sequence-specific DNA binding"));

stuff <- data.frame(revigo.data);
names(stuff) <- revigo.names;

stuff$abslog10pvalue <- as.numeric( as.character(stuff$abslog10pvalue) );
stuff$freqInDbPercent <- as.numeric( as.character(stuff$freqInDbPercent) );
stuff$uniqueness <- as.numeric( as.character(stuff$uniqueness) );
stuff$dispensability <- as.numeric( as.character(stuff$dispensability) );

# by default, outputs to a PDF file
cairo_pdf( file="../../figures/revigo_treemap_hiz_MF_cgn.pdf", width = 6, height = 4); # width and height are in inches

# check the tmPlot command documentation for all possible parameters - there are a lot more
tmPlot(
	stuff,
	index = c("representative","description"),
	vSize = "abslog10pvalue",
	type = "categorical",
	vColor = "representative",
	title = "REVIGO Gene Ontology treemap",
	inflate.labels = FALSE,      # set this to TRUE for space-filling group labels - good for posters
	lowerbound.cex.labels = 0,   # try to draw as many labels as possible (still, some small squares may not get a label)
	bg.labels = "#CCCCCCAA",     # define background color of group labels
												       # "#CCCCCC00" is fully transparent, "#CCCCCCAA" is semi-transparent grey, NA is opaque
	position.legend = "none"
)

dev.off()
