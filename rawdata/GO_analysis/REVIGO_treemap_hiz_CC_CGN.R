

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
revigo.data <- rbind(c("GO:0000788","nuclear nucleosome",0.004,14.8539,0.682,0.000,"nuclear nucleosome"),
c("GO:0060198","clathrin-sculpted vesicle",0.000,1.6992,0.827,0.140,"nuclear nucleosome"),
c("GO:0044815","DNA packaging complex",0.216,8.3979,0.868,0.208,"nuclear nucleosome"),
c("GO:0032994","protein-lipid complex",0.011,1.3472,0.888,0.221,"nuclear nucleosome"),
c("GO:0032982","myosin filament",0.003,6.2518,0.660,0.229,"nuclear nucleosome"),
c("GO:0031463","Cul3-RING ubiquitin ligase complex",0.034,2.1192,0.806,0.239,"nuclear nucleosome"),
c("GO:0032580","Golgi cisterna membrane",0.039,2.4535,0.705,0.245,"nuclear nucleosome"),
c("GO:0043235","receptor complex",0.115,4.4685,0.873,0.262,"nuclear nucleosome"),
c("GO:0032993","protein-DNA complex",0.418,8.2924,0.862,0.291,"nuclear nucleosome"),
c("GO:0005852","eukaryotic translation initiation factor 3 complex",0.117,2.9208,0.744,0.302,"nuclear nucleosome"),
c("GO:0005856","cytoskeleton",1.542,1.9401,0.692,0.326,"nuclear nucleosome"),
c("GO:1902495","transmembrane transporter complex",0.864,7.7212,0.804,0.330,"nuclear nucleosome"),
c("GO:1990351","transporter complex",0.885,7.4949,0.855,0.357,"nuclear nucleosome"),
c("GO:0042555","MCM complex",0.057,1.9300,0.859,0.361,"nuclear nucleosome"),
c("GO:0043228","non-membrane-bounded organelle",8.408,1.8422,0.785,0.374,"nuclear nucleosome"),
c("GO:0005576","extracellular region",2.375,3.2218,0.986,0.000,"extracellular region"),
c("GO:0016020","membrane",61.592,5.3565,0.995,0.000,"membrane"),
c("GO:0030054","cell junction",0.445,3.0862,0.986,0.000,"cell junction"),
c("GO:0031012","extracellular matrix",0.275,20.8539,0.844,0.000,"extracellular matrix"),
c("GO:0032991","macromolecular complex",14.008,2.5654,0.988,0.000,"macromolecular complex"),
c("GO:0045202","synapse",0.299,3.6778,0.986,0.000,"synapse"),
c("GO:0099080","supramolecular complex",0.540,7.7696,0.986,0.000,"supramolecular complex"),
c("GO:0099512","supramolecular fiber",0.538,7.6198,0.802,0.000,"supramolecular fiber"),
c("GO:0044459","plasma membrane part",2.405,10.8239,0.833,0.044,"plasma membrane part"),
c("GO:0071944","cell periphery",11.583,7.6576,0.901,0.105,"plasma membrane part"),
c("GO:0016021","integral component of membrane",55.868,4.0269,0.956,0.110,"plasma membrane part"),
c("GO:0044425","membrane part",57.394,3.9586,0.963,0.327,"plasma membrane part"),
c("GO:0042383","sarcolemma",0.024,1.4901,0.879,0.334,"plasma membrane part"),
c("GO:0098589","membrane region",0.121,2.6676,0.919,0.045,"membrane region"),
c("GO:0032154","cleavage furrow",0.012,1.5253,0.928,0.048,"cleavage furrow"),
c("GO:0098862","cluster of actin-based cell projections",0.034,1.4161,0.933,0.053,"cluster of actin-based cell projections"),
c("GO:0043209","myelin sheath",0.049,1.3644,0.932,0.054,"myelin sheath"),
c("GO:0044297","cell body",0.087,1.5807,0.930,0.057,"cell body"),
c("GO:0043005","neuron projection",0.190,3.9208,0.827,0.062,"neuron projection"),
c("GO:0009986","cell surface",0.241,2.4449,0.925,0.063,"cell surface"),
c("GO:0097458","neuron part",0.320,3.5229,0.924,0.065,"neuron part"),
c("GO:0016528","sarcoplasm",0.019,2.8386,0.889,0.071,"sarcoplasm"),
c("GO:0042995","cell projection",1.050,2.3107,0.918,0.074,"cell projection"));

stuff <- data.frame(revigo.data);
names(stuff) <- revigo.names;

stuff$abslog10pvalue <- as.numeric( as.character(stuff$abslog10pvalue) );
stuff$freqInDbPercent <- as.numeric( as.character(stuff$freqInDbPercent) );
stuff$uniqueness <- as.numeric( as.character(stuff$uniqueness) );
stuff$dispensability <- as.numeric( as.character(stuff$dispensability) );

# by default, outputs to a PDF file
cairo_pdf( file="../../figures/revigo_treemap_hiz_CC_cgn.pdf", width = 6, height = 4); # width and height are in inches

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
