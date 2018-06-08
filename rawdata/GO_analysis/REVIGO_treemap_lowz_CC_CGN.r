

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
revigo.data <- rbind(c("GO:0030140","trans-Golgi network transport vesicle",0.037,3.5686,0.461,0.000,"trans-Golgi network transport vesicle"),
c("GO:0098798","mitochondrial protein complex",0.289,3.5376,0.314,0.320,"trans-Golgi network transport vesicle"),
c("GO:0005773","vacuole",0.455,1.3061,0.539,0.391,"trans-Golgi network transport vesicle"),
c("GO:0031301","integral component of organelle membrane",0.214,2.0788,0.352,0.357,"trans-Golgi network transport vesicle"),
c("GO:0042629","mast cell granule",0.003,1.7744,0.633,0.199,"trans-Golgi network transport vesicle"),
c("GO:0034045","pre-autophagosomal structure membrane",0.018,1.6165,0.644,0.218,"trans-Golgi network transport vesicle"),
c("GO:0005739","mitochondrion",2.156,1.7607,0.489,0.373,"trans-Golgi network transport vesicle"),
c("GO:0097346","INO80-type complex",0.056,1.8459,0.454,0.386,"trans-Golgi network transport vesicle"),
c("GO:0001891","phagocytic cup",0.005,1.3249,0.892,0.034,"phagocytic cup"),
c("GO:0031975","envelope",2.324,1.3161,0.866,0.053,"envelope"),
c("GO:0033202","DNA helicase complex",0.093,1.6165,0.712,0.089,"DNA helicase complex"),
c("GO:0034358","plasma lipoprotein particle",0.010,1.5809,0.775,0.208,"DNA helicase complex"),
c("GO:0032994","protein-lipid complex",0.011,1.4923,0.839,0.210,"DNA helicase complex"));

stuff <- data.frame(revigo.data);
names(stuff) <- revigo.names;

stuff$abslog10pvalue <- as.numeric( as.character(stuff$abslog10pvalue) );
stuff$freqInDbPercent <- as.numeric( as.character(stuff$freqInDbPercent) );
stuff$uniqueness <- as.numeric( as.character(stuff$uniqueness) );
stuff$dispensability <- as.numeric( as.character(stuff$dispensability) );

# by default, outputs to a PDF file
cairo_pdf( file="../../figures/revigo_treemap_lowz_CC_cgn.pdf", width=6, height=4 ) # width and height are in inches

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
