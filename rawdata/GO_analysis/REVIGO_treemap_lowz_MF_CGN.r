

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
revigo.data <- rbind(c("GO:0008565","protein transporter activity",0.537,1.9884,0.850,0.000,"protein transporter activity"),
c("GO:1904680","peptide transmembrane transporter activity",0.040,1.4532,0.860,0.388,"protein transporter activity"),
c("GO:0072349","modified amino acid transmembrane transporter activity",0.110,1.3501,0.859,0.390,"protein transporter activity"),
c("GO:0009055","electron carrier activity",0.607,1.5427,0.955,0.000,"electron carrier activity"),
c("GO:0019003","GDP binding",0.341,3.5686,0.774,0.000,"GDP binding"),
c("GO:0033038","bitter taste receptor activity",0.133,4.0223,0.939,0.000,"bitter taste receptor activity"),
c("GO:0008527","taste receptor activity",0.173,3.1612,0.939,0.370,"bitter taste receptor activity"),
c("GO:0042056","chemoattractant activity",0.168,2.3420,0.955,0.000,"chemoattractant activity"),
c("GO:0055102","lipase inhibitor activity",0.104,1.4848,0.955,0.000,"lipase inhibitor activity"),
c("GO:0061650","ubiquitin-like protein conjugating enzyme activity",0.173,1.9935,0.889,0.000,"ubiquitin-like protein conjugating enzyme activity"),
c("GO:0097200","cysteine-type endopeptidase activity involved in execution phase of apoptosis",0.075,1.7333,0.893,0.102,"ubiquitin-like protein conjugating enzyme activity"),
c("GO:0050664","oxidoreductase activity, acting on NAD(P)H, oxygen as acceptor",0.092,1.6430,0.928,0.104,"ubiquitin-like protein conjugating enzyme activity"),
c("GO:0052689","carboxylic ester hydrolase activity",0.815,1.4669,0.911,0.265,"ubiquitin-like protein conjugating enzyme activity"),
c("GO:0005126","cytokine receptor binding",1.543,2.4881,0.809,0.006,"cytokine receptor binding"),
c("GO:0050662","coenzyme binding",1.052,1.5378,0.886,0.006,"coenzyme binding"),
c("GO:0048037","cofactor binding",1.554,1.4008,0.949,0.007,"cofactor binding"),
c("GO:0036312","phosphatidylinositol 3-kinase regulatory subunit binding",0.058,2.0670,0.939,0.028,"phosphatidylinositol 3-kinase regulatory subunit binding"),
c("GO:0019864","IgG binding",0.064,2.0670,0.920,0.029,"IgG binding"),
c("GO:0005070","SH3/SH2 adaptor activity",0.329,1.8681,0.901,0.033,"SH3/SH2 adaptor activity"),
c("GO:0070491","repressing transcription factor binding",0.329,1.8289,0.936,0.033,"repressing transcription factor binding"),
c("GO:0051082","unfolded protein binding",0.670,1.3296,0.935,0.036,"unfolded protein binding"));

stuff <- data.frame(revigo.data);
names(stuff) <- revigo.names;

stuff$abslog10pvalue <- as.numeric( as.character(stuff$abslog10pvalue) );
stuff$freqInDbPercent <- as.numeric( as.character(stuff$freqInDbPercent) );
stuff$uniqueness <- as.numeric( as.character(stuff$uniqueness) );
stuff$dispensability <- as.numeric( as.character(stuff$dispensability) );

# by default, outputs to a PDF file
cairo_pdf( file="../../figures/revigo_treemap_lowz_MF_cgn.pdf", width=6, height=4 ) # width and height are in inches

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
