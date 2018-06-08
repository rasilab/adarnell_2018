

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
revigo.data <- rbind(c("GO:0006813","potassium ion transport",0.476,7.4202,0.941,0.000,"potassium ion transport"),
c("GO:0032879","regulation of localization",0.726,5.1308,0.831,0.268,"potassium ion transport"),
c("GO:0048519","negative regulation of biological process",1.984,5.0757,0.836,0.299,"potassium ion transport"),
c("GO:0048518","positive regulation of biological process",1.744,3.6990,0.838,0.332,"potassium ion transport"),
c("GO:0065008","regulation of biological quality",3.395,5.0000,0.834,0.345,"potassium ion transport"),
c("GO:0006811","ion transport",5.344,4.0177,0.955,0.350,"potassium ion transport"),
c("GO:0007154","cell communication",7.219,3.9586,0.960,0.000,"cell communication"),
c("GO:0007155","cell adhesion",0.544,4.6990,0.980,0.000,"cell adhesion"),
c("GO:0007610","behavior",0.170,4.3098,0.994,0.000,"behavior"),
c("GO:0009611","response to wounding",0.127,5.4559,0.940,0.000,"response to wounding"),
c("GO:0071495","cellular response to endogenous stimulus",0.402,5.3279,0.935,0.303,"response to wounding"),
c("GO:0070848","response to growth factor",0.150,4.0458,0.929,0.307,"response to wounding"),
c("GO:0009719","response to endogenous stimulus",0.526,4.3768,0.936,0.343,"response to wounding"),
c("GO:0022610","biological adhesion",0.550,4.6198,0.994,0.000,"biological adhesion"),
c("GO:0023052","signaling",6.765,4.5528,0.995,0.000,"signaling"),
c("GO:0032501","multicellular organismal process",2.373,12.6383,0.994,0.000,"multicellular organismal process"),
c("GO:0032502","developmental process",2.812,13.7959,0.994,0.000,"developmental process"),
c("GO:0040007","growth",0.317,4.4437,0.994,0.000,"growth"),
c("GO:0040011","locomotion",0.997,4.4202,0.994,0.000,"locomotion"),
c("GO:0044699","single-organism process",46.569,4.0555,0.997,0.000,"single-organism process"),
c("GO:0044707","single-multicellular organism process",1.829,19.4318,0.540,0.000,"single-multicellular organism process"),
c("GO:0044763","single-organism cellular process",27.536,7.1487,0.868,0.154,"single-multicellular organism process"),
c("GO:0048856","anatomical structure development",2.540,15.1871,0.529,0.000,"anatomical structure development"),
c("GO:0050896","response to stimulus",12.210,3.3565,0.995,0.000,"response to stimulus"),
c("GO:0051703","intraspecies interaction between organisms",0.017,2.9508,0.991,0.000,"intraspecies interaction between organisms"),
c("GO:0065007","biological regulation",20.498,3.6778,0.995,0.000,"biological regulation"),
c("GO:0071840","cellular component organization or biogenesis",8.568,4.0809,0.995,0.000,"cellular component organization or biogenesis"),
c("GO:0000183","chromatin silencing at rDNA",0.016,14.2147,0.637,0.064,"chromatin silencing at rDNA"),
c("GO:0045103","intermediate filament-based process",0.010,6.0000,0.891,0.104,"chromatin silencing at rDNA"),
c("GO:0007267","cell-cell signaling",0.407,6.6778,0.819,0.131,"chromatin silencing at rDNA"),
c("GO:0030048","actin filament-based movement",0.021,5.9208,0.862,0.134,"chromatin silencing at rDNA"),
c("GO:2000179","positive regulation of neural precursor cell proliferation",0.009,3.4949,0.829,0.149,"chromatin silencing at rDNA"),
c("GO:0050878","regulation of body fluid levels",0.074,5.7212,0.858,0.158,"chromatin silencing at rDNA"),
c("GO:0030029","actin filament-based process",0.398,3.0555,0.862,0.169,"chromatin silencing at rDNA"),
c("GO:0006928","movement of cell or subcellular component",0.973,6.2676,0.852,0.184,"chromatin silencing at rDNA"),
c("GO:0032881","regulation of polysaccharide metabolic process",0.013,2.6968,0.869,0.214,"chromatin silencing at rDNA"),
c("GO:0040029","regulation of gene expression, epigenetic",0.130,5.3010,0.845,0.257,"chromatin silencing at rDNA"),
c("GO:0043113","receptor clustering",0.010,3.4949,0.785,0.271,"chromatin silencing at rDNA"),
c("GO:0050808","synapse organization",0.070,7.3279,0.782,0.304,"chromatin silencing at rDNA"),
c("GO:0071880","adenylate cyclase-activating adrenergic receptor signaling pathway",0.003,5.6576,0.803,0.316,"chromatin silencing at rDNA"),
c("GO:0030198","extracellular matrix organization",0.060,4.0132,0.781,0.332,"chromatin silencing at rDNA"),
c("GO:0043062","extracellular structure organization",0.061,4.0000,0.784,0.332,"chromatin silencing at rDNA"),
c("GO:0032200","telomere organization",0.133,6.8239,0.822,0.351,"chromatin silencing at rDNA"),
c("GO:0006357","regulation of transcription from RNA polymerase II promoter",1.273,2.7959,0.787,0.360,"chromatin silencing at rDNA"),
c("GO:0031344","regulation of cell projection organization",0.123,3.6383,0.655,0.366,"chromatin silencing at rDNA"),
c("GO:0051290","protein heterotetramerization",0.004,11.8539,0.839,0.388,"chromatin silencing at rDNA"),
c("GO:0044708","single-organism behavior",0.111,4.1805,0.911,0.076,"single-organism behavior"),
c("GO:0006029","proteoglycan metabolic process",0.028,5.4202,0.932,0.083,"proteoglycan metabolism"),
c("GO:1903510","mucopolysaccharide metabolic process",0.016,4.8539,0.952,0.306,"proteoglycan metabolism"));

stuff <- data.frame(revigo.data);
names(stuff) <- revigo.names;

stuff$abslog10pvalue <- as.numeric( as.character(stuff$abslog10pvalue) );
stuff$freqInDbPercent <- as.numeric( as.character(stuff$freqInDbPercent) );
stuff$uniqueness <- as.numeric( as.character(stuff$uniqueness) );
stuff$dispensability <- as.numeric( as.character(stuff$dispensability) );

# by default, outputs to a PDF file
cairo_pdf( file="../../figures/revigo_treemap_hiz_BP_cgn.pdf", width = 6, height = 4); # width and height are in inches

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
