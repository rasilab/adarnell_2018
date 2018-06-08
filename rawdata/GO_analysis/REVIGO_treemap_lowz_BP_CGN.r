

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
revigo.data <- rbind(c("GO:0001580","detection of chemical stimulus involved in sensory perception of bitter taste",0.225,3.2007,0.898,0.000,"detection of chemical stimulus involved in sensory perception of bitter taste"),
c("GO:0070498","interleukin-1-mediated signaling pathway",0.115,1.9722,0.911,0.169,"detection of chemical stimulus involved in sensory perception of bitter taste"),
c("GO:0048010","vascular endothelial growth factor receptor signaling pathway",0.537,1.7027,0.920,0.258,"detection of chemical stimulus involved in sensory perception of bitter taste"),
c("GO:0038093","Fc receptor signaling pathway",1.431,1.4142,0.800,0.327,"detection of chemical stimulus involved in sensory perception of bitter taste"),
c("GO:0002521","leukocyte differentiation",2.689,1.3915,0.789,0.383,"detection of chemical stimulus involved in sensory perception of bitter taste"),
c("GO:0001906","cell killing",0.629,1.3757,0.989,0.000,"cell killing"),
c("GO:0007006","mitochondrial membrane organization",0.600,3.2366,0.844,0.000,"mitochondrial membrane organization"),
c("GO:0010256","endomembrane system organization",3.353,1.4860,0.939,0.156,"mitochondrial membrane organization"),
c("GO:0061024","membrane organization",7.299,1.5229,0.935,0.220,"mitochondrial membrane organization"),
c("GO:0000045","autophagosome assembly",0.410,1.9842,0.924,0.249,"mitochondrial membrane organization"),
c("GO:0007033","vacuole organization",0.964,1.9590,0.923,0.271,"mitochondrial membrane organization"),
c("GO:0007005","mitochondrion organization",3.756,1.5449,0.913,0.333,"mitochondrial membrane organization"),
c("GO:0048284","organelle fusion",1.321,1.4052,0.921,0.347,"mitochondrial membrane organization"),
c("GO:0070328","triglyceride homeostasis",0.190,1.8935,0.944,0.000,"triglyceride homeostasis"),
c("GO:0090559","regulation of membrane permeability",0.433,1.8416,0.948,0.218,"triglyceride homeostasis"),
c("GO:1990542","mitochondrial transmembrane transport",0.450,3.4949,0.831,0.000,"mitochondrial transmembrane transport"),
c("GO:0006897","endocytosis",4.212,2.2457,0.839,0.203,"mitochondrial transmembrane transport"),
c("GO:0015696","ammonium transport",0.560,1.7602,0.870,0.209,"mitochondrial transmembrane transport"),
c("GO:0030705","cytoskeleton-dependent intracellular transport",0.825,1.6712,0.836,0.219,"mitochondrial transmembrane transport"),
c("GO:0006839","mitochondrial transport",1.829,2.1433,0.873,0.243,"mitochondrial transmembrane transport"),
c("GO:0017038","protein import",1.875,1.3427,0.812,0.244,"mitochondrial transmembrane transport"),
c("GO:0033036","macromolecule localization",16.624,1.4013,0.866,0.263,"mitochondrial transmembrane transport"),
c("GO:0071705","nitrogen compound transport",4.472,1.6962,0.860,0.278,"mitochondrial transmembrane transport"),
c("GO:0051641","cellular localization",14.905,1.3480,0.868,0.353,"mitochondrial transmembrane transport"),
c("GO:0071702","organic substance transport",15.586,1.3686,0.835,0.358,"mitochondrial transmembrane transport"),
c("GO:0043162","ubiquitin-dependent protein catabolic process via the multivesicular body sorting pathway",0.104,2.1409,0.954,0.002,"ubiquitin-dependent protein catabolism via the multivesicular body sorting pathway"),
c("GO:0006501","C-terminal protein lipidation",0.421,1.5382,0.944,0.145,"ubiquitin-dependent protein catabolism via the multivesicular body sorting pathway"),
c("GO:0006360","transcription from RNA polymerase I promoter",0.358,1.3317,0.965,0.157,"ubiquitin-dependent protein catabolism via the multivesicular body sorting pathway"),
c("GO:0006356","regulation of transcription from RNA polymerase I promoter",0.185,1.3260,0.943,0.199,"ubiquitin-dependent protein catabolism via the multivesicular body sorting pathway"),
c("GO:0071545","inositol phosphate catabolic process",0.063,1.9555,0.918,0.278,"ubiquitin-dependent protein catabolism via the multivesicular body sorting pathway"),
c("GO:0006457","protein folding",1.396,2.6308,0.987,0.002,"protein folding"),
c("GO:0006914","autophagy",2.643,1.9003,0.987,0.002,"autophagy"),
c("GO:0050999","regulation of nitric-oxide synthase activity",0.260,1.8502,0.940,0.024,"regulation of nitric-oxide synthase activity"),
c("GO:1902188","positive regulation of viral release from host cell",0.035,1.5724,0.946,0.024,"positive regulation of viral release from host cell"),
c("GO:0050918","positive chemotaxis",0.323,1.5650,0.944,0.312,"positive regulation of viral release from host cell"),
c("GO:0051132","NK T cell activation",0.058,2.0794,0.851,0.028,"NK T cell activation"),
c("GO:0046838","phosphorylated carbohydrate dephosphorylation",0.058,2.0794,0.944,0.029,"phosphorylated carbohydrate dephosphorylation"),
c("GO:0006637","acyl-CoA metabolic process",0.565,1.9876,0.891,0.216,"phosphorylated carbohydrate dephosphorylation"),
c("GO:0033127","regulation of histone phosphorylation",0.075,1.7455,0.880,0.220,"phosphorylated carbohydrate dephosphorylation"),
c("GO:0002068","glandular epithelial cell development",0.138,2.4425,0.876,0.030,"glandular epithelial cell development"),
c("GO:0031018","endocrine pancreas development",0.260,1.8959,0.905,0.126,"glandular epithelial cell development"),
c("GO:0060711","labyrinthine layer development",0.271,1.7625,0.919,0.230,"glandular epithelial cell development"),
c("GO:0010623","programmed cell death involved in cell development",0.075,1.7455,0.901,0.284,"glandular epithelial cell development"),
c("GO:0055025","positive regulation of cardiac muscle tissue development",0.231,1.5471,0.848,0.312,"glandular epithelial cell development"),
c("GO:0042574","retinal metabolic process",0.069,1.8450,0.955,0.030,"retinal metabolism"),
c("GO:0006099","tricarboxylic acid cycle",0.179,1.3710,0.950,0.156,"retinal metabolism"),
c("GO:0042157","lipoprotein metabolic process",1.235,1.3896,0.974,0.052,"lipoprotein metabolism"),
c("GO:0046425","regulation of JAK-STAT cascade",0.837,1.9183,0.843,0.066,"regulation of JAK-STAT cascade"),
c("GO:1900103","positive regulation of endoplasmic reticulum unfolded protein response",0.075,1.7455,0.866,0.258,"regulation of JAK-STAT cascade"),
c("GO:0097296","activation of cysteine-type endopeptidase activity involved in apoptotic signaling pathway",0.075,1.7455,0.830,0.317,"regulation of JAK-STAT cascade"),
c("GO:0060192","negative regulation of lipase activity",0.087,1.5724,0.957,0.327,"regulation of JAK-STAT cascade"),
c("GO:0097696","STAT cascade",0.987,2.0685,0.901,0.341,"regulation of JAK-STAT cascade"),
c("GO:0006953","acute-phase response",0.271,1.8055,0.946,0.067,"acute-phase response"),
c("GO:0032660","regulation of interleukin-17 production",0.127,2.6180,0.886,0.082,"regulation of interleukin-17 production"));

stuff <- data.frame(revigo.data);
names(stuff) <- revigo.names;

stuff$abslog10pvalue <- as.numeric( as.character(stuff$abslog10pvalue) );
stuff$freqInDbPercent <- as.numeric( as.character(stuff$freqInDbPercent) );
stuff$uniqueness <- as.numeric( as.character(stuff$uniqueness) );
stuff$dispensability <- as.numeric( as.character(stuff$dispensability) );

# by default, outputs to a PDF file
cairo_pdf( file="../../figures/revigo_treemap_lowz_BP_cgn.pdf", width = 6, height = 4); # width and height are in inches

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

