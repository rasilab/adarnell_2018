

# A plotting R script produced by the REVIGO server at http://revigo.irb.hr/
# If you found REVIGO useful in your work, please cite the following reference:
# Supek F et al. "REVIGO summarizes and visualizes long lists of Gene Ontology
# terms" PLoS ONE 2011. doi:10.1371/journal.pone.0021800


# --------------------------------------------------------------------------
# If you don't have the ggplot2 package installed, uncomment the following line:
# install.packages( "ggplot2" );
library( ggplot2 );
# --------------------------------------------------------------------------
# If you don't have the scales package installed, uncomment the following line:
# install.packages( "scales" );
library( scales );
library ( ggrepel );
# --------------------------------------------------------------------------
# Here is your data from REVIGO. Scroll down for plot configuration options.

revigo.names <- c("term_ID","description","frequency_%","plot_X","plot_Y","plot_size","log10_p_value","uniqueness","dispensability");
revigo.data <- rbind(c("GO:0005654","nucleoplasm", 1.372, 6.055,-3.800, 5.130,-9.8539,0.419,0.000),
c("GO:0031974","membrane-enclosed lumen", 2.742, 4.104, 4.914, 5.431,-6.7959,0.976,0.000),
c("GO:0032991","macromolecular complex",14.008,-1.143, 2.697, 6.139,-2.2291,0.979,0.000),
c("GO:0043226","organelle",20.788,-6.348, 0.880, 6.311,-5.7212,0.981,0.000),
c("GO:1902494","catalytic complex", 3.734,-4.270,-4.952, 5.565,-3.0706,0.770,0.000),
c("GO:0030496","midbody", 0.040,-3.267, 4.429, 3.590,-1.3292,0.902,0.051),
c("GO:0098858","actin-based cell projection", 0.041,-0.532, 5.446, 3.611,-2.5058,0.833,0.051),
c("GO:0031252","cell leading edge", 0.086,-4.779, 2.543, 3.929,-1.6664,0.898,0.054),
c("GO:0030427","site of polarized growth", 0.091, 1.800, 5.031, 3.950,-1.9208,0.897,0.055),
c("GO:0005829","cytosol", 2.553, 3.906,-4.836, 5.400,-5.5086,0.731,0.167),
c("GO:0044424","intracellular part",35.654, 5.759,-0.593, 6.545,-9.6990,0.742,0.222),
c("GO:0042629","mast cell granule", 0.003, 4.641,-6.970, 2.526,-2.4330,0.716,0.259),
c("GO:0008305","integrin complex", 0.031,-3.192,-6.652, 3.489,-1.7249,0.748,0.298),
c("GO:0098636","protein complex involved in cell adhesion", 0.032,-4.853,-3.375, 3.501,-1.5696,0.832,0.298),
c("GO:0005856","cytoskeleton", 1.542, 3.728,-2.585, 5.181,-3.4949,0.577,0.305),
c("GO:0017053","transcriptional repressor complex", 0.042,-5.223,-4.464, 3.618,-1.6617,0.830,0.306),
c("GO:0005622","intracellular",41.175, 6.635, 1.094, 6.608,-9.2076,0.821,0.321),
c("GO:0005732","small nucleolar ribonucleoprotein complex", 0.072,-2.811,-4.228, 3.852,-1.3194,0.704,0.321),
c("GO:0033391","chromatoid body", 0.003,-0.049,-5.400, 2.516,-1.7307,0.634,0.322),
c("GO:1902493","acetyltransferase complex", 0.152,-4.543,-5.755, 4.175,-1.5174,0.770,0.346),
c("GO:0044444","cytoplasmic part",12.659, 4.979,-2.608, 6.095,-3.1805,0.716,0.347),
c("GO:0048471","perinuclear region of cytoplasm", 0.135, 3.377,-7.198, 4.123,-1.7926,0.786,0.351),
c("GO:0030140","trans-Golgi network transport vesicle", 0.037, 4.898,-5.582, 3.563,-1.9020,0.598,0.361),
c("GO:0044448","cell cortex part", 0.204, 2.343,-6.840, 4.303,-1.5643,0.755,0.366),
c("GO:0044464","cell part",52.385, 6.047, 1.420, 6.712,-1.3686,0.824,0.378));

one.data <- data.frame(revigo.data);
names(one.data) <- revigo.names;
one.data <- one.data [(one.data$plot_X != "null" & one.data$plot_Y != "null"), ];
one.data$plot_X <- as.numeric( as.character(one.data$plot_X) );
one.data$plot_Y <- as.numeric( as.character(one.data$plot_Y) );
one.data$plot_size <- as.numeric( as.character(one.data$plot_size) );
one.data$log10_p_value <- as.numeric( as.character(one.data$log10_p_value) );
one.data$frequency <- as.numeric( as.character(one.data$frequency) );
one.data$uniqueness <- as.numeric( as.character(one.data$uniqueness) );
one.data$dispensability <- as.numeric( as.character(one.data$dispensability) );
#head(one.data);


# --------------------------------------------------------------------------
# Names of the axes, sizes of the numbers and letters, names of the columns,
# etc. can be changed below

p1 <- ggplot( data = one.data );
p1 <- p1 + geom_point( aes( plot_X, plot_Y, colour = log10_p_value, size = plot_size), alpha = I(1) ) + scale_size_area();
p1 <- p1 + scale_colour_gradientn( colours = c("blue", "white"), limits = c( min(one.data$log10_p_value), max(one.data$log10_p_value)) );
p1 <- p1 + geom_point( aes(plot_X, plot_Y, size = plot_size), shape = 21, fill = "transparent", colour = I (alpha ("black", 1) )) + scale_size_area();
p1 <- p1 + scale_size( range=c(2, 12)) + theme_bw(); # + scale_fill_gradientn(colours = heat_hcl(7), limits = c(-300, 0) );
ex <- one.data [ one.data$dispensability < 0.15, ]; 
p1 <- p1 + geom_text_repel( data = ex, aes(plot_X, plot_Y, label = description), colour = I(alpha("black", 0.85)), size = 2.5 );
p1 <- p1 + labs (y = "semantic space x", x = "semantic space y");
p1 <- p1 + theme(legend.key = element_blank(),legend.text = element_text(size=8), axis.title = element_text(size=8), axis.text = element_text(size=8), legend.title = element_text(size=8)) ;
one.x_range = max(one.data$plot_X) - min(one.data$plot_X);
one.y_range = max(one.data$plot_Y) - min(one.data$plot_Y);
p1 <- p1 + xlim(min(one.data$plot_X)-one.x_range/5,max(one.data$plot_X)+one.x_range/5);
p1 <- p1 + ylim(min(one.data$plot_Y)-one.y_range/12,max(one.data$plot_Y)+one.y_range/20);



# --------------------------------------------------------------------------
# Output the plot to screen

#p1;

# Uncomment the line below to also save the plot to a file.
# The file type depends on the extension (default=pdf).

ggsave("../../figures/revigo_lowz_CC_CGN.pdf", device = cairo_pdf, width = 6, height = 4, units = "in", dpi = 300);
