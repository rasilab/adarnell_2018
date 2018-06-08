

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
library( ggrepel );


# --------------------------------------------------------------------------
# Here is your data from REVIGO. Scroll down for plot configuration options.

revigo.names <- c("term_ID","description","frequency_%","plot_X","plot_Y","plot_size","log10_p_value","uniqueness","dispensability");
revigo.data <- rbind(c("GO:0001227","transcriptional repressor activity, RNA polymerase II transcription regulatory region sequence-specific binding", 0.039,-1.173, 3.468, 3.743,-1.7020,0.973,0.000),
c("GO:0004871","signal transducer activity", 2.778, 6.450, 0.159, 5.592,-5.3565,0.987,0.000),
c("GO:0004872","receptor activity", 2.694,-2.210, 1.196, 5.579,-7.8239,0.695,0.000),
c("GO:0005085","guanyl-nucleotide exchange factor activity", 0.244, 4.571,-2.013, 4.537,-2.8601,0.961,0.000),
c("GO:0005198","structural molecule activity", 3.268,-4.870, 5.387, 5.663,-4.9208,0.987,0.000),
c("GO:0005200","structural constituent of cytoskeleton", 0.079, 4.788,-4.335, 4.046,-4.2596,0.972,0.000),
c("GO:0005215","transporter activity", 8.494, 5.602, 1.960, 6.078,-1.9115,0.988,0.000),
c("GO:0005261","cation channel activity", 0.299, 3.121, 0.024, 4.625,-11.2441,0.628,0.000),
c("GO:0005509","calcium ion binding", 0.967, 4.445, 4.574, 5.134,-7.6990,0.937,0.000),
c("GO:0030898","actin-dependent ATPase activity", 0.002, 1.577,-6.306, 2.389,-4.9586,0.975,0.000),
c("GO:0060089","molecular transducer activity", 2.707,-0.406,-0.021, 5.581,-7.8239,0.987,0.000),
c("GO:0008146","sulfotransferase activity", 0.093,-2.347,-6.643, 4.118,-4.5086,0.976,0.016),
c("GO:0042165","neurotransmitter binding", 0.003, 1.247, 0.693, 2.616,-3.2218,0.964,0.034),
c("GO:1901338","catecholamine binding", 0.004,-1.757,-1.979, 2.700,-2.6840,0.959,0.035),
c("GO:0071813","lipoprotein particle binding", 0.004, 6.680,-1.142, 2.744,-2.2118,0.938,0.035),
c("GO:0001540","beta-amyloid binding", 0.007,-2.713, 5.087, 3.006,-2.9626,0.954,0.037),
c("GO:0060090","binding, bridging", 0.052, 4.950, 0.254, 3.866,-1.7433,0.959,0.043),
c("GO:0005516","calmodulin binding", 0.072,-7.137,-0.241, 4.008,-6.2076,0.825,0.044),
c("GO:0005539","glycosaminoglycan binding", 0.107,-6.011,-5.398, 4.178,-3.0362,0.952,0.045),
c("GO:0033218","amide binding", 0.413, 1.484,-2.195, 4.764,-1.3684,0.954,0.051),
c("GO:0044877","macromolecular complex binding", 0.740,-2.552,-4.114, 5.018,-2.6073,0.952,0.054),
c("GO:0043565","sequence-specific DNA binding", 2.222,-1.753, 7.289, 5.495,-2.4908,0.921,0.074),
c("GO:0043167","ion binding",33.492, 2.171, 6.701, 6.673,-3.5229,0.937,0.106),
c("GO:0008296","3'-5'-exodeoxyribonuclease activity", 0.006, 0.287,-5.552, 2.915,-1.9307,0.974,0.113),
c("GO:0005543","phospholipid binding", 0.270, 1.699, 4.336, 4.579,-1.4516,0.940,0.138),
c("GO:0005540","hyaluronic acid binding", 0.011,-5.341,-5.287, 3.186,-2.4776,0.949,0.182),
c("GO:0008143","poly(A) binding", 0.008,-4.022, 7.079, 3.024,-1.7001,0.950,0.184),
c("GO:0003676","nucleic acid binding",21.226, 0.364, 7.059, 6.475,-1.6799,0.930,0.197),
c("GO:0016782","transferase activity, transferring sulfur-containing groups", 0.427,-1.699,-7.135, 4.779,-3.8861,0.975,0.198),
c("GO:0042166","acetylcholine binding", 0.003, 4.435, 3.705, 2.571,-2.0453,0.954,0.207),
c("GO:0004842","ubiquitin-protein transferase activity", 0.352,-3.205,-6.636, 4.695,-2.6716,0.975,0.222),
c("GO:0019787","ubiquitin-like protein transferase activity", 0.378,-2.588,-7.242, 4.726,-2.3288,0.975,0.223),
c("GO:0043169","cation binding",15.599, 3.083, 5.500, 6.342,-4.0506,0.930,0.224),
c("GO:0070405","ammonium ion binding", 0.018, 5.425, 3.807, 3.401,-2.1630,0.950,0.236),
c("GO:0000146","microfilament motor activity", 0.004, 2.230,-6.105, 2.755,-4.1487,0.974,0.238),
c("GO:0097159","organic cyclic compound binding",41.137, 1.625, 6.710, 6.763,-1.5290,0.936,0.262),
c("GO:0003774","motor activity", 0.399, 1.998,-6.429, 4.750,-3.6990,0.972,0.321),
c("GO:0034185","apolipoprotein binding", 0.002,-6.485, 1.007, 2.500,-2.1445,0.854,0.333),
c("GO:0031690","adrenergic receptor binding", 0.002,-6.629,-0.475, 2.520,-2.9208,0.837,0.333),
c("GO:0050431","transforming growth factor beta binding", 0.003,-7.015, 1.880, 2.667,-2.0453,0.845,0.340),
c("GO:0030506","ankyrin binding", 0.003,-7.501,-0.618, 2.691,-5.9586,0.843,0.341),
c("GO:0043548","phosphatidylinositol 3-kinase binding", 0.005,-7.019, 1.407, 2.874,-1.6803,0.848,0.349),
c("GO:0043236","laminin binding", 0.006,-7.504, 1.361, 2.890,-1.8557,0.847,0.350),
c("GO:0097110","scaffold protein binding", 0.007,-7.692, 0.172, 3.000,-1.4484,0.845,0.355),
c("GO:0017147","Wnt-protein binding", 0.009,-7.162, 0.178, 3.084,-3.2291,0.844,0.360),
c("GO:0043394","proteoglycan binding", 0.011,-6.903,-1.310, 3.195,-1.5262,0.837,0.365),
c("GO:0097493","structural molecule activity conferring elasticity", 0.002, 4.238,-4.256, 2.450,-2.9208,0.974,0.368),
c("GO:1990837","sequence-specific double-stranded DNA binding", 0.279,-2.584, 7.284, 4.595,-1.8642,0.918,0.375),
c("GO:0008307","structural constituent of muscle", 0.006, 4.854,-4.009, 2.912,-2.1643,0.973,0.392),
c("GO:0019838","growth factor binding", 0.034,-6.892, 0.567, 3.677,-1.9570,0.832,0.392),
c("GO:0030674","protein binding, bridging", 0.045,-7.408, 0.744, 3.801,-1.4806,0.830,0.399),
c("GO:0003690","double-stranded DNA binding", 0.508,-2.137, 7.459, 4.855,-1.4740,0.929,0.399));

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
p1 <- p1 + scale_colour_gradientn( colours = c("red", "white"), limits = c( min(one.data$log10_p_value), max(one.data$log10_p_value)) );
p1 <- p1 + geom_point( aes(plot_X, plot_Y, size = plot_size), shape = 21, fill = "transparent", colour = I (alpha ("black", 1) )) + scale_size_area();
p1 <- p1 + scale_size( range=c(1, 6)) + theme_bw(); # + scale_fill_gradientn(colours = heat_hcl(7), limits = c(-300, 0) );
ex <- one.data [ one.data$dispensability < 0.025, ]; 
p1 <- p1 + geom_text_repel( data = ex, aes(plot_X, plot_Y, label = description), colour = I(alpha("black", 1)), size = 2.5 );
p1 <- p1 + labs (y = "GO term semantic space, x", x = "GO term semantic space, y");
p1 <- p1 + theme(legend.key = element_blank(), legend.text = element_text(size=8), axis.title = element_text(size=8), axis.text = element_text(size=8), legend.title = element_text(size=8)) ;
one.x_range = max(one.data$plot_X) - min(one.data$plot_X);
one.y_range = max(one.data$plot_Y) - min(one.data$plot_Y);
p1 <- p1 + xlim(min(one.data$plot_X)-one.x_range/10,max(one.data$plot_X)+one.x_range/10);
p1 <- p1 + ylim(min(one.data$plot_Y)-one.y_range/20,max(one.data$plot_Y)+one.y_range/20);



# --------------------------------------------------------------------------
# Output the plot to screen

#p1;

# Uncomment the line below to also save the plot to a file.
# The file type depends on the extension (default=pdf).

ggsave("../../figures/revigo_hiz_MF_CGN.pdf", device = cairo_pdf, width = 4, height = 3, units = "in", dpi = 300)