

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
revigo.data <- rbind(c("GO:0005976","polysaccharide metabolic process", 0.906, 0.769, 3.136, 5.065,-2.3507,0.946,0.000),
c("GO:0007610","behavior", 0.170,-2.867,-0.915, 4.338,-4.9208,0.994,0.000),
c("GO:0032501","multicellular organismal process", 2.373, 0.044,-0.653, 5.483,-6.6383,0.994,0.000),
c("GO:0032502","developmental process", 2.812,-0.634,-3.229, 5.557,-4.7212,0.994,0.000),
c("GO:0034728","nucleosome organization", 0.129,-2.953, 5.794, 4.219,-21.4949,0.839,0.000),
c("GO:0043269","regulation of ion transport", 0.244, 6.083, 2.518, 4.496,-6.7959,0.749,0.000),
c("GO:0043588","skin development", 0.044,-1.757,-6.038, 3.753,-11.3768,0.657,0.000),
c("GO:0051703","intraspecies interaction between organisms", 0.017, 5.125,-7.284, 3.346,-2.9547,0.979,0.000),
c("GO:0045103","intermediate filament-based process", 0.010,-5.143,-3.468, 3.094,-6.0000,0.910,0.048),
c("GO:0072089","stem cell proliferation", 0.029,-3.628, 1.912, 3.566,-3.0132,0.939,0.051),
c("GO:0044763","single-organism cellular process",27.536,-4.073,-5.522, 6.548,-3.2147,0.894,0.091),
c("GO:0034435","cholesterol esterification", 0.003,-5.683, 0.978, 2.590,-2.3516,0.892,0.095),
c("GO:0030204","chondroitin sulfate metabolic process", 0.006,-0.366, 5.125, 2.876,-3.7447,0.860,0.098),
c("GO:0016567","protein ubiquitination", 0.523, 1.407, 6.433, 4.827,-2.2848,0.934,0.113),
c("GO:0006112","energy reserve metabolic process", 0.168,-4.973,-1.824, 4.334,-2.8697,0.883,0.120),
c("GO:0010939","regulation of necrotic cell death", 0.004, 5.288,-2.513, 2.747,-2.7825,0.805,0.165),
c("GO:0015874","norepinephrine transport", 0.001, 3.910, 6.605, 2.212,-3.4437,0.896,0.165),
c("GO:0050878","regulation of body fluid levels", 0.074, 7.396,-2.845, 3.976,-4.4089,0.806,0.191),
c("GO:0007187","G-protein coupled receptor signaling pathway, coupled to cyclic nucleotide second messenger", 0.046, 4.381,-3.753, 3.771,-5.8539,0.758,0.193),
c("GO:0032026","response to magnesium ion", 0.002, 0.652,-7.201, 2.389,-2.0283,0.938,0.199),
c("GO:0002227","innate immune response in mucosa", 0.002, 2.542,-7.086, 2.401,-5.3010,0.907,0.199),
c("GO:0040029","regulation of gene expression, epigenetic", 0.130, 7.208,-0.045, 4.224,-6.6778,0.809,0.209),
c("GO:0042053","regulation of dopamine metabolic process", 0.003, 7.758, 1.632, 2.588,-2.0283,0.833,0.219),
c("GO:0071514","genetic imprinting", 0.006, 5.852, 0.533, 2.881,-2.1124,0.834,0.231),
c("GO:0048519","negative regulation of biological process", 1.984, 7.148,-1.593, 5.406,-2.5850,0.801,0.266),
c("GO:0032879","regulation of localization", 0.726, 6.214, 1.905, 4.969,-2.2848,0.788,0.299),
c("GO:0006044","N-acetylglucosamine metabolic process", 0.055,-1.879, 3.746, 3.850,-2.4802,0.959,0.300),
c("GO:0051172","negative regulation of nitrogen compound metabolic process", 0.792, 6.628,-0.495, 5.007,-6.6576,0.705,0.326),
c("GO:0048518","positive regulation of biological process", 1.744, 7.696,-0.932, 5.350,-2.0467,0.803,0.332),
c("GO:0009611","response to wounding", 0.127, 2.034,-7.125, 4.212,-4.2840,0.931,0.336),
c("GO:0065008","regulation of biological quality", 3.395, 7.122,-1.234, 5.639,-2.4473,0.798,0.345),
c("GO:0043279","response to alkaloid", 0.017, 1.181,-7.029, 3.350,-2.0186,0.917,0.346),
c("GO:0006811","ion transport", 5.344, 5.176, 5.757, 5.836,-6.2007,0.932,0.350),
c("GO:0019433","triglyceride catabolic process", 0.007,-5.730, 0.340, 2.975,-2.3002,0.893,0.353),
c("GO:0010628","positive regulation of gene expression", 0.653, 7.194, 0.269, 4.923,-3.4559,0.766,0.379),
c("GO:0007267","cell-cell signaling", 0.407,-2.910,-3.899, 4.718,-4.8539,0.843,0.381),
c("GO:0006323","DNA packaging", 0.227,-3.558, 5.242, 4.465,-19.5376,0.869,0.382),
c("GO:1904837","beta-catenin-TCF complex assembly", 0.000,-2.297, 6.299, 1.799,-4.3372,0.893,0.387),
c("GO:0015844","monoamine transport", 0.014, 4.388, 6.269, 3.247,-3.0044,0.897,0.393));

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
ex <- one.data [ one.data$dispensability < 0.02, ]; 
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

ggsave("../../figures/revigo_hiz_BP_CGN.pdf", device = cairo_pdf, width = 5, height = 4, units = "in", dpi = 300)