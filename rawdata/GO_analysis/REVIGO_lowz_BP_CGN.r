

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
revigo.data <- rbind(c("GO:0006996","organelle organization", 3.595,-2.604, 5.896, 5.664,-6.5376,0.769,0.000),
c("GO:0008152","metabolic process",75.387, 2.575, 3.561, 6.986,-1.5394,0.998,0.000),
c("GO:0030705","cytoskeleton-dependent intracellular transport", 0.056, 0.893, 2.107, 3.857,-2.4168,0.891,0.000),
c("GO:0032259","methylation", 3.103, 1.647,-2.509, 5.600,-1.9523,0.980,0.000),
c("GO:0043547","positive regulation of GTPase activity", 0.470, 6.026,-4.033, 4.780,-5.9208,0.787,0.000),
c("GO:0071840","cellular component organization or biogenesis", 8.568, 6.751, 4.114, 6.041,-4.7212,0.994,0.000),
c("GO:0006807","nitrogen compound metabolic process",38.744, 3.003, 0.612, 6.696,-2.5482,0.974,0.031),
c("GO:0007017","microtubule-based process", 0.658, 2.764, 5.449, 4.927,-3.1249,0.872,0.039),
c("GO:0006396","RNA processing", 3.210, 0.339,-6.714, 5.615,-3.3565,0.814,0.048),
c("GO:0045730","respiratory burst", 0.004,-4.208, 0.907, 2.732,-2.0846,0.942,0.054),
c("GO:0046838","phosphorylated carbohydrate dephosphorylation", 0.064,-4.195, 3.822, 3.915,-1.9481,0.912,0.055),
c("GO:0008340","determination of adult lifespan", 0.020,-0.912,-0.389, 3.413,-1.7169,0.820,0.060),
c("GO:0016337","single organismal cell-cell adhesion", 0.148, 6.081, 3.402, 4.280,-1.6788,0.927,0.071),
c("GO:1901360","organic cyclic compound metabolic process",30.324,-4.107,-2.127, 6.590,-2.2857,0.942,0.097),
c("GO:0008214","protein dealkylation", 0.019,-4.025,-5.959, 3.378,-2.5817,0.875,0.117),
c("GO:0046483","heterocycle metabolic process",29.664,-6.770,-0.163, 6.580,-3.2518,0.896,0.128),
c("GO:0032878","regulation of establishment or maintenance of cell polarity", 0.013, 6.109,-0.620, 3.217,-1.7352,0.792,0.176),
c("GO:0044237","cellular metabolic process",53.061,-6.475, 1.629, 6.833,-1.9234,0.917,0.176),
c("GO:0030029","actin filament-based process", 0.398, 3.622, 5.022, 4.708,-1.4693,0.876,0.177),
c("GO:0048524","positive regulation of viral process", 0.016, 6.918,-4.596, 3.321,-2.2000,0.771,0.179),
c("GO:0006349","regulation of gene expression by genetic imprinting", 0.004, 4.407,-6.176, 2.689,-1.3750,0.795,0.192),
c("GO:0007049","cell cycle", 1.885, 1.943, 5.498, 5.384,-1.3845,0.860,0.208),
c("GO:0008156","negative regulation of DNA replication", 0.032, 3.610,-5.704, 3.608,-2.8416,0.715,0.210),
c("GO:0043170","macromolecule metabolic process",39.491,-5.244,-2.078, 6.705,-2.0996,0.939,0.211),
c("GO:0051056","regulation of small GTPase mediated signal transduction", 0.218, 6.713,-2.542, 4.447,-3.8239,0.708,0.219),
c("GO:0010467","gene expression",19.671,-0.967,-7.447, 6.402,-5.0315,0.894,0.222),
c("GO:0006725","cellular aromatic compound metabolic process",29.628,-6.441, 0.399, 6.580,-2.9355,0.896,0.245),
c("GO:1902170","cellular response to reactive nitrogen species", 0.003, 7.725, 0.490, 2.539,-1.8897,0.913,0.255),
c("GO:0046834","lipid phosphorylation", 0.177, 0.235, 4.777, 4.357,-1.3982,0.850,0.258),
c("GO:0065009","regulation of molecular function", 1.726, 5.983,-4.495, 5.345,-3.6778,0.797,0.268),
c("GO:0098732","macromolecule deacylation", 0.085,-4.539,-6.343, 4.040,-1.6000,0.895,0.268),
c("GO:0071108","protein K48-linked deubiquitination", 0.015,-4.586,-5.331, 3.276,-2.2388,0.870,0.271),
c("GO:0006450","regulation of translational fidelity", 0.303, 2.601,-5.228, 4.590,-1.3750,0.689,0.288),
c("GO:0035601","protein deacylation", 0.082,-3.393,-6.334, 4.020,-1.6273,0.858,0.301),
c("GO:0043162","ubiquitin-dependent protein catabolic process via the multivesicular body sorting pathway", 0.014,-5.597,-3.834, 3.263,-1.9751,0.892,0.317),
c("GO:0048010","vascular endothelial growth factor receptor signaling pathway", 0.013, 7.184,-2.016, 3.207,-3.2676,0.761,0.336),
c("GO:0008380","RNA splicing", 0.413, 1.403,-7.564, 4.725,-2.9788,0.835,0.346),
c("GO:0048583","regulation of response to stimulus", 1.120, 6.462,-3.015, 5.158,-1.3504,0.774,0.348),
c("GO:0071545","inositol phosphate catabolic process", 0.062,-2.531, 3.361, 3.897,-1.8256,0.853,0.361),
c("GO:0043173","nucleotide salvage", 0.126,-1.085, 1.651, 4.207,-1.5303,0.799,0.372),
c("GO:0016071","mRNA metabolic process", 0.798, 0.653,-7.699, 5.010,-1.7144,0.843,0.373),
c("GO:0007028","cytoplasm organization", 0.006,-1.601, 5.864, 2.900,-3.0000,0.820,0.374),
c("GO:0006306","DNA methylation", 0.190,-1.871,-6.595, 4.387,-2.5735,0.798,0.380),
c("GO:0043954","cellular component maintenance", 0.011,-0.886, 6.066, 3.134,-1.3693,0.814,0.390),
c("GO:0090304","nucleic acid metabolic process",21.449,-0.973,-5.301, 6.440,-3.9208,0.801,0.391),
c("GO:0097190","apoptotic signaling pathway", 0.117, 5.615,-1.108, 4.177,-1.8008,0.715,0.394));

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
p1 <- p1 + geom_point( aes(plot_X, plot_Y, size = plot_size), shape = 21, fill = "transparent", colour = I (alpha ("black", 0.85) )) + scale_size_area();
p1 <- p1 + scale_size( range=c(0.5, 3)) + theme_bw(); # + scale_fill_gradientn(colours = heat_hcl(7), limits = c(-300, 0) );
ex <- one.data [ one.data$dispensability < 0.05, ]; 
p1 <- p1 + geom_text_repel( data = ex, aes(plot_X, plot_Y, label = description), colour = I(alpha("black", 1)), size = 2 );
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

ggsave("../../figures/revigo_lowz_BP_cgn.pdf", device = cairo_pdf, width = 4, height = 3, units = "in", dpi = 300);
