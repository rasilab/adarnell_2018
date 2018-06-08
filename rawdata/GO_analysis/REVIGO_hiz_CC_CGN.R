

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
revigo.data <- rbind(c("GO:0000786","nucleosome", 0.191, 4.174,-4.195, 4.274,-30.0000,0.555,0.000),
c("GO:0005576","extracellular region", 2.375,-4.044, 3.201, 5.369,-3.5686,0.974,0.000),
c("GO:0031012","extracellular matrix", 0.275,-2.777,-0.834, 4.433,-10.5528,0.805,0.000),
c("GO:0032991","macromolecular complex",14.008,-0.788, 2.186, 6.139,-1.4087,0.977,0.000),
c("GO:0098589","membrane region", 0.121,-5.700,-2.155, 4.077,-2.1972,0.908,0.000),
c("GO:0099080","supramolecular complex", 0.540, 2.873, 5.850, 4.726,-8.0000,0.974,0.000),
c("GO:0043209","myelin sheath", 0.049,-4.507,-5.988, 3.683,-2.7645,0.926,0.044),
c("GO:0043005","neuron projection", 0.190, 5.399, 5.453, 4.271,-2.2464,0.868,0.049),
c("GO:0031226","intrinsic component of plasma membrane", 1.299,-5.275, 0.848, 5.106,-8.1549,0.794,0.058),
c("GO:0071944","cell periphery",11.583,-3.555,-4.148, 6.057,-2.1278,0.899,0.095),
c("GO:0016021","integral component of membrane",55.868,-1.067, 5.266, 6.740,-1.9905,0.936,0.104),
c("GO:0060198","clathrin-sculpted vesicle", 0.000, 3.409,-8.281, 1.176,-1.7129,0.812,0.131),
c("GO:0032994","protein-lipid complex", 0.011, 5.446, 0.222, 3.025,-1.3655,0.867,0.220),
c("GO:0031463","Cul3-RING ubiquitin ligase complex", 0.034, 7.062, 0.085, 3.525,-1.6772,0.796,0.237),
c("GO:0043235","receptor complex", 0.115, 7.361,-1.077, 4.054,-2.5482,0.850,0.259),
c("GO:0005852","eukaryotic translation initiation factor 3 complex", 0.117, 6.327,-3.688, 4.060,-2.0620,0.747,0.260),
c("GO:0044815","DNA packaging complex", 0.216, 5.628,-1.244, 4.327,-30.0000,0.844,0.273),
c("GO:0034702","ion channel complex", 0.200, 7.493,-2.254, 4.293,-5.9586,0.774,0.274),
c("GO:0043228","non-membrane-bounded organelle", 8.408, 0.371,-7.839, 5.918,-2.8041,0.768,0.285),
c("GO:0032993","protein-DNA complex", 0.418, 5.971,-2.036, 4.614,-26.9208,0.838,0.291),
c("GO:0034399","nuclear periphery", 0.066, 1.192,-7.485, 3.811,-1.4292,0.726,0.324),
c("GO:0044425","membrane part",57.394,-0.269, 5.512, 6.752,-1.8225,0.946,0.327),
c("GO:1990351","transporter complex", 0.885, 6.679,-2.304, 4.940,-4.6198,0.830,0.331),
c("GO:0005882","intermediate filament", 0.090, 0.073,-6.877, 3.949,-17.5850,0.588,0.351));

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
p1 <- p1 + geom_point( aes(plot_X, plot_Y, size = plot_size), shape = 21, fill = "transparent", colour = I (alpha ("black", 0.85) )) + scale_size_area();
p1 <- p1 + scale_size( range=c(0.5, 3)) + theme_bw(); # + scale_fill_gradientn(colours = heat_hcl(7), limits = c(-300, 0) );
ex <- one.data [ one.data$dispensability < 0.25, ]; 
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

ggsave("../../figures/revigo_hiz_CC_CGN.pdf", device = cairo_pdf, width = 4, height = 3, units = "in", dpi = 300)