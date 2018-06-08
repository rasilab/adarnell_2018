

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
revigo.data <- rbind(c("GO:0000988","transcription factor activity, protein binding", 0.843, 5.554, 1.893, 5.074,-4.0915,0.983,0.000),
c("GO:0003712","transcription cofactor activity", 0.250,-5.596, 3.422, 4.547,-3.3979,0.940,0.000),
c("GO:0004702","signal transducer, downstream of receptor, with serine/threonine kinase activity", 0.081, 4.227, 3.942, 4.060,-2.8996,0.893,0.000),
c("GO:0005041","low-density lipoprotein receptor activity", 0.002,-0.627,-0.171, 2.393,-1.6123,0.970,0.000),
c("GO:0005085","guanyl-nucleotide exchange factor activity", 0.244,-0.952, 2.415, 4.537,-4.6576,0.896,0.000),
c("GO:0005488","binding",55.656,-6.898, 3.476, 6.894,-3.6576,0.992,0.000),
c("GO:0070063","RNA polymerase binding", 0.069,-6.133,-4.398, 3.986,-5.0177,0.748,0.000),
c("GO:0072349","modified amino acid transmembrane transporter activity", 0.019, 1.047, 3.497, 3.431,-1.9666,0.978,0.000),
c("GO:0098772","molecular function regulator", 1.115, 1.622, 0.421, 5.196,-3.1249,0.983,0.000),
c("GO:0032452","histone demethylase activity", 0.009,-2.633, 0.972, 3.108,-2.1487,0.972,0.017),
c("GO:0050664","oxidoreductase activity, acting on NAD(P)H, oxygen as acceptor", 0.011,-5.302, 4.908, 3.199,-1.5237,0.972,0.018),
c("GO:0032451","demethylase activity", 0.015, 4.244, 1.041, 3.322,-1.5405,0.971,0.018),
c("GO:0036459","thiol-dependent ubiquitinyl hydrolase activity", 0.184, 0.850, 5.667, 4.413,-2.0550,0.898,0.021),
c("GO:0005515","protein binding", 4.410, 1.056,-4.305, 5.793,-3.0269,0.914,0.050),
c("GO:0003682","chromatin binding", 0.220,-4.028, 3.271, 4.492,-1.8758,0.934,0.056),
c("GO:0044212","transcription regulatory region DNA binding", 0.275, 6.125,-2.575, 4.587,-2.8827,0.864,0.058),
c("GO:0000287","magnesium ion binding", 1.785, 6.190,-0.117, 5.400,-1.7749,0.909,0.073),
c("GO:0050662","coenzyme binding", 3.468,-0.136,-3.710, 5.689,-2.4776,0.897,0.080),
c("GO:1901363","heterocyclic compound binding",41.115, 1.551,-7.504, 6.763,-8.1549,0.883,0.130),
c("GO:0048037","cofactor binding", 5.599,-1.928,-8.050, 5.897,-2.0414,0.912,0.138),
c("GO:0005543","phospholipid binding", 0.270, 3.702,-1.101, 4.579,-1.7126,0.915,0.146),
c("GO:0016279","protein-lysine N-methyltransferase activity", 0.087,-2.643, 5.784, 4.087,-2.5622,0.862,0.174),
c("GO:0097367","carbohydrate derivative binding",17.252,-0.472,-8.681, 6.385,-2.1433,0.898,0.196),
c("GO:0001067","regulatory region nucleic acid binding", 0.282, 5.637,-3.015, 4.598,-2.8447,0.899,0.207),
c("GO:0003676","nucleic acid binding",21.226, 3.176,-6.562, 6.475,-5.2518,0.858,0.213),
c("GO:0036094","small molecule binding",21.337, 0.429,-7.533, 6.478,-3.0605,0.895,0.214),
c("GO:0000049","tRNA binding", 0.670, 5.963,-3.622, 4.974,-1.4642,0.907,0.224),
c("GO:0016887","ATPase activity", 4.560,-0.101, 5.657, 5.808,-1.9957,0.865,0.248),
c("GO:0101005","ubiquitinyl hydrolase activity", 0.184,-0.333, 6.092, 4.413,-1.9952,0.927,0.249),
c("GO:0002161","aminoacyl-tRNA editing activity", 0.232,-1.311, 5.350, 4.515,-1.8190,0.926,0.255),
c("GO:0043168","anion binding",20.942, 2.746,-5.537, 6.470,-3.4089,0.877,0.260),
c("GO:0043167","ion binding",33.492, 0.879,-8.231, 6.673,-2.0535,0.887,0.262),
c("GO:0005548","phospholipid transporter activity", 0.056, 2.204, 3.221, 3.901,-1.5029,0.978,0.281),
c("GO:0097159","organic cyclic compound binding",41.137, 1.316,-7.805, 6.763,-7.3565,0.883,0.292),
c("GO:0003723","RNA binding", 5.283, 5.050,-4.605, 5.871,-3.3188,0.865,0.307),
c("GO:0036312","phosphatidylinositol 3-kinase regulatory subunit binding", 0.002,-6.863,-1.795, 2.405,-1.9408,0.817,0.328),
c("GO:0034185","apolipoprotein binding", 0.002,-5.936,-2.490, 2.500,-1.3687,0.815,0.332),
c("GO:0000166","nucleotide binding",20.185, 3.656,-6.380, 6.454,-4.4437,0.804,0.332),
c("GO:1901265","nucleoside phosphate binding",20.185, 3.323,-6.888, 6.454,-4.4318,0.859,0.332),
c("GO:0031994","insulin-like growth factor I binding", 0.002,-7.096,-2.170, 2.526,-1.7102,0.815,0.333),
c("GO:0043560","insulin receptor substrate binding", 0.002,-5.549,-3.683, 2.529,-1.8190,0.815,0.333),
c("GO:0051879","Hsp90 protein binding", 0.006,-5.831,-4.774, 2.935,-1.7510,0.806,0.351),
c("GO:0046875","ephrin receptor binding", 0.007,-6.344,-2.751, 3.003,-1.4697,0.797,0.355),
c("GO:0050321","tau-protein kinase activity", 0.001, 3.942, 4.301, 2.318,-1.7102,0.920,0.359),
c("GO:0019894","kinesin binding", 0.010,-6.348,-3.250, 3.157,-2.7496,0.801,0.362),
c("GO:0045309","protein phosphorylated amino acid binding", 0.013,-6.845,-3.635, 3.247,-2.2388,0.790,0.367),
c("GO:0070491","repressing transcription factor binding", 0.014,-7.179,-2.867, 3.282,-2.1046,0.784,0.369),
c("GO:0035064","methylated histone binding", 0.015,-6.195,-3.876, 3.331,-2.2449,0.797,0.371),
c("GO:0051219","phosphoprotein binding", 0.020,-6.903,-4.151, 3.447,-1.8435,0.794,0.378),
c("GO:0008013","beta-catenin binding", 0.021,-6.574,-4.248, 3.473,-1.4203,0.793,0.379),
c("GO:0008022","protein C-terminus binding", 0.035,-7.084,-3.327, 3.687,-1.8611,0.788,0.391));

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
p1 <- p1 + ylim(min(one.data$plot_Y)-one.y_range/20,max(one.data$plot_Y)+one.y_range/20);



# --------------------------------------------------------------------------
# Output the plot to screen

#p1;

# Uncomment the line below to also save the plot to a file.
# The file type depends on the extension (default=pdf).

ggsave("../../figures/revigo_lowz_MF_CGN.pdf", device = cairo_pdf, width = 6, height = 4, units = "in", dpi = 300);
