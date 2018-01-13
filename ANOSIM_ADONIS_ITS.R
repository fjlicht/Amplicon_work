library(ggplot2)
library(vegan)
library(dplyr)
library(scales)
library(grid)
library(reshape2)
library(phyloseq)

psY = readRDS("~/Dropbox/16S_OREI/16S_dada2_output/ITS/Year1_2_3.rds")
psY
psY = readRDS("~/Dropbox/16S_OREI/16S_dada2_output/ITS/T3_4_5.rds")

psY = readRDS("~/Dropbox/16S_OREI/16S_dada2_output/ITS/Year4.rds")


sample_data(psY)$Location_Year

set.seed(1000)
sessionInfo()
psY


get_taxa_unique(psY, "Kingdom")
get_taxa_unique(psY, "Class")
get_taxa_unique(psY, "Genus")

pslog <- transform_sample_counts(psY, function(x) log(1 + x))
sample_data(pslog)$Location_Year <- cut(sample_data(pslog)$Year_N)

rel_abundance<- t(apply(otu_table(psY), 1, function(x) x/sum(x)))
qplot(rel_abundance[, 12], geom = "histogram") + xlab("Relative Abundance in Year 1-3")

library("phyloseq"); packageVersion("phyloseq")
library("ggplot2"); packageVersion("ggplot2")
library("plyr"); packageVersion("plyr")
library("dplyr")
theme_set(theme_bw())
ITS = psY
wh0 = genefilter_sample(ITS, filterfun_sample(function(x) x > 5), A=0.5*nsamples(ITS))
ITS1 = prune_taxa(wh0, ITS)

ITS1 = transform_sample_counts(ITS1, function(x) 1E6 * x/sum(x))
phylum.sum = tapply(taxa_sums(ITS1), tax_table(ITS1)[, "Phylum"], sum, na.rm=TRUE)
top5phyla = names(sort(phylum.sum, TRUE))[1:5]
ITS1 = prune_taxa((tax_table(ITS1)[, "Phylum"] %in% top5phyla), ITS1)

GP.ord <- ordinate(ITS1, "NMDS", "bray")

p1 = plot_ordination(ITS1, GP.ord,type="species", color="Phylum", title="Taxonomy")
plot_ordination('list')
print(p1)
p1 = plot_ordination(ITS1, GP.ord,type="sites", color="Location", title="Taxonomy")
print(p1)
p1 = plot_ordination(ITS1, GP.ord,type="biplot", color="Location", title="Taxonomy")
print(p1)
p1 = plot_ordination(ITS1, GP.ord,type="split", color="Location", title="Taxonomy")
print(p1)

GP.ord <- ordinate(ITS1, "NMDS", "bray")
p1 = plot_ordination(ITS1, GP.ord, type="taxonomy", color="Phylum", title="Taxonomy")
print(p1)


p2 = plot_ordination(ITS1, GP.ord, type="sites", color="Location", shape="Year") 
p2 + geom_polygon(aes(fill=Location)) + geom_point(size=5) + ggtitle("Treatment")

##biplot
p3 = plot_ordination(ITS1, GP.ord, type="biplot", color="Location", shape="Year_N", title="biplot")
# Some stuff to modify the automatic shape scale
GP1.shape.names = get_taxa_unique(ITS1, "Phylum")
GP1.shape <- 15:(15 + length(GP1.shape.names) - 1)
names(GP1.shape) <- GP1.shape.names
GP1.shape["Year_N"] <- 16
p3 + scale_shape_manual(values=GP1.shape)


#####
#http://deneflab.github.io/MicrobeMiseq/demos/mothur_2_phyloseq.html
######
library(vegan)
set.seed(1)

# Calculate bray curtis distance matrix
ITS_bray <- phyloseq::distance(psY, method = "bray")

# make a data frame from the sample_data
sampledf <- data.frame(sample_data(psY))

# Adonis test
adonis(ITS_bray ~ Year_N, data = sampledf)
adonis(ITS_bray ~ Location_Year, data = sampledf)
adonis(ITS_bray ~ Location_Year + Treatment, data = sampledf)
adonis(ITS_bray ~ Location_Year + Treatment, data = sampledf)

# Homogeneity of dispersion test
beta <- betadisper(ITS_bray, sampledf$Location_Year)
permutest(beta)


#################################################
# Remove data points with missing metadata
ITS_not_na <- psY %>%
  subset_samples(
    !is.na(Average_FB) & 
      !is.na(CEC) &
      !is.na(pH) & 
      !is.na(Organic_Matter) & 
      !is.na(obs)
  )

bray_not_na <- phyloseq::distance(physeq = ITS_not_na, method = "bray")
sample_data(bray_not_na)$Location_Year

###remove control
# Remove data points with missing metadata
ITS_no_control <- ITS_not_na %>%
  subset_samples(drop(!Treatment=="T_8")
  )
psY <- ITS_no_control

bray_no_control <- phyloseq::distance(physeq = ITS_no_control, method = "bray")

sample_data(psY)$Location_Year

# CAP ordinate
cap_ord <- ordinate(
  physeq = ITS_no_control, 
  method = "CAP",
  distance = bray_no_control,
  formula = ~ pH + Organic_Matter + Winter_precip + CEC
)

# CAP plot
cap_plot <- plot_ordination(
  physeq = ITS_no_control, 
  ordination = cap_ord, 
  color = "Treatment", 
  axes = c(1,2)
) + 
  aes(shape = Location_Year) + 
  geom_point(aes(colour = Treatment), alpha = 0.4, size = 4) + 
  geom_point(colour = "grey90", size = 1.5) + 
  scale_color_manual(values = c("#a65628", "red", "#ffae19", "#4daf4a", 
                                "#1919ff", "darkorchid3", "magenta","green", "orange","violet","blue")
  )


# Now add the environmental variables as arrows
arrowmat <- vegan::scores(cap_ord, display = "bp")

# Add labels, make a data.frame
arrowdf <- data.frame(labels = rownames(arrowmat), arrowmat)

# Define the arrow aesthetic mapping
arrow_map <- aes(xend = CAP1, 
                 yend = CAP2, 
                 x = 0, 
                 y = 0, 
                 shape = NULL, 
                 color = NULL, 
                 label = labels)

label_map <- aes(x = 1.3 * CAP1, 
                 y = 1.3 * CAP2, 
                 shape = NULL, 
                 color = NULL, 
                 label = labels)

arrowhead = arrow(length = unit(0.02, "npc"))

# Make a new graphic
cap_plot + 
  geom_segment(
    mapping = arrow_map, 
    size = .5, 
    data = arrowdf, 
    color = "gray", 
    arrow = arrowhead
  ) + 
  geom_text(
    mapping = label_map, 
    size = 4,  
    data = arrowdf, 
    show.legend = FALSE
  )
##
anova(cap_ord)

##########################
##Alpha Diversity
min_lib <- min(sample_sums(psY))
# Initialize matrices to store richness and evenness estimates
nsamp = nsamples(psY)
trials = 100

richness <- matrix(nrow = nsamp, ncol = trials)
row.names(richness) <- sample_names(psY)

evenness <- matrix(nrow = nsamp, ncol = trials)
row.names(evenness) <- sample_names(psY)

# It is always important to set a seed when you subsample so your result is replicable 
set.seed(3)

for (i in 1:100) {
  # Subsample
  r <- rarefy_even_depth(psY, sample.size = min_lib, verbose = FALSE, replace = TRUE)
  
  # Calculate richness
  rich <- as.numeric(as.matrix(estimate_richness(r, measures = "Observed")))
  richness[ ,i] <- rich
  
  # Calculate evenness
  even <- as.numeric(as.matrix(estimate_richness(r, measures = "InvSimpson")))
  evenness[ ,i] <- even
}
# Create a new dataframe to hold the means and standard deviations of richness estimates
SAMPLEID_S <- row.names(richness)
mean <- apply(richness, 1, mean)
sd <- apply(richness, 1, sd)
measure <- rep("Richness", nsamp)
rich_stats <- data.frame(SAMPLEID_S, mean, sd, measure)

# Create a new dataframe to hold the means and standard deviations of evenness estimates
SAMPLEID_S <- row.names(evenness)
mean <- apply(evenness, 1, mean)
sd <- apply(evenness, 1, sd)
measure <- rep("Inverse Simpson", nsamp)
even_stats <- data.frame(SAMPLEID_S, mean, sd, measure)
#Now we will combine our estimates for richness and evenness into one dataframe
alpha <- rbind(rich_stats, even_stats)
#add the sample metadata into this dataframe using the merge() command
s <- data.frame(sample_data(psY))
alphadiv <- merge(alpha, s, by = "SAMPLEID_S") 
#reorder some factors in this dataset before plotting them
alphadiv <- order_dates(alphadiv)
#plot the two alpha diversity measures in a timeseries using a facet
p<-ggplot(alphadiv, aes(x = Treatment, y = mean, color = Year_N, group = Year_N, shape = Location)) +
  geom_point(size = 2) + 
  geom_line(size = 0.8) +
  facet_wrap(measure ~ Location, ncol = 2, scales = "free") +
  scale_color_brewer(type="div", palette = 'Set1') +
  scale_x_discrete(
    breaks = c("T_1", "T_2", "T_3","T_4", "T_5", "T_6","T_7", "T_8"),
    labels = c("T_1", "T_2", "T_3","T_4", "T_5", "T_6","T_7", "T_8"), 
    drop = FALSE
  ) +
  theme(
    axis.title.x = element_blank(),
    axis.title.y = element_blank()
  )
p
#T3,4,5
pd <- position_dodge(0.1) # move them .05 to the left and right
p<-ggplot(alphadiv, aes(x = Year_N, y = mean, color = Treatment, group = Treatment)) +
  geom_point(size = 2, shape=21, fill="white") + 
  geom_line(size = 0.8) +
  geom_errorbar(aes(ymin = mean-sd, ymax = mean+sd),colour="black",position= pd) +
  facet_wrap(measure ~ Location, ncol = 2, scales = "free") +
  scale_color_brewer(type="div", palette = 'Set1') +
  ylab("Alpha Diversity") +
  scale_colour_hue(name="Treatments",    # Legend label, use darker colors
                   breaks=c("T_3", "T_4", "T_5"),
                   labels=c("Treatment 3", "Treatment 4","Treatment 5"),
                   l=40) +                    # Use darker colors, lightness=40
  ggtitle("The Richness and Inverse Shannon: Alpha Diversity of ITS") +
  expand_limits() +                        # Expand y range
  scale_y_continuous() +         # Set tick every 4
  theme_bw(
  )
p
###No Control Treatments 1-7 Year 1-3
#plot the two alpha diversity measures in a timeseries using a facet
p<-ggplot(alphadiv, aes(x = Year_N, y = mean, color = Treatment, group = Treatment)) +
  geom_point(size = 2, shape=21, fill="white") + 
  geom_line(size = 0.8) +
  geom_errorbar(aes(ymin = mean-sd, ymax = mean+sd),colour="black",position= pd) +
  facet_wrap(measure ~ Location, ncol = 2, scales = "free") +
  scale_color_brewer(type="div", palette = 'Set1') +
  ylab("Alpha Diversity") +
  scale_colour_hue(name="Treatments",    # Legend label, use darker colors
                   breaks=c("T_1", "T_2", "T_3","T_4", "T_5", "T_6","T_7"),
                   labels=c("Treatment 1","Treatment 2","Treatment 3", "Treatment 4","Treatment 5","Treatment 6","Treatment 7"),
                   l=40) +                    # Use darker colors, lightness=40
  ggtitle("The Richness and Inverse Shannon: Alpha Diversity of ITS") +
  expand_limits() +                        # Expand y range
  scale_y_continuous() +
  theme_bw(
  )
p
#####http://joey711.github.io/phyloseq-demo/Restroom-Biogeography.html
surface_group = get_variable(restroomR, "SURFACE")
surface_ano = anosim(distance(restroomR, "bray"), surface_group)
surface_ano$signif
surface_ano$statistic

restroomRg = subset_samples(restroomR, !GENDER == "None")
gender_group = get_variable(restroomRg, "GENDER")
gender_ano = anosim(distance(restroomRg, "bray"), gender_group)
gender_ano$signif

df = as(sample_data(restroomR), "data.frame")
d = distance(restroomR, "bray")
restroomadonis = adonis(d ~ SURFACE + GENDER + BUILDING + FLOOR, df)
restroomadonis

plot(restroomadonis$aov.tab)
