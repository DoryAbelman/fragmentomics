# file: git_06-plotting.R
# author: Derek Wong, Ph.D
# date: July 22nd, 2021
# edited by Dory Abelman from Feb 4-11, 2022 


# manipulated id to match M4 study naming scheme

id_edited <- id 

# if the id_edited contains an M4 tag, process it like an M4 sample, else don't change:
if (grepl("M4-", id_edited, fixed = TRUE)){
  id_edited <- sub("M4-", "", id) # remove start (M4-)
  id_edited <- sub("\\DNA[^DNA]*$", "", id_edited) # remove everything after 'DNA'
  # Attach sample type
  id_edited <- sub("-O-", "_Bone_marrow_cells", id_edited)
  id_edited <- sub("-P-", "_Blood_cell_free_DNA", id_edited)
  id_edited <- sub("-N-", "_Bone_marrow_cfDNA", id_edited)
  id_edited <- sub("-B-", "_Blood_buffy_coat", id_edited)
}

## Load in healthy controls
load(healthy)

## Set themes and plot layouts - slightly modified theme parameters to match desired style
mytheme <- theme_classic(base_size=12) + theme(
  axis.text.x = element_blank(),
  axis.ticks.x=element_blank(),
  strip.text.x = element_text(size=11),
  strip.text.y = element_text(size=12),
  axis.title.x = element_text(face="bold", size=17),
  axis.title.y = element_text(size=15),
  axis.text.y = element_text(size=15),
  plot.title = element_text(size=15),
  legend.position = "right",
  legend.title = element_text(size=10),
  legend.text = element_text(size=10),
  panel.grid.major = element_blank(),
  panel.grid.minor = element_blank(),
  strip.background=element_rect(fill="white", color="white"),
  panel.spacing.x=unit(0.1, "lines"))

armlevels <- c("1p","1q","2p","2q","3p","3q","4p","4q","5p","5q","6p","6q",
               "7p","7q","8p","8q", "9p", "9q","10p","10q","11p","11q","12p",
               "12q","13q","14q","15q","16p","16q","17p","17q","18p","18q",
               "19p", "19q","20p","20q","21q","22q")
df.fr3$arm <- factor(df.fr3$arm, levels=armlevels)
healthy_median$arm <- factor(healthy_median$arm, levels=armlevels)

arm <- df.fr3 %>% group_by(arm) %>%
  dplyr::summarize(n=n()) %>%
  mutate(arm = as.character(arm))
small.arms <- setNames(c("", "10q", "", "12q", "", "16",
                         "", "17q", "", "18q",
                         "", "", "", "",
                         "", ""),
                       c("10p", "10q", "12p", "12q", "16p", "16q",
                         "17p", "17q", "18p", "18q",
                         "19p", "19q", "20p", "20q",
                         "21q", "22q"))
arm.labels <- setNames(arm$arm, arm$arm)
arm.labels[names(small.arms)] <- small.arms

# Setting up plots - added colour parameter for legend

colours <- c("Normal" = "black", "Sample" = "red")

# Modified plots to include legend, less white space, and titles + subtitles

## Plot Fragmentation profile
f1 <- ggplot(df.fr3, aes(x=bin, y=ratio_centered, group=sample_id, color="Sample")) + 
  geom_line(size=0.75, alpha=0.9)
f1 <- f1 + geom_line(data=healthy_median, aes(x=bin, y=median_ratio_centered, color = "Normal"), size=0.75, alpha=0.6)
f1 <- f1 + labs(title=id_edited, subtitle="Fraction of small DNA fragments (100–150 bp) to larger DNA fragments (151–220 bp)", x="Chromosome position", y="Fragmentation profile\n", color="Legend")
f1 <- f1 + facet_grid(~arm, switch="x",space="free_x", scales="free_x", labeller=labeller(arm=arm.labels))
f1 <- f1 + coord_cartesian(xlim = NULL, ylim=c(-.7,.9), expand = TRUE)
f1 <- f1 + scale_color_manual(values=colours) + mytheme
ggsave(file.path(outdir, paste0(id, "_fragment.pdf")), f1, width=17, height=5, units="in")

## Plot short fragment coverage profile
c1 <- ggplot(df.fr3, aes(x=bin, y=coverage_centered, group=sample_id, color="Sample")) + 
  geom_line(size=0.75, alpha=0.9)
c1 <- c1 + geom_line(data=healthy_median, aes(x=bin, y=median_coverage_centered, color = "Normal"), size=0.75, alpha=0.6)
c1 <- c1 + labs(title=id_edited, subtitle="Coverage profile of short (≤150 bp) fragments", x="Chromosome position", y="Coverage profile\n", color="Legend")
c1 <- c1 + facet_grid(~arm, switch="x",space="free_x", scales="free_x", labeller=labeller(arm=arm.labels))
c1 <- c1 + coord_cartesian(xlim = NULL, ylim=c(-.7,.9), expand = TRUE)
c1 <- c1 + scale_color_manual(values=colours) + mytheme
ggsave(file.path(outdir, paste0(id, "_coverage.pdf")), c1, width=17, height=5, units="in")

## Plot combined profile
b1 <- ggplot(df.fr3, aes(x=bin, y=combined_centered, group=sample_id, color="Sample")) + 
  geom_line(size=0.75, alpha=0.9)
b1 <- b1 + geom_line(data=healthy_median, aes(x=bin, y=median_combined_centered, color = "Normal"), size=0.75, alpha=0.6)
b1 <- b1 + labs(title=id_edited, subtitle="Adjusted ratio of short to long fragments based on GC correction", x="Chromosome position", y="Combined profile\n", color="Legend")
b1 <- b1 + facet_grid(~arm, switch="x",space="free_x", scales="free_x", labeller=labeller(arm=arm.labels))
b1 <- b1 + coord_cartesian(xlim = NULL, ylim=c(-.7,.9), expand = TRUE)
b1 <- b1 + scale_color_manual(values=colours) + mytheme
ggsave(file.path(outdir, paste0(id, "_combined.pdf")), b1, width=17, height=5, units="in")

# Print exit command
print("Plotting completed") 
