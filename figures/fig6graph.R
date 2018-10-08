
#figure 6 tables of chemistry components: 

library(xlsx)
library(ggplot2)
library(reshape2)
library(plyr)
library(dplyr)
library(readr)

library(cowplot)

data <- read.csv2("final_chemistry_for_AN2_pub.csv", header = T, sep = "\t")  

data
sds <- data %>% select(contains("sd"))

sds <- cbind(data$Species, sds)
sds <- melt(sds)

chemicals <- data %>% select(-contains("sd"))

chemicals <- melt(chemicals)

data_pl <- cbind(chemicals, sds[,3])

head(data_pl)
names(data_pl) <- c("Species", "variable", "value", "sd")

data_pl <- dplyr::filter(data_pl, !grepl("Peonidin",variable))
data_pl <- dplyr::filter(data_pl, !grepl("Cyanidin",variable))

data_pl$group <- c(rep("anthocyanins", 12), rep("flavonols", 12))

data_pl$variable <- factor(data_pl$variable, levels=c("Delphinidin", "Petunidin", "Malvidin", "Kaempferol", "Quercetin", "Myricetin"))


windowsFonts()

graphs = list()


#subset according to species: 
species = c("P.axillaris","P.inflata","P.secreta")

e = 0

for (i in species){
sset <-  data_pl %>% subset(Species == i)
print(sset)
e <- e + 1

graph <- ggplot(sset, aes(x = group, y = value, fill = variable)) +
  geom_histogram(colour = "black" , stat = "identity", position = "dodge") +
  
  geom_errorbar(aes(ymin = value - sd, ymax = value + sd), width = 0.2, position = position_dodge(0.9))+
  #facet_grid(. ~ Species, scales = "free") + 
  scale_fill_manual(values=c("darkmagenta", "darkorchid", "darkorchid4", "gray80", "gray50", "grey30")) +
  #fix scale 
  scale_y_continuous(breaks = seq(0, 2, by = 0.25)) + 
  ylim(0, 1.2)+
  theme(text=element_text(family="Arial", size = 20),
                          axis.text=element_text(size=20), 
          legend.position="none")+
  xlab("") + ylab("µg flavonoids per mg FW")
  

graphs[[i]] = graph

print(graph)
  
}


#get legend: 
sset <-  data_pl %>% subset(Species == "P.axillaris")
t <- ggplot(sset, aes(x = group, y = value, fill = variable)) +
  geom_histogram(colour = "black" , stat = "identity", position = "dodge") +
  
  geom_errorbar(aes(ymin = value - sd, ymax = value + sd), width = 0.2, position = position_dodge(0.9))+
  #facet_grid(. ~ Species, scales = "free") + 
  scale_fill_manual(values=c("darkmagenta", "darkorchid", "darkorchid4", "gray80", "gray50", "grey30")) +
  #fix scale 
  scale_y_continuous(breaks = seq(0, 2, by = 0.25)) + 
  ylim(0, 1.2)+
  theme(text=element_text(family="Arial", size = 20),
        axis.text=element_text(size=20))+ 
  xlab("") + ylab("µg flavonoids per mg FW")


library(gridExtra)
library(grid)
# extract Legend 
g_legend <- function(a.gplot){ 
  tmp <- ggplot_gtable(ggplot_build(a.gplot)) 
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box") 
  legend <- tmp$grobs[[leg]] 
  return(legend)} 


t <- g_legend(t) 

grid.newpage()
grid.draw(t)
ggsave("legend_fig6.tiff", t, device = "tiff", 
       dpi = 320)



graphs[["P.axillaris"]]

graphs[[1]]
ggsave("inflata.FLAVO20.tiff", graphs[[2]], device = "tiff", 
       dpi = 320)

ggsave("axillaris.FLAVO20.tiff", graphs[[1]], device = "tiff", 
       dpi = 320)

ggsave("secreta.FLAVO20.tiff", graphs[[3]], device = "tiff", 
       dpi = 320)




#complete graph 
ggplot(data_pl, aes(x = Species, y = value, fill = variable)) +
  geom_histogram(stat = "identity", position = "dodge") +
  
  geom_errorbar(aes(ymin = value - sd, ymax = value + sd), width = 0.2, position = position_dodge(0.9))+
  facet_grid(. ~ Species, scales = "free") + 
  scale_fill_manual(values=c("darkmagenta", "darkorchid", "darkorchid4", "gray80", "gray50", "grey30"))


ggplot(data_pl, aes(x = Species, y = value, fill = variable)) +
  geom_histogram(stat = "identity", position = "dodge") +
  
  geom_errorbar(aes(ymin = value - sd, ymax = value + sd), width = 0.2, position = position_dodge(0.9))+
    facet_grid(. ~ Species, scales = "free") + 
    scale_fill_manual(values=c("darkmagenta", "darkorchid", "darkorchid4", "gray80", "gray50", "grey30"))

