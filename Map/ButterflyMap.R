## Author: Hannah E. Aichelman
## This script creates a map of the western United States where Limenitis butterflies of interest were collected for this study.

#### Set-Up ####
#install.packages('usmap')

library(usmap)
library(ggplot2)

#### Create Map ####
# example of creating a map of the entire US:
plot_usmap(regions = "counties") + 
  labs(title = "US Counties",
       subtitle = "This is a blank map of the counties of the United States.") + 
  theme(panel.background = element_rect(color = "black", fill = "lightblue"))


# plot only western states of interest for butterfly ranges:
map = plot_usmap(include = c("CA", "ID", "NV", "OR", "WA"), fill = c("#74c476", "#fed976", "#fb6a4a", "#fed976", "#fed976"), alpha = 0.7) +
  #labs(title = "Western US States") +
  theme_bw()
map
ggsave(map, file = "/Users/hannahaichelman/Dropbox/BU_Postdoc/ButterflyGenomes/Github/Butterfly_Mimicry/Map/RangeMap.pdf", width=4, height=5, units=c("in"), useDingbats=FALSE)

