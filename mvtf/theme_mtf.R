library(tidyverse)

theme_mtf = theme_gray() + 
  theme(plot.title=element_text(family='sans',face='bold',size=10,hjust=0.5),
        axis.title.y=element_text(family='sans',face='bold',size=8),
        axis.text.y=element_text(family='sans',size=7), 
        axis.ticks.x=element_blank(),
        legend.text=element_text(family='sans',size=8),
        legend.title=element_text(family='sans',face='bold',size=8),
        legend.key.size=unit(0.4,'cm'),
        axis.text.x=element_text(family='sans',face='italic',size=7),
        axis.title.x=element_text(family='sans',face='bold',size=8),
        strip.text=element_text(family='sans',face='bold',size=8))
