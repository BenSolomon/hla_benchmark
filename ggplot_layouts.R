library(tidyverse)
library(ggh4x)


gg_accuracy <- function(df){
  df %>% 
    ggplot(aes(x = field, y = accuracy))+
    stat_summary(fun = mean, geom = "point", position = position_dodge(width=0.9)) +
    stat_summary(fun.data = mean_se, geom = "errorbar", position=position_dodge(width=0.9)) +
    theme_bw() +
    facet_nested(.~locus+genotyper, scales = "free_y") +
    coord_cartesian(ylim = c(0,1))+
    scale_y_continuous(n.breaks = 6)+
    labs(y = "Accuracy", x = "")+
    theme(axis.text.x = element_text(angle = 90))+
    theme(panel.spacing = unit(0.1,"line"))
}


gg_frequency <- function(df){
  df %>% 
    ggplot(aes(x = field, y = frequency))+
    stat_summary(fun = mean, geom = "point", position = position_dodge(width=0.9)) +
    stat_summary(fun.data = mean_se, geom = "errorbar", position=position_dodge(width=0.9)) +
    theme_bw() +
    facet_nested(.~locus+genotyper, scales = "free_y") +
    coord_cartesian(ylim = c(0,1))+
    scale_y_continuous(n.breaks = 6)+
    labs(y = "Rate of allele prediction", x = "")+
    theme(axis.text.x = element_text(angle = 90))+
    theme(panel.spacing = unit(0.2,"line"))
}


strip_y <- function(df){
  df + theme(axis.title.y=element_blank(),
             axis.text.y=element_blank(),
             axis.ticks.y=element_blank())
}
