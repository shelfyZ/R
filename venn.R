#-------------1.加载包------------------
pacman::p_load(tidyverse,ggplot2,RColorBrewer,ggvenn,ggsci,cowplot,clusterProfiler,
               org.At.tair.db,patchwork,ggpubr)

#-------------2.整理数据------------------
mvslup<-  read.table('mvslup.txt', header = T) %>% pull('Gene')
fvslup <- read.table('fvslup.txt', header = T) %>% pull('Gene')
mvsfup <- read.table('mvsfup.txt', header = T) %>% pull('Gene')
fvsmup <- read.table('fvsmup.txt', header = T) %>% pull('Gene')

common <- list('M vs L' = mvslup , 'F vs L' = fvslup,'M vs F' = mvsfup,'F vs M' = fvsmup)


p1 <- ggvenn(common,c('M vs L','F vs L'),fill_alpha = 0.7,set_name_size = 5) +
  scale_fill_npg()

p2 <- ggvenn(common,c('M vs L','M vs F'),fill_alpha = 0.7,set_name_size = 5) +
  scale_fill_npg()

p3 <- ggvenn(common,c('F vs L','F vs M'),fill_alpha = 0.7,set_name_size = 5)+
      scale_fill_npg()



plot_grid(p1,p2,p3,labels = LETTERS)


(p1|p2|p3)+
  plot_layout(heights = c(1,1.2))+
  plot_annotation(tag_levels = 'A')


##########################################################
#
#########################################################
p <-venn.diagram(x=list(ve = vege_expr2, re = repro_expr1),
                 filename = '2sets.tif',
                 col = 'black',
                 #fill = brewer_pal(palette = 'set2')(2))
                 fill = pal_npg(palette = c('nrc'),alpha=1)(2))