if(!requireNamespace('fpc',quietly = TRUE))
  install.packages('fpc')

library(fpc) # 加载fpc包，用于统计不同p值的轮廓系数

rm(list=ls())
setwd('D:/cluster')

#读入数据
expr=read.delim('新建文本文档.txt',header=T,row.names=1,sep="\t")
expr_s=scale(log2(expr+1))

#=======================选K值：1.使用轮廓系数=======================
K <- 2:10  #设置遍历的K值范围，从2到10
round <- 2  # 实际分析设置为20以上，迭代计算取均值，避免局部最优
#计算每个K分类的轮廓系数
bestK <- sapply(K,function(i){
  print(paste("K=",i))
  mean(sapply(1:round,function(r){
    print(paste("Round",r))
    result <- kmeans(expr, i)
    stats <- cluster.stats(dist(expr), result$cluster)
    stats$avg.silwidth
  }))
})

#直观展示不同K值对应轮廓系数
plot(K,bestK,type='l',main='轮廓系数与K的关系', ylab='轮廓系数')   

# 自动挑选轮廓系数最大值对应的K
Selected_K=K[match(max(bestK),bestK)]

#==================选K值：2.手动赋值===================================
#也可以手动设置K值，直接赋值给Selected_K变量：
Selected_K=6


#=================Kmeans聚类=================================
set.seed(107)     #设置随机种子，固定随机数，方便结果重现，此处数字可以任意设置
result=kmeans(expr_s,Selected_K)  #使用选择的K值进行kmeans聚类
table(result$cluster)


#================设置用于展示的表达矩阵=====================
for_plot=log2(expr+1)  #使用log2转换缩小显示范围
for_plot$Group=as.factor(result$cluster)    #添加K均值聚类的分类结果
for_plot=for_plot[order(for_plot$Group),]   #按照分类对表达矩阵排序



#================将K均值分类结果输出到文件=================
write.table(for_plot,file='kmeans-Result.txt',sep="\t",quote=F)


#================可视化1.热图绘制======================
#加载画热图函数包，绘制热图与分类结果图
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager", update=F)
if (!requireNamespace("ComplexHeatmap", quietly = TRUE))
  BiocManager::install("ComplexHeatmap", update=F)
library(ComplexHeatmap)


#设置指示条
ha = rowAnnotation(df = for_plot$Group,width = unit(5, "mm"),col=list(Group=setNames(rainbow(Selected_K),levels(for_plot$Group))))

DrawExp=for_plot
DrawExp$Group=NULL

#设置热图
hexp=Heatmap(as.matrix(DrawExp),na_col="grey",cluster_rows=F,cluster_columns = F,color=cc(100)  #此处的cluster为层次聚类，因此不用
             #, width = unit(10, "mm") # 每个cell的大小，如果不用，便针对输出大小进行自适应设置
              ,column_title ="Gene\nExpression"  #图标题
              ,column_title_gp = gpar(fontsize = 10)  
              ,heatmap_legend_param = list(title = "log2(FPKM+1)") # 图注
              ,show_column_names=T, # 是否显示列名
              show_row_names=FALSE)#是否显示行名
)

# 输出png文件，可以适当调节高和宽
png('kmeans.heat.png',h=3000,w=2000,res=300,units = 'px')
draw(ha+hexp, gap = unit(c(10, 15), "mm"))
dev.off()

#========================可视化2.折线图绘制=================
if (!requireNamespace("dplyr", quietly = TRUE))
  install.packages("dplyr", update=F)
if (!requireNamespace("reshape2", quietly = TRUE))
  BiocManager::install("reshape2", update=F)
if (!requireNamespace("RColorBrewer", quietly = TRUE))
  BiocManager::install("RColorBrewer", update=F)
if (!requireNamespace("ggplot2", quietly = TRUE))
  BiocManager::install("ggplot2", update=F)

library(dplyr)
library(reshape2)
library(RColorBrewer)
library(ggplot2)

#设置绘图用的颜色方案
usedColors=colorRampPalette(brewer.pal(7,'Set1'))(Selected_K)
#也可以使用下面方法自定义
#usedColors=c('red','blue','yellow','pink','black','purple')


#构建绘图的数据结构，包含表达量、基因ID以及每个基因的分类
result = for_plot %>% 
  select(-Group)%>%
  mutate(ID=rownames(for_plot),clust = paste0("cluster_", as.character(for_plot$Group))) %>%
  melt(id.var=c('ID','clust'))


#给每条线赋予不同分组的颜色
lineColors=usedColors[as.factor(result$clust)]

#可指定样品显示顺序，可以手动填写每个样品名，如下所示：
#order=c("S1","S2", "S3", "S4", "S5", "S6")
#for_plot$variable<- factor(for_plot$variable, levels = order)


#构建ggplot对象，并设置x轴为样品列，y轴为表达量列
#添加geom_line图层，设置每条线的透明度为0.1
g<-ggplot(result,aes(x =  variable , y = value , group = ID)) + geom_line(alpha = 0.1 , aes(col =lineColors))


#对每个cluster添加一个中位数线
g<-g + stat_summary(fun.y=median, colour="black", geom="line",group=1)


#设置外观
g<-g + theme_bw()


#按照cluster进行分面显示
g<-g + facet_wrap(~clust)


#添加坐标线以及主标题
g<-g + labs(x = "Samples",y = "log2(FPKM+1)",title = "K-means Clustering")


#微调字体等信息
g<-g+theme(
  legend.position = "none" ,
  plot.title=element_text(size=15,hjust=.5),
  axis.text.x = element_text(color="black", size=10, angle=90),
)

#输出到图片
png('kmeans.linePlot.png',w=3000,h=2000,res=300,units="px")
g
dev.off()
