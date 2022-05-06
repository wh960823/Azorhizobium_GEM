#Root nodule VS CK
library(edgeR)
data <-read.table(file.choose(), header=TRUE, sep="\t")
d <- data[,2:7]
rownames(d) <- data[,1]
group <- c(rep("CK",3),rep("RN",3))
my_list <- DGEList(counts = d, group=group)
TMM <- calcNormFactors(my_list, method="TMM")
design <- model.matrix(~group)
dge <- estimateDisp(TMM, design, robust = TRUE)
fit <- glmFit(dge, design, robust = TRUE)
lrt <- topTags(glmLRT(fit), n = nrow(my_list$counts))
write.table(glmLRT(fit)$table, 'control_treat.glmLRT1_RN.txt', sep = '\t', col.names = NA, quote = FALSE)

# Stem nodule vs CK
data <-read.table(file.choose(), header=TRUE, sep="\t")
d <- data[,2:7]
rownames(d) <- data[,1]
group <- c(rep("CK",3),rep("SN",3))
my_list <- DGEList(counts = d, group=group)
TMM <- calcNormFactors(my_list, method="TMM")
design <- model.matrix(~group)
dge <- estimateDisp(TMM, design, robust = TRUE)
fit <- glmFit(dge, design, robust = TRUE)
lrt <- topTags(glmLRT(fit), n = nrow(my_list$counts))
write.table(lrt, 'control_treat.glmLRT_SN.txt', sep = '\t', col.names = NA, quote = FALSE)

# 8mM VS CK
data <-read.table(file.choose(), header=TRUE, sep="\t")
d <- data[,2:7]
rownames(d) <- data[,1]
group <- c(rep("CK",3),rep("SN",3))
my_list <- DGEList(counts = d, group=group)
TMM <- calcNormFactors(my_list, method="TMM")
design <- model.matrix(~group)
dge <- estimateDisp(TMM, design, robust = TRUE)
fit <- glmFit(dge, design, robust = TRUE)
lrt <- topTags(glmLRT(fit), n = nrow(my_list$counts))
write.table(lrt, 'control_treat.glmLRT_8mM.txt', sep = '\t', col.names = NA, quote = FALSE)

# 4mM VS CK
data <-read.table(file.choose(), header=TRUE, sep="\t")
d <- data[,2:7]
rownames(d) <- data[,1]
group <- c(rep("CK",3),rep("SN",3))
my_list <- DGEList(counts = d, group=group)
TMM <- calcNormFactors(my_list, method="TMM")
design <- model.matrix(~group)
dge <- estimateDisp(TMM, design, robust = TRUE)
fit <- glmFit(dge, design, robust = TRUE)
lrt <- topTags(glmLRT(fit), n = nrow(my_list$counts))
write.table(lrt, 'control_treat.glmLRT_4mM.txt', sep = '\t', col.names = NA, quote = FALSE)

# 0.01mM VS CK
data <-read.table(file.choose(), header=TRUE, sep="\t")
d <- data[,2:7]
rownames(d) <- data[,1]
group <- c(rep("CK",3),rep("CK1",3))
my_list <- DGEList(counts = d, group=group)
TMM <- calcNormFactors(my_list, method="TMM")
design <- model.matrix(~group)
dge <- estimateDisp(TMM, design, robust = TRUE)
fit <- glmFit(dge, design, robust = TRUE)
lrt <- topTags(glmLRT(fit), n = nrow(my_list$counts))
write.table(lrt, 'control_treat.glmLRT_001mM.txt', sep = '\t', col.names = NA, quote = FALSE)

# WGCNA analysis
library(WGCNA)
library(stringr)
trans_data <- read.table(file.choose(),head=T) 
Expr <-as.data.frame(t(trans_data[,1:12]))
dim(Expr)
MEs = net$MEs
MEs = orderMEs(MEs)
modNames = substring(names(MEs),3)
module = "red"
moduleGenes = moduleColors == module
column = match(module, modNames)
sizeGrWindow(7, 7);
par(mfrow = c(1,1));
verboseScatterplot(abs(geneModuleMembership[moduleGenes, column]),3
                   abs(geneTraitSignificance[moduleGenes, 1]),
                   xlab = paste("Module Membership in", module, "module"),
                   ylab = "Gene significance for body weight",
                   main = paste("Module membership vs. gene significance\n"),
                   cex.main = 1.2, cex.lab = 1.2, cex.axis = 1.2, col = module)
geneModuleMembership = as.data.frame(cor(Expr, MEs, use = "p"));
MMPvalue = as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership), nSamples));
names(geneModuleMembership) = paste("MM", modNames, sep="");
names(MMPvalue) = paste("p.MM", modNames, sep="");
geneTraitSignificance = as.data.frame(cor(Expr, weight, use = "p"));
GSPvalue = as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance), nSamples));
sizeGrWindow(7, 7);
par(mfrow = c(1,1));
verboseScatterplot(abs(geneModuleMembership[moduleGenes, column]),
                   3 weight_g length_cmab_fat other_fat total_fat X100xfat_weight Trigly Total_Chol HDL_CholUCFFA Glucose LDL_plus_VLDL MCP_1_phys Insulin_ug_l Glucose_Insulin Leptin_pg_ml Adiponectin Aortic.lesions Aneurysm Aortic_cal_M Aortic_cal_L CoronaryArtery_Cal Myocardial_cal BMD_all_limbs BMD_femurs_only
                   abs(geneTraitSignificance[moduleGenes, 1]),
                   xlab = paste("Module Membership in", module, "module"),
                   ylab = "Gene significance for body weight",
                   main = paste("Module membership vs. gene significance\n"),
                   cex.main = 1.2, cex.lab = 1.2, cex.axis = 1.2, col = module)

# PCA ploting
pca_data=t(as.matrix(rpkm_matrix))
pca_data1 = t(as.matrix(rpkm_matrix[,1:12]))
pca_data2 = t(as.matrix(rpkm_matrix[,13:18]))
pca_group=factor(c(rep('CK',3),rep('0.01mM',3), rep('4mM', 3), rep('8mM', 3),rep('RN', 3), rep('SN', 3)))
pca_group=factor(c(rep('CK',3),rep('0.01mM',3), rep('4mM', 3), rep('8mM', 3)))
pca_group=factor(c(rep('RN', 3), rep('SN', 3)))
pca1=prcomp(pca_data1,center=TRUE,retx=T)
p2 <-ggord(pca1, grp_in = pca_group, arrow=0,vec_ext = 0,txt=NULL,ellipse=T,ellipse_pro=0.95,hull=T,
           cols=c("#F8766D", "#B79F00", "#00BA38", "#00BFC4"), xlims=c(-100,80), ylims = c(-60,60))

#draw of fit line 
my_data <- read.table(file.choose(), header=TRUE, sep="\t")
colnames(my_data) <- c("locus_tag","RN/CK","SN/CK")
library(ggpubr)
ggscatter(my_data,x="RN/CK",y="SN/CK",add = "reg.line",conf.int = T,ggtheme=theme_classic(),ellipse.level=0.95,palette="npg",size=2,fullrange=T,font.label = c(15, "bold"))+
  stat_cor(label.x = -15, label.y = 15) + xlim(-15,15)+ylim(-15,15)

# plot the expression pattern of co-expression modules
my_data <- read.table(file.choose(), header = F)
df <- as.data.frame(my_data)
group_v = factor("CK", "0.01mM", "4mM", "8mM")
ggplot(data = df, mapping = aes(x = df$V2,y = df$V3,group=group_v))+geom_line(aes(group=df$V1), lwd =1.5)+geom_point(size=2)

# one way ANOVA
my_data <- read.table(file.choose(), header = F)
my_data <- as.data.frame(my_data)[1:4,]
group = as.factor(c("CK","0.01mM","4mM","8mM"))
qqPlot(lm(my_data$V3~as.factor(my_data$V2), data = my_data), simulate = TRUE, main = 'QQ Plot', labels = FALSE)

my_data <- read.table(file.choose(), header = T, sep = "\t")
expr <- my_data[,3:20]
rownames(expr) = my_data[,1]
pheatmap((expr),scale="row", cluster_cols = F,cluster_rows=T,color = colorRampPalette(c("MediumBlue","white","red"))(256),
         cellwidth = 15, cellheight =10)

# Boxplot
a1 <- read.table(file.choose(),sep="\t",header =T)
a1$Treatment <- factor(a1$Treatment, levels = c('CK','0.01mM','4mM','8mM','RN','SN'), ordered=TRUE)
p1 <- ggplot(data = a1,aes(x=a1$Treatment,y=a1$RPKM,colour=factor(a1$locus_tag)))+
  geom_boxplot(lwd=1,fill='white')+facet_wrap(~a1$gene.name,nrow = 10,scales='free')+guides(fill=F)
p2 <- p1+theme_bw()
p3 <- p2+ scale_color_aaas()

#
my_data <- read.table(file.choose(),sep="\t", header=F)
for_hclust <- my_data[,2:4]
rownames(for_hclust) <- my_data[,1]
my_dist <- dist(for_hclust, method = "euclidean", diag = FALSE, upper = FALSE)
a <- hclust(my_dist,method="average")
b <- cutree(a, k = 10)

#plot expression data in a pathway
library(reshape2)
library(ggplot2)
my_data <- read.table(file.choose(), header = T, sep="\t")
rownames(my_data) <- my_data[,1]
my_data1 <- melt(my_data[,1:5], id = "gene_ID", variable.name="treatment",value.name="logFC")
ggplot(data = my_data1, aes(x = treatment, y = logFC, group = gene_ID)) + 
  geom_point(size = 2, color = "blue") + 
  geom_line(size = 1.5, color = "blue", alpha = .6) + 
  labs(x = "Treatment", y = "logFC") + theme_bw()

#draw bar plot
my_data <- read.table(file.choose(),sep="\t",header=T)
ggplot(data = my_data, aes(x = Pathway, y =Count)) + geom_col(color = 'blue', fill='white', lwd=1.5)+coord_flip()+theme_bw()

#
library(NbClust)
library(ape)
my_data <- read.table(file.choose(), header=TRUE, sep ="\t")
for_clust <- as.data.frame(t(my_data))
for_clust_scaled <- scale(for_clust)
d<-dist(for_clust_scaled)
clust_fit <- hclust(d, method="average")
tree = as.phylo(clust_fit)
write.tree(tree, file="tree.nwk") 
ggtree(clust_fit,layout="rectangular", size=0.8, col="deepskyblue3",branch.length = 300)+geom_tiplab(size=3, color="purple4")+xlim(NA,150)
plot(clust_fit,hang=-1,cex=.8,main="average")
nc <- NbClust(for_clust_scaled, distance = "euclidean", method = "average", min.nc = 2,max.nc = 15)
for_clust1 <- for_clust[1:4,]
fviz_nbclust(for_clust1,kmeans,method = "wss") + geom_vline(xintercept = 4,linetype=2)
km_result <- kmeans(for_clust1,2,nstart = 24)
dd <- cbind(for_clust1,cluster=km_result$cluster)
fviz_cluster(km_result,data=for_clust1,palette=c("#2E9FDF","#00AFBB",
                                                "#E7B800","#FC4E07"),
             ellipse.type ="euclid",star.plot=TRUE,repel=TRUE,ggtheme = theme_minimal())

#draw bubble plot
library("ggsci")
library(ggplot2)
my_file <- read.table(file.choose(),header=T,sep="\t")
ggplot(data=my_file,aes(x=odd.ratio,y=Metabolite.Set,size=Hits,color=-log10(FDR))) +
  geom_point() + theme_bw() +
  scale_colour_gradient(low="lightblue",high="darkblue") + 
  labs(x="Odd Ratio",y="Metabolite Category",title="Enriched Pathways",
       colour=expression(-log[10]("FDR")),size="Metabolite number") +
  theme_bw()+scale_fill_npg() + theme(plot.title = element_text(hjust = 0.5))


#draw valcano
my_data <- read.table(file.choose(),sep="\t",header=T)
ggplot(data=my_data,aes(x=GC.content,y=Abundance,size=Genome.size,color=class)) +
  geom_point(alpha=0.5) + theme_bw()+scale_color_aaas()

# calculate correlation
library(ggcorrplot)
my_data <- read.table(file.choose(),sep="\t",header=T)
cor_my_data <- round(cor(my_data), 3)
p_my_data <- cor_pmat(my_data)
ggcorrplot(cor_my_data,hc.order = T,  
           ggtheme = ggplot2::theme_void(base_size = 15), 
           colors = c("CornflowerBlue","white","Salmon"), 
           lab = T,lab_size = 5,    
           tl.cex = 15,             
           p.mat = p_my_data,         
           sig.level = 0.01,        
           pch = 3,                 
           pch.cex = 10,
           show.diag =FALSE,
           type="upper")            


# draw valcano plot
my_data <- read.table(file.choose(),sep="\t",header=T)
themel <-theme_bw()+theme(#legend.position='top',#图例位置
  legend.title = element_blank(),#图例标题去掉
  #text= element_text(family=''),#字体
  #panel.border = ekement_blank(),#图形边界
  #strip.background = element_rect(fill='white'),#分面图小标题背景
  strip.text = element_text(size =12,face = 'bold',colour ='black'),#分面图小标题字体
  plot.title = element_text(hjust=0.5,size=16,vjust=0.5),#标题位置
  panel.grid.major=element_blank(),#网格线
  panel.grid.minor=element_blank(),#次级网格线
  #legned.title=element_text(size=10,colour='black',vjust=-0.5),#图例标题
  legend.text=element_text(size=10,colour='black',face='bold'),#图例文字
  #legend.background= element_rect(),#图例背景
  axis.text=element_text(size=10,colour='black',face='bold'),#坐标轴文字
  #axis.text.x=element_text(angle=45,hjust=0.9),
  #axis.line=element_line(size=0.8,colour='black'),#坐标轴颜色
  axis.title=element_text(size=12,colour='black',face='bold')#坐标轴标签
  #panel.spacing=unit(10,'mm')
)#画布大小
ggplot(data=my_data,aes(x=my_data$VIP.0.01mM,y=-log(my_data$FDR.0.01mM)))+geom_point(color = "orange",size=exp(my_data$FC.0.01mM))+themel + geom_vline(xintercept = 1,color="blue",lwd=1.2)


#draw heatmap
library(ComplexHeatmap)
my_data <- read.table(file.choose(),sep="\t",header=T)
Heatmap(my_data[,3:6],col =colorRampPalette(c("navy", "white", "firebrick3"))(50),
        )
