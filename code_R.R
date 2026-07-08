
###################################### codes for figures ################################################

###### load packages ######
library(readxl)
library(ggplot2)
library(ggrepel)
library(VennDiagram)
library(tidyverse)
library(pheatmap)
library(ComplexHeatmap)
library(circlize)
library(reshape2)
library(mgcv)
library(forestplot)
library(igraph)
library(ggraph)

###### Fig.2 ######
###fig.2a
data<-read_excel("Source Data Fig.2.xlsx", sheet = "fig2a-d")
data$lnhr<-log(data$hr)
data$lnll<-log(data$ll)
data$lnul<-log(data$ul)
data$q<-p.adjust(data$p,method="fdr")
data$logq<--log10(data$q)
data$lnhrb<-log(data$hrb)
data$lnllb<-log(data$llb)
data$lnulb<-log(data$ulb)
data$qb<-p.adjust(data$pb,method="fdr")
data$logqb<--log10(data$qb)

ggplot(data,aes(x=lnhrb,y=logqb))+
  geom_point(data=subset(data,qb > 0.05),aes(size=abs(lnhrb)),pch=21,color="black",fill="grey80",show.legend = FALSE)+
  geom_point(data=subset(data,lnhrb > 0 & qb < 0.05),aes(size=abs(lnhrb)),pch=21,color="black",fill="#c25759",show.legend = FALSE)+
  geom_point(data=subset(data,lnhrb < 0 & qb < 0.05),aes(size=abs(lnhrb)),pch=21,color="black",fill="#599cb4",show.legend = FALSE)+
  scale_size_continuous(range = c(1, 6)) +  
  geom_text_repel(data = subset(data, qb < 0.05),aes(label=EntrezGeneSymbol),size=3,direction="both",force=1)+
  scale_x_continuous(limits=c(-0.5,0.5),breaks=c(-0.5,-0.4,-0.3,-0.2,-0.1,0,0.1,0.2,0.3,0.4,0.5))+
  scale_y_continuous(limits=c(0,5),breaks=c(0,1,2,3,4,5))+
  theme_bw() +
  theme(panel.grid=element_blank(),legend.box.background = element_rect(color="black"),
        panel.background = element_blank())+ 
  theme(axis.text=element_text(size=8))+
  geom_hline(aes(yintercept=1.30103),color="darkgrey",linetype="dashed",linewidth=0.5)+
  xlab("Prospective effect")+ylab("-log(q)")+theme(axis.title=element_text(size=8))

###fig.2b
ggplot(data,aes(x=lnhr,y=logq))+
  geom_point(data=subset(data,q > 0.05),aes(size=abs(lnhr)),pch=21,color="black",fill="grey80",show.legend = FALSE)+
  geom_point(data=subset(data,lnhr > 0 & q < 0.05),aes(size=abs(lnhr)),pch=21,color="black",fill="#c25759",show.legend = FALSE)+
  geom_point(data=subset(data,lnhr < 0 & q < 0.05),aes(size=abs(lnhr)),pch=21,color="black",fill="#599cb4",show.legend = FALSE)+
  scale_size_continuous(range = c(1, 6)) +  
  geom_text_repel(data = subset(data, q < 0.05),aes(label=EntrezGeneSymbol),size=3,direction="both",force=1)+
  scale_x_continuous(limits=c(-0.5,0.5),breaks=c(-0.5,-0.4,-0.3,-0.2,-0.1,0,0.1,0.2,0.3,0.4,0.5))+
  scale_y_continuous(limits=c(0,5),breaks=c(0,1,2,3,4,5))+
  theme_bw() +
  theme(panel.grid=element_blank(),legend.box.background = element_rect(color="black"),
        panel.background = element_blank())+ 
  theme(axis.text=element_text(size=8))+
  geom_hline(aes(yintercept=1.30103),color="darkgrey",linetype="dashed",linewidth=0.5)+
  xlab("Time-varying effect")+ylab("-log(q)")+theme(axis.title=element_text(size=8))

###fig.2c
ggplot(data,aes(x=lnhr,y=lnhrb))+
  geom_point(data=subset(data,qb > 0.05&q>0.05),size=4,pch=21,color="black",fill="grey80")+
  geom_point(data=subset(data,qb > 0.05&q<0.05&lnhr>0),size=4,pch=21,color="black",fill="#e69191")+
  geom_point(data=subset(data,qb > 0.05&q<0.05&lnhr<0),size=4,pch=21,color="black",fill="#92b5ca")+
  geom_point(data=subset(data,qb < 0.05&q>0.05),size=4,pch=21,color="black",fill="#7c44af")+
  geom_point(data=subset(data,lnhrb > 0 & qb < 0.05&q<0.05),size=4,pch=21,color="black",fill="#c25759")+
  geom_point(data=subset(data,lnhrb < 0 & qb < 0.05&q<0.05),size=4,pch=21,color="black",fill="#599cb4")+
  geom_text_repel(data = subset(data, qb < 0.05|q<0.05),aes(label=EntrezGeneSymbol),size=4,direction="both",force=1)+
  scale_x_continuous(limits=c(-0.5,0.5),breaks=seq(-0.5,0.5,0.1))+
  scale_y_continuous(limits=c(-0.5,0.5),breaks=seq(-0.5,0.5,0.1))+
  theme_bw() +
  theme(panel.grid=element_blank(),legend.box.background = element_rect(color="black"),panel.background = element_blank())+
  theme(axis.text=element_text(size=8))+
  geom_hline(aes(yintercept=0),color="darkgrey",linetype="dashed",size=0.5)+geom_vline(aes(xintercept=0),color="darkgrey",linetype="dashed",size=0.5)+ ##设置参考线
  xlab("Time-varying effect")+ylab("Prospective effect")+theme(axis.title=element_text(size=8)) 


###fig.2d
set1 <- data$EntrezGeneSymbol[data$q < 0.05]
set2 <- data$EntrezGeneSymbol[data$qb < 0.05]

venn.plot <- draw.pairwise.venn(
  area1 = length(set1),
  area2 = length(set2),
  cross.area = length(intersect(set1, set2)),
  category = c("Time-varying", "Prospective"),
  fill = c("#c25759", "#599cb4"),
  alpha = 0.5,
  cex = 1.5,
  cat.cex = 1.5,
  lwd = 2)

grid.draw(venn.plot)
dev.off()


###fig.2e
df_sensi<-read_excel("Source Data Fig.2.xlsx", sheet = "fig2e")
df_sensi<-as.data.frame(df_sensi)

df_sensi <- df_sensi %>%  mutate(
  qa = p.adjust(pa, method = "fdr"),
  qb = p.adjust(pb, method = "fdr"),
  qc = p.adjust(pc, method = "fdr"),
  qd = p.adjust(pd, method = "fdr"),
  qe = p.adjust(pe, method = "fdr"))

df_sensi <- df_sensi %>%  mutate(
  hra = ifelse(qa > 0.05, NA, hra),
  hrb = ifelse(qb > 0.05, NA, hrb),
  hrc = ifelse(qc > 0.05, NA, hrc),
  hrd = ifelse(qd > 0.05, NA, hrd),
  hre = ifelse(qe > 0.05, NA, hre))

df_sensi <- df_sensi %>% mutate(
  lnhra = log(hra),
  lnhrb = log(hrb),
  lnhrc = log(hrc),
  lnhrd = log(hrd),
  lnhre = log(hre))

row.names(df_sensi)<-df_sensi$protein
mat <- df_sensi %>% select(lnhra, lnhrb, lnhrc, lnhrd, lnhre) %>% as.matrix()

row_anno <- data.frame(sig = df_sensi$sig)
row.names(row_anno) <- df_sensi$protein
ann_colors <- list(sig = c("pos" = "#c25759","neg" = "#599cb4"))

pheatmap(
  mat,
  color = colorRampPalette(c("#599cb4", "white", "#c25759"))(100),
  breaks = seq(-0.6, 0.6, length.out = 101),
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  annotation_row = row_anno,
  annotation_colors = ann_colors,
  fontsize_row = 10,
  fontsize_col = 10,
  border_color = "grey80",
  na_col = "white")


###fig.2f
df <- read_excel("Source Data Fig.2.xlsx", sheet = "fig2f")
df <- as.data.frame(df)
protein_names <- df$protein

immune_mat <- df[, c("innate", "adaptive", "cytokine")]
immune_mat <- as.matrix(immune_mat)
immune_label_mat <- matrix( "0",nrow = nrow(immune_mat),ncol = ncol(immune_mat))
rownames(immune_label_mat) <- protein_names
colnames(immune_label_mat) <- c("Innate", "Adaptive", "Signaling")
immune_label_mat[immune_mat[, "innate"] == 1, "Innate"] <- "Innate"
immune_label_mat[immune_mat[, "adaptive"] == 1, "Adaptive"] <- "Adaptive"
immune_label_mat[immune_mat[, "cytokine"] == 1, "Signaling"] <- "Signaling"

ptype_mat <- df[, c("ig", "c", "o")]
ptype_mat <- as.matrix(ptype_mat)
ptype_label_mat <- matrix("0",nrow = nrow(ptype_mat), ncol = ncol(ptype_mat))
rownames(ptype_label_mat) <- protein_names
colnames(ptype_label_mat) <- c("Immunoglobulins", "Complements", "Others")
ptype_label_mat[ptype_mat[, "ig"] == 1, "Immunoglobulins"] <- "Immunoglobulins"
ptype_label_mat[ptype_mat[, "c"] == 1, "Complements"] <- "Complements"
ptype_label_mat[ptype_mat[, "o"] == 1, "Others"] <- "Others"
col_immune <- c("0" = "white","Innate" = "#b3c6bf", "Adaptive" = "#efcdae", "Signaling" = "#a8b9cc")
col_ptype <- c("0" = "white", "Immunoglobulins" = "#dda29f","Complements" = "#78b9bb","Others" = "#a98dac")

ht1 <- Heatmap(
  immune_label_mat,
  name = "Immune system",
  col = col_immune,
  cluster_rows = FALSE,
  cluster_columns = FALSE,
  show_column_names = FALSE,
  show_row_names = TRUE,
  row_names_side = "left",
  row_names_gp = gpar(fontsize = 12),
  column_title = "Immune system",
  column_title_gp = gpar(fontsize = 13, fontface = "plain"),
  rect_gp = gpar(col = "black", lwd = 0.5),
  heatmap_legend_param = list(
    title = "Immune system",
    at = c("Innate", "Adaptive", "Signaling"),
    labels = c("Innate", "Adaptive", "Signaling")))

ht2 <- Heatmap(
  ptype_label_mat,
  name = "Protein type",
  col = col_ptype,
  cluster_rows = FALSE,
  cluster_columns = FALSE,
  show_column_names = FALSE,
  show_row_names = FALSE,
  column_title = "Protein type",
  column_title_gp = gpar(fontsize = 13, fontface = "plain"),
  rect_gp = gpar(col = "black", lwd = 0.5),
  heatmap_legend_param = list(
    title = "Protein type",
    at = c("Immunoglobulins", "Complements", "Others"),
    labels = c("Immunoglobulins", "Complements", "Others")))

draw( ht1 + ht2,  gap = unit(4, "mm"),  heatmap_legend_side = "right",  merge_legend = FALSE)


###fig2g
df <- read_excel("Source Data Fig.2.xlsx", sheet = "fig2g")
df <- as.data.frame(df)
df$lnhrcancer <- log(df$hrcancer)
df$lnllcancer <- log(df$llcancer)
df$lnulcancer <- log(df$ulcancer)
df$lnhrcvd <- log(df$hrcvd)
df$lnllcvd <- log(df$llcvd)
df$lnulcvd <- log(df$ulcvd)
df$q <- p.adjust(df$p, method = "fdr")
df$qcancer <- p.adjust(df$pcancer, method = "fdr")
df$qcvd <- p.adjust(df$pcvd, method = "fdr")

df_long <- df %>%
  mutate(id = row_number()) %>%
  pivot_longer(cols = c(hr, ll, ul, hrcancer, llcancer, ulcancer, hrcvd, llcvd, ulcvd), names_to = "name", values_to = "value") %>%
  mutate(outcome = case_when(grepl("cancer", name) ~ "Cancer", grepl("cvd", name) ~ "CVD", TRUE ~ "All-cause"), metric = case_when(grepl("^hr", name) ~ "hr", grepl("^ll", name) ~ "ll", grepl("^ul", name) ~ "ul")) %>%
  select(-name) %>%
  pivot_wider(names_from = metric, values_from = value)

df_long <- df_long %>%
  mutate(p_use = case_when(outcome == "All-cause" ~ p, outcome == "Cancer" ~ pcancer, outcome == "CVD" ~ pcvd), q_use = case_when(outcome == "All-cause" ~ q, outcome == "Cancer" ~ qcancer, outcome == "CVD" ~ qcvd)) %>%
  mutate(point_group = case_when(p_use >= 0.05 ~ "Not significant", p_use < 0.05 & q_use >= 0.05 & hr > 1 ~ "Positive, P < 0.05", p_use < 0.05 & q_use >= 0.05 & hr < 1 ~ "Negative, P < 0.05", q_use < 0.05 & hr > 1 ~ "Positive, FDR < 0.05", q_use < 0.05 & hr < 1 ~ "Negative, FDR < 0.05")) %>%
  mutate(outcome = factor(outcome, levels = c("All-cause", "Cancer", "CVD"))) %>%
  arrange(id) %>%
  mutate(protein = factor(protein, levels = rev(unique(protein))))

ggplot(df_long, aes(x = protein, y = hr)) +
  geom_errorbar(aes(ymin = ll, ymax = ul, color = point_group), width = 0.3) +
  geom_point(aes(color = point_group, fill = point_group), shape = 21, size = 3, stroke = 0.8) +
  geom_hline(yintercept = 1, linetype = 2) +
  scale_y_log10(breaks = c(0.5, 1, 2), limits = c(0.4, 2.6)) +
  coord_flip() +
  facet_wrap(~ outcome, nrow = 1) +
  theme_bw(base_size = 8) +
  theme(panel.grid.major = element_line(linetype = "dashed", color = "grey85"), panel.grid.minor = element_line(linetype = "dashed", color = "grey85")) +
  scale_color_manual(name = "Significance", values = c("Not significant" = "grey70", "Positive, P < 0.05" = "#c25759", "Negative, P < 0.05" = "#599cb4", "Positive, FDR < 0.05" = "#c25759", "Negative, FDR < 0.05" = "#599cb4")) +
  scale_fill_manual(name = "Significance", values = c("Not significant" = "grey70", "Positive, P < 0.05" = "white", "Negative, P < 0.05" = "white", "Positive, FDR < 0.05" = "#c25759", "Negative, FDR < 0.05" = "#599cb4")) +
  labs(y = "Hazard Ratio (95% CI)", x = NULL)


###### Fig.3 ######
###fig3a
data <- read_excel("Source Data Fig.3.xlsx", sheet = "fig3a")
data<-as.data.frame(data)
data$q<-p.adjust(data$p,method="fdr")
data$qtime<-p.adjust(data$ptime,method="fdr")
data$qgroup<-p.adjust(data$pgroup,method="fdr")
mat <- data %>% select(protein, coeftime, coefgroup, coef)
mat <- mat %>% column_to_rownames("protein")
annotation_row <- data %>% select(protein, type1,type2) %>% column_to_rownames("protein")
ann_colors <- list(type1 = c("c"="#ccebc5", "ig"="#fde0ef","o"="#bebada"), type2 = c("pos" = "#c25759","neg" = "#599cb4"))

pheatmap(
  mat,
  scale = "none",
  cluster_rows = TRUE,
  cluster_cols = FALSE,
  clustering_method = "median",
  annotation_row = annotation_row,
  annotation_colors = ann_colors,
  fontsize_row = 10,
  border_color = "grey",
  color = colorRampPalette(c("#01665e", "white", "#b35806"))(100),
  breaks = seq(-0.4, 0.4, length.out = 101))

###fig3b
df <- read_excel("Source Data Fig.3.xlsx", sheet = "fig3b")
protein_vars <- c("stdlogp10643","stdlogp01700", "stdlogp01701", "stdlogp01834", "stdlogp01714", "stdlogp02747","stdlogp06312", "stdlogp04003",
                  "stdlogp08603", "stdlogp27918","stdlogp19652", "stdlogq9hc29", "stdlogp35579", "stdlogq92835", "stdlogp04004", "stdlogp10909")

for (pro in protein_vars) {
  df_summary <- df %>% group_by(time, failureall) %>% summarise(mean_protein = mean(.data[[pro]], na.rm = TRUE),.groups="drop")
  df_summary$time <- factor(df_summary$time, levels = c(1, 2, 3))
  p<- ggplot(df_summary, aes(x = time,y = mean_protein,group = failureall,color = as.factor(failureall))) +
    geom_line() +
    geom_point(shape = 15, size = 2) +
    labs( title = pro,x = "Time",y = "Protein Z-Score",color = "Failure Status") +
    scale_color_manual(values = c("#4c72b0", "#C44E52"),labels = c("Survivors", "Deceased")) +
    scale_x_discrete(expand = c(0.05, 0.05)) +
    scale_y_continuous(limits = c(-1, 1)) +
    theme_classic() +
    theme(legend.position = "none",
          axis.line = element_line(color = "black"),
          axis.ticks = element_line(linewidth = 0.5))
  print(p)}


###fig3c
df <- read_excel("Source Data Fig.3.xlsx", sheet = "fig3c")
df <-as.data.frame(df)

node_survivor <- tibble(node=unique(c(df$source,df$target))) %>%mutate(label=node,type3=case_when(grepl("^IG",node) ~ "ig", node %in% c("C4BPA","CFH") ~ "c",TRUE ~ "o"))
g_survivor <- graph_from_data_frame(d=df %>% select(source,target), vertices=node_survivor %>% select(node,label,type3),directed=FALSE)
V(g_survivor)$degree <- degree(g_survivor)
E(g_survivor)$beta <- df$beta
E(g_survivor)$sign <- ifelse(df$beta>0,"positive","negative")
type_colors <- c("c"="#1f78b4","ig"="#fb9a06","o"="#33a02c")
edge_colors <- c("positive"="#FF69B4","negative"="#009E9E")
degree_min <- min(V(g_survivor)$degree,na.rm=TRUE)
degree_max <- max(V(g_survivor)$degree,na.rm=TRUE)

ggraph(g_survivor,layout="kk") +
  geom_edge_link(alpha=0.6,aes(color=sign)) +
  geom_node_point(aes(color=type3,size=degree)) +
  geom_node_text(aes(label=label),repel=TRUE,size=3) +
  scale_color_manual(values=type_colors) +
  scale_size_continuous(limits=c(degree_min,degree_max),range=c(3,8), breaks=c(1,2,3)) +
  scale_edge_color_manual(values=edge_colors) +
  theme_void() +
  theme(legend.position="right")

###fig3d
df <- read_excel("Source Data Fig.3.xlsx", sheet = "fig3d")
df <-as.data.frame(df)

node_decedent <- tibble(node=unique(c(df$source,df$target)))
node_decedent <- node_decedent %>% mutate(label=node)
node_decedent <- node_decedent %>% mutate(type3=case_when(grepl("^IG",node) ~ "ig", node %in% c("C1QC","C4BPA","CFH","CFP","C7","VTN","CLU") ~ "c", TRUE ~ "o"))
g_decedent <- graph_from_data_frame(d=df %>% select(source,target),vertices=node_decedent %>% select(node,label,type3),directed=FALSE)
V(g_decedent)$degree <- degree(g_decedent)
E(g_decedent)$beta <- df$beta
E(g_decedent)$sign <- ifelse(df$beta>0,"positive","negative")
type_colors <- c("c"="#1f78b4","ig"="#fb9a06","o"="#33a02c")
edge_colors <- c("positive"="#FF69B4","negative"="#009E9E")
degree_min <- min(V(g_decedent)$degree,na.rm=TRUE)
degree_max <- max(V(g_decedent)$degree,na.rm=TRUE)

ggraph(g_decedent,layout="kk") + 
  geom_edge_link(alpha=0.6,aes(color=sign)) + 
  geom_node_point(aes(color=type3,size=degree)) + 
  geom_node_text(aes(label=label),repel=TRUE,size=3) + 
  scale_color_manual(values=type_colors) + 
  scale_size_continuous(limits=c(degree_min,degree_max),range=c(3,8),breaks=c(2,4,6)) + 
  scale_edge_color_manual(values=edge_colors) + 
  theme_void() + 
  theme(legend.position="right")

###fig3e,f
df <- read_excel("Source Data Fig.3.xlsx", sheet = "fig3e-h")
df <-as.data.frame(df)
df_long <- df %>% pivot_longer(cols=c(dec,sur),names_to="group",values_to="value") %>%mutate(group=factor(group,levels=c("sur","dec"),labels=c("Survivors","Decedents")))
df_top <- df_long %>%  filter(var %in% c("n_edges","mean degree","density")) %>%mutate(var=factor(var,levels=c("n_edges","mean degree","density"),labels=c("No. of edges","Mean degree","Density")))
df_edge <- df_long %>%  filter(var %in% c("n_edges_ig_ig","n_edges_c_c","n_edges_c_ig")) %>%mutate(var=factor(var,levels=c("n_edges_ig_ig","n_edges_c_c","n_edges_c_ig"),labels=c("Ig–Ig edges","C–C edges","C–Ig edges")))
group_cols <- c("Survivors"="#4E79A7","Decedents"="#E15759")

group_cols <- c("Survivors"="#7FA7C9","Decedents"="#E89A93")
group_cols <- c("Survivors"="#9BBBD6","Decedents"="#F0B2AB")
group_cols <- c("Survivors"="#6F9EC4","Decedents"="#D9827A")

ggplot(df_top,aes(x=group,y=value,fill=group)) +
  geom_col(width=0.65,color="black",linewidth=0.25) +
  geom_text(aes(label=value),vjust=-0.4,size=3.5) +
  facet_wrap(~var,nrow=1,scales="free_y") +
  scale_fill_manual(values=group_cols) +
  scale_y_continuous(expand=expansion(mult=c(0,0.15))) +
  labs(x=NULL,y=NULL) +
  theme_classic(base_size=12) +
  theme(legend.position="none",
    strip.background=element_blank(),
    strip.text=element_text(size=12,face="bold"),
    axis.text.x=element_text(size=10),
    axis.text.y=element_text(size=10))

ggplot(df_edge,aes(x=var,y=value,fill=group)) +
  geom_col(position=position_dodge(width=0.7),width=0.6,color="black",linewidth=0.25) +
  geom_text(aes(label=value),
            position=position_dodge(width=0.7),
            vjust=-0.4,size=3.5) +
  scale_fill_manual(values=group_cols) +
  scale_y_continuous(expand=expansion(mult=c(0,0.15))) +
  labs(x=NULL,y="No. of edges",fill=NULL) +
  theme_classic(base_size=12) +
  theme(legend.position="top",
    axis.text.x=element_text(size=11),
    axis.text.y=element_text(size=10))


###### Fig.4 ######
###fig4a
df_num <- read_excel("Source Data Fig.4.xlsx", sheet = "fig4a,c")
ggplot(df_num[1:4,], aes(x=factor(var, levels=var), y=num)) +
  geom_bar(stat="identity", fill="saddlebrown", width=0.8) +
  theme_bw() +
  labs(x=NULL, y=NULL) +
  scale_y_continuous(limits=c(0,10),breaks=seq(0,10,by=5)) +
  theme(axis.text.x=element_text(angle=45,hjust=1,color="black"),
        axis.text.y=element_text(color="black"),
        panel.border=element_rect(colour="black", fill=NA, linewidth=1),
        panel.grid=element_blank())

###fig4b
df_bubble <- read_excel("Source Data Fig.4.xlsx", sheet = "fig4b,d")
df_bubble <- as.data.frame(df_bubble)
tissues <- c("Lymphoid tissue", "Intestine", "Liver")
df_long <- melt(df_bubble, id.vars="protein", measure.vars=tissues)
df_long$value[is.na(df_long$value)] <- 0  # 填充NA值为0

ggplot(df_long, aes(x=factor(protein, levels=df_bubble$protein), y=variable)) +
  geom_point(aes(color=ifelse(value == 1, "saddlebrown", NA), size=5, alpha=ifelse(value == 1, 1, 0))) +
  scale_color_identity() +
  theme_bw() +
  labs(x=NULL, y=NULL) +
  theme(panel.grid.major=element_line(linetype="dashed", color="grey"),
        panel.grid.minor=element_blank(),
        axis.text.x=element_text(angle=90, hjust=1),
        text=element_text(color="black"), 
        legend.position="none",
        panel.border=element_rect(colour="black", fill=NA, linewidth=1))

###fig4c
df_num <- read_excel("Source Data Fig.4.xlsx", sheet = "fig4a,c")
ggplot(df_num[5:7,], aes(x=factor(var, levels=var), y=num)) +
  geom_bar(stat="identity", fill="forestgreen", width=0.8) +
  theme_bw() +
  labs(x=NULL, y=NULL) +
  scale_y_continuous(limits=c(0,10),breaks=seq(0,10,by=5)) +
  theme(axis.text.x=element_text(angle=45,hjust=1,color="black"),
        axis.text.y=element_text(color="black"),
        panel.border=element_rect(colour="black", fill=NA, linewidth=1),
        panel.grid=element_blank()) 

###fig4d
df_bubble <- read_excel("Source Data Fig.4.xlsx", sheet = "fig4b,d")
df_bubble <- as.data.frame(df_bubble)
tissues <- c("Plasma cells", "Hepatocytes", "Kupffer cells")
df_long <- melt(df_bubble, id.vars="protein", measure.vars=tissues)
df_long$value[is.na(df_long$value)] <- 0  # 填充NA值为0

ggplot(df_long, aes(x=factor(protein, levels=df_bubble$protein), y=variable)) +
  geom_point(aes(color=ifelse(value == 1, "forestgreen", NA), size=5, alpha=ifelse(value == 1, 1, 0))) +
  scale_color_identity() +
  theme_bw() +
  labs(x=NULL, y=NULL) +
  theme(panel.grid.major=element_line(linetype="dashed", color="grey"),
        panel.grid.minor=element_blank(),
        axis.text.x=element_text(angle=90, hjust=1),
        text=element_text(color="black"), 
        legend.position="none",
        panel.border=element_rect(colour="black", fill=NA, linewidth=1))

###fig4e
df<-read_excel("Source Data Fig.4.xlsx", sheet = "fig4e")
df<-as.data.frame(df)
df_plot <- df %>%  mutate(pathways = factor(pathways, levels = pathways[order(logq)]),fill_group = case_when(
  is.na(zscore) | zscore == 0 ~ "Z = 0 / NA",
  zscore <= -2               ~ "Z ≤ -2",
  zscore < 0 & zscore > -2   ~ "-2 < Z < 0",
  zscore > 2                ~ "Z > 2",
  zscore > 0 & zscore <= 2   ~ "0 < Z ≤ 2"))
fill_colors <- c( "Z ≤ -2"= "#00CED1","-2 < Z < 0"= "#99F0EB","Z = 0 / NA"  = "grey90", "0 < Z ≤ 2"   = "#fff7bc",  "Z > 2" = "#FFA500")

ggplot(df_plot, aes(x = logq, y = pathways, fill = fill_group)) +
  geom_col(width = 0.7, color = "black", linewidth = 0.2) +
  scale_fill_manual(values = fill_colors, name = "Z-score") +
  theme_bw() +
  theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank()) +
  labs( x = "log(p value)",y = "pathways"  )

mat <- df
row.names(mat)<-df[,1]
mat<-mat[,6:19]
my_colors <- c("#4F94CD", "white", "#CD5555")
my_breaks <- c(-1.5, -0.5, 0.5, 1.5)

pheatmap(mat, color = my_colors,
         breaks = my_breaks,
         cluster_rows = FALSE,
         cluster_cols = T,
         clustering_method = "ward.D",  
         border_color = "grey80",
         legend=FALSE,
         legend_breaks = c(-1, 0, 1),
         annotation_colors = ann_colors)

###fig4g
data<-read_excel("Source Data Fig.4.xlsx", sheet = "fig4f-h")
data<-as.data.frame(data)
data$code_id<-as.factor(data$code_id)
data_age <- data.frame(age = seq(min(data$age, na.rm=TRUE),max(data$age, na.rm=TRUE),length.out = 200))

fit_pathway1 <- gamm(pathway_score1 ~ s(age, k = 4), random = list(code_id = ~1), data = data)
summary(fit_pathway1$gam)
fit_pathway2 <- gamm(pathway_score2 ~ s(age, k = 4), random = list(code_id = ~1), data = data)
summary(fit_pathway2$gam)

pred_pathway1 <- predict(fit_pathway1$gam, data_age, se.fit = TRUE)
data_age$pathway1 <- pred_pathway1$fit
data_age$upper_pathway1 <-pred_pathway1$fit + 1.96 * pred_pathway1$se.fit
data_age$lower_pathway1 <-pred_pathway1$fit - 1.96 * pred_pathway1$se.fit

pred_pathway2 <- predict(fit_pathway2$gam, data_age, se.fit = TRUE)
data_age$pathway2 <- pred_pathway2$fit
data_age$upper_pathway2 <-pred_pathway2$fit + 1.96 * pred_pathway2$se.fit
data_age$lower_pathway2 <-pred_pathway2$fit - 1.96 * pred_pathway2$se.fit

ggplot() +
  geom_line(data = data_age, aes(x = age, y = pathway1), size = 1.2, color = "#1F77B4") +  # 蓝色线
  geom_ribbon(data = data_age, aes(x = age, ymin = lower_pathway1, ymax = upper_pathway1), alpha = 0.2, fill = "#1F77B4") +  # 蓝色置信区间
  geom_line(data = data_age, aes(x = age, y = pathway2), size = 1.2, color = "#FF7F0E") +  # 橙色线
  geom_ribbon(data = data_age, aes(x = age, ymin = lower_pathway2, ymax = upper_pathway2), alpha = 0.2, fill = "#FF7F0E") +  # 橙色置信区间
  scale_y_continuous(limits = c(0, 2)) +
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  labs(x = "Age", y = "Pathway Score")

###fig4h
ggplot(data, aes(x = factor(time), y = pathway_score1, fill = factor(death))) +
  geom_boxplot(outliers = FALSE) +
  labs(x = "Time", y = "score_pathway1", fill = "Death") +
  scale_fill_manual(values = c("0" = "#5D8B8E", "1" = "#FFA07A"))+
  scale_y_continuous(limits=c(-8,8))+
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) 

###fig4i
ggplot(data, aes(x = factor(time), y = pathway_score2, fill = factor(death))) +
  geom_boxplot(outliers = FALSE) +
  labs(x = "Time", y = "score_pathway2", fill = "Death") +
  scale_fill_manual(values = c("0" = "#5D8B8E", "1" = "#FFA07A"))+
  scale_y_continuous(limits=c(-8,8))+
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) 



###### Fig.5 ######
###fig5a
df_coef<-read_excel("Source Data Fig.5.xlsx", sheet = "fig5a")
df_coef<-as.data.frame(df_coef)
df_coef$protein <- factor(df_coef$protein, levels = rev(df_coef$protein))

ggplot(df_coef, aes(y = protein)) +
  geom_tile(aes(x = -0.06, fill = type3), width = 0.035,height = 0.7, show.legend = TRUE) +
  geom_col( aes(x = abs(lnhr), fill = lnhr > 0),width = 0.7) +
  scale_fill_manual(values = c("TRUE"  = "#C44E52","FALSE" = "#4C72B0",acti = "#fdc086", regu = "#a6cee3"),
                    breaks = c("TRUE", "FALSE", "ig", "c", "o"),
                    labels = c("Positive","Negative", "Immunoglobulins","Complements","Others")) +
  scale_x_continuous( limits = c(-0.1, 0.5), breaks = seq(0, 0.5, by = 0.1)) +
  labs(x = "Weight (|ln HR|)",y = NULL) +
  theme_bw() +
  theme(panel.grid.major.y = element_blank(),panel.grid.minor = element_blank(),axis.ticks.y = element_blank())

###fig5b
df<-read_excel("Source Data Fig.5.xlsx", sheet = "fig5b-e")
df<-as.data.frame(df)
df$time <- factor(df$time, levels = c(1, 2, 3))
df<-df %>% mutate(death = ifelse(is.na(death), 0, death))

ggplot(df, aes(x = stdscore_regu, y = stdscore_acti)) +
  geom_point(data = df[df$death == 0, ], aes(shape=time),color= "grey90", fill = "grey90",alpha = 0.6, size = 1.5) +
  geom_point(data = df[df$death == 1, ], aes(shape=time),color= "#c44e52",fill= "#c44e52",alpha = 0.6,size = 1.5) +
  geom_density_2d(data=df[df$death==0, ],aes(x = stdscore_regu,y=stdscore_acti),color= "#1f77b4", alpha=0.7, linewidth=0.8,breaks=c(0.05,0.1)) +
  geom_density_2d(data=df[df$death==1, ],aes(x = stdscore_regu,y=stdscore_acti),color= "red3", alpha=0.85, linewidth=0.95, breaks=c(0.05,0.1)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", linewidth = 0.6,color = "black") +
  theme_bw() +
  theme(panel.grid = element_blank()) +
  coord_cartesian(xlim = c(-5, 5),ylim = c(-5, 5)) +
  scale_shape_manual(values = c("1" = 15, "2" = 16, "3" = 17),name = "Time") +
  labs(x = "Standardized Regulation Score",y = "Standardized Activation Score")

###fig5c
ggplot(df, aes(x = time, y = stdscore_imba, fill = factor(death))) +
  geom_boxplot(width = 0.7, alpha = 0.85,outliers=FALSE) +
  scale_fill_manual(values = c("0" = "#4c72b0", "1" = "#c44e52")) +
  labs(title = "CDS", x = "Time", y = "stdscore_imba") +
  coord_cartesian(ylim = c(-5, 5)) +
  theme_bw()+
  theme(panel.grid = element_blank())

###fig5d
ggplot(df, aes(x = time, y = stdscore_acti, fill = factor(death))) +
  geom_boxplot(width = 0.7, alpha = 0.85,outliers=FALSE) +
  scale_fill_manual(values = c("0" = "#4c72b0", "1" = "#c44e52")) +
  labs(title = "CAS", x = "Time", y = "stdscore_imba") +
  coord_cartesian(ylim = c(-5, 5)) +
  theme_bw()+
  theme(panel.grid = element_blank())

###fig5e
ggplot(df, aes(x = time, y = stdscore_regu, fill = factor(death))) +
  geom_boxplot(width = 0.7, alpha = 0.85,outliers=FALSE) +
  scale_fill_manual(values = c("0" = "#4c72b0", "1" = "#c44e52")) +
  labs(title = "CRS", x = "Time", y = "stdscore_imba") +
  coord_cartesian(ylim = c(-5, 5)) +
  theme_bw()+
  theme(panel.grid = element_blank())

###fig5f
df_forest<-read_excel("Source Data Fig.5.xlsx", sheet = "fig5f,g")
df_forest<-as.data.frame(df_forest)

df_forest %>% filter(cohort=="discovery") %>%  ggplot(aes(y = score, x = log(hr))) + 
  geom_point(aes(fill = score),  shape = 22,size = 5, stroke = 1.2) +
  geom_errorbar(aes(xmin = log(ll), xmax = log(ul)),width = 0.3,linewidth = 0.9) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey50", linewidth = 0.7) + 
  scale_fill_manual(values = c("CDS" = "#d73027","CAS" = "#fdae61","CRS" = "#1a9850")) +
  scale_y_discrete(limits = c("CRS", "CAS", "CDS")) +
  scale_x_continuous(limits = c(log(0.49), log(2.1)),breaks = log(c(0.5, 1, 2)),labels = c("0.5", "1", "2")) +
  labs(title = "Discovery cohort",x = "Hazard Ratio (95% CI)",y = "") +
  theme_bw() +
  theme(panel.grid = element_blank(),
        axis.line = element_line(color = "black"),
        axis.ticks = element_line(color = "black"),
        legend.position = "none")

###fig5g
df_forest %>% filter(cohort=="validation") %>%  ggplot(aes(y = score, x = log(hr))) + 
  geom_point(aes(fill = score),  shape = 22,size = 5, stroke = 1.2) +
  geom_errorbar(aes(xmin = log(ll), xmax = log(ul)),width = 0.3,linewidth = 0.9) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey50", linewidth = 0.7) + 
  scale_fill_manual(values = c("CDS" = "#d73027","CAS" = "#fdae61","CRS" = "#1a9850")) +
  scale_y_discrete(limits = c("CRS", "CAS", "CDS")) +
  scale_x_continuous(limits = c(log(0.49), log(2.1)),breaks = log(c(0.5, 1, 2)),labels = c("0.5", "1", "2")) +
  labs(title = "Validation cohort",x = "Hazard Ratio (95% CI)",y = "") +
  theme_bw() +
  theme(panel.grid = element_blank(),
        axis.line = element_line(color = "black"),
        axis.ticks = element_line(color = "black"),
        legend.position = "none")

###fig5h
df_imba<-read_excel("Source Data Fig.5.xlsx", sheet = "fig5h")
df_imba<-as.data.frame(df_imba)
ggplot(df_imba, aes(age_center, mean_score)) +
  geom_line(linewidth = 1.4, color = "#F4A582") +
  geom_ribbon(aes(ymin = lower, ymax = upper),fill = "#F4A582", alpha = 0.25) +
  geom_vline(xintercept = 64.4,linetype = "dashed",linewidth = 1,color = "grey50") +
  geom_rug(data = df[df$age >= 45 & df$age <= 80,],aes(x = age),inherit.aes = FALSE,sides = "b",alpha = 0.2,color="#00693433")+
  scale_x_continuous(limits = c(45, 80), breaks = c(45,50, 60, 70, 80)) +
  scale_y_continuous(limits=c(-0.5,1),breaks = c(-0.5, 0, 0.5,1)) +
  labs(x = "Age (years)", y = "Imbalance score")+
  theme_bw() +
  theme(panel.grid = element_blank(),
        axis.line = element_line(color = "black"),
        axis.ticks = element_line(color = "black"),
        axis.text = element_text(color = "black"),
        axis.title = element_text(color = "black")) 

###fig5i
df_acti_regu<-read_excel("Source Data Fig.5.xlsx", sheet = "fig5i")
df_acti_regu<-as.data.frame(df_acti_regu)
ggplot(df_acti_regu, aes(age_center, mean_score, color = score, fill = score)) +
  geom_line(linewidth = 1) +
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.2, color = NA) +
  geom_vline(xintercept = 64.2, linetype = "dashed", linewidth = 1, color = "grey50") +
  geom_vline(xintercept = 64.7, linetype = "dashed", linewidth = 1, color = "grey50") +
  geom_rug(data = df[df$age >= 45 & df$age <= 80,],aes(x = age),inherit.aes = FALSE,sides = "b",alpha = 0.2,color="#00693433")+
  scale_color_manual(values = c("stdscore_acti" = "#E6AB02", "stdscore_regu" = "#92C5DE")) +
  scale_fill_manual(values = c("stdscore_acti" = "#E6AB02", "stdscore_regu" = "#92C5DE")) +
  scale_x_continuous(limits = c(45, 80), breaks = c(45,50, 60, 70, 80)) +
  scale_y_continuous(limits=c(-0.5,1),breaks = c(-0.5, 0, 0.5,1)) +
  labs(x = "Age (years)", y = "Imbalance score")+
  theme_bw() +
  theme(panel.grid = element_blank(),
        axis.line = element_line(color = "black"),
        axis.ticks = element_line(color = "black"),
        axis.text = element_text(color = "black"),
        axis.title = element_text(color = "black"))

###fig5j
hr_imba<-read_excel("Source Data Fig.5.xlsx", sheet = "fig5j")
hr_imba<-as.data.frame(hr_imba)
ggplot(hr_imba, aes(age, HR)) +
  geom_line(linewidth = 1.2,color="#f4a582") +
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.2,fill="#f4a582") +
  geom_hline(yintercept = 1, linetype = 2) +
  geom_vline(xintercept = 65.5,linetype = "dashed",linewidth = 1,color = "grey50") +
  geom_rug(data = df[df$age >= 45 & df$age <= 80,],aes(x = age),inherit.aes = FALSE,sides = "b",alpha = 0.2,color="#00693433")+
  scale_y_continuous(limits=c(0,5),breaks = c(0,1,2,3,4,5)) +
  theme_bw()+
  theme(panel.grid = element_blank(),
        axis.line = element_line(color = "black"),
        axis.ticks = element_line(color = "black"),
        axis.text = element_text(color = "black"),
        axis.title = element_text(color = "black"),
        plot.title = element_text(color = "black"))

###fig5k
hr_acti_regu<-read_excel("Source Data Fig.5.xlsx", sheet = "fig5k")
hr_acti_regu<-as.data.frame(hr_acti_regu)
ggplot(hr_acti_regu, aes(age, HR, color = score, fill = score)) +
  geom_line(linewidth = 1) +
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.15, color = NA) +
  geom_hline(yintercept = 1, linetype = 2) +
  geom_vline(xintercept = 62.5,linetype = "dashed",linewidth = 1,color = "grey50") +
  geom_vline(xintercept = 63.5,linetype = "dashed",linewidth = 1,color = "grey50") +
  geom_rug(data = df[df$age >= 45 & df$age <= 80,],aes(x = age),inherit.aes = FALSE,sides = "b",alpha = 0.2,color="#00693433")+
  scale_color_manual(values = c("stdscore_acti" = "#e6ab02","stdscore_regu" = "#92c5de")) +
  scale_fill_manual(values = c("stdscore_acti" = "#e6ab02","stdscore_regu" = "#92c5de")) +
  scale_y_continuous(limits=c(0,5),breaks = c(0,1,2,3,4,5)) +
  theme_bw()+
  theme(panel.grid = element_blank(),
        axis.line = element_line(color = "black"),
        axis.ticks = element_line(color = "black"),
        axis.text = element_text(color = "black"),
        axis.title = element_text(color = "black"),
        plot.title = element_text(color = "black"))

###fig5l
df_auc<-read_excel("Source Data Fig.5.xlsx", sheet = "fig5l")
df_auc<-as.data.frame(df_auc)

ggplot(df_auc, aes(x = time, y = AUC, color = marker, group = marker)) +
  geom_line(linewidth = 1.2) +
  geom_point(size = 2,shape=16) +
  scale_color_manual(values = c("stdscore_imba" = "#B2182B","stdscore_acti" = "#F4A582","stdscore_regu" = "#92C5DE")) +
  scale_x_continuous(breaks = unique(df_auc$time)) +
  labs(x = "Time", y = "Time-dependent AUC", color = "Score") +
  scale_y_continuous(limits = c(0.4, 1), breaks = seq(0.2, 1, 0.2)) +
  theme_bw() +
  theme(
    panel.grid = element_blank(),
    axis.line = element_line(color = "black"),
    axis.ticks = element_line(color = "black"),
    axis.text = element_text(color = "black"),
    axis.title = element_text(color = "black"))

###fig5m
df_landmark<-read_excel("Source Data Fig.5.xlsx", sheet = "fig5m")
df_landmark<-as.data.frame(df_landmark)

ggplot()+ 
  geom_segment(data = df_landmark, aes(x = x_start, xend = x_end,y = AUC, yend = AUC,color = marker),linewidth = 1.2) +
  geom_segment(data = df_landmark, aes(x = x, xend = x,y = AUC_left, yend = AUC_right,color = marker),linewidth = 1) +
  scale_color_manual(values = c("stdscore_imba" = "#B2182B","stdscore_acti" = "#F4A582","stdscore_regu" = "#92C5DE")) +
  scale_x_continuous(limits = c(2, 16),breaks = seq(2, 16, by = 2)) +
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.2)) +
  labs(x = "Follow-up time (years)", y = "Landmark AUC", color = "Score") +
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())



###### Fig.6 ######

###fig6a
df<-read_excel("Source Data Fig.6.xlsx", sheet = "fig6a")
df<-as.data.frame(df)
df$Score <- factor(df$Score, levels = c("CImbS", "CActS", "CRegS"))
df$y <- 4 - as.numeric(df$Score)
label_df <- df %>%  mutate(label = paste0(sprintf("%.2f", coef)," (", sprintf("%.2f", ll),", ", sprintf("%.2f", ul),"), p=", p) )

ggplot(label_df, aes(x = coef, y = y)) +
  geom_errorbar(aes(xmin = ll, xmax = ul), width = 0.15, linewidth = 0.5) +
  geom_point(aes(fill = factor(y)), shape = 22, size = 2.8, color = "black", stroke = 0.4) +
  geom_vline(xintercept = 0, linetype = "dashed", linewidth = 0.4) +
  geom_text(aes(label = label), hjust = -0.05, size = 3) +
  scale_y_continuous(breaks = c(3, 2, 1), labels = c("CImbS", "CActS", "CRegS")) +
  scale_fill_manual(values = c("3" = "#d73027","2" = "#fdae61","1" = "#1a9850"),guide = "none") +
  labs(x = "Coefficient (95% CI)", y = NULL) +
  theme_classic(base_size = 12) +
  theme(axis.line.y = element_blank(),axis.ticks.y = element_blank()) +
  coord_cartesian(xlim = c(min(df$ll, na.rm = TRUE) - 0.02, max(df$ul, na.rm = TRUE) + 0.08),clip = "off"  )

###fig6b
df <- read_excel("Source Data Fig.6.xlsx", sheet = "fig6b")
df <- as.data.frame(df)
df$`Nutrient replaced` <- factor(df$`Nutrient replaced`,levels = c("Fat", "Protein", "Fat and protein"))
df$group <- factor(df$group, levels = c("CImbS", "CRegS"))

ggplot(df, aes(x = coef, y = `Nutrient replaced`)) +
  geom_vline(xintercept = 0, linetype = "dashed", linewidth = 0.4, color = "black") +
  geom_errorbar(aes(xmin = ll, xmax = ul), width = 0.15, color = "black") +
  geom_point(aes(fill = group), shape = 22, size = 2.8, color = "black", stroke = 0.5) +
  facet_grid(. ~ group) +
  scale_fill_manual(values = c("CImbS" = "#d73027", "CRegS" = "#1a9850"), guide = "none") +
  labs(x = "Coefficient (95% CI)", y = NULL) +
  theme_classic(base_size = 12) +
  theme(strip.background = element_blank(), strip.text = element_text(face = "bold"), axis.line.y = element_blank(), axis.ticks.y = element_blank(), panel.border = element_rect(color = "black", fill = NA, linewidth = 0.6)) +
  coord_cartesian(xlim = c(min(df$ll, na.rm = TRUE) - 0.02, max(df$ul, na.rm = TRUE) + 0.03), clip = "off")

###fig6c
df <- read_excel("Source Data Fig.6.xlsx", sheet = "fig6c")
df <- as.data.frame(df)
df$Protein <- factor(df$Protein, levels = c("CFH", "C4BPA", "VTN", "CLU"))
df$y <- 5 - as.numeric(df$Protein)
df$estimate_text <- sprintf("%.2f (%.2f, %.2f)", df$coef, df$ll, df$ul)
df$p_text <- ifelse(df$p < 0.001, "<0.001", sprintf("%.2f", df$p))

ggplot(df, aes(x = coef, y = y)) +
  geom_errorbar(aes(xmin = ll, xmax = ul), width = 0.15, linewidth = 0.5) +
  geom_point(shape = 21, size = 3, fill = "#00a0e9", color = "black", stroke = 0.5) +
  geom_vline(xintercept = 0, linetype = "dashed", linewidth = 0.4) +
  geom_text(aes(x = 0.125, label = estimate_text),hjust = 0, size = 3.5) +geom_text(aes(x = 0.255, label = p_text),hjust = 0.5, size = 3.5) +
  annotate("text", x = 0.125, y = 4.65,label = "Coefficient (95% CI)",hjust = 0, size = 3.5)+annotate("text", x = 0.255, y = 4.65,label = "P",hjust = 0.5, size = 3.5) +
  scale_y_continuous(breaks = c(4, 3, 2, 1),labels = c("CFH", "C4BPA", "VTN", "CLU"), limits = c(0.5, 4.8)) +
  scale_x_continuous(breaks = c(-0.1, 0, 0.1)) +
  labs(x = "Coefficient (95% CI)", y = NULL) +
  theme_classic(base_size = 12) +
  theme(axis.line.y = element_blank(),axis.ticks.y = element_blank(),plot.margin = margin(5.5, 180, 5.5, 5.5)) +
  coord_cartesian(xlim = c(-0.1, 0.1), clip = "off")


###fig6d
df <- read_excel("Source Data Fig.6.xlsx", sheet = "fig6d")
df <- as.data.frame(df)

ggplot() +
  geom_ribbon(data = df, aes(x = eacarbo, ymin = lci, ymax = uci), alpha = 0.25) +
  geom_line(data = df, aes(x = eacarbo, y = fit), linewidth = 1) +
  scale_x_continuous(limits = c(50, 400), breaks = seq(50, 400, 50), labels = function(x) sprintf("%.0f", x)) +
  coord_cartesian(ylim = c(-1, 1)) +
  theme_classic() +
  labs(x = "Energy-adjusted carbohydrate intake", y = "VTN")

###fig6e
df <- read_excel("Source Data Fig.6.xlsx", sheet = "fig6e")
df <- as.data.frame(df)

ggplot() +
  geom_ribbon(data = df, aes(x = eacarbo, ymin = lci, ymax = uci), alpha = 0.25) +
  geom_line(data = df, aes(x = eacarbo, y = fit), linewidth = 1) +
  scale_x_continuous(limits = c(50, 400), breaks = seq(50, 400, 50), labels = function(x) sprintf("%.0f", x)) +
  coord_cartesian(ylim = c(-1, 1)) +
  theme_classic() +
  labs(x = "Energy-adjusted carbohydrate intake", y = "CLU")

###fig6f
df <- read_excel("Source Data Fig.6.xlsx", sheet = "fig6f", col_names = FALSE)
df <- as.data.frame(df)
df <- df[rowSums(!is.na(df)) > 0, 1:5]
names(df) <- c("Exposure", "T1_label", "T2_label", "T3_label", "P_interaction")

label_data <- as.matrix(df)
label_data[is.na(label_data)] <- ""

extract_num <- function(x, k) {
  m <- regexec("\\s*([0-9.]+)\\s*\\(([0-9.]+),\\s*([0-9.]+)\\)", as.character(x))
  res <- regmatches(as.character(x), m)
  out <- rep(NA_real_, length(x))
  ok <- lengths(res) == 4
  out[ok] <- as.numeric(vapply(res[ok], function(z) z[k + 1], character(1)))
  return(out)}

HR <- data.frame(
  T1_beta = extract_num(df$T1_label, 1),
  T1_lci = extract_num(df$T1_label, 2),
  T1_uci = extract_num(df$T1_label, 3),
  T2_beta = extract_num(df$T2_label, 1),
  T2_lci = extract_num(df$T2_label, 2),
  T2_uci = extract_num(df$T2_label, 3),
  T3_beta = extract_num(df$T3_label, 1),
  T3_lci = extract_num(df$T3_label, 2),
  T3_uci = extract_num(df$T3_label, 3))

mean_matrix <- as.matrix(HR[, c("T1_beta", "T2_beta", "T3_beta")])
lower_matrix <- as.matrix(HR[, c("T1_lci", "T2_lci", "T3_lci")])
upper_matrix <- as.matrix(HR[, c("T1_uci", "T2_uci", "T3_uci")])

xticks <- c(0.5, 1, 5)
attr(xticks, "labels") <- sprintf("%.1f", xticks)

forestplot(
  labeltext = label_data,
  mean = mean_matrix,
  lower = lower_matrix,
  upper = upper_matrix,
  zero = 1,
  lwd.zero = 1,
  graph.pos = 2,
  boxsize = 0.1,
  lineheight = unit(18, "mm"),
  colgap = unit(5, "mm"),
  align = "c",
  xlog = TRUE,
  fn.ci_norm = fpDrawNormalCI,
  xticks = xticks,
  graphwidth = unit(90, "mm"),
  hrzl_lines = list("2" = gpar(lty = 1, columns = c(3:5), col = "black"),"3" = gpar(lty = 1, columns = c(1:6), col = "black")),
  txt_gp = fpTxtGp(
    label = gpar(cex = 1),
    ticks = gpar(cex = 1),
    xlab = gpar(cex = 1),
    title = gpar(cex = 1),
    legend = gpar(cex = 0.85),
    legend.title = gpar(cex = 0.1)),
  legend = c("T1", "T2", "T3"),
  legend_args = fpLegend(pos = list("topright"), r = unit(NA, "snpc")),
  col = fpColors(box = c("#7570b3", "#1b9e77", "#d95f02"),line = c("grey40", "grey40", "grey40")),
  ci.vertices = TRUE,
  ci.vertices.height = 0.05,
  ci.vertices.tot = 0.1) |>
  fp_add_lines(h_4 = gpar(lty = 2), h_5 = gpar(lty = 2),h_6 = gpar(lty = 2)) |>
  fp_set_zebra_style("#EEE")


###### Fig.7 ######
###fig7b
df_long <- read_excel("Source Data Fig.7.xlsx", sheet = "fig7b-d")
df_long <- as.data.frame(df_long)
df_long <- df_long %>% mutate(week=case_when(group==1 & time==1 ~ 1, group==1 & time==2 ~ 2, group==2 & time==1 ~ 3, group==2 & time==2 ~ 4))
summary_df <- df_long %>% summarise(across(c(stdscore_imba,stdscore_acti,stdscore_regu),list(mean=~mean(.x,na.rm=TRUE),se=~sd(.x,na.rm=TRUE)/sqrt(sum(!is.na(.x))),n=~sum(!is.na(.x))),.names="{.col}_{.fn}"),.by=week)
summary_imba <- summary_df %>% transmute(week,mean=stdscore_imba_mean,se=stdscore_imba_se,n=stdscore_imba_n,ci=qt(0.975,n-1)*se)
summary_acti <- summary_df %>% transmute(week,mean=stdscore_acti_mean,se=stdscore_acti_se,n=stdscore_acti_n,ci=qt(0.975,n-1)*se)
summary_regu <- summary_df %>% transmute(week,mean=stdscore_regu_mean,se=stdscore_regu_se,n=stdscore_regu_n,ci=qt(0.975,n-1)*se)

ggplot(df_long,aes(x=week,y=stdscore_imba,group=code_id)) +
  geom_line(color="grey90",linewidth=0.2) +
  geom_point(color="grey90",size=4,shape=16) +
  geom_errorbar(data=summary_imba,aes(x=week,ymin=mean-ci,ymax=mean+ci),inherit.aes=FALSE,width=0.12,color="#cb3127",linewidth=0.8) +
  geom_line(data=summary_imba,aes(x=week,y=mean,group=1),inherit.aes=FALSE,color="#cb3127",linewidth=1.1) +
  geom_point(data=summary_imba,aes(x=week,y=mean),inherit.aes=FALSE,color="#cb3127",size=4,shape=16) +
  scale_x_continuous(breaks=1:4,labels=c("HFLC before","HFLC after","LFHC before","LFHC after")) +
  scale_y_continuous(limits=c(-3,3)) +
  labs(x=NULL,y="CImbS") +
  theme_classic() +
  theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank(),
        legend.position="none",panel.border=element_blank(),
        axis.line=element_line(color="black"),
        axis.text=element_text(color="black"),axis.title=element_text(color="black"),
        axis.ticks=element_line(color="black"))

###fig7c
ggplot(df_long,aes(x=week,y=stdscore_acti,group=code_id)) +
  geom_line(color="grey90",linewidth=0.6) +
  geom_point(color="grey90",size=4,shape=16) +
  geom_errorbar(data=summary_acti,aes(x=week,ymin=mean-ci,ymax=mean+ci),inherit.aes=FALSE,width=0.12,color="#f4a961",linewidth=0.8) +
  geom_line(data=summary_acti,aes(x=week,y=mean,group=1),inherit.aes=FALSE,color="#f4a961",linewidth=1.1) +
  geom_point(data=summary_acti,aes(x=week,y=mean),inherit.aes=FALSE,color="#f4a961",size=4,shape=16) +
  scale_x_continuous(breaks=1:4,labels=c("HFLC before","HFLC after","LFHC before","LFHC after")) +
  scale_y_continuous(limits=c(-3,3)) +
  labs(x=NULL,y="CActS") +
  theme_classic() +
  theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank(),
        legend.position="none",panel.border=element_blank(),
        axis.line=element_line(color="black"),
        axis.text=element_text(color="black"),axis.title=element_text(color="black"),
        axis.ticks=element_line(color="black"))

###fig7d
ggplot(df_long,aes(x=week,y=stdscore_regu,group=code_id)) +
  geom_line(color="grey90",linewidth=0.6) +
  geom_point(color="grey90",size=4,shape=16) +
  geom_errorbar(data=summary_regu,aes(x=week,ymin=mean-ci,ymax=mean+ci),inherit.aes=FALSE,width=0.12,color="#1d924e",linewidth=0.8) +
  geom_line(data=summary_regu,aes(x=week,y=mean,group=1),inherit.aes=FALSE,color="#1d924e",linewidth=1.1) +
  geom_point(data=summary_regu,aes(x=week,y=mean),inherit.aes=FALSE,color="#1d924e",size=4,shape=16) +
  scale_x_continuous(breaks=1:4,labels=c("HFLC before","HFLC after","LFHC before","LFHC after")) +
  scale_y_continuous(limits=c(-3,3)) +
  labs(x=NULL,y="CActS") +
  theme_classic() +
  theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank(),
        legend.position="none",panel.border=element_blank(),
        axis.line=element_line(color="black"),
        axis.text=element_text(color="black"),axis.title=element_text(color="black"),
        axis.ticks=element_line(color="black"))

###fig7e
df_wide <- read_excel("Source Data Fig.7.xlsx", sheet = "fig7e-h")
plot_imba <- df_wide %>% select(dg1_stdscore_imba,dg2_stdscore_imba) %>% pivot_longer(cols=everything(),names_to="diet",values_to="value") %>%
  mutate(diet=case_when(diet=="dg1_stdscore_imba" ~ "HFLC", diet=="dg2_stdscore_imba" ~ "LFHC"),diet=factor(diet,levels=c("HFLC","LFHC")))

ggplot(plot_imba,aes(x=diet,y=value,fill=diet,color=diet)) +
  geom_violin(width=0.5,alpha=0.10,trim=FALSE,linewidth=0.3,color="grey80") +
  geom_boxplot(width=0.25,outlier.shape=NA,alpha=0.75,linewidth=0.5,color="black") +
  geom_jitter(width=0.1,size=1.4,alpha=0.65,shape=16) +
  scale_fill_manual(values=c(HFLC="#1F77B4",LFHC="#E68E15")) +
  scale_color_manual(values=c(HFLC="#1F77B4",LFHC="#E68E15")) +
  scale_y_continuous(limits=c(-5.5,5),breaks = c(-5,-2.5,0,2.5,5)) +
  labs(x=NULL,y=expression(Delta*" CImbS")) +
  theme_classic() +
  theme(legend.position="none",
        axis.line=element_line(color="black",linewidth=0.5),
        axis.text=element_text(color="black"),
        axis.title=element_text(color="black"),
        axis.ticks=element_line(color="black"))

###fig7f
plot_acti <- df_wide %>% select(dg1_stdscore_acti,dg2_stdscore_acti) %>% pivot_longer(cols=everything(),names_to="diet",values_to="value") %>%
  mutate(diet=case_when(diet=="dg1_stdscore_acti" ~ "HFLC", diet=="dg2_stdscore_acti" ~ "LFHC"),diet=factor(diet,levels=c("HFLC","LFHC")))

ggplot(plot_acti,aes(x=diet,y=value,fill=diet,color=diet)) +
  geom_violin(width=0.5,alpha=0.10,trim=FALSE,linewidth=0.3,color="grey80") +
  geom_boxplot(width=0.25,outlier.shape=NA,alpha=0.75,linewidth=0.5,color="black") +
  geom_jitter(width=0.1,size=1.4,alpha=0.65,shape=16) +
  scale_fill_manual(values=c(HFLC="#1F77B4",LFHC="#E68E15")) +
  scale_color_manual(values=c(HFLC="#1F77B4",LFHC="#E68E15")) +
  scale_y_continuous(limits=c(-5,5)) +
  labs(x=NULL,y=expression(Delta*" CActS")) +
  theme_classic() +
  theme(legend.position="none",
        axis.line=element_line(color="black",linewidth=0.5),
        axis.text=element_text(color="black"),
        axis.title=element_text(color="black"),
        axis.ticks=element_line(color="black"))

###fig7g
plot_regu <- df_wide %>% select(dg1_stdscore_regu,dg2_stdscore_regu) %>% pivot_longer(cols=everything(),names_to="diet",values_to="value") %>%
  mutate(diet=case_when(diet=="dg1_stdscore_regu" ~ "HFLC", diet=="dg2_stdscore_regu" ~ "LFHC"),diet=factor(diet,levels=c("HFLC","LFHC")))

ggplot(plot_regu,aes(x=diet,y=value,fill=diet,color=diet)) +
  geom_violin(width=0.5,alpha=0.10,trim=FALSE,linewidth=0.3,color="grey80") +
  geom_boxplot(width=0.25,outlier.shape=NA,alpha=0.75,linewidth=0.5,color="black") +
  geom_jitter(width=0.1,size=1.4,alpha=0.65,shape=16) +
  scale_fill_manual(values=c(HFLC="#1F77B4",LFHC="#E68E15")) +
  scale_color_manual(values=c(HFLC="#1F77B4",LFHC="#E68E15")) +
  scale_y_continuous(limits=c(-5,5)) +
  labs(x=NULL,y=expression(Delta*" CRegS")) +
  theme_classic() +
  theme(legend.position="none",
        axis.line=element_line(color="black",linewidth=0.5),
        axis.text=element_text(color="black"),
        axis.title=element_text(color="black"),
        axis.ticks=element_line(color="black"))

###fig7h
protein_list <- c("stdlogp08603","stdlogp04003","stdlogp04004","stdlogp10909")
protein_lab <- c(stdlogp08603="CFH",stdlogp04003="C4BPA",stdlogp04004="VTN",stdlogp10909="CLU")

plot_df <- df_wide %>%  mutate(id=row_number()) %>%
  select(id,all_of(c(paste0("dg1_",protein_list),paste0("dg2_",protein_list)))) %>%
  pivot_longer(cols=-id,names_to=c("diet_code","protein"),names_pattern="^(dg[12])_(.*)$",values_to="value") %>%
  mutate(diet=case_when(diet_code=="dg1" ~ "HFLC",diet_code=="dg2" ~ "LFHC"), diet=factor(diet,levels=c("HFLC","LFHC")),
         protein_label=factor(protein_lab[protein],levels=protein_lab[protein_list]))

ggplot(plot_df,aes(x=diet,y=value)) +
  geom_hline(yintercept=0,linetype="dashed",color="grey50",linewidth=0.5) +
  geom_boxplot(aes(fill=diet),width=0.55,outlier.shape=NA,alpha=0.45,linewidth=0.5,color="black") +
  geom_jitter(aes(color=diet),width=0.12,size=1.3,alpha=0.75,shape=16) +
  facet_wrap(~protein_label,nrow=1,scales="fixed") +
  scale_fill_manual(values=c(HFLC="#1F77B4",LFHC="#E68E15")) +
  scale_color_manual(values=c(HFLC="#1F77B4",LFHC="#E68E15")) +
  scale_y_continuous(limits=c(-5,5)) +
  labs(x=NULL,y=expression("Changes in protein level")) +
  theme_classic() +
  theme(legend.position="none",
        strip.background=element_blank(),
        strip.text=element_text(color="black",face="bold",size=12),
        axis.line=element_line(color="black",linewidth=0.5),
        axis.text=element_text(color="black"),
        axis.title=element_text(color="black"),
        axis.ticks=element_line(color="black"))





###################################### codes for data analyses ################################################

#################### co-variation network analysis #######################
##### load package
library(tidyverse)
library(lmerTest)
library(igraph)
library(ggraph)
library(ggridges)

df <- read.csv("pro_longi.csv",header=TRUE)
df$code_id <- as.factor(df$code_id)
df$time <- as.factor(df$time)
df$sexadj <- as.factor(df$sexadj)

protein_list <- c("stdlogp10643","stdlogp01700","stdlogp01701","stdlogp01834","stdlogp01714","stdlogp02747","stdlogp06312","stdlogp04003",
                  "stdlogp08603","stdlogp27918","stdlogp19652","stdlogq9hc29","stdlogp35579","stdlogq92835","stdlogp04004","stdlogp10909")

dat <- df %>% group_by(code_id) %>% mutate(across(all_of(protein_list),list(mean=~mean(.x,na.rm=TRUE),w=~.x-mean(.x,na.rm=TRUE)),.names="{.col}_{.fn}")) %>% ungroup()
dat <- dat %>% group_by(code_id) %>% mutate(failureall=max(failureall,na.rm=TRUE)) %>% ungroup()

dat <- dat %>% group_by(code_id) %>% filter(n() >= 2) %>% ungroup()

dat0 <- dat %>% filter(failureall==0)
dat1 <- dat %>% filter(failureall==1)


q_cut <- 0.1
beta_cut <- 0.10
B <- 1000

n_decedent <- length(unique(dat1$code_id))
survivor_ids <- unique(dat0$code_id)

cat("Number of decedent code_id with repeated measurements:",n_decedent,"\n")
cat("Number of survivor code_id with repeated measurements:",length(survivor_ids),"\n")
cat("Number of decedent observations:",nrow(dat1),"\n")
cat("Number of survivor observations:",nrow(dat0),"\n")


####### 1. Match survivor sampling by visit structure and age strata 

n_age_group <- 4

get_id_visit_age <- function(data) {
  data %>%
    group_by(code_id) %>%
    summarise(
      n_visit=n(),
      age_id=min(ageadj,na.rm=TRUE),
      .groups="drop"
    )
}

id_visit_age_decedent <- get_id_visit_age(dat1)
id_visit_age_survivor <- get_id_visit_age(dat0)

age_breaks <- quantile(id_visit_age_decedent$age_id,probs=seq(0,1,length.out=n_age_group+1),na.rm=TRUE)
age_breaks <- unique(as.numeric(age_breaks))

if (length(age_breaks) < 3) {
  age_breaks <- pretty(id_visit_age_decedent$age_id,n=n_age_group)
}

age_breaks[1] <- -Inf
age_breaks[length(age_breaks)] <- Inf

id_visit_age_decedent <- id_visit_age_decedent %>%
  mutate(age_group=cut(age_id,breaks=age_breaks,include.lowest=TRUE,right=FALSE))

id_visit_age_survivor <- id_visit_age_survivor %>%
  mutate(age_group=cut(age_id,breaks=age_breaks,include.lowest=TRUE,right=FALSE))

visit_age_dist_decedent <- id_visit_age_decedent %>%
  count(n_visit,age_group,name="n_id") %>%
  arrange(n_visit,age_group)

print(visit_age_dist_decedent)

check_survivor_pool <- id_visit_age_survivor %>%
  count(n_visit,age_group,name="n_survivor_pool") %>%
  right_join(visit_age_dist_decedent,by=c("n_visit","age_group")) %>%
  mutate(n_survivor_pool=replace_na(n_survivor_pool,0),
         enough=n_survivor_pool>=n_id)

print(check_survivor_pool)

if (any(!check_survivor_pool$enough)) {
  stop("Some n_visit-age_group strata do not have enough survivor IDs. Reduce n_age_group, for example from 4 to 3.")
}

sample_survivor_ids_matched_visit_age <- function(id_visit_age_survivor,visit_age_dist_decedent) {
  sampled <- list()
  for (k in seq_len(nrow(visit_age_dist_decedent))) {
    v <- visit_age_dist_decedent$n_visit[k]
    a <- visit_age_dist_decedent$age_group[k]
    n_need <- visit_age_dist_decedent$n_id[k]
    pool <- id_visit_age_survivor %>%
      filter(n_visit==v,age_group==a) %>%
      pull(code_id)
    if (length(pool) < n_need) stop(paste0("Not enough survivor IDs with n_visit = ",v," and age_group = ",a))
    sampled[[k]] <- sample(pool,size=n_need,replace=FALSE)
  }
  unlist(sampled)
}

run_pairwise_lmer <- function(data,protein_list) {
  results <- list()
  for (i in protein_list) {
    for (j in protein_list) {
      if (i==j) next
      w_var <- paste0(j,"_w")
      mean_var <- paste0(j,"_mean")
      form <- as.formula(paste0(i," ~ ",w_var," + ",mean_var," + ageadj + sexadj + factor(time) + (1 | code_id)"))
      vars_needed <- c(i,w_var,mean_var,"ageadj","sexadj","time","code_id")
      N <- sum(complete.cases(data[,vars_needed]))
      model <- tryCatch(lmer(form,data=data,REML=FALSE),error=function(e) NULL)
      if (is.null(model)) {
        results[[paste(i,j,sep="_")]] <- data.frame(outcome=i,predictor=j,term=NA,estimate=NA,se=NA,p=NA,N=N)
      } else {
        coef_summary <- summary(model)$coef
        results[[paste(i,j,sep="_")]] <- data.frame(outcome=i,predictor=j,term=rownames(coef_summary),estimate=coef_summary[,"Estimate"],se=coef_summary[,"Std. Error"],p=coef_summary[,"Pr(>|t|)"],N=N)
      }
    }
  }
  bind_rows(results)
}

extract_within <- function(res) {
  res %>% filter(!is.na(term),grepl("_w$",term)) %>% select(outcome,predictor,estimate,se,p,N) %>%
    rename(within_estimate=estimate,within_se=se,within_p=p) %>% mutate(q=p.adjust(within_p,method="fdr"))
}

make_directed_edges <- function(res_within,q_cut=0.05,beta_cut=0.10) {
  res_sig <- res_within %>% filter(q<q_cut,abs(within_estimate)>beta_cut)
  if (nrow(res_sig)==0) {
    return(data.frame(outcome=character(),predictor=character(),within_estimate=numeric(),within_se=numeric(),within_p=numeric(),N=integer(),q=numeric(),pair=character()))
  }
  res_sig %>% mutate(pair=as.character(ifelse(outcome<predictor,paste0(outcome,"__",predictor),paste0(predictor,"__",outcome))))
}

merge_directed_edges <- function(res_sig) {
  if (nrow(res_sig)==0) {
    return(data.frame(pair=character(),source=character(),target=character(),outcome=character(),predictor=character(),n_edges=integer(),beta=numeric(),se=numeric(),min_p=numeric(),min_q=numeric()))
  }
  res_sig %>% group_by(pair) %>%
    summarise(outcome=paste(outcome,collapse=";"),predictor=paste(predictor,collapse=";"),betas=list(within_estimate),ses=list(within_se),ps=list(within_p),qs=list(q),.groups="drop") %>%
    rowwise() %>%
    mutate(n_edges=length(betas),
           beta=if(n_edges==1) betas[[1]] else sum((1/unlist(ses)^2)*unlist(betas))/sum(1/unlist(ses)^2),
           se=if(n_edges==1) ses[[1]] else sqrt(1/sum(1/unlist(ses)^2)),
           min_p=min(unlist(ps),na.rm=TRUE),
           min_q=min(unlist(qs),na.rm=TRUE)) %>%
    ungroup() %>%
    separate(pair,into=c("source","target"),sep="__",remove=FALSE) %>%
    select(pair,source,target,outcome,predictor,n_edges,beta,se,min_p,min_q)
}

calc_edge_type_metrics <- function(edge_df,node_df) {
  if (nrow(edge_df)==0) {
    return(data.frame(n_edges_c_ig=0,n_edges_c_c=0,n_edges_ig_ig=0,prop_edges_c_ig=NA_real_,prop_edges_c_c=NA_real_,prop_edges_ig_ig=NA_real_))
  }
  node_type <- node_df %>% select(node,type3)
  edge_type <- edge_df %>%
    left_join(node_type,by=c("source"="node")) %>% rename(source_type=type3) %>%
    left_join(node_type,by=c("target"="node")) %>% rename(target_type=type3)
  total_edges <- nrow(edge_type)
  n_c_ig <- sum((edge_type$source_type=="c" & edge_type$target_type=="ig") | (edge_type$source_type=="ig" & edge_type$target_type=="c"),na.rm=TRUE)
  n_c_c <- sum(edge_type$source_type=="c" & edge_type$target_type=="c",na.rm=TRUE)
  n_ig_ig <- sum(edge_type$source_type=="ig" & edge_type$target_type=="ig",na.rm=TRUE)
  data.frame(n_edges_c_ig=n_c_ig,n_edges_c_c=n_c_c,n_edges_ig_ig=n_ig_ig,
             prop_edges_c_ig=n_c_ig/total_edges,prop_edges_c_c=n_c_c/total_edges,prop_edges_ig_ig=n_ig_ig/total_edges)
}

calc_network_metrics <- function(edge_df,node_df,group_name=NA,iteration=NA) {
  edge_type_metrics <- calc_edge_type_metrics(edge_df,node_df)
  g <- graph_from_data_frame(d=edge_df %>% select(source,target),vertices=node_df %>% select(node,label,type3),directed=FALSE)
  deg <- degree(g)
  if (ecount(g)==0) {
    net_metrics <- data.frame(group=group_name,iteration=iteration,n_nodes=sum(deg>0),n_edges=0,mean_degree=mean(deg),density=0,
                              global_efficiency=NA_real_,mean_betweenness=NA_real_,max_degree=max(deg))
    return(bind_cols(net_metrics,edge_type_metrics))
  }
  d <- distances(g)
  inv_d <- 1/d
  inv_d[is.infinite(inv_d)] <- 0
  diag(inv_d) <- 0
  n <- vcount(g)
  net_metrics <- data.frame(group=group_name,iteration=iteration,n_nodes=sum(deg>0),n_edges=ecount(g),mean_degree=mean(deg),density=edge_density(g,loops=FALSE),
                            global_efficiency=sum(inv_d)/(n*(n-1)),mean_betweenness=mean(betweenness(g,directed=FALSE,normalized=TRUE)),max_degree=max(deg))
  bind_cols(net_metrics,edge_type_metrics)
}

empirical_compare <- function(metric_name,dec_value,surv_df) {
  x <- surv_df[[metric_name]]
  x <- x[!is.na(x)]
  n_valid <- length(x)
  if (n_valid==0 | is.na(dec_value)) {
    return(data.frame(metric=metric_name,decedent_value=dec_value,
                      survivor_median=NA_real_,survivor_mean=NA_real_,survivor_sd=NA_real_,
                      survivor_q25=NA_real_,survivor_q75=NA_real_,
                      survivor_p2.5=NA_real_,survivor_p97.5=NA_real_,
                      empirical_p_lower=NA_real_,empirical_p_upper=NA_real_,
                      z_vs_survivor=NA_real_,n_valid=n_valid))
  }
  p_lower <- (sum(x<=dec_value)+1)/(n_valid+1)
  p_upper <- (sum(x>=dec_value)+1)/(n_valid+1)
  z_value <- ifelse(sd(x,na.rm=TRUE)==0,NA_real_,(dec_value-mean(x,na.rm=TRUE))/sd(x,na.rm=TRUE))
  data.frame(metric=metric_name,decedent_value=dec_value,
             survivor_median=median(x,na.rm=TRUE),survivor_mean=mean(x,na.rm=TRUE),survivor_sd=sd(x,na.rm=TRUE),
             survivor_q25=quantile(x,0.25,na.rm=TRUE),survivor_q75=quantile(x,0.75,na.rm=TRUE),
             survivor_p2.5=quantile(x,0.025,na.rm=TRUE),survivor_p97.5=quantile(x,0.975,na.rm=TRUE),
             empirical_p_lower=p_lower,
             empirical_p_upper=p_upper,
             z_vs_survivor=z_value,
             n_valid=n_valid)
}

node <- read.csv("node_table.csv",header=TRUE)
node$node <- trimws(node$node)
node$type3 <- trimws(tolower(node$type3))

res1 <- run_pairwise_lmer(dat1,protein_list)
res1_within <- extract_within(res1)
res1_sig <- make_directed_edges(res1_within,q_cut=q_cut,beta_cut=beta_cut)
edge1 <- merge_directed_edges(res1_sig)
metrics1 <- calc_network_metrics(edge1,node,group_name="decedent_full",iteration=1)

surv_edge_list <- vector("list",B)
surv_metric_list <- vector("list",B)
surv_directed_list <- vector("list",B)

for (b in 1:B) {
  sampled_ids <- sample_survivor_ids_matched_visit_age(id_visit_age_survivor,visit_age_dist_decedent)
  boot_data <- dat0 %>% filter(code_id %in% sampled_ids)
  boot_res <- run_pairwise_lmer(boot_data,protein_list)
  boot_within <- extract_within(boot_res) %>%
    mutate(iteration=b,n_sampled_code_id=length(unique(boot_data$code_id)),n_obs=nrow(boot_data))
  boot_sig <- make_directed_edges(boot_within,q_cut=q_cut,beta_cut=beta_cut) %>% mutate(iteration=b)
  boot_edge <- merge_directed_edges(boot_sig) %>% mutate(iteration=b)
  boot_metric <- calc_network_metrics(boot_edge,node,group_name="survivor_size_obs_matched",iteration=b) %>%
    mutate(n_sampled_code_id=length(unique(boot_data$code_id)),n_obs=nrow(boot_data))
  surv_directed_list[[b]] <- boot_sig
  surv_edge_list[[b]] <- boot_edge
  surv_metric_list[[b]] <- boot_metric
  cat("Finished survivor network",b,"of",B,"\n")
}

surv_directed_all <- bind_rows(surv_directed_list)
surv_edges_all <- bind_rows(surv_edge_list)
surv_metrics <- bind_rows(surv_metric_list)

edge0_freq <- surv_edges_all %>% group_by(pair,source,target) %>%
  summarise(freq=n()/B,mean_beta=mean(beta,na.rm=TRUE),median_beta=median(beta,na.rm=TRUE),mean_n_edges=mean(n_edges,na.rm=TRUE),.groups="drop") %>%
  arrange(desc(freq))
edge0_rep <- edge0_freq %>% filter(freq>=0.50) %>% mutate(beta=mean_beta,n_edges=round(mean_n_edges))
metrics0_rep <- calc_network_metrics(edge0_rep,node,group_name="survivor_representative_freq50",iteration=NA)

metric_names <- c("n_nodes","n_edges","mean_degree","density","global_efficiency","mean_betweenness","max_degree",
                  "n_edges_c_ig","n_edges_c_c","n_edges_ig_ig","prop_edges_c_ig","prop_edges_c_c","prop_edges_ig_ig")

metric_compare <- bind_rows(lapply(metric_names,function(m) empirical_compare(m,metrics1[[m]][1],surv_metrics)))

metric_compare <- metric_compare %>%
  mutate(empirical_p_two_sided=pmin(1,2*pmin(empirical_p_lower,empirical_p_upper)),
         direction=case_when(
           decedent_value > survivor_median ~ "higher_in_decedents",
           decedent_value < survivor_median ~ "lower_in_decedents",
           TRUE ~ "similar"
         ),
         empirical_p_observed_tail=case_when(
           decedent_value > survivor_median ~ empirical_p_upper,
           decedent_value < survivor_median ~ empirical_p_lower,
           TRUE ~ 1
         ),
         empirical_p_directional=empirical_p_observed_tail)


############# slide_window analysis with longitudinal modelling ##################
library(dplyr)
library(purrr)
library(lme4)
library(segmented)

df_gnhs<-read.csv("gnhs.csv",header=TRUE)
df_gnhs$time <- factor(df_gnhs$time, levels = c(1, 2, 3))
df_gnhs<-df_gnhs %>% mutate(death = ifelse(is.na(death), 0, death))

fit_slide <- function(center_age, data, varname, window_width) {
  df_sub <- data %>% filter(ageadj >= center_age - window_width / 2,ageadj <  center_age + window_width / 2)
  n_id  <- n_distinct(df_sub$id)
  n_obs <- nrow(df_sub)
  if (n_id < 50 | n_obs < 50) return(NULL)
  use_mixed <- n_obs > n_id
  if (use_mixed) {
    fit <- lmer(as.formula(paste0(varname, " ~ 1 + (1 | id)")), data = df_sub)
    est <- fixef(fit)[1]
    se  <- sqrt(vcov(fit)[1, 1])
  } else {
    fit <- lm(as.formula(paste0(varname, " ~ 1")), data = df_sub)
    est <- coef(fit)[1]
    se  <- summary(fit)$coefficients[1, 2]}
  tibble(age_center = center_age, mean_score = est, lower= est - 1.96 * se, upper = est + 1.96 * se,score = varname)}

window_width <- 10
step_size    <- 1
age_seq <- seq(45, 80, by = step_size)
res_imba <- map_dfr(age_seq, fit_slide,data = df_gnhs,varname = "stdscore_imba", window_width = window_width)
res_acti <- map_dfr(age_seq, fit_slide, data = df_gnhs, varname = "stdscore_acti",window_width = window_width)
res_regu <- map_dfr(age_seq, fit_slide,data = df_gnhs,varname = "stdscore_regu",window_width = window_width)

fit_segmented <- function(df, score_name){
  lm_fit <- lm(as.formula(paste0(score_name, " ~ age_center")), data = df)
  seg_fit <- segmented(lm_fit, seg.Z = ~ age_center, npsi = 1)
  s <- summary(seg_fit)
  bp <- s$psi[1, "Est."]
  beta1 <- coef(seg_fit)["age_center"]
  beta_change <- s$coefficients["U1.age_center", "Estimate"]
  beta2 <- beta1 + beta_change
  se   <- s$coefficients["U1.age_center", "Std. Error"]
  z <- beta_change / se
  log_p <- log(2) + pnorm(abs(z), lower.tail = FALSE, log.p = TRUE)
  p <- exp(log_p)
  data.frame(score = score_name,breakpoint = bp,slope_pre = beta1,slope_post = beta2,slope_change = beta_change,se = se,z = z,p_value = p)
}

res_imba_seg <- fit_segmented(res_imba, "mean_score")
res_acti_seg <- fit_segmented(res_acti, "mean_score")
res_regu_seg <- fit_segmented(res_regu, "mean_score")
res_all_seg <- bind_rows(res_imba_seg, res_acti_seg,res_regu_seg)
res_all_seg


####################### age_specific HR ############################
library(survival)
library(splines)
library(ggplot2)

df_gnhs<-read.csv("gnhs.csv",header=TRUE)
df_gnhs$time <- factor(df_gnhs$time, levels = c(1, 2, 3))
df_gnhs<-df_gnhs %>% mutate(death = ifelse(is.na(death), 0, death))
df_gnhs$age_entry <- df_gnhs$age_baseline + df_gnhs$entry
df_gnhs$age_exit  <- df_gnhs$age_baseline + df_gnhs$exit

get_hr_curve <- function(varname, data) {
  fit <- coxph(as.formula(paste0("Surv(age_entry, age_exit, failureall) ~ ",varname, " + tt(", varname, ") + sexadj")),
               data = data,tt = function(x, t, ...) {x * ns(t, df = 3)})
  age_seq <- seq(45, 80, by = 0.5)
  B <- ns(age_seq, df = 3)
  coef_all <- coef(fit)
  idx <- grep(varname, names(coef_all))
  beta <- coef_all[idx]
  X <- cbind(1, B)
  logHR <- X %*% beta
  HR <- exp(logHR)
  vc <- vcov(fit)[idx, idx]
  se <- sqrt(diag(X %*% vc %*% t(X)))
  data.frame(age = age_seq, HR = as.vector(HR),lower = exp(logHR - 1.96 * se),upper = exp(logHR + 1.96 * se),score = varname)
}

df_imba <- get_hr_curve("stdscore_imba", df_gnhs)
df_acti <- get_hr_curve("stdscore_acti", df_gnhs)
df_regu <- get_hr_curve("stdscore_regu", df_gnhs)

get_hr_extrema <- function(df){ max_idx <- which.max(df$HR)
min_idx <- which.min(df$HR)
data.frame(age_max = df$age[max_idx], HR_max  = df$HR[max_idx], age_min = df$age[min_idx],HR_min  = df$HR[min_idx])}

ext_imba <- get_hr_extrema(df_imba)
ext_acti <- get_hr_extrema(df_acti)
ext_regu <- get_hr_extrema(df_regu)
ext_all <- rbind(cbind(score="imba", ext_imba), cbind(score="acti", ext_acti),cbind(score="regu", ext_regu))
ext_all


####################### landmark AUC ############################

library(dplyr)
library(tidyr)
library(purrr)
library(ggplot2)
library(pROC)
library(timeROC)

df_long <- df_gnhs
markers   <- c("stdscore_imba", "stdscore_acti", "stdscore_regu")
landmarks <- seq(2, 14, by = 2)
horizon   <- 2

#landmark_auc_df <- expand.grid(marker = markers,landmark = landmarks, stringsAsFactors = FALSE) %>% rowwise() %>% 
mutate(datL = list(make_landmark_data_interval(df=df_long,tL=landmark,horizon=horizon,marker=marker)),n=nrow(datL), events=sum(datL$event_L), 
       AUC = if (length(unique(datL$event_L)) > 1) {roc_obj <- pROC::roc(response = datL$event_L, predictor = datL$score,quiet = TRUE)
       as.numeric(pROC::auc(roc_obj))} else {NA_real_}) %>%ungroup() %>% select(marker, landmark, horizon, n, events, AUC)

landmark_auc_df <- expand.grid(marker = markers,landmark = landmarks, stringsAsFactors = FALSE) %>% rowwise() %>% 
  mutate(datL = list(make_landmark_data_interval(df=df_long,tL=landmark,horizon=horizon,marker=marker)),n=nrow(datL), events=sum(datL$event_L),
         AUC = if (length(unique(datL$event_L)) > 1) {roc_obj <-pROC::roc(response = datL$event_L, predictor = datL$score,quiet = TRUE)
         as.numeric(pROC::auc(roc_obj))} else {NA_real_}) %>% ungroup()

compare_auc <- function(df, landmark_time){df_sub <- df %>% filter(landmark == landmark_time)
if(nrow(df_sub) < 3) return(NULL)
get_roc <- function(marker_name){ row_idx <- which(df_sub$marker == marker_name)
if(length(row_idx)==0) return(NULL)
dat <- df_sub$datL[[row_idx]]
if(length(unique(dat$event_L)) < 2) return(NULL)
pROC::roc(dat$event_L, dat$score, quiet=TRUE) }
roc_imba <- get_roc("stdscore_imba")
roc_acti <- get_roc("stdscore_acti")
roc_regu <- get_roc("stdscore_regu")
tibble(landmark = landmark_time,
       p_imba_vs_acti = if(!is.null(roc_imba)&!is.null(roc_acti)) pROC::roc.test(roc_imba, roc_acti, method="delong")$p.value else NA,
       p_imba_vs_regu = if(!is.null(roc_imba)&!is.null(roc_regu)) pROC::roc.test(roc_imba, roc_regu, method="delong")$p.value else NA,
       p_acti_vs_regu = if(!is.null(roc_acti)&!is.null(roc_regu)) pROC::roc.test(roc_acti, roc_regu, method="delong")$p.value else NA
)}
auc_compare_res <- map_dfr(landmarks, ~compare_auc(landmark_auc_df, .x))


landmark_age_map <- tibble(landmark = landmarks) %>% rowwise() %>%
  mutate(risk_set = list(df_long %>% filter(entry < landmark, exit > landmark)),
         age_mean = mean(risk_set$age_baseline + landmark, na.rm = TRUE),
         age_sd   = sd(risk_set$age_baseline + landmark, na.rm = TRUE),
         age_med  = median(risk_set$age_baseline + landmark, na.rm = TRUE),
         age_q1   = quantile(risk_set$age_baseline + landmark, 0.25, na.rm = TRUE),
         age_q3   = quantile(risk_set$age_baseline + landmark, 0.75, na.rm = TRUE),
         n = nrow(risk_set)) %>% ungroup() %>% select(landmark, n, age_mean, age_sd, age_med, age_q1, age_q3)
landmark_auc_with_age <- landmark_auc_df %>%left_join(landmark_age_map, by = "landmark")

df_plot <- landmark_auc_with_age %>%mutate(x_start = landmark,x_end = landmark + horizon)
df_vline <- df_plot %>% select(marker,landmark,AUC) %>% mutate(landmark_right=landmark+horizon) %>% rename(landmark_left=landmark, AUC_left = AUC) %>%
  left_join(df_plot %>% select(marker, landmark, AUC) %>% rename(landmark_right = landmark, AUC_right = AUC),by = c("marker", "landmark_right")) %>%
  filter(!is.na(AUC_right)) %>% mutate(x = landmark_right)


####################### td AUC ############################
library(dplyr)
library(timeROC)
library(tidyr)

df_long <- df_long %>% mutate(time = as.numeric(as.character(time)))
df_base <- df_long %>% group_by(id) %>% arrange(time, .by_group = TRUE) %>% summarise(
  entry = min(entry, na.rm = TRUE),
  exit  = max(exit, na.rm = TRUE),
  failureall = max(failureall, na.rm = TRUE),
  stdscore_imba = stdscore_imba[time == min(time)][1],
  stdscore_acti = stdscore_acti[time == min(time)][1],
  stdscore_regu = stdscore_regu[time == min(time)][1]) %>% ungroup() %>%
  mutate(surtime = exit - entry )

df_base <- df_base %>%  mutate(stdscore_regu = -stdscore_regu)

markers <- c("stdscore_imba", "stdscore_acti", "stdscore_regu")
auc_list <- lapply(markers, function(m) {
  roc <- timeROC(
    T = df_base$surtime,
    delta = df_base$failureall,
    marker = df_base[[m]],
    cause = 1,
    times = 2:10,
    iid = TRUE)
  data.frame(marker = m,time = 2:10,AUC = as.numeric(roc$AUC))})

auc_df <- bind_rows(auc_list)

