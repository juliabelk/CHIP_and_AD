suppressMessages({
    library(dplyr)
    library(ggplot2)
    library(RColorBrewer)
    library(ggrepel)
})

brain_vaf <- list(
    "ACT2_ns_DNMT"=0.0036,
    "ACT3_ns_DNMT"=0.035,
    "ACT4_ns_SF3B"=0.0042,
    "ACT5_ns_TET2"=0.0069,
    "ACT6_ns_TET2_oc"=0.0052,
    "ACT6_ns_TET2_ce"=0.001567036,
    "ACT6_ns_TET2_p"=0.0166,
    "ACT6_s_TET2_p"=0.0243,
    "ACT7_ns_ASXL_ASXL"=0.0092,
    "ACT7_ns_ASXL_TET2"=0.0052
)
blood_vaf <- list(
    "ACT2_ns_DNMT"=0.14,
    "ACT4_ns_SF3B"=0.18,
    "ACT7_ns_ASXL_ASXL"=0.327,
    "ACT3_ns_DNMT"=0.385,
    "ACT5_ns_TET2"=0.111,
    "ACT6_ns_TET2_oc"=0.143
)
pct_mut_mg <- list()

counts_df <- read.delim("cell_count_by_cluster_and_sample.csv",sep=",")
counts_df$total  <- rowSums(counts_df)
counts_df$non_mg <- rowSums(counts_df[,c("C7.Monocyte","C8.cDC","C9.T.cell","C10.B.cell")])
counts_df$mg_pct <- counts_df$C6.Microglia / counts_df$total

print("35")
print(counts_df)

calculate_adj_vaf_mut_mg <- function(df, row_id, ky=NA) { # vaf, 
    
    if (is.na(ky)) ky <- row_id
    
    vaf <- brain_vaf[[ky]]
    
    mut_mg <- round( df[row_id,"total"] * vaf*2) - df[row_id,"non_mg"]
    #VAF estimate assuming non-microglia mutant cells are 100% VAF
    vaf_adj <- mut_mg / df[row_id,"total"] / 2
    set.seed(1234)
    x = rbinom(n=10^6, size=200000, prob=vaf_adj)
    y = rbinom(n=10^6, size=df[row_id,"total"],  prob=df[row_id,"mg_pct"])
    z = 2*(x/200000)/(y/df[row_id,"total"])
    #hist(z, prob=T, br=100, col="skyblue2", main=NULL, xlab="Proportion")
    
    q <- quantile(z, probs = c(.025, .5, .975))
    print(row_id)
    print(sprintf("Adjusted VAF: %s ---- Fraction mutant mg: %s (%s - %s)",
            round(vaf_adj,4),round(q[2],4),round(q[1],4),round(q[3],4)))
    return(q)
}

pct_mut_mg[["ACT6_ns_TET2_ce"]]   <- calculate_adj_vaf_mut_mg(counts_df,"ACT6_ns_TET2_ce")
pct_mut_mg[["ACT6_ns_TET2_p"]]    <- calculate_adj_vaf_mut_mg(counts_df,"ACT6_ns_TET2_p")
pct_mut_mg[["ACT6_s_TET2_p"]]     <- calculate_adj_vaf_mut_mg(counts_df,"ACT6_s_TET2_p")
pct_mut_mg[["ACT2_ns_DNMT"]]      <- calculate_adj_vaf_mut_mg(counts_df,"ACT2_ns_DNMT")
pct_mut_mg[["ACT4_ns_SF3B"]]      <- calculate_adj_vaf_mut_mg(counts_df,"ACT4_ns_SF3B")
pct_mut_mg[["ACT7_ns_ASXL_ASXL"]] <- calculate_adj_vaf_mut_mg(counts_df,"ACT7_ns_ASXL","ACT7_ns_ASXL_ASXL")
pct_mut_mg[["ACT7_ns_ASXL_TET2"]] <- calculate_adj_vaf_mut_mg(counts_df,"ACT7_ns_ASXL","ACT7_ns_ASXL_TET2")
pct_mut_mg[["ACT3_ns_DNMT"]]      <- calculate_adj_vaf_mut_mg(counts_df,"ACT3_ns_DNMT")
pct_mut_mg[["ACT5_ns_TET2"]]      <- calculate_adj_vaf_mut_mg(counts_df,"ACT5_ns_TET2")
pct_mut_mg[["ACT6_ns_TET2_oc"]]   <- calculate_adj_vaf_mut_mg(counts_df,"ACT6_ns_TET2_oc")


df <- as.data.frame(matrix(nrow=0,ncol=5))
for (i in names(pct_mut_mg)) {
    #if (i %in% names(blood_vaf))
    df <- rbind(df,c(i,pct_mut_mg[[i]],brain_vaf[[i]]))
}
colnames(df) <- c("id","mg_q1","mg_q2","mg_q3","brain_vaf")#,"brain_vaf")
df$blood_vaf <- unlist(lapply(df$id,function(id) {
    if (id %in% names(blood_vaf)) blood_vaf[[id]]
    else NA
}))
df$donor <- unlist(lapply(df$id,function(id) {
    tmp <- strsplit(id,"_")[[1]]
    paste0(tmp[3],"_",tmp[1])
} ))
for (i in 2:5) df[[i]] <- as.numeric(df[[i]])

plt1 <- ggplot(df,aes(x=id,y=mg_q2,fill=donor)) +
    geom_bar(stat="identity",color="black") +
    geom_errorbar(aes(ymin=mg_q1,ymax=mg_q3),width=0.2)+
    theme_classic() + 
    #scale_fill_manual(values=pal) +
    theme(text = element_text(size=24,face="bold")) +
    guides(colour = guide_legend(override.aes = list(size=10))) +
    theme(axis.ticks = element_line(colour = "black", size = 1), axis.ticks.length=unit(5, "pt")) +
    theme(panel.grid.major = element_blank()) +
    theme(aspect.ratio=0.5) +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
    ggtitle("Percent mutant microglia")

#plt
ggsave("pct_mut_mg.eps",plt1,width=12,height=8)

blood_df <- df[which(!is.na(df$blood_vaf)),]

#blood_df

rval <- cor.test(blood_df$blood_vaf*2, blood_df$mg_q2)#, method=c("pearson", "kendall", "spearman"))
rval

plt2 <- ggplot(blood_df,aes(x=blood_vaf*2,y=mg_q2,fill=donor,label=id)) +
    geom_point(size=5,shape=21,color="black") +
    #geom_bar(stat="identity",color="black") +
    #geom_errorbar(aes(ymin=mg_q1,ymax=mg_q3),width=0.2)+
    theme_bw() + 
    #scale_fill_manual(values=pal) +
    theme(text = element_text(size=24,face="bold")) +
    guides(colour = guide_legend(override.aes = list(size=10))) +
    theme(axis.ticks = element_line(colour = "black", size = 1), axis.ticks.length=unit(5, "pt")) +
    theme(panel.grid.minor = element_blank()) +
    scale_x_continuous(limits=c(0,1),breaks=seq(0,1,0.2)) +
    scale_y_continuous(limits=c(0,1),breaks=seq(0,1,0.2)) +
    geom_text_repel(size=6,box.padding=1.5,max.overlaps=Inf) +
    ggtitle(paste0("R = ",rval$estimate,"\np = ",rval$p.value)) +
    theme(aspect.ratio=1)


#plt
ggsave("blood_mg_pct_mut_corr.eps",plt2,width=12,height=8)





