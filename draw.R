require(ggplot2); require(scales); 
require(data.table)
require(reshape2)

d=(read.csv("./CSV_Files/res_NewErrRates.csv", sep=",", header=F))
names(d) <- c("E", "DR", "X", "Diameter", "PD", "N", "ErrLen", "NumErrSeqDiv", "Rep", "FP0", "FN0", "TP0", "TN0", "FP", "FN", "TP", "TN")

d = d[d$DR!="concatenation",]
d = d[d$N > 19 | !grepl("16S",d$E),]
nlabels = c("1","2%","5%","10%","20%")
#levels(d$DR)[levels(d$DR)=="concatenation"] <- "concat"

d$n=with(d,as.factor(round(100/((NumErrSeqDiv==0&grepl("ErrLen$",E))*20+(NumErrSeqDiv!=N|grepl("ErrLen$",E))*NumErrSeqDiv+(NumErrSeqDiv==N&!grepl("ErrLen$",E))*100))))
d$nb=d$n
levels(d$n) <- c(levels(d$n)[1],paste(levels(d$n)[0:-1],"%",sep=""),"~5%")
d[grepl("General$",d$E),"n"]="~5%"

d$ErrLen = (d$ErrLen==0)*8+d$ErrLen
d$ErrLenT = paste(d$ErrLen, intToUtf8(215), ifelse(grepl("small",d$E),"11","11"),sep="")
#d[grepl("small",d$E),]$ErrLenT = paste(d[grepl("small",d$E),]$ErrLen, intToUtf8(215), "4",sep="")
d[grepl("General$",d$E),"ErrLenT"]="~50"
d$ErrLenT = factor(d$ErrLenT,levels=c("2×11","2×7","3×11", "3×7","4×11", "4×7","8×11", "8×7","16×11","16×7", "32×11", "32×7","64×11","~50" ))

d$SL = with(d,(TN+FP+FN+TP)/as.numeric(as.character(N)))


# For Hackett, make sure sure error length is at most half of the sequence length
d = d[d$ErrLen!=64 | !grepl("Hack",d$E),]
d = d[!d$ErrLen==32 | !grepl("Hack",d$E) | d$SL>=704,]


# Aggregate Sum function
summ_roc <- function(d2,form) {
  ad2 = dcast(d2, form ,fun.aggregate=sum,value.var = c("FP"))
  ad2=cbind(dcast(d2, form ,fun.aggregate=sum,value.var = c("FP")),
            dcast(d2, form ,fun.aggregate=sum,value.var = c("FN"))[,length(ad2)],
            dcast(d2, form ,fun.aggregate=sum,value.var = c("TP"))[,length(ad2)],
            dcast(d2, form ,fun.aggregate=sum,value.var = c("TN"))[,length(ad2)],
            dcast(d2, form ,fun.aggregate=sum,value.var = c("FN0"))[,length(ad2)],
            dcast(d2, form ,fun.aggregate=sum,value.var = c("TP0"))[,length(ad2)],
            dcast(d2, form ,fun.aggregate=sum,value.var = c("TN0"))[,length(ad2)],
            dcast(d2, form ,fun.aggregate=sum,value.var = c("FP0"))[,length(ad2)]
  )
  names(ad2)[(length(names(ad2))-7):(length(names(ad2)))]=c("FP","FN","TP","TN", "FN0", "TP0", "TN0", "FP0")
  ad2
}

# 16S.B: K - Recall vs Diameter
ggplot(aes(x=Diameter,y=TP/(TP+FN), 
           color=as.factor(NumErrSeqDiv)),data=d[d$E=="16S.B_K",])+geom_point(alpha=0.5)+
  theme_classic()+geom_smooth()+scale_y_continuous("Recall")+
  scale_shape(name="")+scale_color_brewer(palette = "Paired",name="k")
ggsave("Figures/Union_Figures/Recall-k.pdf",width = 6,height = 6)

options(digits = 2)
d2=summ_roc(d[d$E=="16S.B_K" & d$N > 19 & d$ErrLen!=0,], NumErrSeqDiv+cut(Diameter, breaks = c(0, 0.1, 0.2, 0.5, 0.8, 1), right = F)~.)
A = data.frame(x=d2$FP/(d2$FP+d2$TN),y=d2$TP/(d2$TP+d2$FN), NumErrSeqDiv=d2$NumErrSeqDiv, DR=d2$`cut(Diameter, breaks = c(0, 0.1, 0.2, 0.5, 0.8, 1), right = F)`)
B = data.frame(x=d2$FP0/(d2$FP0+d2$TN0),y=as.vector(matrix(1.1,nrow=nrow(d2))), NumErrSeqDiv=d2$NumErrSeqDiv, DR=d2$`cut(Diameter, breaks = c(0, 0.1, 0.2, 0.5, 0.8, 1), right = F)`)
ggplot(data=A, aes(x, y, group=as.factor(NumErrSeqDiv),color=as.factor(NumErrSeqDiv), shape=as.factor(DR))) + geom_point(alpha=1)+
  theme_light()+theme(legend.position = "right")+geom_line()+
  scale_shape(name="Diameter")+scale_color_brewer(name=expression(k),palette = "Paired")+
  scale_x_continuous(name="FPR",labels=percent)+
  scale_y_continuous("Recall",labels=percent,breaks = c(0.2,0.4,0.6,0.8,1))+
  ggtitle("16S with varying k, fixed error length 88 and 5% error frequency")
ggsave("Figures/Union_Figures/16S_ks_fixed_length_ROC.pdf", width=6, height=6)


# General Results: Recall vs Diameter
# Uses normal distribution to draw the error lengths and number of erroneuous sequences
ggplot(aes(x=Diameter, y=TP/(TP+FN)), data=d[d$E=="16S.B_General" & d$N > 10 & d$ErrLen==8,]) + theme_classic() + geom_point() + geom_smooth() + scale_y_continuous("Recall") + ggtitle("16S.B: Recall vs Diameter")
ggsave("Figures/General_Figures/16SB_General_Recall.pdf", width=6, height=6)
ggplot(aes(x=Diameter, y=TP/(TP+FN)), data=d[d$E=="Hackett_General" & d$N > 10 & d$ErrLen==8,]) + theme_classic() + geom_point() + geom_smooth() + scale_y_continuous("Recall") + ggtitle("Hackett: Recall vs Diameter")
ggsave("Figures/General_Figures/Hackett_General_Recall.pdf", width=6, height=6)
ggplot(aes(x=Diameter, y=TP/(TP+FN)), data=d[d$E=="small-10-aa_Res" & d$N > 10 & d$ErrLen==8,]) + theme_classic() + geom_point() + geom_smooth() + scale_y_continuous("Recall") + ggtitle("small-10-aa: Recall vs Diameter")
ggsave("Figures/General_Figures/small10aa_General_Recall.pdf", width=6, height=6)






# # ROC all 
# options(digits = 2)
# d$E2=sub("_.*","",d$E)
# d2=summ_roc(d[grepl("Err",d$E) & (d$N > 19|d$N==0)  ,], E2+n+ErrLenT+cut(Diameter, breaks = c(0, 0.1, 0.2, 0.5, 0.8, 1), right = F)~.)
# A = data.frame(x=d2$FP/(d2$FP+d2$TN),y=d2$TP/(d2$TP+d2$FN),E =d2$E2, ErrLenT=d2$ErrLenT,  n=d2$n, DR=d2$`cut(Diameter, breaks = c(0, 0.1, 0.2, 0.5, 0.8, 1), right = F)`)
# B = data.frame(x=d2$FP0/(d2$FP0+d2$TN0),y=as.vector(matrix(1,nrow=nrow(d2))), E =d2$E2, n=d2$n, ErrLenT=d2$ErrLenT, DR=d2$`cut(Diameter, breaks = c(0, 0.1, 0.2, 0.5, 0.8, 1), right = F)`)
# ggplot(data=A, aes(x, y, shape=interaction(n,ErrLenT,sep=", "),color=as.factor(DR))) + 
#   geom_point(alpha=1)+facet_wrap(~sub("_.*","",E),ncol=1,scales="free_x")+
#   geom_path(aes(group=interaction(DR,n)),linetype=2)+geom_path(aes(group=interaction(DR,ErrLenT)),linetype=1)+
#   theme_bw()+theme(legend.position = "right",legend.text.align = 1)+
#   scale_shape_manual(name="Err Freq%, Len",values=c(1,2,5,15,17,8,19,18,3,4,6,10,20:40))+
#   scale_color_brewer(name="Diameter",palette = "Dark2")+
#   scale_x_continuous(name="FPR",labels=percent)+
#   scale_y_continuous("Recall",labels=percent)+
#   geom_linerange(aes(x=x,ymin=0.995,ymax=1.005,color=as.factor(DR)),data=B,linetype=1,size=1)
# ggsave("Figures/ErrParam_Figures/All_both_ROC.pdf", width=6, height=10)



# 16S.B - Varying both


d$DR2 = cut(d$Diameter, breaks = c(0, 0.1, 0.2, 0.5, 0.8, 1), right = F)
# ROC for 16S.B with varying error lengths and fixed percentage of erroneous sequences
options(digits = 2)
d2=summ_roc(d[d$E %in% c( "16S.B_ErrLen","16S.B_NumErrAlns", "16S.B_General") & d$N > 19 ,], n+ErrLenT+DR2+E~.)
A = data.frame(x=d2$FP/(d2$FP+d2$TN),y=d2$TP/(d2$TP+d2$FN), ErrLenT=d2$ErrLenT,  n=d2$n, DR=d2$DR2,E=d2$E)
B = data.frame(x=d2$FP0/(d2$FP0+d2$TN0),y=as.vector(matrix(1,nrow=nrow(d2))), n=d2$n, ErrLenT=d2$ErrLenT, DR=d2$DR2)
ggplot(data=A, aes(x, y, shape=interaction(ErrLenT,n,sep=", "),color=as.factor(DR))) + 
  geom_point(alpha=1)+
  geom_path(aes(group=interaction(DR,n)),data=A[A$n!="~5%",],linetype=1)+
  geom_path(aes(group=interaction(DR,ErrLenT)),data=A[A$n!="~5%",],linetype=2)+
  theme_bw()+theme(legend.position = "right",legend.text.align = 1)+
  scale_shape_manual(name="Length, Freq",values=c(15,17,1,2,5,8,9,7,6,19,18,3,90))+
  scale_color_brewer(name="Diameter",palette = "Dark2")+
  scale_x_continuous(name="FPR",labels=percent)+
  scale_y_continuous("Recall",labels=percent)+
  geom_linerange(aes(x=x,ymin=0.995,ymax=1.005,color=as.factor(DR)),data=B,linetype=1,size=1)
ggsave("Figures/ErrParam_Figures/16SB_ErrLenNumErr_ROC.pdf", width=6.6, height=5)

ggplot(data=A, aes(x, y, shape=interaction(ErrLenT,n,sep=", "),color=as.factor(DR))) + 
  geom_point(alpha=1)+
  geom_path(aes(group=DR),data=A[A$n!="~5%",],linetype=1)+
  theme_classic()+theme(legend.position = c(0.7,0.2),legend.text.align = 1,legend.direction = "horizontal")+
  scale_shape_manual(name="Err Len, Freq",values=c(15,17,1,2,5,8,9,7,6,19,18,3))+
  scale_color_brewer(name="Diameter",palette = "Dark2")+
  scale_x_continuous(name="FPR",labels=percent)+
  scale_y_continuous("Recall",labels=percent)+
  facet_wrap(~(E=="16S.B_NumErrAlns"),labeller = function(x) list(c("Changing Error Length", "Changing Error Frequency")))+
  geom_linerange(aes(x=x,ymin=0.995,ymax=1.005,color=as.factor(DR)),data=B,linetype=1,size=1)
ggsave("Figures/ErrParam_Figures/16SB_ErrLenNumErr_ROC_faceted.pdf", width=10, height=5)



d2=summ_roc(d[d$E %in% c( "16S.B_ErrLen") & d$N > 19 ,], ErrLen+cut(Diameter, breaks = c(0, 0.1, 0.2, 0.5, 0.8, 1), right = F)~.)
A = data.frame(x=d2$FP/(d2$FP+d2$TN),y=d2$TP/(d2$TP+d2$FN), ErrLen=d2$ErrLen, DR=d2$`cut(Diameter, breaks = c(0, 0.1, 0.2, 0.5, 0.8, 1), right = F)`)
B = data.frame(x=d2$FP0/(d2$FP0+d2$TN0),y=as.vector(matrix(1,nrow=nrow(d2))), ErrLen=d2$ErrLen, DR=d2$`cut(Diameter, breaks = c(0, 0.1, 0.2, 0.5, 0.8, 1), right = F)`)
ggplot(data=A, aes(x, y, shape=as.factor(ErrLen), group=as.factor(DR),color=as.factor(DR))) +
  geom_point(alpha=1)+geom_line()+
  theme_bw()+theme(legend.position = "right")+
  geom_linerange(aes(x=x,ymin=0.995,ymax=1.005,color=as.factor(DR)),data=B,linetype=1,size=1)+
  scale_shape_manual(name="Error Length",values=c(1,2,5,6,15,17,19),labels = function(x) (paste(x, intToUtf8(215), "11")))+
  scale_color_brewer(name="Diameter",palette = "Dark2")+
  scale_x_continuous(name="FPR",labels=percent)+
  scale_y_continuous("Recall",labels=percent)
ggsave("Figures/ErrParam_Figures/16SB_ErrLen_ROC.pdf", width=6, height=4.5)

# ROC for 16S.B with varying percentages of erroneous sequences and fixed error lengths
d2=d[d$E=="16S.B_NumErrAlns",]
d2$n=with(d2,as.factor(100/((NumErrSeqDiv!=N)*NumErrSeqDiv+(NumErrSeqDiv==N)*100)))
ggplot(aes(x=FP/(FP+TN),y=TP/(TP+FN), color=as.factor(n) ),data=summ_roc(d2,n~.))+
  geom_point(alpha=1)+
  theme_light()+theme(legend.position = c(.85,.25))+
  scale_color_brewer(name="n",labels=nlabels, palette="Paired")+
  scale_x_continuous(name="FPR",labels=percent)+scale_y_continuous("Recall")+
  ggtitle("16S.B with Varying Number of Erroneous Sequences: ROC")
ggsave("Figures/ErrParam_Figures/16SB_NumErrAlns_ROC.pdf", width=6, height=6)



fit = lm((TP/(TP+FN))~ErrLenT*n*Diameter*N,d[d$E %in% c( "16S.B_ErrLen","16S.B_NumErrAlns") & d$N > 19,])
af <- anova(fit)
afss <- af$"Sum Sq"
require(Hmisc)
options(digits=3)
latex(cbind(round(af[,1:4],2),Pvalue=round(af[,5:5],6),PctExp=round(afss/sum(afss)*100,1)),file = "16S-anova.tex")


ggplot(aes(x=DR,y=FN/(FN+TN),group=interaction(ErrLenT,n,sep=", "),color=interaction(ErrLenT,n,sep=", "),linetype="After filtering"),
    data=d[d$E %in% c( "16S.B_ErrLen","16S.B_NumErrAlns") & d$N > 19,])+
  #geom_boxplot(outlier.alpha = .5, outlier.size = 0.4)+#geom_point(alpha=0.5,size=1)+
  stat_summary(position = position_dodge(width=0.3),geom="linerange",size=0.8)+
  stat_summary(position = position_dodge(width=0.3),geom="line")+
  theme_classic()+
  theme(legend.position = "bottom",legend.direction = "horizontal", legend.text.align = 1)+
  scale_y_continuous("Percent error",labels=percent)+
  scale_shape(name="")+scale_x_discrete(name="")+
  facet_wrap(~E,labeller = function(x) list(E=c("Changing Error Length","Changing Error Frequency")),scales="free")+
  scale_linetype_manual(name="",values=c(1,3))+
  stat_summary(aes(y=(FN+TP)/(FN+FP+TP+TN),linetype="Before filtering"),position = position_dodge(width=0.3),alpha=0.99,geom="line")+
  scale_color_brewer(palette = "Paired",name="Error Length")
ggsave("Figures/ErrParam_Figures/16S.B_ErrLenNumErr_percenterror.pdf",width = 9,height = 4.5)



ggplot(aes(color=DR2,yend=FN/(FN+TN),y=(FN+TP)/(FN+FP+TP+TN),x=interaction(ErrLenT,n,sep=", "),xend=interaction(ErrLenT,n,sep=", ")),
       data=data.table::dcast(setDT(d[d$E %in% c( "16S.B_ErrLen","16S.B_NumErrAlns") & d$N > 19,]),ErrLenT+n+DR2+E~.,fun.aggregate = mean,value.var=c("FP","TP","TN","FN")))+
  #geom_boxplot(outlier.alpha = .5, outlier.size = 0.4)+#geom_point(alpha=0.5,size=1)+
  geom_segment(position = position_dodge2(width=0.8),size=0.8,arrow = arrow(length=unit(0.2,"cm")))+
  #stat_summary(position = position_dodge(width=0.3),geom="line")+
  theme_classic()+
  theme(legend.position = c(.27,.1),legend.direction = "horizontal", legend.text.align = 1)+
  scale_y_log10("Percent error",labels=percent)+
  scale_shape(name="")+scale_x_discrete(name="Error Len, Freq")+
  facet_wrap(~E,labeller = function(x) list(E=c("Changing Error Length","Changing Error Frequency")),scales="free_x",shrink = T)+
  scale_color_brewer(palette = "Dark2",name="Diameter")#+geom_hline(yintercept = 0.0003)
ggsave("Figures/ErrParam_Figures/16S.B_ErrLenNumErr_percenterror_arrow_log.pdf",width = 10,height =5)



bs = read.csv('CSV_Files/res_treecmp.csv',head=F)
names(bs) <- c("AlignmentName","DiameterRangeGene","X","Diameter","PD","N","ErrLen",
               "NumErrSeqDiv","Rep","ms_err","pd_err","rf_err","rfw_err","ms_res","pd_res","rf_res","rfw_res")

bs2 = cbind(melt(bs[c(1:9,10:13)],id.vars = 1:9), melt(bs[c(1:9,14:17)],id.vars = 1:9)[,10:11])
names(bs2)[c(12,13)]=c("v","after")
names(bs2)

bh = cbind(dcast(DiameterRangeGene+N+Diameter~.,data=bs[bs$AlignmentName=="HackettGenes",],value.var = "rf_res",fun.aggregate = mean),
           dcast(DiameterRangeGene~.,data=bs[bs$AlignmentName=="HackettGenes",],value.var = "rf_err",fun.aggregate = mean)[,2])
names(bh) = c("DiameterRangeGene","N","Diameter", "rf_res","rf_err")

bh = cbind(dcast(DiameterRangeGene+N+Diameter~.,data=bs[bs$AlignmentName=="HackettGenes",],value.var = "rfw_res",fun.aggregate = mean),
           dcast(DiameterRangeGene~.,data=bs[bs$AlignmentName=="HackettGenes",],value.var = "rfw_err",fun.aggregate = mean)[,2])
names(bh) = c("DiameterRangeGene","N","Diameter", "rf_res","rf_err")

ggplot(aes(color=N,x=DiameterRangeGene,y=-(rf_err-rf_res)/(2*N-6)),data=bs[bs$AlignmentName=="16S.B",])+geom_jitter(width = 0.2,alpha=0.5)+
  #geom_point(aes(y=rf_res/(2*N-6)),color="red")+
  theme_classic()+geom_boxplot(fill="transparent",color="red",outlier.alpha = 0,size=0.8)+
  scale_color_gradient(low = "#55BBFF",high = "#112244")+
  scale_y_continuous("Change in Normalized RF (after-before)",labels=percent)+
  scale_x_discrete(name="Diameter")+
  geom_hline(yintercept = 0)+
  theme(legend.position = c(.2,.9),legend.direction = "horizontal", legend.text.align = 1)+
  ggsave("Figures/ErrParam_Figures/16S-RF.pdf",width = 5.2,height = 5)


ggplot(aes(color=N,x=DiameterRangeGene,y=-(rfw_err-rfw_res)),data=bs[bs$AlignmentName=="16S.B",])+geom_jitter(width = 0.2,alpha=0.5)+
  #geom_point(aes(y=rf_res/(2*N-6)),color="red")+
  theme_classic()+geom_boxplot(fill="transparent",color="red",outlier.alpha = 0,size=0.8)+
  scale_color_gradient(low = "#55BBFF",high = "#112244")+
  scale_y_continuous("Change in  WRF (after-before)")+
  scale_x_discrete(name="Diameter")+
  geom_hline(yintercept = 0)+
  theme(legend.position = c(.2,.1),legend.direction = "horizontal", legend.text.align = 1)+
  ggsave("Figures/ErrParam_Figures/16S-RFw.pdf",width = 5.2,height = 5)

ggplot(aes(color=rfw_err,x=DiameterRangeGene,y=(rfw_err-rfw_res)/rfw_err),data=bs[bs$AlignmentName=="16S.B",])+geom_jitter(width = 0.2,alpha=0.5)+
  #geom_point(aes(y=rf_res/(2*N-6)),color="red")+
  theme_classic()+geom_boxplot(fill="transparent",color="red",outlier.alpha = 0,size=0.8)+
  scale_color_gradient(low = "#55BBFF",high = "#112244",name="Error before")+
  scale_y_continuous("Relative reduction in  WRF",labels=percent)+
  scale_x_discrete(name="Diameter")+
  geom_hline(yintercept = 0)+
  theme(legend.position = c(.3,.1),legend.direction = "horizontal", legend.text.align = 0)+
  ggsave("Figures/ErrParam_Figures/16S-RFw-change.pdf",width = 5.2,height = 5)





ggplot(aes(color=N,x=DiameterRangeGene,y=(value-after)),data=bs2[bs2$AlignmentName=="16S.B",])+geom_point()+
  #geom_point(aes(y=rf_res/(2*N-6)),color="red")+
  theme_bw()+stat_summary(fill="transparent",color="red")+
  facet_wrap(.~v,scales = "free")


ggplot(aes(yend=rf_res/(2*N-6),y=rf_err/(2*N-6),x=reorder(DiameterRangeGene,rf_err/(2*N-6)),xend=reorder(DiameterRangeGene,rf_err),
           shape=rf_res>rf_err,color=Diameter),
       data=bh)+
  #geom_boxplot(outlier.alpha = .5, outlier.size = 0.4)+#geom_point(alpha=0.5,size=1)+
  geom_segment(size=0.8,arrow = arrow(length=unit(0.2,"cm")))+
  #stat_summary(position = position_dodge(width=0.3),geom="line")+
  theme_classic()+
  theme(legend.position = c(.27,.891),legend.direction = "horizontal", legend.text.align = 1,axis.text.x = element_text(angle=90))+
  scale_y_continuous(name="Normalized RF error",labels=percent)+
  scale_shape(name="")+
  scale_x_discrete(name="")+
  scale_color_gradient(low = "#55BBFF",high = "#112244",name="Diameter")+
  ggsave("Figures/ErrParam_Figures/Hacket-RF.pdf",width = 5.2,height = 5)


ggplot(aes(yend=rf_res,y=rf_err,x=reorder(DiameterRangeGene,rf_err/(2*N-6)),xend=reorder(DiameterRangeGene,rf_err),
           shape=rf_res>rf_err,color=Diameter),
       data=bh)+
  #geom_boxplot(outlier.alpha = .5, outlier.size = 0.4)+#geom_point(alpha=0.5,size=1)+
  geom_segment(size=0.8,arrow = arrow(length=unit(0.2,"cm")))+
  #stat_summary(position = position_dodge(width=0.3),geom="line")+
  theme_classic()+
  theme(legend.position = c(.27,.891),legend.direction = "horizontal", legend.text.align = 1,axis.text.x = element_text(angle=90))+
  scale_y_continuous(name="WRF error")+
  scale_shape(name="")+
  scale_x_discrete(name="")+
  scale_color_gradient(low = "#55BBFF",high = "#112244",name="Diameter")+
  ggsave("Figures/ErrParam_Figures/Hacket-wRF.pdf",width = 5.2,height = 5)

ggplot(aes(yend=rf_res/(2*N-6),y=rf_err/(2*N-6),x=reorder(Rep,rf_err),xend=reorder(Rep,rf_err),
           shape=rf_res>rf_err,color=AlignmentName),
       data=bs[grepl("small-10-aa-RV100-BBA0039",bs$AlignmentName),])+
  #geom_boxplot(outlier.alpha = .5, outlier.size = 0.4)+#geom_point(alpha=0.5,size=1)+
  geom_segment(position = position_dodge(width=0.4),size=0.8,arrow = arrow(length=unit(0.1,"cm")))+
  #stat_summary(position = position_dodge(width=0.3),geom="line")+
  theme_classic()+
  theme(legend.position = c(.27,.891),legend.direction = "horizontal", legend.text.align = 1,axis.text.x = element_text(angle=90))+
  scale_y_continuous(name="Normalized RF error",labels=percent)+
  scale_shape(name="")+
  scale_x_discrete(name="")+
  scale_color_brewer(name="",palette = "Dark2",labels=c("TAPER","DivA"))+
  ggsave("Figures/ErrParam_Figures/AA-RF.pdf",width = 5.2,height = 5)

ggplot(aes(yend=rfw_res,y=rfw_err,x=reorder(Rep,rf_err/(2*N-6)),xend=reorder(Rep,rf_err),
           shape=rf_res>rf_err),
       data=bs[bs$AlignmentName=="small-10-aa-RV100-BBA0039",])+
  #geom_boxplot(outlier.alpha = .5, outlier.size = 0.4)+#geom_point(alpha=0.5,size=1)+
  geom_segment(position = position_dodge2(width=0.8),size=0.8,arrow = arrow(length=unit(0.2,"cm")),color="#3377AA")+
  #stat_summary(position = position_dodge(width=0.3),geom="line")+
  theme_classic()+
  theme(legend.position = c(.27,.891),legend.direction = "horizontal", legend.text.align = 1,axis.text.x = element_text(angle=90))+
  scale_y_continuous(name="WRF error")+
  scale_shape(name="")+
  scale_x_discrete(name="")+
  ggsave("Figures/ErrParam_Figures/AA-wRF.pdf",width = 5.2,height = 5)

  #facet_wrap(~E,labeller = function(x) list(E=c("Changing Error Length","Changing Error Frequency")),scales="free_x",shrink = T)+
  #scale_color_brewer(palette = "Dark2",name="Diameter",guide="none")#+geom_hline(yintercept = 0.0003)

ggplot(aes(x=DiameterRangeGene,y=(rf_err-rf_res)/(2*N-6)),data=bs)+geom_point()+
  #geom_point(aes(y=rf_res/(2*N-6)),color="red")+
  theme_bw()




ggplot(aes(color=DR2,yend=FN/(FN+TN),y=(FN+TP)/(FN+FP+TP+TN),x=interaction(ErrLenT,n,sep=", "),xend=interaction(ErrLenT,n,sep=", ")),
       data=data.table::dcast(setDT(d[d$E %in% c( "16S.B_ErrLen","16S.B_NumErrAlns") & d$N > 19,]),ErrLenT+n+DR2+E~.,fun.aggregate = mean,value.var=c("FP","TP","TN","FN")))+
  #geom_boxplot(outlier.alpha = .5, outlier.size = 0.4)+#geom_point(alpha=0.5,size=1)+
  geom_segment(position = position_dodge2(width=0.8),size=0.8,arrow = arrow(length=unit(0.2,"cm")))+
  #stat_summary(position = position_dodge(width=0.3),geom="line")+
  theme_classic()+
  theme(legend.position = c(.27,.1),legend.direction = "horizontal", legend.text.align = 1)+
  scale_y_sqrt("Percent error",labels=percent)+
  scale_shape(name="")+scale_x_discrete(name="Error Len, Freq")+
  facet_wrap(~E,labeller = function(x) list(E=c("Changing Error Length","Changing Error Frequency")),scales="free_x",shrink = T)+
  scale_color_brewer(palette = "Dark2",name="Diameter")#+geom_hline(yintercept = 0.0003)
ggsave("Figures/ErrParam_Figures/16S.B_ErrLenNumErr_percenterror_arrow_sqrt.pdf",width = 10,height =5)

ggplot(aes(color=DR2,yend=FN/(FN+TN),y=(FN+TP)/(FN+FP+TP+TN),x=interaction(ErrLenT,n,sep=", "),xend=interaction(ErrLenT,n,sep=", ")),
       data=data.table::dcast(setDT(d[d$E %in% c( "16S.B_ErrLen","16S.B_NumErrAlns") & d$N > 19,]),ErrLenT+n+DR2+E~.,fun.aggregate = mean,value.var=c("FP","TP","TN","FN")))+
  #geom_boxplot(outlier.alpha = .5, outlier.size = 0.4)+#geom_point(alpha=0.5,size=1)+
  geom_segment(position = position_dodge2(width=0.8),size=0.8,arrow = arrow(length=unit(0.2,"cm")))+
  #stat_summary(position = position_dodge(width=0.3),geom="line")+
  theme_classic()+
  theme(legend.position = c(.233,.87),legend.direction = "horizontal", legend.text.align = 1)+
  scale_y_continuous("Percent error",labels=percent)+
  scale_shape(name="")+scale_x_discrete(name="Error Len, Freq")+
  facet_wrap(~E,labeller = function(x) list(E=c("Changing Error Length","Changing Error Frequency")),scales="free_x",shrink = T)+
  scale_color_brewer(palette = "Dark2",name="D")#+geom_hline(yintercept = 0.0003)
ggsave("Figures/ErrParam_Figures/16S.B_ErrLenNumErr_percenterror_arrow.pdf",width = 10,height =5)


ggplot(aes(color=DR2, y=(FN+TN)/(FN+FP+TP+TN),x=interaction(ErrLenT,n,sep=", ")),
       data=d[d$E %in% c( "16S.B_ErrLen","16S.B_NumErrAlns") & d$N > 19,])+
  stat_summary(position = position_dodge(width=0.7),size=0.4)+
  theme_classic()+
  theme(legend.position = "bottom",legend.direction = "horizontal", legend.text.align = 1)+
  scale_y_continuous("Alignment size reduction",labels=percent)+
  scale_shape(name="")+scale_x_discrete(name="Error Length / Freq")+
  facet_wrap(~E,labeller = function(x) list(E=c("Changing Error Length","Changing Error Frequency")),ncol=2,scales="free_x")+
  scale_color_brewer(palette = "Dark2",name="Diameter")
ggsave("Figures/ErrParam_Figures/16SB_ErrLenNumErrAlns_sizechange.pdf",width = 10,height =5)


ggplot(aes(x=Diameter,y=TP/(TP+FN),color=interaction(ErrLenT,n,sep=", ")),
       data=d[d$E %in% c( "16S.B_ErrLen","16S.B_NumErrAlns") & d$N > 19,])+
  geom_point(alpha=0.4,size=.5)+
  theme_classic()+theme(legend.position = c(.75,.15),legend.direction = "horizontal", legend.text.align = 1)+
  geom_smooth()+scale_y_continuous("Recall")+
  scale_shape(name="")+facet_wrap(~E,labeller = function(x) list(E=c("Changing Error Length","Changing Error Frequency")))+
  scale_color_brewer(palette = "Paired",name="")
ggsave("Figures/ErrParam_Figures/16S.B_ErrLenNumErr_Recall.pdf",width = 9,height = 4.5)

ggplot(aes(x=N,y=TP/(TP+FN),color=interaction(ErrLenT,n,sep=", ")),data=d[d$E %in% c( "16S.B_ErrLen","16S.B_NumErrAlns") & d$N > 19,])+
  geom_point(alpha=0.4,size=.5)+
  theme_classic()+theme(legend.position = c(.75,.15),legend.direction = "horizontal", legend.text.align = 1)+
  geom_smooth(method="lm")+scale_y_continuous("Recall")+
  scale_shape(name="")+facet_wrap(~E,labeller = function(x) list(E=c("Changing Error Length","Changing Error Frequency")))+
  scale_color_brewer(palette = "Paired",name="")
ggsave("Figures/ErrParam_Figures/16S.B_ErrLenNumErr_Recall_N.pdf",width = 9,height = 4.5)

ggplot(aes(x=N,y=FP/(FP+TN),color=interaction(ErrLenT,n,sep=", ")),data=d[d$E %in% c( "16S.B_ErrLen","16S.B_NumErrAlns") & d$N > 19,])+
  geom_point(alpha=0.4,size=.5)+
  theme_classic()+theme(legend.position = "bottom",legend.direction = "horizontal", legend.text.align = 1)+
  geom_smooth(se=F)+scale_y_continuous("FPR")+
  scale_shape(name="")+facet_wrap(~E,labeller = function(x) list(E=c("Changing Error Length","Changing Error Frequency")))+
  scale_color_brewer(palette = "Paired",name="")
ggsave("Figures/ErrParam_Figures/16S.B_ErrLenNumErr_FPR_N.pdf",width = 9,height = 5)


ggplot(aes(x=Diameter,y=FP/(TP+FP),color=interaction(ErrLenT,n,sep=", ")),data=d[d$E %in% c( "16S.B_ErrLen","16S.B_NumErrAlns") & d$N > 19,])+
  geom_point(alpha=0.4,size=.5)+
  theme_classic()+theme(legend.position ="bottom",legend.direction = "horizontal", legend.text.align = 1)+
  geom_smooth()+scale_y_continuous("FDR")+
  scale_shape(name="")+facet_wrap(~E,labeller = function(x) list(E=c("Changing Error Length","Changing Error Frequency")))+
  scale_color_brewer(palette = "Paired",name="")
ggsave("Figures/ErrParam_Figures/16S.B_ErrLenNumErr_FDR.pdf",width = 9,height = 4.5)


ggplot(aes(x=Diameter,y=FP/(FP+TN),color=interaction(ErrLenT,n,sep=", ")),data=d[d$E %in% c( "16S.B_ErrLen","16S.B_NumErrAlns") & d$N > 19,])+
  geom_point(alpha=0.4,size=.5)+
  theme_classic()+theme(legend.position ="bottom",legend.direction = "horizontal", legend.text.align = 1)+
  geom_smooth()+scale_y_continuous("FPR",labels=percent)+
  scale_shape(name="")+facet_wrap(~E,labeller = function(x) list(E=c("Changing Error Length","Changing Error Frequency")))+
  scale_color_brewer(palette = "Paired",name="")
ggsave("Figures/ErrParam_Figures/16S.B_ErrLenNumErr_FPR.pdf",width = 9,height = 4.5)


# Recall vs Diameter

# 16S.B - Varying error lengths and fixed percentage of erroneous sequences
ggplot(aes(x=Diameter,y=TP/(TP+FN),color=ErrLenT),data=d[d$E=="16S.B_ErrLen" & d$N > 19,])+
  geom_point(alpha=0.5)+theme_classic()+geom_smooth()+scale_y_continuous("Recall")+
  scale_shape(name="")+scale_color_brewer(palette = "Paired",name="error len")+
  ggtitle("16S.B with Varying Error Lengths: Recall vs Diameter")
ggsave("Figures/ErrParam_Figures/16S.B_ErrLen_Recall.pdf",width = 6,height = 6)


# 16S.B - Varying percentage of erroneous sequences and fixed error lengths
ggplot(aes(x=Diameter,y=TP/(TP+FN), group= as.factor(100/((NumErrSeqDiv!=N)*NumErrSeqDiv+(NumErrSeqDiv==N)*100)),
  color=as.factor(100/((NumErrSeqDiv!=N)*NumErrSeqDiv+(NumErrSeqDiv==N)*100)), shape=cut((FP/(FP+TN)),breaks=c(-1,0,0.001,0.1,1))),
  data=d[d$E=="16S.B_NumErrAlns" ,])+geom_point(alpha=0.5)+
  theme_classic()+geom_smooth(se=F)+scale_shape_manual(name="FPR",values=c(1,16,4))+
  scale_color_brewer(palette = "Paired",name="n",labels=nlabels)+
  scale_y_continuous(name="Recall")+coord_cartesian(ylim=c(0.35,1))+
  ggtitle("16S.B with Varying Number of Erroneous Sequences: Recall vs Diameter")
ggsave("Figures/ErrParam_Figures/16S.B_NumErrAlns_Recall.pdf",width = 6,height = 6)






############################## Hacket


ggplot(aes(x=reorder(paste(DR,round(Diameter,3),round(SL,0),as.numeric(as.character(N)),sep="\n"),TP/(TP+FN)),y=TP/(TP+FN),color=interaction(ErrLenT,n,sep=", ")),
       data=d[d$E %in% c("Hackett_Genes_ErrLen","Hackett_Genes_NumErrAlns") & d$ErrLen<64,])+
  stat_summary(position = position_dodge(width=0.6))+
  #geom_point(alpha=0.5,size=1)+
  theme_bw()+theme(legend.position = "bottom",legend.direction = "horizontal", legend.text.align = 1)+
  scale_y_continuous("Recall",labels=percent)+
  scale_shape(name="")+scale_x_discrete(name="Gene")+facet_wrap(~E,labeller = function(x) list(E=c("Changing Error Length","Changing Error Frequency")),ncol=1)+
  scale_color_brewer(palette = "Paired",name="Error Frequency")#+stat_summary(aes(label=round(..y..,2)),geom="text",position = "jitter")
ggsave("Figures/ErrParam_Figures/Hackett_NumErrErrLen_Recall.pdf",width = 9,height = 8)

ggplot(aes(x=reorder(paste(DR,round(Diameter,3),round(SL,0),as.numeric(as.character(N)),sep="\n"),TP/(TP+FN)),y=FP/(TP+FP),color=interaction(ErrLenT,n,sep=", ")),
       data=d[d$E %in% c("Hackett_Genes_ErrLen","Hackett_Genes_NumErrAlns") & d$ErrLen<64,])+
  stat_summary(position = position_dodge(width=0.6))+
  #geom_point(alpha=0.5,size=1)+
  theme_bw()+theme(legend.position = "bottom",legend.direction = "horizontal", legend.text.align = 1)+
  scale_y_continuous("FDR",labels=percent)+
  scale_shape(name="")+scale_x_discrete(name="Gene")+facet_wrap(~E,labeller = function(x) list(E=c("Changing Error Length","Changing Error Frequency")),ncol=1)+
  scale_color_brewer(palette = "Paired",name="Error Frequency")
ggsave("Figures/ErrParam_Figures/Hackett_NumErrErrLen_FDR.pdf",width = 9,height = 8)

ggplot(aes(x=reorder(paste(DR,round(Diameter,3),round(SL,0),as.numeric(as.character(N)),sep="\n"),TP/(TP+FN)),y=FP/(TN+FP),color=interaction(ErrLenT,n,sep=", ")),
       data=d[d$E %in% c("Hackett_Genes_ErrLen","Hackett_Genes_NumErrAlns") & d$ErrLen<64,])+
  stat_summary(position = position_dodge(width=0.6))+
  #geom_point(alpha=0.5,size=1)+
  theme_bw()+theme(legend.position = "bottom",legend.direction = "horizontal", legend.text.align = 1)+
  scale_y_continuous("FPR",labels=percent)+
  scale_shape(name="")+scale_x_discrete(name="Gene")+facet_wrap(~E,labeller = function(x) list(E=c("Changing Error Length","Changing Error Frequency")),ncol=1,scales="free_y")+
  scale_color_brewer(palette = "Paired",name="Error Frequency")
ggsave("Figures/ErrParam_Figures/Hackett_NumErrErrLen_FPR.pdf",width = 9,height = 8)


ggplot(aes(x=reorder(paste(DR,round(Diameter,3),round(SL,0),as.numeric(as.character(N)),sep="\n"),TP/(TP+FN)),y=TP/(TP+FN),color=n),data=d[d$E =="Hackett_Genes_NumErrAlns",])+ # %in% c( "Hackett_ErrLen","Hackett_NumErrAlns","Hackett_General") ,])+
  stat_summary(position = position_dodge(width=0.6))+
  #geom_point(alpha=0.5,size=1)+
  theme_bw()+theme(legend.position = "bottom",legend.direction = "horizontal", legend.text.align = 1)+
  scale_y_continuous("Recall",labels=percent)+
  scale_shape(name="")+scale_x_discrete(name="Gene")+
  scale_color_brewer(palette = "Paired",name="Error Frequency")
ggsave("Figures/ErrParam_Figures/Hackett_NumErr_Recall.pdf",width = 9,height = 6)



ggplot(aes(x=Diameter,y=TP/(TP+FN),color=interaction(ErrLenT,n,sep=", ")),
       data=d[d$E %in% c( "Hackett_Genes_ErrLen","Hackett_Genes_NumErrAlns"),])+
  stat_summary(position = position_dodge(width=0.01),alpha=0.99)+
  #geom_point(alpha=0.5,size=1)+
  theme_bw()+theme(legend.position = "bottom",legend.direction = "horizontal", legend.text.align = 1)+
  geom_smooth(se=F,method="lm")+scale_y_continuous("Recall",labels=percent)+
  scale_shape(name="")+#scale_x_discrete(name="Gene")+
  scale_color_brewer(palette = "Paired",name="Error Len, Freq")+
  #facet_wrap(~E,labeller = function(x) list(E=c("Changing Error Length","Changing Error Frequency")),ncol=1)+
  geom_text(aes(label=DR,y=rep(c(0.3,0.3,0.3,0.3,0.3),4)[1:19]),data=d[d$E =="Hackett_Genes_NumErrAlns"  & d$n=="2%" &d$Rep==1,],
            position = position_jitter(width = 0,height = 0.0),color="black")
ggsave("Figures/ErrParam_Figures/Hackett_NumErrErrLen_Recall_vs_Diameter_2.pdf",width = 9,height = 4.5)



ggplot(aes(x=SL,y=TP/(TP+FN),color=interaction(ErrLenT,n,sep=", ")),
       data=d[d$E %in% c( "Hackett_Genes_ErrLen","Hackett_Genes_NumErrAlns") & d$DR != "concat",])+
  stat_summary(position = position_dodge(width=0.01),alpha=0.8)+
  #geom_point(alpha=0.5,size=1)+
  theme_bw()+theme(legend.position = "bottom",legend.direction = "horizontal", legend.text.align = 1)+
  geom_smooth(se=F,method="lm")+scale_y_continuous("Recall",labels=percent)+
  scale_shape(name="")+scale_x_continuous(name="Sequence Length")+
  scale_color_brewer(palette = "Paired",name="Error Len, Freq")+
  #facet_wrap(~E,labeller = function(x) list(E=c("Changing Error Length","Changing Error Frequency")),ncol=1)+
  geom_text(aes(label=DR,y=rep(c(0.3,0.3,0.3,0.3,0.3),4)[1:19]),data=d[d$E =="Hackett_Genes_NumErrAlns"& d$DR != "concat"  & d$n=="2%" &d$Rep==1 ,],
            position = position_jitter(width = 0,height = 0.0),color="black")
ggsave("Figures/ErrParam_Figures/Hackett_NumErrErrLen_Recall_vs_SL_2.pdf",width = 9,height = 4.5)


ggplot(aes(x=N,y=TP/(TP+FN),color=interaction(ErrLenT,n,sep=", "),group=interaction(DR,ErrLenT,n,sep=", ")),
       data=d[d$E %in% c( "Hackett_Genes_ErrLen","Hackett_Genes_NumErrAlns") & d$DR != "concat",])+
  stat_summary(position = position_dodge(width=0.01),alpha=0.75)+
  #geom_point(alpha=0.5,size=1)+
  theme_bw()+theme(legend.position = "bottom",legend.direction = "horizontal", legend.text.align = 1)+
  geom_smooth(aes(group=interaction(ErrLenT,n,sep=", ")),se=F,method="lm")+scale_y_continuous("Recall",labels=percent)+
  scale_x_continuous(name="Sequence count")+
  scale_color_brewer(palette = "Paired",name="Error Len, Freq")+
  #facet_wrap(~E,labeller = function(x) list(E=c("Changing Error Length","Changing Error Frequency")),ncol=1)+
  geom_text(aes(label=DR,y=rep(c(0.3,0.3,0.3,0.3,0.3),4)[1:19]),data=d[d$E =="Hackett_Genes_NumErrAlns"& d$DR != "concat"  & d$n=="2%" &d$Rep==1 ,],
            position = position_jitter(width = 0,height = 0.0),color="black")
ggsave("Figures/ErrParam_Figures/Hackett_NumErrErrLen_Recall_vs_N_2.pdf",width = 9,height = 4.5)


fit = lm((TP/(TP+FN))~ErrLenT*n*Diameter*SL*N,d[d$E %in% c( "Hackett_Genes_ErrLen","Hackett_Genes_NumErrAlns") ,])
af <- anova(fit)
afss <- af$"Sum Sq"
require(Hmisc)
options(digits=3)
latex(cbind(round(af[,1:4],2),Pvalue=round(af[,5:5],6),PctExp=round(afss/sum(afss)*100,1)),file = "anova-hackett.tex")

ggplot(aes(x=reorder(paste(DR,round(Diameter,3),round(SL,0),as.numeric(as.character(N)),sep="\n"),
                     -(FN+TP)/(FN+FP+TP+TN)#-FN/(FN+TN)
                     #SL*as.numeric(as.character(N))
                     #Diameter
                     ),y=FN/(FN+TN),group=ErrLenT,color=ErrLenT,linetype="After filtering"),data=d[d$E =="Hackett_Genes_ErrLen" &d$DR != "concatenation"&d$ErrLen<64,])+ # %in% c( "Hackett_ErrLen","Hackett_NumErrAlns","Hackett_General") ,])+
  #geom_boxplot(outlier.alpha = .5, outlier.size = 0.4)+#geom_point(alpha=0.5,size=1)+
  stat_summary(position = position_dodge(width=0.3),geom="linerange",size=0.8)+
  stat_summary(position = position_dodge(width=0.3),geom="line")+
  theme_classic()+theme(legend.position = "bottom",legend.direction = "horizontal", legend.text.align = 1)+
  scale_y_continuous("Percent error",labels=percent)+
  scale_shape(name="")+scale_x_discrete(name="Gene")+
  scale_linetype_manual(name="",values=c(1,3))+
  stat_summary(aes(y=(FN+TP)/(FN+FP+TP+TN),linetype="Before filtering"),position = position_dodge(width=0.3),alpha=0.9,geom="line")+
  scale_color_brewer(palette = "Dark2",name="Error Length")
ggsave("Figures/ErrParam_Figures/Hackett_NumErr_percenterror.pdf",width = 9,height = 6.5)


ggplot(aes(x=reorder(paste(DR,round(Diameter,3),round(SL,0),as.numeric(as.character(N)),sep="\n"),
                     SL),
           xend=reorder(paste(DR,round(Diameter,3),round(SL,0),as.numeric(as.character(N)),sep="\n"),
                        SL),
           yend=FN/(FN+TN),y=(FN+TP)/(FN+FP+TP+TN),color=interaction(ErrLenT,n,sep=", ")),
       data=data.table::dcast(setDT(d[ (d$E =="Hackett_Genes_ErrLen" | d$E=="Hackett_Genes_NumErrAlns") &d$ErrLen<64& d$DR!="concat",]),ErrLenT+n+SL+N+Diameter+DR+E~.,fun.aggregate = mean,value.var=c("FP","TP","TN","FN")))+
  geom_segment(position = position_dodge(width=0.8),size=0.8,arrow = arrow(length=unit(0.2,"cm")))+
  theme_classic()+
  theme(legend.position = "bottom",legend.direction = "horizontal", legend.text.align = 1)+
  scale_y_log10("Percent error",labels=percent)+
  scale_shape(name="")+scale_x_discrete(name="")+
  facet_wrap(~E,labeller = function(x) list(E=c("Changing Error Length","Changing Error Frequency")),ncol=1)+
  scale_color_brewer(palette = "Paired",name="Error Length / Freq")
ggsave("Figures/ErrParam_Figures/Hackett_ErrLenNumErrAlns_percenterror_arrow_log.pdf",width = 9,height =9)


ggplot(aes(x=reorder(paste(DR,round(Diameter,3),round(SL,0),as.numeric(as.character(N)),sep="\n"),
                     SL),
           y=(FN+TN)/(FN+FP+TP+TN),color=interaction(ErrLenT,n,sep=", ")),
       data=d[ (d$E =="Hackett_Genes_ErrLen" | d$E=="Hackett_Genes_NumErrAlns") &d$ErrLen<64,])+
  stat_summary(position = position_dodge(width=0.8))+
  theme_classic()+
  theme(legend.position = "bottom",legend.direction = "horizontal", legend.text.align = 1)+
  scale_y_continuous("Alignment size reduction",labels=percent)+
  scale_shape(name="")+scale_x_discrete(name="")+
  facet_wrap(~E,labeller = function(x) list(E=c("Changing Error Length","Changing Error Frequency")),ncol=1)+
  scale_color_brewer(palette = "Paired",name="Error Length / Freq")
ggsave("Figures/ErrParam_Figures/Hackett_ErrLenNumErrAlns_sizechange.pdf",width = 9,height =9)

ggplot(aes(x=reorder(paste(DR,round(Diameter,3),round(SL,0),as.numeric(as.character(N)),sep="\n"),
                     -(FN+TP)/(FN+FP+TP+TN)),
           xend=reorder(paste(DR,round(Diameter,3),round(SL,0),as.numeric(as.character(N)),sep="\n"),
                        -(FN+TP)/(FN+FP+TP+TN)),
           yend=FN/(FN+TN),y=(FN+TP)/(FN+FP+TP+TN),color=interaction(ErrLenT,n,sep=", ")),
       data=data.table::dcast(setDT(d[ (d$E =="Hackett_Genes_ErrLen" | d$E=="Hackett_Genes_NumErrAlns") &d$ErrLen<64,]),ErrLenT+n+SL+N+Diameter+DR+E~.,fun.aggregate = mean,value.var=c("FP","TP","TN","FN")))+
  geom_segment(position = position_dodge(width=0.8),size=0.8,arrow = arrow(length=unit(0.2,"cm")))+
  theme_classic()+
  theme(legend.position = "bottom",legend.direction = "horizontal", legend.text.align = 1)+
  scale_y_continuous("Percent error",labels=percent)+
  scale_shape(name="")+scale_x_discrete(name="")+
  facet_wrap(~E,labeller = function(x) list(E=c("Changing Error Length","Changing Error Frequency")),ncol=1)+
  scale_color_brewer(palette = "Paired",name="Error Length / Freq")
ggsave("Figures/ErrParam_Figures/Hackett_ErrLenNumErrAlns_percenterror_arrow.pdf",width = 9,height =9)


ggplot(aes(x=reorder(paste(DR,round(Diameter,3),round(SL,0),as.numeric(as.character(N)),sep="\n"), -(FN+TP)/(FN+FP+TP+TN)),y=(FN)/(FN+TN),group=ErrLenT,color=ErrLenT,linetype="After filtering"),data=d[d$E =="Hackett_Genes_ErrLen" &d$DR != "concatenation"&d$ErrLen<64,])+ # %in% c( "Hackett_ErrLen","Hackett_NumErrAlns","Hackett_General") ,])+
  #geom_boxplot(outlier.alpha = .5, outlier.size = 0.4)+#geom_point(alpha=0.5,size=1)+
  stat_summary(position = position_dodge(width=0.3),geom="linerange",size=0.2)+
  stat_summary(position = position_dodge(width=0.3),geom="line")+
  theme_classic()+theme(legend.position = "bottom",legend.direction = "horizontal", legend.text.align = 1)+
  scale_y_sqrt("Percent error",labels=percent,breaks=c(.1,.2,.3,0.5,1,2,3,4)*0.01)+
  scale_shape(name="")+scale_x_discrete(name="Gene")+
  scale_linetype_manual(name="",values=c(1,3))+
  stat_summary(aes(y=(FN+TP)/(FN+FP+TP+TN),linetype="Before filtering"),position = position_dodge(width=0.3),alpha=0.9,geom="line")+
  scale_color_brewer(palette = "Dark2",name="Error Length")
ggsave("Figures/ErrParam_Figures/Hackett_NumErr_percenterror_sqrt.pdf",width = 8,height = 6.5)

ggplot(aes(x=reorder(paste(DR,round(Diameter,3),round(SL,0),as.numeric(as.character(N)),sep="\n"),TP/(TP+FN)),y=TP/(TP+FN),color=ErrLenT),data=d[d$E =="Hackett_Genes_ErrLen" &d$ErrLen<64,])+ # %in% c( "Hackett_ErrLen","Hackett_NumErrAlns","Hackett_General") ,])+
  #geom_boxplot(outlier.alpha = .5, outlier.size = 0.4)+#geom_point(alpha=0.5,size=1)+
  stat_summary(position = position_dodge(width=0.7))+
  theme_bw()+theme(legend.position = "bottom",legend.direction = "horizontal", legend.text.align = 1)+
  scale_y_continuous("Recall",labels=percent)+
  scale_shape(name="")+scale_x_discrete(name="Gene")+
  scale_color_brewer(palette = "Paired",name="Error Length")
  #scale_fill_brewer(palette = "Spectral",name="Error Length")
ggsave("Figures/ErrParam_Figures/Hackett_ErrLen_Recall.pdf",width = 9,height = 6)

ggplot(aes(x=Diameter,y=TP/(TP+FN),color=ErrLenT),data=d[d$E =="Hackett_Genes_ErrLen" &d$DR != "concatenation"&d$ErrLen<64,])+ # %in% c( "Hackett_ErrLen","Hackett_NumErrAlns","Hackett_General") ,])+
  stat_summary(position = position_dodge(width=0.01))+
  #geom_point(alpha=0.5,size=1)+
  theme_bw()+theme(legend.position = "bottom",legend.direction = "horizontal", legend.text.align = 1)+
  geom_smooth(se=F,method="lm")+scale_y_continuous("Recall",labels=percent)+
  scale_shape(name="")+#scale_x_discrete(name="Gene")+
  scale_color_brewer(palette = "Paired",name="Error Length")+
  geom_text(aes(label=DR,y=0.3),data=d[d$E =="Hackett_Genes_NumErrAlns" &d$DR != "concatenation" & d$n=="2%" &d$Rep==1,],
            position = position_jitter(width = 0,height = 0.07))
ggsave("Figures/ErrParam_Figures/Hackett_ErrLen_Recall_vs_Diameter.pdf",width = 9,height = 5)


ggplot(aes(x=SL,y=TP/(TP+FN),color=ErrLenT),data=d[d$E =="Hackett_Genes_ErrLen" &d$DR != "concat"&d$ErrLen<64,])+ # %in% c( "Hackett_ErrLen","Hackett_NumErrAlns","Hackett_General") ,])+
  stat_summary(position = position_dodge(width=0.01))+
  #geom_point(alpha=0.5,size=1)+
  theme_bw()+theme(legend.position = "bottom",legend.direction = "horizontal", legend.text.align = 1)+
  geom_smooth(se=F,method="lm")+scale_y_continuous("Recall",labels=percent)+
  scale_shape(name="")+#scale_x_discrete(name="Gene")+
  scale_color_brewer(palette = "Paired",name="Error Length")+
  geom_text(aes(label=DR,y=0.3),data=d[d$E =="Hackett_Genes_NumErrAlns" &d$DR != "concat" & d$n=="2%" &d$Rep==1,],
            position = position_jitter(width = 0,height = 0.07))
ggsave("Figures/ErrParam_Figures/Hackett_ErrLen_Recall_vs_SL.pdf",width = 9,height = 5)

ggplot(aes(x=as.numeric(as.character(N)),y=TP/(TP+FN),color=ErrLenT),data=d[d$E =="Hackett_Genes_ErrLen" &d$DR != "concatenation"&d$ErrLen<64,])+ # %in% c( "Hackett_ErrLen","Hackett_NumErrAlns","Hackett_General") ,])+
  stat_summary(position = position_dodge(width=0.01))+
  #geom_point(alpha=0.5,size=1)+
  theme_bw()+theme(legend.position = "bottom",legend.direction = "horizontal", legend.text.align = 1)+
  geom_smooth(se=F,method="lm")+scale_y_continuous("Recall",labels=percent)+
  scale_shape(name="")+scale_x_discrete(name="Number of sequences")+
  scale_color_brewer(palette = "Paired",name="Error Length")+
  geom_text(aes(label=DR,y=0.3),data=d[d$E =="Hackett_Genes_NumErrAlns" &d$DR != "concatenation" & d$n=="2%" &d$Rep==1,],
            position = position_jitter(width = 0,height = 0.07))
ggsave("Figures/ErrParam_Figures/Hackett_ErrLen_Recall_vs_N.pdf",width = 9,height = 5)


# ROC for Hackett with varying error lengths and fixed percentage of erroneous sequences
options(digits = 2)
d2=summ_roc(d[d$E %in% c( "Hackett_Genes_ErrLen","Hackett_Genes_General","Hackett_Genes_NumErrAlns") ,], ErrLenT+n~.)
A = data.frame(x=d2$FP/(d2$FP+d2$TN),y=d2$TP/(d2$TP+d2$FN), n=d2$n, ErrLen=d2$ErrLenT)
B = data.frame(x=d2$FP0/(d2$FP0+d2$TN0),y=as.vector(matrix(1.1,nrow=nrow(d2))), n=d2$n, ErrLen=d2$ErrLenT)
ggplot(data=A, aes(x, y, color=ErrLen, shape=n)) + geom_point(alpha=1)+
  theme_light()+theme(legend.position = "right")+
  scale_shape(name="Error Frequency")+scale_color_brewer(name="Error Length",palette = "Paired")+
  scale_x_continuous(name="FPR",labels=percent)+
  scale_y_continuous("Recall",labels=percent)+
ggsave("Figures/ErrParam_Figures/Hackett_ErrLenFreq_ROC.pdf", width=6, height=6)

options(digits = 2)
d2=summ_roc(d[d$E %in% c( "Hackett_Genes_ErrLen","Hackett_Genes_NumErrAlns") ,], ErrLenT+n+DR~.)
A = data.frame(x=d2$FP/(d2$FP+d2$TN),y=d2$TP/(d2$TP+d2$FN), n=d2$n, ErrLen=d2$ErrLenT, DR=d2$DR)
ggplot(data=A, aes(x, y, size=n, shape=ErrLen,color=reorder(DR,y/x))) + geom_point(alpha=.99)+
  theme_light()+theme(legend.position = "right",axis.text = element_text(size=15),axis.title = element_text(size=16))+
  scale_shape_manual(name="Error Length",values=c(6,5,4,1,3,2,19,17,18,16,15,14,90))+scale_color_discrete(name="Gene")+
  scale_x_continuous(name="FPR",labels=percent)+
  scale_size_manual(name="Error Frequency",values=sqrt(c(1,4,10,20,40)))+
  scale_y_continuous("Recall",labels=percent)
ggsave("Figures/ErrParam_Figures/Hackett_ErrLenFreq_ROC_genes.pdf", width=9.5, height=9.5)

# # ROC for Hackett with varying percentages of erroneous sequences and fixed error lengths
# d2=d[d$E=="Hackett_NumErrAlns",]
# d2$n=with(d2,as.factor(100/((NumErrSeqDiv!=N)*NumErrSeqDiv+(NumErrSeqDiv==N)*100)))
# ggplot(aes(x=FP/(FP+TN),y=TP/(TP+FN), color=as.factor(n) ),data=summ_roc(d2,n~.))+
#   geom_point(alpha=1)+
#   theme_light()+theme(legend.position = c(.85,.25))+
#   scale_color_brewer(name="n",labels=nlabels, palette="Paired")+
#   scale_x_continuous(name="FPR",labels=percent)+scale_y_continuous("Recall")+
#   ggtitle("Hackett with Varying Number of Erroneous Sequences: ROC")
# ggsave("Figures/ErrParam_Figures/Hackett_NumErrAlns_ROC.pdf", width=6, height=6)

# # Hackett - Varying error lengths and fixed percentage of erroneous sequences
# ggplot(aes(x=Diameter,y=TP/(TP+FN),color=as.factor(ErrLen)),data=d[d$E=="Hackett_ErrLen" & d$N > 19,])+
#   geom_point(alpha=0.5)+theme_classic()+geom_smooth(se=F)+scale_y_continuous("Recall")+
#   scale_shape(name="")+scale_color_brewer(palette = "Paired",name="error len")+
#   ggtitle("Hackett with Varying Error Lengths: Recall vs Diameter")
# ggsave("Figures/ErrParam_Figures/Hackett_ErrLen_Recall.pdf",width = 6,height = 6)
# 
# # Hackett - Varying percentage of erroneous sequences and fixed error lengths
# ggplot(aes(x=Diameter,y=TP/(TP+FN), group= as.factor(100/((NumErrSeqDiv!=N)*NumErrSeqDiv+(NumErrSeqDiv==N)*100)),
#            color=as.factor(100/((NumErrSeqDiv!=N)*NumErrSeqDiv+(NumErrSeqDiv==N)*100)), shape=cut((FP/(FP+TN)),breaks=c(-1,0,0.001,0.1,1))),
#        data=d[d$E=="Hackett_NumErrAlns" ,])+geom_point(alpha=0.5)+
#   theme_classic()+geom_smooth(se=F)+scale_shape_manual(name="FPR",values=c(1,16,4))+
#   scale_color_brewer(palette = "Paired",name="n",labels=nlabels)+
#   scale_y_continuous(name="Recall")+coord_cartesian(ylim=c(0.35,1))+
#   ggtitle("Hackett with Varying Number of Erroneous Sequences: Recall vs Diameter")
# ggsave("Figures/ErrParam_Figures/Hackett_NumErrAlns_Recall.pdf",width = 6,height = 6)



ggplot(aes(x=n,y=TP/(TP+FN),color=ErrLenT),data=d[d$E %in% c( "small-10-aa_ErrLen","small-10-aa_NumErrAlns") ,])+
  geom_boxplot()+
  theme_classic()+theme(legend.position = c(.85,.25),legend.direction = "vertical", legend.text.align = 1)+
  geom_smooth(se=F,method="lm")+scale_y_continuous("Recall",labels=percent)+
  scale_shape(name="")+scale_x_discrete(name="Error Frequency")+
  scale_color_brewer(palette = "Paired",name="")+
  geom_boxplot(aes(x="5%",fill="DivA",color="8×11"),data=t)+
  scale_fill_discrete(name="")+
ggsave("Figures/ErrParam_Figures/small-10-aa_ErrLenNumErr_Recall.pdf",width = 5,height = 4.5)

ggplot(aes(x=n,y=FP/(TP+FP),color=ErrLenT),data=d[d$E %in% c( "small-10-aa_ErrLen","small-10-aa_NumErrAlns"),])+
  geom_boxplot()+
  theme_classic()+theme(legend.position = c(.85,.85),legend.direction = "vertical", legend.text.align = 1)+
  geom_smooth(se=F,method="lm")+scale_y_continuous("FDR",labels=percent)+scale_x_discrete(name="Error Frequency")+
  scale_shape(name="")+#facet_wrap(~E,labeller = function(x) list(E=c("Changing Error Length","Changing Error Frequency")))+
  scale_color_brewer(palette = "Paired",name="")+
  geom_boxplot(aes(x="5%",fill="DivA",color="8×11"),data=t)
ggsave("Figures/ErrParam_Figures/small-10-aa_ErrLenNumErr_FDR.pdf",width = 5,height = 4.5)

ggplot(aes(x=n,y=FP/(TN+FP),color=ErrLenT),data=d[d$E %in% c( "small-10-aa_ErrLen","small-10-aa_NumErrAlns"),])+
  geom_boxplot()+
  theme_classic()+theme(legend.position =  c(.85,.85),legend.direction = "vertical", legend.text.align = 1)+
  geom_smooth(se=F,method="lm")+scale_y_continuous("FPR",labels=percent)+
  scale_shape(name="")+#facet_wrap(~E,labeller = function(x) list(E=c("Changing Error Length","Changing Error Frequency")))+
  scale_color_brewer(palette = "Paired",name="")+scale_x_discrete(name="Error Frequency")+
  geom_boxplot(aes(x="5%",fill="DivA",color="8×11"),data=t)
ggsave("Figures/ErrParam_Figures/small-10-aa_ErrLenNumErr_FPR.pdf",width = 5,height = 4.5)


ggplot(aes(x=interaction(ErrLenT,n,sep=","),color="TAPER",xend=interaction(ErrLenT,n,sep=","),yend=FN/(FN+TN),y=(FN+TP)/(FN+FP+TP+TN)), #,color=interaction(ErrLenT,n,sep=", ")),
       data=data.table::dcast(setDT(d[d$E %in% c( "small-10-aa_ErrLen","small-10-aa_NumErrAlns"),]),ErrLenT+n+E~.,fun.aggregate = mean,value.var=c("FP","TP","TN","FN")))+
  #geom_boxplot(outlier.alpha = .5, outlier.size = 0.4)+#geom_point(alpha=0.5,size=1)+
  geom_segment(position = position_dodge(width=0.8),size=0.8,arrow = arrow(length=unit(0.2,"cm")))+
  #stat_summary(position = position_dodge(width=0.3),geom="line")+
  theme_classic()+
  theme(legend.position = "bottom",legend.direction = "horizontal", legend.text.align = 1)+
  scale_y_log10("Percent error",labels=percent)+
  scale_shape(name="")+scale_x_discrete(name="Error Len, Freq.")+
  facet_wrap(~E,labeller = function(x) list(E=c("Changing Error Length","Changing Error Frequency")),scales="free_x")+
  scale_color_brewer(palette = "Dark2",name="Method")+
  geom_segment(aes(color="DivA",x="8×11,5%",xend="8×11,5%"),
               data=data.table::dcast(setDT(t),ErrLen~.,fun.aggregate = mean,value.var=c("FP","TP","TN","FN")),
               position = position_nudge(x=0.2),size=0.8,arrow = arrow(length=unit(0.2,"cm")))+
ggsave("Figures/ErrParam_Figures/small-10-aa_ErrLenNumErr_percenterror_arrow_log.pdf",width = 7.5,height =4)

ggplot(aes(x=interaction(ErrLenT,n,sep=", "),color="TAPER",
           y=(FN+TN)/(FN+FP+TP+TN)), #,color=interaction(ErrLenT,n,sep=", ")),
       data=d[d$E %in% c( "small-10-aa_ErrLen","small-10-aa_NumErrAlns"),])+
  stat_summary(position = position_dodge(width=0.8))+
  theme_classic()+
  theme(legend.position = "bottom",legend.direction = "horizontal", legend.text.align = 1)+
  scale_y_continuous("Alignment size reduction",labels=percent)+
  scale_shape(name="")+scale_x_discrete(name="")+
  facet_wrap(~E,labeller = function(x) list(E=c("Changing Error Length","Changing Error Frequency")),scales="free_x")+
  scale_color_brewer(palette = "Dark2",name="Method")+
  stat_summary(aes(color="DivA",x="8×11, 5%"),data=t)+
ggsave("Figures/ErrParam_Figures/small-10-aa_ErrLenNumErrAlns_sizechange.pdf",width = 7.5,height =4)

# # small-10-aa - Varying error lengths and fixed percentage of erroneous sequences
# ggplot(aes(x=n,y=TP/(TP+FN),color=as.factor(ErrLen)),data=d[d$E=="small-10-aa_ErrLen" ,])+
#   geom_point(alpha=0.5)+theme_classic()+geom_smooth()+scale_y_continuous("Recall")+
#   scale_shape(name="")+scale_color_brewer(palette = "Paired",name="error len")+
#   ggtitle("small-10-aa with Varying Error Lengths: Recall vs Diameter")
# ggsave("Figures/ErrParam_Figures/small10aa_ErrLen_Recall.pdf",width = 6,height = 6)
# 
# # small-10-aa - Varying percentage of erroneous sequences and fixed error lengths
# ggplot(aes(x=Diameter,y=TP/(TP+FN), group= as.factor(100/((NumErrSeqDiv!=N)*NumErrSeqDiv+(NumErrSeqDiv==N)*100)),
#            color=as.factor(100/((NumErrSeqDiv!=N)*NumErrSeqDiv+(NumErrSeqDiv==N)*100)), shape=cut((FP/(FP+TN)),breaks=c(-1,0,0.001,0.1,1))),
#        data=d[d$E=="small-10-aa_NumErrAlns" ,])+geom_point(alpha=0.5)+
#   theme_classic()+geom_smooth(se=F)+scale_shape_manual(name="FPR",values=c(1,16,4,2))+
#   scale_color_brewer(palette = "Paired",name="n",labels=nlabels)+
#   scale_y_continuous(name="Recall")+coord_cartesian(ylim=c(0.35,1))+
#   ggtitle("small-10-aa with Varying Number of Erroneous Sequences: Recall vs Diameter")
# ggsave("Figures/ErrParam_Figures/small10aa_NumErrAlns_Recall.pdf",width = 6,height = 6)


# ROC for small-10-aa with varying error lengths and fixed percentage of erroneous sequences
options(digits = 2)
d2=summ_roc(d[d$E %in% c( "small-10-aa_ErrLen","small-10-aa_NumErrAlns"),], ErrLenT+n~.)
A = data.frame(x=d2$FP/(d2$FP+d2$TN),y=d2$TP/(d2$TP+d2$FN), n=d2$n, ErrLen=d2$ErrLenT)
B = data.frame(x=d2$FP0/(d2$FP0+d2$TN0),y=as.vector(matrix(1.1,nrow=nrow(d2))), n=d2$n, ErrLen=d2$ErrLenT)
ggplot(data=A, aes(x, y, color=as.factor(ErrLen))) + geom_point(alpha=.991)+
  theme_light()+theme(legend.position = "none")+ #,legend.direction = "horizontal")+
  scale_shape(name=element_blank())+scale_color_brewer(name=element_blank(),palette = "Paired")+
  scale_x_continuous(name="FPR",labels=percent)+
  scale_y_continuous("Recall",labels=percent,breaks = c(0.2,0.4,0.6,0.8,1))+
  geom_text(aes(label=paste(n,ErrLen,sep=", ")),nudge_y = 0.012,size=2.75)+
  coord_cartesian(xlim=c(0.0124,0.01425))
  #ggtitle("small-10-aa with Varying Error Lengths: ROC")
ggsave("Figures/ErrParam_Figures/small10aa_ErrLenNumErr_ROC.pdf", width=4, height=4)


ggplot(data=A, aes(x, y, size=n, shape=ErrLen),color="black") + geom_point(alpha=.99)+
  theme_light()+theme(legend.position = "right",axis.text = element_text(size=15),axis.title = element_text(size=16))+
  scale_shape_manual(name="Error Length",values=c(6,5,4,1,3,2,19,17,18,16,15,14,90))+scale_color_discrete(name="Gene")+
  scale_x_continuous(name="FPR",labels=percent)+
  scale_size_manual(name="Error Frequency",values=sqrt(c(1,4,10,20,40)))+
  scale_y_continuous("Recall",labels=percent)
ggsave("Figures/ErrParam_Figures/small10aa_ErrLenNumErr_ROC_2.pdf", width=6, height=5)

# # ROC for small-10-aa with varying percentages of erroneous sequences and fixed error lengths
# d2=d[d$E=="small-10-aa_NumErrAlns",]
# d2$n=with(d2,as.factor(100/((NumErrSeqDiv!=N)*NumErrSeqDiv+(NumErrSeqDiv==N)*100)))
# ggplot(aes(x=FP/(FP+TN),y=TP/(TP+FN), color=as.factor(n) ),data=summ_roc(d2,n~.))+
#   geom_point(alpha=1)+
#   theme_light()+theme(legend.position = c(.85,.25))+
#   scale_color_brewer(name="n",labels=nlabels, palette="Paired")+
#   scale_x_continuous(name="FPR",labels=percent)+scale_y_continuous("Recall")+
#   ggtitle("small-10-aa with Varying Number of Erroneous Sequences: ROC")
# ggsave("Figures/ErrParam_Figures/small10aa_NumErrAlns_ROC.pdf", width=6, height=6)
# 


# ROC for General 16S.B
options(digits = 5)
d2=summ_roc(d[d$E=="16S.B_General" & d$N > 10 & d$ErrLen==8,], ErrLen+cut(Diameter, breaks = c(0, 0.1, 0.2, 0.5, 0.8, 1), right = F)~.)
A = data.frame(x=d2$FP/(d2$FP+d2$TN),y=d2$TP/(d2$TP+d2$FN), ErrLen=d2$ErrLen, DR=d2$`cut(Diameter, breaks = c(0, 0.1, 0.2, 0.5, 0.8, 1), right = F)`)
B = data.frame(x=d2$FP0/(d2$FP0+d2$TN0),y=as.vector(matrix(1.1,nrow=nrow(d2))), ErrLen=d2$ErrLen, DR=d2$`cut(Diameter, breaks = c(0, 0.1, 0.2, 0.5, 0.8, 1), right = F)`)
ggplot(data=A, aes(x, y, color=as.factor(DR))) + geom_point(alpha=1)+
  theme_light()+theme(legend.position = "right")+geom_point(data=B)+
  scale_color_brewer(name="Diameter",palette = "Paired",labels = function(x) (paste(x)))+
  scale_x_continuous(name="FPR",labels=percent)+
  scale_y_continuous("Recall",labels=percent,breaks = c(0.2,0.4,0.6,0.8,0.9,0.95,1,1.1))+
  ggtitle("16S.B: ROC")
ggsave("Figures/General_Figures/16SB_General_ROC.pdf", width=6, height=6)

# ROC for General Hackett
options(digits = 5)
d2=summ_roc(d[d$E=="Hackett_General" & d$N > 10,], ErrLen+cut(Diameter, breaks = c(0, 0.1, 0.2, 0.5, 0.8, 1), right = F)~.)
A = data.frame(x=d2$FP/(d2$FP+d2$TN),y=d2$TP/(d2$TP+d2$FN), ErrLen=d2$ErrLen, DR=d2$`cut(Diameter, breaks = c(0, 0.1, 0.2, 0.5, 0.8, 1), right = F)`)
B = data.frame(x=d2$FP0/(d2$FP0+d2$TN0),y=as.vector(matrix(1.1,nrow=nrow(d2))), ErrLen=d2$ErrLen, DR=d2$`cut(Diameter, breaks = c(0, 0.1, 0.2, 0.5, 0.8, 1), right = F)`)
ggplot(data=A, aes(x, y, color=as.factor(DR))) + geom_point(alpha=1)+
  theme_light()+theme(legend.position = "right")+geom_point(data=B)+
  scale_color_brewer(name="Diameter",palette = "Paired",labels = function(x) (paste(x)))+
  scale_x_continuous(name="FPR",labels=percent)+
  scale_y_continuous("Recall",labels=percent,breaks = c(0.2,0.4,0.6,0.8,0.9,0.95,1,1.1))+
  ggtitle("Hackett: ROC")
ggsave("Figures/General_Figures/Hackett_General_ROC.pdf", width=6, height=6)

# ROC for General small-10-aa
options(digits = 5)
d2=summ_roc(d[d$E=="small-10-aa-RV100-BBA0039" & d$N > 19,], ErrLen+cut(Diameter, breaks = c(0, 0.1, 0.2, 0.5, 0.8, 1), right = F)~.)
A = data.frame(x=d2$FP/(d2$FP+d2$TN),y=d2$TP/(d2$TP+d2$FN), ErrLen=d2$ErrLen, DR=d2$`cut(Diameter, breaks = c(0, 0.1, 0.2, 0.5, 0.8, 1), right = F)`)
B = data.frame(x=d2$FP0/(d2$FP0+d2$TN0),y=as.vector(matrix(1.1,nrow=nrow(d2))), ErrLen=d2$ErrLen, DR=d2$`cut(Diameter, breaks = c(0, 0.1, 0.2, 0.5, 0.8, 1), right = F)`)
ggplot(data=A, aes(x, y, color=as.factor(DR))) + geom_point(alpha=1)+
  theme_light()+theme(legend.position = "right")+geom_point(data=B)+
  scale_color_brewer(name="Diameter",palette = "Paired",labels = function(x) (paste(x)))+
  scale_x_continuous(name="FPR",labels=percent)+
  scale_y_continuous("Recall",labels=percent,breaks = c(0.2,0.4,0.6,0.8,0.9,0.95,1,1.1))+ggtitle("small-10-aa: ROC")
ggsave("Figures/General_Figures/small10aa_General_ROC.pdf", width=6, height=6)





# Plots figures for union of running the correction algorithm on different k values
d=(read.csv("./CSV_Files/variedUnion.csv", sep=",", header=F))
d$E=paste(d$V1,revalue(d$V2, c("1k"="")),sep="")
names(d) <- c("X", "ks", "DR", "Diameter", "PD", "N", "ErrLen", "Rep", "FP0", "FN0", "TP0", "TN0", "FP", "FN", "TP", "TN","E")

# Recall vs Diameter
ggplot(aes(x=Diameter,y=TP/(TP+FN),color=as.factor(ErrLen)),data=d[d$E=="16S.B_ErrLen" & d$N > 10 & d$ErrLen!=0,])+
  geom_point(alpha=0.1)+theme_classic()+geom_smooth(se=F)+scale_y_continuous("Recall")+
  scale_shape(name="")+scale_color_brewer(palette = "Paired",name="error len", labels = function(x) (paste(x, intToUtf8(215), "11")))+
  ggtitle("16S.B without Union: Recall vs Diameter")
ggsave("Figures/Union_Figures/16SB_NoUnion_Recall.pdf", width=6, height=6)

ggplot(aes(x=Diameter,y=TP/(TP+FN),color=as.factor(ErrLen)),data=d[d$E=="16S.B_UnionErrLen2k" & d$N > 10 & d$ErrLen!=0,])+
  geom_point(alpha=0.1)+theme_classic()+geom_smooth(se=F)+scale_y_continuous("Recall")+
  scale_shape(name="")+scale_color_brewer(palette = "Paired",name="error len", labels = function(x) (paste(x, intToUtf8(215), "11")))+
  ggtitle("16S.B using Union of 2k: Recall vs Diameter")
ggsave("Figures/Union_Figures/16SB_2k_Recall.pdf", width=6, height=6)

ggplot(aes(x=Diameter,y=TP/(TP+FN),color=as.factor(ErrLen)),data=d[d$E=="16S.B_UnionErrLen3k" & d$N > 10 & d$ErrLen!=0,])+
  geom_point(alpha=0.1)+theme_classic()+geom_smooth(se=F)+scale_y_continuous("Recall")+
  scale_shape(name="")+scale_color_brewer(palette = "Paired",name="error len", labels = function(x) (paste(x, intToUtf8(215), "11")))+
  ggtitle("16S.B using Union of 3k: Recall vs Diameter")
ggsave("Figures/Union_Figures/16SB_3k_Recall.pdf", width=6, height=6)

ggplot(aes(x=Diameter,y=TP/(TP+FN),color=as.factor(ErrLen)),data=d[d$E=="16S.B_UnionErrLen4k" & d$N > 10 & d$ErrLen!=0,])+
  geom_point(alpha=0.1)+theme_classic()+geom_smooth(se=F)+scale_y_continuous("Recall")+
  scale_shape(name="")+scale_color_brewer(palette = "Paired",name="error len", labels = function(x) (paste(x, intToUtf8(215), "11")))+
  ggtitle("16S.B using Union of 4k: Recall vs Diameter")
ggsave("Figures/Union_Figures/16SB_4k_Recall.pdf", width=6, height=6)

# False Discovery Rate vs Diameter
ggplot(aes(x=Diameter,y=FP/(TP+FP),color=as.factor(ErrLen)),data=d[d$E=="16S.B_ErrLen" & d$N > 10 & d$ErrLen!=0,])+
  geom_point(alpha=0.1)+theme_classic()+geom_smooth(se=F)+scale_y_continuous("FDR")+
  scale_shape(name="")+scale_color_brewer(palette = "Paired",name="error len", labels = function(x) (paste(x, intToUtf8(215), "11")))+
  ggtitle("16S.B without Union: FDR vs Diameter")
ggsave("Figures/Union_Figures/16SB_NoUnion_FDR.pdf", width=6, height=6)

ggplot(aes(x=Diameter,y=FP/(TP+FP),color=as.factor(ErrLen)),data=d[d$E=="16S.B_UnionErrLen2k" & d$N > 10 & d$ErrLen!=0,])+
  geom_point(alpha=0.1)+theme_classic()+geom_smooth(se=F)+scale_y_continuous("FDR")+
  scale_shape(name="")+scale_color_brewer(palette = "Paired",name="error len", labels = function(x) (paste(x, intToUtf8(215), "11")))+
  ggtitle("16S.B using Union of 2k: FDR vs Diameter")
ggsave("Figures/Union_Figures/16SB_2k_FDR.pdf", width=6, height=6)

ggplot(aes(x=Diameter,y=FP/(TP+FP),color=as.factor(ErrLen)),data=d[d$E=="16S.B_UnionErrLen3k" & d$N > 10 & d$ErrLen!=0,])+
  geom_point(alpha=0.1)+theme_classic()+geom_smooth(se=F)+scale_y_continuous("FDR")+
  scale_shape(name="")+scale_color_brewer(palette = "Paired",name="error len", labels = function(x) (paste(x, intToUtf8(215), "11")))+
  ggtitle("16S.B using Union of 3k: FDR vs Diameter")
ggsave("Figures/Union_Figures/16SB_3k_FDR.pdf", width=6, height=6)

ggplot(aes(x=Diameter,y=FP/(TP+FP),color=as.factor(ErrLen)),data=d[d$E=="16S.B_UnionErrLen4k" & d$N > 10 & d$ErrLen!=0,])+
  geom_point(alpha=0.1)+theme_classic()+geom_smooth(se=F)+scale_y_continuous("FDR")+
  scale_shape(name="")+scale_color_brewer(palette = "Paired",name="error len", labels = function(x) (paste(x, intToUtf8(215), "11")))+
  ggtitle("16S.B using Union of 4k: FDR vs Diameter")
ggsave("Figures/Union_Figures/16SB_4k_FDR.pdf", width=6, height=6)


ggplot(data=d[d$N > 19 &d$ks %in% c("3k") ,])+
  geom_point(aes(x=Diameter,y=FP/(TN+FP),linetype="After error",color="After error"),alpha=0.4,size=.6)+
  geom_smooth(aes(x=Diameter,y=FP/(TN+FP),linetype="After error",color="After error"),se=F)+
  geom_point(aes(x=Diameter,y=FP0/(TN0+FP0),linetype="Before error",color="Before error"),alpha=0.4,size=.6)+
  geom_smooth(aes(x=Diameter,y=FP0/(TN0+FP0),linetype="Before error",color="Before error"),se=F)+
  scale_shape(name="",solid = T)+scale_color_brewer(palette = "Dark2",name="")+
  facet_wrap(~ErrLen,nrow=2,labeller = function(x) {list(ErrLen=paste(x$ErrLen, intToUtf8(215), "11"))})+
  theme_classic()+theme(legend.position=c(.89,.15))+
  scale_linetype_manual(name="",values=c(1,1))+
  scale_x_continuous(breaks=c(1/4,1/2,3/4,1))+scale_y_log10("FPR")#+coord_cartesian(ylim=c(0,0.003))
ggsave("Figures/Union_Figures/16SB_before_FPR.pdf", width=8, height=4)  

ggplot(aes(x=Diameter,y=FP/(TN+FP),color=as.factor(ks)),data=d[d$N > 19 ,])+
  geom_point(alpha=0.3,size=0.7)+
  geom_smooth(se=F)+
  facet_wrap(~ErrLen,nrow=2,labeller = function(x) {list(ErrLen=paste(x$ErrLen, intToUtf8(215), "11"))})+
  scale_shape(name="")+scale_color_brewer(palette = "Dark2",name="", labels = function(x) (paste(x, "setting")))+
  theme_classic()+theme(legend.position = c(.89,.15))+
  scale_x_continuous(breaks=c(1/4,1/2,3/4,1))+scale_y_log10("FPR")
ggsave("Figures/Union_Figures/16SB_allk_FPR.pdf", width=8, height=4)  


ggplot(aes(x=Diameter,y=FP/(TP+FP),color=as.factor(ks)),data=d[d$N > 19 ,])+
  geom_point(alpha=0.3,size=0.7)+
  geom_smooth(se=F)+
  facet_wrap(~ErrLen,nrow=2,labeller = function(x) {list(ErrLen=paste(x$ErrLen, intToUtf8(215), "11"))})+
  scale_shape(name="")+scale_color_brewer(palette = "Dark2",name="", labels = function(x) (paste(x, "setting")))+
  theme_classic()+theme(legend.position = c(.89,.15))+
  scale_x_continuous(breaks=c(1/4,1/2,3/4,1))+scale_y_continuous("FDR")
ggsave("Figures/Union_Figures/16SB_allk_FDR.pdf", width=8, height=4)  

ggplot(aes(x=Diameter,y=TP/(TP+FN),color=as.factor(ks)),data=d[d$N > 19 ,])+
  geom_point(alpha=0.3,size=0.7)+
  geom_smooth(se=F)+
  facet_wrap(~ErrLen,nrow=2,labeller = function(x) {list(ErrLen=paste(x$ErrLen, intToUtf8(215), "11"))})+
  scale_shape(name="")+scale_color_brewer(palette = "Dark2",name="", labels = function(x) (paste(x, "setting")))+
  theme_classic()+theme(legend.position = c(.89,.15))+
  scale_x_continuous(breaks=c(1/4,1/2,3/4,1))+scale_y_continuous("Recall")
ggsave("Figures/Union_Figures/16SB_allk_Recall.pdf", width=8, height=4)  

# ROC all

options(digits = 2)
d2=summ_roc(d[d$N > 19,], ks+ErrLen+cut(Diameter, breaks = c(0, 0.1, 0.2, 0.5, 0.8, 1), right = F)~.)
A = data.frame(x=d2$FP/(d2$FP+d2$TN),y=d2$TP/(d2$TP+d2$FN), ErrLen=d2$ErrLen, ks=d2$ks, DR=d2$`cut(Diameter, breaks = c(0, 0.1, 0.2, 0.5, 0.8, 1), right = F)`)
B = data.frame(x=d2$FP0/(d2$FP0+d2$TN0),y=as.vector(matrix(1.1,nrow=nrow(d2))), ErrLen=d2$ErrLen, ks=d2$ks, DR=d2$`cut(Diameter, breaks = c(0, 0.1, 0.2, 0.5, 0.8, 1), right = F)`)
ggplot(data=A, aes(x, y, color=as.factor(ks), shape=as.factor(DR))) + geom_point(alpha=0.99)+geom_line(aes(group=ks))+
  theme_bw()+theme(legend.position = c(.87,.15),legend.box.just = "top",legend.box = "horizontal")+
  #geom_point(data=B)+
  #geom_line(aes(group=ks),data=B,linetype=2)+
  scale_shape(name="Diameter")+scale_color_brewer(name="Setting",palette = "Dark2",label = function(x) substr(gsub(pattern = ""," ",x),2,4) )+
  scale_x_continuous(name="FPR",labels=percent)+
  facet_wrap(~ErrLen,nrow=2,scales="free_y", labeller = function(x) {list(ErrLen=paste(x$ErrLen, intToUtf8(215), "11"))})+
  scale_y_continuous("Recall",labels=percent)#coord_cartesian(xlim=c(0, 0.0015), ylim=c(0,1))
  #ggtitle("16S.B: ROC")
ggsave("Figures/Union_Figures/16SB_allKs_ROC.pdf", width=10, height=5)


# ROC
options(digits = 2)
d2=summ_roc(d[d$E=="16S.B_ErrLen" & d$N > 19,], ErrLen+cut(Diameter, breaks = c(0, 0.1, 0.2, 0.5, 0.8, 1), right = F)~.)
A = data.frame(x=d2$FP/(d2$FP+d2$TN),y=d2$TP/(d2$TP+d2$FN), ErrLen=d2$ErrLen, DR=d2$`cut(Diameter, breaks = c(0, 0.1, 0.2, 0.5, 0.8, 1), right = F)`)
B = data.frame(x=d2$FP0/(d2$FP0+d2$TN0),y=as.vector(matrix(1.1,nrow=nrow(d2))), ErrLen=d2$ErrLen, DR=d2$`cut(Diameter, breaks = c(0, 0.1, 0.2, 0.5, 0.8, 1), right = F)`)
ggplot(data=A, aes(x, y, color=as.factor(ErrLen), shape=as.factor(DR))) + geom_point(alpha=1)+
  theme_light()+theme(legend.position = "right")+geom_point(data=B)+
  scale_shape(name="Diameter")+scale_color_brewer(name="Error Length",palette = "Paired",labels = function(x) (paste(x, intToUtf8(215), "11")))+
  scale_x_continuous(name="FPR",labels=percent)+
  scale_y_continuous("Recall",labels=percent,breaks = c(0.2,0.4,0.6,0.8,1, 1.1))+#coord_cartesian(xlim=c(0, 0.0015), ylim=c(0,1))
  ggtitle("16S.B without Union: ROC")
ggsave("Figures/Union_Figures/16SB_NoUnion_ROC.pdf", width=6, height=6)

options(digits = 2)
d2=summ_roc(d[d$E=="16S.B_UnionErrLen2k" & d$N > 19,], ErrLen+cut(Diameter, breaks = c(0, 0.1, 0.2, 0.5, 0.8, 1), right = F)~.)
A = data.frame(x=d2$FP/(d2$FP+d2$TN),y=d2$TP/(d2$TP+d2$FN), ErrLen=d2$ErrLen, DR=d2$`cut(Diameter, breaks = c(0, 0.1, 0.2, 0.5, 0.8, 1), right = F)`)
B = data.frame(x=d2$FP0/(d2$FP0+d2$TN0),y=as.vector(matrix(1.1,nrow=nrow(d2))), ErrLen=d2$ErrLen, DR=d2$`cut(Diameter, breaks = c(0, 0.1, 0.2, 0.5, 0.8, 1), right = F)`)
ggplot(data=A, aes(x, y, color=as.factor(ErrLen), shape=as.factor(DR))) + geom_point(alpha=1)+
  theme_light()+theme(legend.position = "right")+geom_point(data=B)+
  scale_shape(name="Diameter")+scale_color_brewer(name="Error Length",palette = "Paired",labels = function(x) (paste(x, intToUtf8(215), "11")))+
  scale_x_continuous(name="FPR",labels=percent)+
  scale_y_continuous("Recall",labels=percent,breaks = c(0.2,0.4,0.6,0.8,1, 1.1))+#coord_cartesian(xlim=c(0, 0.0015), ylim=c(0,1))
  ggtitle("16S.B using Union of 2k: ROC")
ggsave("Figures/Union_Figures/16SB_2k_ROC.pdf", width=6, height=6)

options(digits = 2)
d2=summ_roc(d[d$E=="16S.B_UnionErrLen3k" & d$N > 19,], ErrLen+cut(Diameter, breaks = c(0, 0.1, 0.2, 0.5, 0.8, 1), right = F)~.)
A = data.frame(x=d2$FP/(d2$FP+d2$TN),y=d2$TP/(d2$TP+d2$FN), ErrLen=d2$ErrLen, DR=d2$`cut(Diameter, breaks = c(0, 0.1, 0.2, 0.5, 0.8, 1), right = F)`)
B = data.frame(x=d2$FP0/(d2$FP0+d2$TN0),y=as.vector(matrix(1.1,nrow=nrow(d2))), ErrLen=d2$ErrLen, DR=d2$`cut(Diameter, breaks = c(0, 0.1, 0.2, 0.5, 0.8, 1), right = F)`)
ggplot(data=A, aes(x, y, color=as.factor(ErrLen), shape=as.factor(DR))) + geom_point(alpha=1)+
  theme_light()+theme(legend.position = "right")+geom_point(data=B)+
  scale_shape(name="Diameter")+scale_color_brewer(name="Error Length",palette = "Paired",labels = function(x) (paste(x, intToUtf8(215), "11")))+
  scale_x_continuous(name="FPR",labels=percent)+
  scale_y_continuous("Recall",labels=percent,breaks = c(0.2,0.4,0.6,0.8,1, 1.1))+#coord_cartesian(xlim=c(0, 0.0015), ylim=c(0,1))
  ggtitle("16S.B using Union of 3k: ROC")
ggsave("Figures/Union_Figures/16SB_3k_ROC.pdf", width=6, height=6)

options(digits = 2)
d2=summ_roc(d[d$E=="16S.B_UnionErrLen4k" & d$N > 19,], ErrLen+cut(Diameter, breaks = c(0, 0.1, 0.2, 0.5, 0.8, 1), right = F)~.)
A = data.frame(x=d2$FP/(d2$FP+d2$TN),y=d2$TP/(d2$TP+d2$FN), ErrLen=d2$ErrLen, DR=d2$`cut(Diameter, breaks = c(0, 0.1, 0.2, 0.5, 0.8, 1), right = F)`)
B = data.frame(x=d2$FP0/(d2$FP0+d2$TN0),y=as.vector(matrix(1.1,nrow=nrow(d2))), ErrLen=d2$ErrLen, DR=d2$`cut(Diameter, breaks = c(0, 0.1, 0.2, 0.5, 0.8, 1), right = F)`)
ggplot(data=A, aes(x, y, color=as.factor(ErrLen), shape=as.factor(DR))) + geom_point(alpha=1)+
  theme_light()+theme(legend.position = "right")+geom_point(data=B)+
  scale_shape(name="Diameter")+scale_color_brewer(name="Error Length",palette = "Paired",labels = function(x) (paste(x, intToUtf8(215), "11")))+
  scale_x_continuous(name="FPR",labels=percent)+
  scale_y_continuous("Recall",labels=percent,breaks = c(0.2,0.4,0.6,0.8,1, 1.1))+#coord_cartesian(xlim=c(0, 0.0015), ylim=c(0,1))
  ggtitle("16S.B using Union of 4k: ROC")
ggsave("Figures/Union_Figures/16SB_4k_ROC.pdf", width=6, height=6)






d=(read.csv("./CSV_Files/res_16SK.csv", sep=",", header=F))
names(d) <- c("E", "DR", "X", "Diameter", "PD", "N", "ErrLen", "K", "Rep", "FP0", "FN0", "TP0", "TN0", "FP", "FN", "TP", "TN")

d = d[d$N > 19, ]
nlabels = c("1","2%","5%","10%","20%")
#levels(d$DR)[levels(d$DR)=="concatenation"] <- "concat"
d$n="~5%"

d$ErrLenT = paste(d$ErrLen, intToUtf8(215), ifelse(grepl("small",d$E),"11","11"),sep="")
d$ErrLenT = factor(d$ErrLenT,levels=c("2×11","3×11", "4×11","8×11","16×11", "32×11", "64×11" ))

d$SL = with(d,(TN+FP+FN+TP)/as.numeric(as.character(N)))

d = d[d$E == "16S.B_K_p0.1_q0.5",]

options(digits = 2)
d2=summ_roc(d[d$N > 19 ,], K+E+ErrLen+cut(Diameter, breaks = c(0, 0.1, 0.2, 0.5, 0.8, 1), right = F)~.)
A = data.frame(x=d2$FP/(d2$FP+d2$TN),y=d2$TP/(d2$TP+d2$FN), E=d2$E, ErrLen=d2$ErrLen, K=d2$K, DR=d2$`cut(Diameter, breaks = c(0, 0.1, 0.2, 0.5, 0.8, 1), right = F)`)
#B = data.frame(x=d2$FP0/(d2$FP0+d2$TN0),y=as.vector(matrix(1.1,nrow=nrow(d2))), ErrLen=d2$ErrLen, K=d2$K, DR=d2$`cut(Diameter, breaks = c(0, 0.1, 0.2, 0.5, 0.8, 1), right = F)`)
ggplot(data=A[A$ErrLen <= 32,], aes(x, y, color=as.factor(K), shape=as.factor(DR))) + 
  geom_point(alpha=0.99)+
  geom_line(aes(group=K))+
  theme_bw()+theme(legend.position = c(.87,.15),legend.box.just = "top",legend.box = "horizontal")+
  #geom_point(data=B)+
  #geom_line(aes(group=K),data=B,linetype=2)+
  scale_shape(name="Diameter")+scale_color_brewer(name="Setting",palette = "Spectral",label = function(x) substr(gsub(pattern = ""," ",x),2,4) )+
  scale_x_continuous(name="FPR",labels=percent)+
  facet_wrap(~ErrLen, labeller = function(x) {list(ErrLen=paste(x$ErrLen, intToUtf8(215), "11"))})+
  scale_y_continuous("Recall",labels=percent)#coord_cartesian(xlim=c(0, 0.0015), ylim=c(0,1))
#ggtitle("16S.B: ROC")
ggsave("Figures/ErrParam_Figures/16SB_K_byerr_faceted.pdf", width=10, height=6)


ggplot(data=A[A$ErrLen <= 64,], aes(x, y, color=as.factor(ErrLen), shape=as.factor(DR))) + 
  geom_point(alpha=0.99)+
  geom_line(aes(group=ErrLen))+
  theme_bw()+theme(legend.position = c(.9,.2))+
  #geom_point(data=B)+
  #geom_line(aes(group=K),data=B,linetype=2)+
  scale_shape(name="Diameter")+
  scale_color_brewer(name="Error Len",palette = "Paired" ,labels = function(x) {paste(x, intToUtf8(215), "11")} )+
  scale_x_continuous(name="FPR",labels=percent,breaks=c(0.0001,0.0004,0.0007))+
  facet_wrap(~K, nrow=2,labeller = label_both )+
  scale_y_continuous("Recall",labels=percent)+#coord_cartesian(xlim=c(0, 0.0015), ylim=c(0,1))
ggsave("Figures/ErrParam_Figures/16SB_K_byK_faceted.pdf", width=10, height=8)


