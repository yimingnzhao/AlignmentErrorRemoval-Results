require(ggplot2); require(scales); require(data.table); require(reshape2)

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

dc = d[d$E %in% c("16S.B-1-cutoffs" , "16S.B_General"),]

d= d[d$E!="16S.B-1-cutoffs",]

t = read.csv('CSV_Files/res_other_methods.csv',header=F)
names(t)=c("AlignmentName" ,"DiameterRange" ,"X" ,"Diameter", "PD", "N", "ErrLen", "NumErrSeqDiv", "Rep", "real_time", "user_time", "sys_time", "FP", "FN", "TP", "TN")
#t$time=as.numeric(sub("m.*","",t$real_time))*60+as.numeric(sub("s$","",sub(".*m","",t$real_time)))
t = t[t$N > 19 | !grepl("16S",t$AlignmentName),]
t$time=as.numeric(sub("user$","",sub("s$","",sub(".*m","",t$real_time)))) + as.numeric(ifelse(grepl("user",t$real_time), 0, sub("[m].*","",t$real_time)))*60
unique(t$AlignmentName)

md= merge(t[t$AlignmentName=="16S.B-Divvier",],rbind(d[d$ErrLenT=="~50",],dc),
           by.y = c("DR","X","Diameter","PD","N","NumErrSeqDiv","Rep"),by.x=c("DiameterRange","X","Diameter","PD","N","NumErrSeqDiv","Rep"))

head(md)
t2=md[,c(18,1:5,9,6:7,10:12,24:27,17,20:23)]
t3=md[,c(8,1:5,9,6:7,10:12,13:16,17,20:23)]
names(t2) = names(t3) = c(names(t),c("FP0","FN0","TP0","TN0"))
t4=rbind(t2,t3)
unique(t4$AlignmentName)
t4$AlignmentName = factor(t4$AlignmentName,labels = c("16S.B","16S.B-1-cutoffs", "16S.B-Divvier"))
t4$DR2 = cut(t4$Diameter, breaks = c(0, 0.1, 0.2, 0.5, 0.8, 1), right = F)
head(t4)

t5 = t[t$AlignmentName %in% c("small-10-aa-RV100-BBA0039-DivA","small-10-aa-RV100-BBA0039-Divvier"),]
t5 = t5[,c(1,13:16,13:16,6,7)]
names(t5)=names(d)[c(1,10:18,20)]
t5$n="5%"
t5$ErrLenT="8×11"
head(t5)

t$DR2 = cut(t$Diameter, breaks = c(0, 0.1, 0.2, 0.5, 0.8, 1), right = F)


tn = read.csv('CSV_Files/res_taper_multi_0.1.7.csv',header=F)
names(tn)=c("AlignmentName" ,"DiameterRange" ,"X" ,"Diameter", "PD", "N", "ErrLen", "NumErrSeqDiv", "Rep", "FP0","FN0","TP0","TN0","FP", "FN", "TP", "TN")
head(tn)
tdm = merge(tn,d[d$ErrLenT=="~50",],by.y = c("DR","X","Diameter","PD","N","NumErrSeqDiv","Rep"),by.x=c("DiameterRange","X","Diameter","PD","N","NumErrSeqDiv","Rep"))
tdm1 = tdm[,1:17]; names(tdm1)[9:17]=c("ErrLen", "FP0","FN0","TP0","TN0","FP", "FN", "TP", "TN")
tdm2 = tdm[,c(1:7,18:27)]; names(tdm2)[8:17]=c("AlignmentName", "ErrLen", "FP0","FN0","TP0","TN0","FP", "FN", "TP", "TN")
tdm=rbind(tdm1,tdm2)
tdm$DR2 = cut(tdm$Diameter, breaks = c(0, 0.1, 0.2, 0.5, 0.8, 1), right = F)






bs = read.csv('CSV_Files/res_treecmp.csv',head=F)
names(bs) <- c("AlignmentName","DiameterRangeGene","X","Diameter","PD","N","ErrLen",
               "NumErrSeqDiv","Rep","ms_err","pd_err","rf_err","rfw_err","ms_res","pd_res","rf_res","rfw_res")

bs2 = cbind(melt(bs[c(1:9,10:13)],id.vars = 1:9), melt(bs[c(1:9,14:17)],id.vars = 1:9)[,10:11])
names(bs2)[c(12,13)]=c("v","after")
names(bs2)
head(bs)
bs$DR2 = cut(bs$Diameter, breaks = c(0, 0.1, 0.2, 0.5, 0.8, 1), right = F)
bs$rf_errN = bs$rf_err/(2*bs$N-6)
levels(bs$AlignmentName)

summary(with(data=bs[bs$AlignmentName%in%c("16S.B"),],(rf_err-rf_res))>0)
summary(with(data=bs[bs$AlignmentName%in%c("16S.B"),],(rf_err-rf_res))==0)
summary(with(data=bs[bs$AlignmentName%in%c("16S.B") & bs$rf_err!=0 ,],(rf_err-rf_res)/rf_err))
summary(with(data=bs[bs$AlignmentName%in%c("16S.B") & bs$rf_err!=0 ,],(rf_err-rf_res)/(2*N-6)))

summary(with(data=bs[bs$AlignmentName%in%c("16S.B"),],(rfw_err-rfw_res))>0)
summary(with(data=bs[bs$AlignmentName%in%c("16S.B"),],(rfw_err-rfw_res))==0)
summary(with(data=bs[bs$AlignmentName%in%c("16S.B") & bs$rfw_err!=0 ,],(rfw_err-rfw_res)/rfw_err))

with(data=bs[bs$AlignmentName%in%c("16S.B") & bs$rfw_err!=0 ,],cor.test((rfw_err-rfw_res)/rfw_err,Diameter))

summary(with(data=bs[bs$AlignmentName%in%c("16S.B-Divvier"),],(rfw_err-rfw_res))==0)

temp = merge(bs[bs$AlignmentName%in%c("16S.B"),],bs[bs$AlignmentName%in%c("16S.B-Divvier"),],
             by.y = c("DiameterRangeGene","X","Diameter","PD","N","NumErrSeqDiv","ms_err","pd_err","rf_err","rfw_err", "ErrLen"),
             by.x=c("DiameterRangeGene","X","Diameter","PD","N","NumErrSeqDiv","ms_err","pd_err","rf_err","rfw_err", "ErrLen"))[,c("AlignmentName.x","DiameterRangeGene","X","Diameter","PD","N","ErrLen","NumErrSeqDiv","Rep.x","ms_err","pd_err","rf_err","rfw_err","ms_res.x","pd_res.x","rf_res.x","rfw_res.x","DR2.x","rf_errN.x")]
names(temp) = names(bs)
temp$AlignmentName="16S.B-both"


bh = cbind(dcast(DiameterRangeGene+N+Diameter~.,data=bs[bs$AlignmentName=="HackettGenes",],value.var = "rf_res",fun.aggregate = mean),
           dcast(DiameterRangeGene~.,data=bs[bs$AlignmentName=="HackettGenes",],value.var = "rf_err",fun.aggregate = mean)[,2])
names(bh) = c("DiameterRangeGene","N","Diameter", "rf_res","rf_err")
bh$m="RF"

bh2 = cbind(dcast(DiameterRangeGene+N+Diameter~.,data=bs[bs$AlignmentName=="HackettGenes",],value.var = "rfw_res",fun.aggregate = mean),
            dcast(DiameterRangeGene~.,data=bs[bs$AlignmentName=="HackettGenes",],value.var = "rfw_err",fun.aggregate = mean)[,2])
names(bh2) = c("DiameterRangeGene","N","Diameter", "rf_res","rf_err")
bh2$m="WRF"

#write.csv(tdm,"oldnew.csv")
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

d2=summ_roc(tdm, AlignmentName+DR2~.)
A = data.frame(x=d2$FP/(d2$FP+d2$TN),y=d2$TP/(d2$TP+d2$FN), DR=d2$DR2,E=d2$AlignmentName)
ggplot(data=A, aes(x, y, shape=E,color=as.factor(DR))) + 
  geom_point(alpha=1)+
  geom_path(aes(group=DR),data=A[A$n!="~5%",],linetype=1)+
  theme_bw()+theme(legend.position = "right",legend.text.align = 1)+
  scale_shape_manual(name="Length, Freq",values=c(15,17,1,2,5,8,9,7,6,19,18,3,90),labels=c("new","old"))+
  scale_color_brewer(name="Diameter",palette = "Dark2")+
  scale_x_continuous(name="FPR",labels=percent)+
  scale_y_continuous("Recall",labels=percent)

#ggsave("old-new.pdf")

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
  theme_classic()+theme(legend.position = c(0.7,0.2),legend.text.align = 1,legend.direction = "horizontal",
                        plot.tag.position = c(0.015, 0.975),plot.tag = element_text(size=15,face = "bold"))+labs(tag = "a)")+
  scale_shape_manual(name="Err Len, Freq",values=c(15,17,1,2,5,8,9,7,6,19,18,3))+
  scale_color_brewer(name="Diameter",palette = "Dark2")+
  scale_x_continuous(name="FPR",labels=percent)+
  scale_y_continuous("Recall",labels=percent)+
  facet_wrap(~(E=="16S.B_NumErrAlns"),labeller = function(x) list(c("Changing Error Length", "Changing Error Frequency")))+
  geom_linerange(aes(x=x,ymin=0.995,ymax=1.005,color=as.factor(DR)),data=B,linetype=1,size=1)+
  #annotate(geom="text",label="a)",x=-0.00033,y=1.024,size=5)+coord_cartesian(clip = 'off',xlim=c(0,0.0016),ylim=c(0.7,1))+
ggsave("Figures/ErrParam_Figures/16SB_ErrLenNumErr_ROC_faceted.pdf", width=10, height=4.6)



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
  theme(legend.position = c(.27,.061),legend.direction = "horizontal", legend.text.align = 1, 
        plot.tag.position = c(0.01, 0.975),plot.tag = element_text(size=15,face = "bold"))+labs(tag = "b)")+
  scale_y_log10("Percent error",labels=percent)+
  scale_shape(name="")+scale_x_discrete(name="Error Len, Freq")+
  facet_wrap(~E,labeller = function(x) list(E=c("Changing Error Length","Changing Error Frequency")),scales="free_x",shrink = T)+
  scale_color_brewer(palette = "Dark2",name="Diameter")+#+
  #coord_cartesian(clip = 'off',ylim=c(0.000003,0.02),xlim=c(1,7))+annotate(geom="text",label="c)",x=-1,y=0.035,size=5)+
ggsave("Figures/ErrParam_Figures/16S.B_ErrLenNumErr_percenterror_arrow_log.pdf",width = 10,height =4.6)


ggplot(aes(color=AlignmentName,yend=FN/(FN+TN),y=(FN+TP)/(FN+FP+TP+TN),x=DR2,xend=DR2),
       data=data.table::dcast(setDT(t4[t4$AlignmentName %in% c("16S.B-Divvier", "16S.B","16S.B-1-cutoffs") &t4$N>19,]),
                              DR2+AlignmentName~.,fun.aggregate = mean,value.var=c("FP","TP","TN","FN")))+
  #geom_boxplot(outlier.alpha = .5, outlier.size = 0.4)+#geom_point(alpha=0.5,size=1)+
  geom_segment(position = position_dodge2(width=0.8),size=0.8,arrow = arrow(length=unit(0.2,"cm")))+
  #stat_summary(position = position_dodge(width=0.3),geom="line")+
  theme_classic()+
  theme(legend.position = c(.82,.1),axis.text.x = element_text(size=7.5),
        plot.tag.position = c(0.015, 0.975),plot.tag = element_text(size=15,face = "bold"))+labs(tag = "d)")+
  scale_y_log10("Percent error",labels=percent)+
  scale_x_discrete(name="Diameter")+
  facet_wrap(~"Changing methods")+
  #scale_x_discrete(name="Error Len, Freq")+
  #facet_wrap(~E,labeller = function(x) list(E=c("Changing Error Length","Changing Error Frequency")),scales="free_x",shrink = T)+
  scale_color_brewer(name="",palette = "Set2" ,labels=c("TAPER","TAPER (no 2D)","Divvier"))+
  ggsave("Figures/ErrParam_Figures/16SB_methods_arrow.pdf", width=2.5*1.2, height=4.8*1.2)

ggplot(aes(color=AlignmentName,yend=FN/(FN+TN),y=(FN+TP)/(FN+FP+TP+TN),x=DR2,xend=DR2),
       data=data.table::dcast(setDT(t[t$AlignmentName %in% c("16S.B-Divvier", "16S.B") &t$N>19,]),
                              DR2+AlignmentName~.,fun.aggregate = mean,value.var=c("FP","TP","TN","FN")))+
  #geom_boxplot(outlier.alpha = .5, outlier.size = 0.4)+#geom_point(alpha=0.5,size=1)+
  geom_segment(position = position_dodge2(width=0.8),size=0.8,arrow = arrow(length=unit(0.2,"cm")))+
  #stat_summary(position = position_dodge(width=0.3),geom="line")+
  theme_classic()+
  theme(legend.position = c(.82,.1),axis.text.x = element_text(size=7.5),
        plot.tag.position = c(0.015, 0.975),plot.tag = element_text(size=15,face = "bold"))+labs(tag = "d)")+
  scale_y_log10("Percent error",labels=percent)+
  scale_x_discrete(name="Diameter")+
  facet_wrap(~"Changing methods")+
  #scale_x_discrete(name="Error Len, Freq")+
  #facet_wrap(~E,labeller = function(x) list(E=c("Changing Error Length","Changing Error Frequency")),scales="free_x",shrink = T)+
  scale_color_brewer(name="",palette = "Set2" ,labels=c("TAPER","Divvier"))

options(digits = 2)
d2=summ_roc(t4[t4$N>19,], DR2+AlignmentName~.)
A = data.frame(x=d2$FP/(d2$FP+d2$TN),y=d2$TP/(d2$TP+d2$FN), DR2=d2$DR2,  AlignmentName=d2$AlignmentName)
ggplot(data=A[A$AlignmentName!="16S.B-1-cutoffs",], aes(x, y, shape=AlignmentName,color=DR2)) + 
  geom_point(alpha=0.99,size=2.5)+
  #geom_path(aes(group=interaction(AlignmentName)))+
  theme_classic()+theme(legend.position = c(0.75,0.32),legend.text.align = 1,axis.title.y = element_text(vjust=-4),
                        plot.tag.position = c(0.075, 0.975),plot.tag = element_text(size=15,face = "bold"))+labs(tag = "b)")+
  scale_shape_manual(name="",values=c(15,17,16,18,3,90),labels=c("TAPER","TAPER-n", "Divvier"))+
  scale_color_brewer(name="",palette = "Dark2")+
  #scale_size_discrete(name="Diameter")+
  scale_x_continuous(name="FPR")+
  scale_y_continuous("Recall",labels=percent)+
  #annotate(geom="text",label="b)",x=-0.008,y=1.02,size=5)+
  coord_cartesian(xlim=c(0,0.03),ylim=c(0.7,1))+
  facet_wrap(~"Changing methods")+
  #geom_line(aes(group=DiameterRange),linetype=3,size=0.3)+
  #geom_linerange(aes(x=x,ymin=0.995,ymax=1.005,color=as.factor(DiameterRange)),data=B,linetype=1,size=1)
ggsave("Figures/ErrParam_Figures/16SB_methods-narrow.pdf", width=2.7*1.1, height=4.8*1.1)

ins=qplot(AlignmentName,time/60,data=t[!is.na(t$DR2),c("AlignmentName","DR2","time")],geom="violin")+
  scale_y_continuous(trans="log10",name="time (mins)")+
  theme_classic()+theme(axis.title.y = element_text(vjust=-3))+
  scale_x_discrete(labels=c("TAPER","Divvier"),name="")

sort(t[!is.na(t$DR2),c("time")])

#merge(A[,c(3,4,1,2)],melt(dcast(DR2~AlignmentName,data=t[!is.na(t$DR2),c("AlignmentName","DR2","time")],fun.aggregate = mean)),by.x=1:2,by.y=1:2)
ggplot(data=A[A$AlignmentName!="16S.B-1-cutoffs",], 
       aes(x, y, shape=AlignmentName,color=DR2)) + 
  geom_point(alpha=0.99,size=2.5)+
  #geom_path(aes(group=interaction(AlignmentName)))+
  theme_classic()+
  theme(legend.position = c(0.67,0.93),legend.direction = "horizontal",
          plot.tag.position = c(0.015, 0.975),plot.tag = element_text(size=16,face = "bold"))+labs(tag = "c)")+
  scale_shape_manual(name="",values=c(15,17,16,18,3,90),labels=c("TAPER","Divvier"))+
  scale_color_brewer(name="",palette = "Dark2", guide="none")+
  #scale_size_discrete(name="Diameter")+
  scale_x_continuous(name="FPR",labels=percent)+
  geom_text(aes(label=DR2),nudge_y = -0.0045,nudge_x=0.0002,size=2.1)+
  #geom_text(aes(label=paste(round(value/60,1),"m")),nudge_y = -0.0035,nudge_x=0.0004,size=2.3)+
  scale_y_continuous("Recall",labels=percent)+
  annotation_custom(ggplotGrob(ins), xmin = 0.008, xmax = 0.0188, ymin = 0.83, ymax = 0.898)+
ggsave("Figures/ErrParam_Figures/16SB_methods.pdf",width = 3.75,height = 3.8)


ggplot(data=A[A$AlignmentName!="16S.B-Divvier",], 
       aes(x, y, shape=AlignmentName,color=DR2)) + 
  geom_point(alpha=0.99,size=2.5)+
  #geom_path(aes(group=interaction(AlignmentName)))+
  theme_classic()+
  theme(legend.position = c(0.73,0.93),legend.direction = "vertical")+
  scale_shape_manual(name="",values=c(15,17,16,18,3,90),labels=c("Default","No Step 4 (p=q=c=1)"))+
  scale_color_brewer(name="",palette = "Dark2", guide="none")+
  #scale_size_discrete(name="Diameter")+
  scale_x_continuous(name="FPR",labels=percent)+
  geom_text(aes(label=DR2),nudge_y = -0.0045,nudge_x=0.00002,size=2.5,data=A[A$AlignmentName=="16S.B-1-cutoffs",])+
  #geom_text(aes(label=paste(round(value/60,1),"m")),nudge_y = -0.0035,nudge_x=0.0004,size=2.3)+
  scale_y_continuous("Recall",labels=percent)+
  ggsave("Figures/ErrParam_Figures/16SB_TAPER_Vars.pdf",width = 3.75,height = 3.8)

#+ ggsave("Figures/ErrParam_Figures/16SB_methods_mins.pdf",width = 1.7,height = 3.8)

 ggplot(data=A, aes(x, y, shape=DiameterRange,color=AlignmentName)) + 
  geom_point(alpha=1)+
  #geom_path(aes(group=interaction(AlignmentName)))+
  theme_bw()+theme(legend.position = "right",legend.text.align = 1)+
  scale_shape_manual(name="Length, Freq",values=c(15,17,1,2,5,8,9,7,6,19,18,3,90))+
  scale_color_brewer(name="Method",palette = "Dark2")+
  scale_x_continuous(name="FPR",labels=percent)+
  scale_y_continuous("Recall",labels=percent)
  #geom_line(aes(group=AlignmentName,linetype=AlignmentName),color="1")
  #geom_linerange(aes(x=x,ymin=0.995,ymax=1.005,color=as.factor(DiameterRange)),data=B,linetype=1,size=1)

ggplot(data=t4, aes(FP/(FP+TN), TP/(TP+FN), shape=DiameterRange,color=AlignmentName)) + 
  geom_point(alpha=1)+
  #geom_path(aes(group=interaction(AlignmentName)))+
  theme_bw()+theme(legend.position = "right",legend.text.align = 1)+
  scale_shape_manual(name="Length, Freq",values=c(15,17,1,2,5,8,9,7,6,19,18,3,90))+
  scale_color_brewer(name="Diameter",palette = "Dark2")+
  scale_x_continuous(name="FPR",labels=percent,trans = "log2")+
  scale_y_continuous("Recall",labels=percent)

ggplot(aes(color=reorder(AlignmentName,rf_res),x=DR2,y=-(rf_err-rf_res)/(2*N-6)),data=rbind(bs[bs$AlignmentName%in%c("16S.B","16S.B-Divvier"),],temp))+
  geom_jitter(alpha=0.5,position = position_jitterdodge())+
  #geom_point(aes(y=rf_res/(2*N-6)),color="red")+
  theme_classic()+
  geom_boxplot(fill="transparent",outlier.alpha = 0,size=0.8)+
  #scale_color_gradient(low = "#55BBFF",high = "#112244")+
  scale_color_manual(name="",values=c("black","grey40","red") ,labels=c("TAPER (all)","TAPER (subset)","Divvier"))+
  scale_y_continuous("Change in Normalized RF (after-before)",labels=percent)+
  scale_x_discrete(name="Diameter")+
  geom_hline(yintercept = 0)+
  theme(legend.position = c(.2,.9),legend.direction = "vertical", legend.text.align = 1)+
  ggsave("Figures/ErrParam_Figures/16S-RF.pdf",width = 5.2,height = 5)


ggplot(aes(color=reorder(AlignmentName,rf_res),x=DR2,y=-(rfw_err-rfw_res)),data=rbind(bs[bs$AlignmentName%in%c("16S.B","16S.B-Divvier"),],temp))+
  geom_jitter(alpha=0.5,position = position_jitterdodge())+
  #facet_wrap(~AlignmentName)+
  #geom_point(aes(y=rf_res/(2*N-6)),color="red")+
  theme_classic()+
  geom_boxplot(fill="transparent",outlier.alpha = 0,size=0.8)+
  #scale_color_gradient(low = "#55BBFF",high = "#112244")+
  scale_color_manual(name="",values=c("#103060","#406090","#D02020") ,labels=c("TAPER (all)","TAPER (subset)","Divvier"))+
  scale_y_continuous("Change in  WRF (after-before)")+
  scale_x_discrete(name="Diameter")+
  geom_hline(yintercept = 0)+
  theme(legend.position = c(.2,.9),legend.direction = "vertical", legend.text.align = 1)+
  ggsave("Figures/ErrParam_Figures/16S-RFw.pdf",width = 5.2,height = 5)



ggplot(aes(fill=reorder(AlignmentName,rf_res),x=DR2,y=ifelse(rf_err==0&rf_res==0,0,(rf_err-rf_res)/rf_err)),
       data=rbind(bs[bs$AlignmentName%in%c("16S.B","16S.B-Divvier"),],temp))+
  geom_hline(yintercept = 0,color="gray70")+
  geom_boxplot(outlier.alpha = 0.7,size=0.5,color="black",outlier.size = 0.2)+
  #geom_violin(draw_quantiles = c(0.25,0.5,0.75))+
  #geom_jitter(alpha=0.5,size=0.4,position = position_jitterdodge())+
  theme_classic()+
  stat_summary(aes(color=reorder(AlignmentName,rf_res)),position = position_dodge(width = 0.75),size=.3)+
  scale_fill_manual(name="",values=c("#4090B9","#60B9F9","#DF4040") ,labels=c("TAPER (all)","TAPER (239 subset)","Divvier (239 subset)"))+
  scale_color_manual(name="",values=c("#204060","#103088","#601010") ,labels=c("TAPER (all)","TAPER (239 subset)","Divvier (239 subset)"))+
  scale_y_continuous("Relative reduction in  RF",labels=percent,lim=c(-1.9,1))+
  scale_x_discrete(name="Diameter")+
  theme(legend.position = c(.25,.2),legend.direction = "vertical", 
        plot.tag.position = c(0.01, 0.975),plot.tag = element_text(size=16,face = "bold"))+labs(tag = "d)")+
  ggsave("Figures/ErrParam_Figures/16S-rel-RF.pdf",width = 4.6*0.9,height = 3.8*0.9)

ggplot(aes(fill=reorder(AlignmentName,rf_res),x=DR2,y=(rfw_err-rfw_res)/rfw_err),
       data=rbind(bs[bs$AlignmentName%in%c("16S.B","16S.B-Divvier"),],temp))+
  geom_hline(yintercept = 0,color="gray70")+
  geom_boxplot(outlier.alpha = 0.7,size=0.5,color="black",outlier.size = 0.2)+
  #geom_violin(draw_quantiles = c(0.25,0.5,0.75))+
  #geom_jitter(alpha=0.5,size=0.4,position = position_jitterdodge())+
  theme_classic()+
  stat_summary(aes(color=reorder(AlignmentName,rf_res)),position = position_dodge(width = 0.75),size=.3)+
  scale_fill_manual(name="",values=c("#4090B9","#60B9F9","#DF4040") ,labels=c("TAPER (all)","TAPER (239 subset)","Divvier (239 subset)"))+
  scale_color_manual(name="",values=c("#204060","#103088","#601010") ,labels=c("TAPER (all)","TAPER (239 subset)","Divvier (239 subset)"))+
  scale_y_continuous("Relative reduction in  RF",labels=percent,lim=c(-1.9,1))+
  scale_x_discrete(name="Diameter")+
  theme(legend.position = c(.25,.2),legend.direction = "vertical", 
        plot.tag.position = c(0.01, 0.975),plot.tag = element_text(size=16,face = "bold"))+labs(tag = "e)")+
  ggsave("Figures/ErrParam_Figures/16S-RFw-change.pdf",width = 4.6*0.9,height = 3.8*0.9)

  ggplot(aes(x=DR2,y=value,color=variable),data=melt(bs[bs$AlignmentName%in%c("16S.B"),c("DR2","N", "rf_errN","rfw_err")],id.vars = c("DR2","N")))+
  geom_jitter(alpha=0.5,position = position_jitterdodge())+
  #geom_density()+
  #geom_point(aes(y=rf_res/(2*N-6)),color="red")+
  theme_classic()+
  geom_boxplot(fill="transparent",outlier.alpha = 0,size=0.8)+
  #scale_color_gradient(low = "#55BBFF",high = "#112244",name="Error before")+
  scale_color_brewer(name="",palette = "Paired" ,labels=c("RF","WRF"))+
  scale_y_continuous("(w)RF error before filtering")+
  scale_x_discrete(name="Diameter")+
  theme(legend.position = c(.2,.9), legend.text.align = 0,
        plot.tag.position = c(0.01, 0.975),plot.tag = element_text(size=12,face = "bold"))+labs(tag = "g)")+
 ggsave("Figures/ErrParam_Figures/16S-RF-error.pdf",width = 2.68*1.1,height = 3.5*1.1)
  
ggplot(aes(color=AlignmentName),
       data=bs[grepl("small-10-a",bs$AlignmentName),])+
  geom_hline(yintercept = 0,color="gray40",linetype=1)+
  geom_violin(aes(y=(rf_err-rf_res)/(rf_err),x="RF"),draw_quantiles = c(0.25,0.5,0.75))+
  geom_violin(aes(y=(rfw_err-rfw_res)/(rfw_err),x="WRF"),draw_quantiles = c(0.25,0.5,0.75))+
  theme_classic()+
  theme(legend.position = c(.29,.94), legend.text.align = 0,
                      plot.tag.position = c(0.01, 0.975),plot.tag = element_text(size=12,face = "bold"))+labs(tag = "f)")+
  scale_y_continuous(name="Relative reduction in RF",labels=percent)+
  scale_color_brewer(palette = "Dark2",name=element_blank(),labels=c("TAPER","DivA"))+
  scale_x_discrete(name=element_blank())
ggsave("Figures/ErrParam_Figures/small-10-aa_RF.pdf",width = 2.2,height =4)

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
       data=bh2)+
  #geom_boxplot(outlier.alpha = .5, outlier.size = 0.4)+#geom_point(alpha=0.5,size=1)+
  geom_segment(size=0.8,arrow = arrow(length=unit(0.2,"cm")))+
  #stat_summary(position = position_dodge(width=0.3),geom="line")+
  theme_classic()+
  theme(legend.position = c(.27,.891),legend.direction = "horizontal", legend.text.align = 1,axis.text.x = element_text(angle=90))+
  scale_y_continuous(name="WRF error")+
  scale_shape(name="")+
  scale_x_discrete(name="")+
  scale_color_gradient(low = "#55BBFF",high = "#112244",name="Diameter")+
  ggsave("Figures/ErrParam_Figures/Hacket-wRF.pdf",width = 5,height = 5)



ggplot(aes(ymax=rf_res/(2*N-6),ymin=rf_err/(2*N-6),
           yend=rf_res/(2*N-6),y=rf_res/(2*N-6),
           x=reorder(Rep,rf_err),xend=reorder(Rep,rf_err),
           shape=rf_res>rf_err,color=AlignmentName),
       data=bs[grepl("small-10-aa-RV100-BBA0039",bs$AlignmentName),])+
  #geom_boxplot(outlier.alpha = .5, outlier.size = 0.4)+#geom_point(alpha=0.5,size=1)+
  geom_linerange(position = position_dodge(width=0.4),size=0.8,arrow = arrow(length=unit(0.1,"cm")))+
  geom_segment(position = position_dodge(width=0.4),size=0.8,arrow = arrow(length=unit(0.1,"cm")))+
  #stat_summary(position = position_dodge(width=0.3),geom="line")+
  theme_classic()+
  theme(legend.position = c(.27,.891),legend.direction = "horizontal", legend.text.align = 1,axis.text.x = element_text(angle=90))+
  scale_y_continuous(name="Normalized RF error",labels=percent)+
  scale_shape(name="")+
  scale_x_discrete(name="")+
  scale_color_brewer(name="",palette = "Dark2",labels=c("TAPER","DivA"))+
  ggsave("Figures/ErrParam_Figures/AA-RF.pdf",width = 6,height = 5)

ggplot(aes(yend=rfw_res,y=rfw_err,x=reorder(Rep,rf_err/(2*N-6)),xend=reorder(Rep,rf_err),
           shape=rf_res>rf_err,color=AlignmentName),
       data=bs[grepl("small-10-aa-RV100-BBA0039",bs$AlignmentName),])+
  #geom_boxplot(outlier.alpha = .5, outlier.size = 0.4)+#geom_point(alpha=0.5,size=1)+
  geom_segment(position = position_dodge2(width=0.8),size=0.8,arrow = arrow(length=unit(0.2,"cm")))+
  #stat_summary(position = position_dodge(width=0.3),geom="line")+
  theme_classic()+
  theme(legend.position = c(.27,.891),legend.direction = "horizontal", legend.text.align = 1,axis.text.x = element_text(angle=90))+
  scale_y_continuous(name="WRF error")+
  scale_shape(name="")+
  scale_x_discrete(name="")+
  scale_color_brewer(name="",palette = "Dark2",labels=c("TAPER","DivA"))+
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


go=d[ (d$E =="Hackett_Genes_ErrLen") &d$ErrLen<64 & d$DR!="concat" &d$ErrLen == 8 &d$n=="5%"&d$Rep==0,"SL"]
hb=ggplot()+
  #geom_boxplot(outlier.alpha = .5, outlier.size = 0.4)+#geom_point(alpha=0.5,size=1)+
  geom_segment(aes(yend=rf_res,y=rf_err,x=reorder(DiameterRangeGene,go),xend=reorder(DiameterRangeGene,go),
                   color="WRF"), size=0.8,arrow = arrow(length=unit(0.2,"cm")),data=bh2[bh2$DiameterRangeGene != "concat",])+
  geom_segment(aes(yend=5*rf_res/(2*N-6),y=5*rf_err/(2*N-6),x=reorder(DiameterRangeGene,go),xend=reorder(DiameterRangeGene,go),
                   color="RF"), size=0.8,arrow = arrow(length=unit(0.2,"cm")),position =position_nudge(x=0.3),data=bh[bh$DiameterRangeGene != "concat",])+
  #stat_summary(position = position_dodge(width=0.3),geom="line")+
  theme_classic()+
  theme(legend.position = c(.2,.91),legend.direction = "horizontal", legend.text.align = 1,#axis.text.x = element_text(angle=90), 
        plot.tag.position = c(0.01, 0.975),plot.tag = element_text(size=15,face = "bold"))+labs(tag = "c)")+
  scale_y_continuous(name="WRF error",
                     sec.axis = sec_axis(~./5, name="RF",labels = percent))+
  scale_x_discrete(name="")+
  scale_color_manual(name="",values = c("#003092","#60D290"))
#scale_color_gradient(low = "#55BBFF",high = "#112244",name="Diameter")+
hb
ggsave("Figures/ErrParam_Figures/Hacket-both.pdf",width = 9,height =3.2)

haf=ggplot(aes(x=reorder(paste(DR,round(Diameter,3),round(SL,0),as.numeric(as.character(N)),sep="\n"),
                     SL),
           xend=reorder(paste(DR,round(Diameter,3),round(SL,0),as.numeric(as.character(N)),sep="\n"),
                        SL),
           yend=FN/(FN+TN),y=(FN+TP)/(FN+FP+TP+TN),color=interaction(ErrLenT,n,sep=", ")),
       data=data.table::dcast(setDT(d[ (d$E =="Hackett_Genes_ErrLen" | d$E=="Hackett_Genes_NumErrAlns") &d$ErrLen<64& d$DR!="concat",]),ErrLenT+n+SL+N+Diameter+DR+E~.,fun.aggregate = mean,value.var=c("FP","TP","TN","FN")))+
  geom_segment(position = position_dodge(width=0.8),size=0.8,arrow = arrow(length=unit(0.2,"cm")))+
  theme_classic()+
  theme(legend.position = "right",legend.direction = "vertical", legend.text.align = 1,
        legend.spacing = unit(0,"pt"),legend.margin = margin(0,0,0,0,"pt"), legend.box.margin = margin(0,0,0,0,"pt"),legend.box.spacing = unit(2,"pt"), 
        plot.tag.position = c(0.01, 0.975),plot.tag = element_text(size=15,face = "bold"))+labs(tag = "b)")+
  scale_y_log10("Percent error",labels=percent)+
  scale_shape(name="")+scale_x_discrete(name="")+
  facet_wrap(~E,labeller = function(x) list(E=c("Changing Error Length","Changing Error Frequency")),ncol=1)+
  scale_color_brewer(palette = "Paired",name="Error Prof.")

haf+ggsave("Figures/ErrParam_Figures/Hackett_ErrLenNumErrAlns_percenterror_arrow_log.pdf",width = 9,height =6.5)

require(gridExtra)

g=arrangeGrob(hr,haf, hb, nrow=2,widths=c(1,1.1),heights=c(2,1),layout_matrix = rbind(c(1, 2), c(1,3)))
ggsave("Figures/ErrParam_Figures/Hackett_all.pdf",width = 16,height =8.8,g)


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
hr = ggplot(data=A, aes(x, y, size=n, shape=ErrLen,color=reorder(DR,y/x))) + geom_point(alpha=.99)+
  theme_light()+
  theme(legend.position = "right",axis.text = element_text(size=15),axis.title = element_text(size=16), legend.spacing = unit(0,"pt"),
        plot.tag.position = c(0.01, 0.975),plot.tag = element_text(size=15,face = "bold"))+labs(tag = "a)")+
  scale_shape_manual(name="Error Length",values=c(6,5,4,1,3,2,19,17,18,16,15,14,90))+scale_color_discrete(name="Gene")+
  scale_x_continuous(name="FPR",labels=percent)+
  scale_size_manual(name="Error Frequency",values=sqrt(c(1,4,10,20,40)))+
  scale_y_continuous("Recall",labels=percent)
hr
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



ggplot(aes(x=n,y=TP/(TP+FN),fill=ErrLenT),data=d[d$E %in% c( "small-10-aa_ErrLen","small-10-aa_NumErrAlns") ,])+
  geom_boxplot(aes(color="TAPER"))+
  theme_classic()+
  theme(legend.position = c(.78,.23),legend.direction = "vertical", legend.text.align = 1, legend.box = "horizontal")+
  geom_smooth(se=F,method="lm")+scale_y_continuous("Recall",labels=percent)+
  scale_shape(name="")+scale_x_discrete(name="Error Frequency")+
  scale_fill_brewer(palette = "Paired",name="")+
  geom_boxplot(aes(x="5%",color="DivA",fill="8×11"), width = 0.2,
               data=t[t$AlignmentName=="small-10-aa-RV100-BBA0039-DivA",])+
  geom_boxplot(aes(x="5%",color="Divvier",fill="8×11"), width = 0.2,
               data=t[t$AlignmentName=="small-10-aa-RV100-BBA0039-Divvier",])+
  scale_color_manual(name="",values=c("red", "orange", "black"))+
ggsave("Figures/ErrParam_Figures/small-10-aa_ErrLenNumErr_Recall.pdf",width = 5,height = 4.5)

ggplot(aes(x=n,y=FP/(TP+FP),fill=ErrLenT),data=d[d$E %in% c( "small-10-aa_ErrLen","small-10-aa_NumErrAlns") ,])+
  geom_boxplot(aes(color="TAPER"))+
  theme_classic()+
  theme(legend.position = c(.78,.83),legend.direction = "vertical", legend.text.align = 1, legend.box = "horizontal")+
  geom_smooth(se=F,method="lm")+scale_y_continuous("FDR",labels=percent)+
  scale_shape(name="")+scale_x_discrete(name="Error Frequency")+
  scale_fill_brewer(palette = "Paired",name="")+
  geom_boxplot(aes(x="5%",color="DivA",fill="8×11"), width = 0.2,
               data=t[t$AlignmentName=="small-10-aa-RV100-BBA0039-DivA",])+
  geom_boxplot(aes(x="5%",color="Divvier",fill="8×11"), width = 0.2,
               data=t[t$AlignmentName=="small-10-aa-RV100-BBA0039-Divvier",])+
  scale_color_manual(name="",values=c("red", "orange", "black"))+
ggsave("Figures/ErrParam_Figures/small-10-aa_ErrLenNumErr_FDR.pdf",width = 5,height = 4.5)

ggplot(aes(x=n,y=FP/(TN+FP),fill=ErrLenT),data=d[d$E %in% c( "small-10-aa_ErrLen","small-10-aa_NumErrAlns") ,])+
  geom_boxplot(aes(color="TAPER"))+
  theme_classic()+
  theme(legend.position = c(.78,.83),legend.direction = "vertical", legend.text.align = 1, legend.box = "horizontal")+
  geom_smooth(se=F,method="lm")+scale_y_continuous("FPR",labels=percent,trans = "log10")+
  scale_shape(name="")+scale_x_discrete(name="Error Frequency")+
  scale_fill_brewer(palette = "Paired",name="")+
  geom_boxplot(aes(x="5%",color="DivA",fill="8×11"), width = 0.4,
               data=t[t$AlignmentName=="small-10-aa-RV100-BBA0039-DivA",])+
  geom_boxplot(aes(x="5%",color="Divvier",fill="8×11"), width = 0.2,
               data=t[t$AlignmentName=="small-10-aa-RV100-BBA0039-Divvier",])+
  scale_color_manual(name="",values=c("red", "orange", "black"))+
ggsave("Figures/ErrParam_Figures/small-10-aa_ErrLenNumErr_FPR.pdf",width = 5,height = 4.5)


ggplot(aes(x=interaction(ErrLenT,n,sep=","),color="aTAPER",xend=interaction(ErrLenT,n,sep=","),yend=FN/(FN+TN),y=(FN+TP)/(FN+FP+TP+TN)), #,color=interaction(ErrLenT,n,sep=", ")),
       data=data.table::dcast(setDT(d[d$E %in% c( "small-10-aa_ErrLen","small-10-aa_NumErrAlns"),]),ErrLenT+n+E~.,fun.aggregate = mean,value.var=c("FP","TP","TN","FN")))+
  #geom_boxplot(outlier.alpha = .5, outlier.size = 0.4)+#geom_point(alpha=0.5,size=1)+
  geom_segment(size=0.8,arrow = arrow(length=unit(0.2,"cm")))+
  #stat_summary(position = position_dodge(width=0.3),geom="line")+
  theme_classic()+
  theme(legend.position = c(.07,.17),legend.direction = "vertical", legend.text.align = 1, 
        plot.tag.position = c(0.01, 0.975),plot.tag = element_text(size=15,face = "bold"))+labs(tag = "e)")+
  scale_y_log10("Percent error",labels=percent)+
  scale_shape(name="")+scale_x_discrete(name="Error Len, Freq.")+
  facet_wrap(~E,labeller = function(x) list(E=c("Changing Error Length","Changing Error Frequency")),scales="free_x")+
  scale_color_brewer(palette = "Dark2",name="",labels=c("TAPER","DivA"))+
  geom_segment(aes(color="DivA",x="8×11,5%",xend="8×11,5%"),
               data=data.table::dcast(setDT(t[t$AlignmentName=="small-10-aa-RV100-BBA0039-DivA",]),ErrLen~.,fun.aggregate = mean,value.var=c("FP","TP","TN","FN")),
               position = position_nudge(x=0.2),size=0.8,arrow = arrow(length=unit(0.2,"cm")))+
  geom_segment(aes(color="Divvier",x="8×11,5%",xend="8×11,5%"),
               data=data.table::dcast(setDT(t[t$AlignmentName=="small-10-aa-RV100-BBA0039-Divvier",]),ErrLen~.,fun.aggregate = mean,value.var=c("FP","TP","TN","FN")),
               position = position_nudge(x=0.4),size=0.8,arrow = arrow(length=unit(0.2,"cm")))+
ggsave("Figures/ErrParam_Figures/small-10-aa_ErrLenNumErr_percenterror_arrow_log.pdf",width = 7.5,height =4)

ggplot(aes(x=interaction(ErrLenT,n,sep=", "),color="TAPER",
           y=(FN+TN)/(FN+FP+TP+TN)), #,color=interaction(ErrLenT,n,sep=", ")),
       data=d[d$E %in% c( "small-10-aa_ErrLen","small-10-aa_NumErrAlns"),])+
  stat_summary(position = position_dodge(width=0.8))+
  theme_classic()+
  theme(legend.position = "bottom",legend.direction = "horizontal", legend.text.align = 1)+
  scale_y_continuous("Portion of alignment retained",labels=percent)+
  scale_shape(name="")+scale_x_discrete(name="")+
  facet_wrap(~E,labeller = function(x) list(E=c("Changing Error Length","Changing Error Frequency")),scales="free_x")+
  scale_color_brewer(palette = "Dark2",name="Method")+
  stat_summary(aes(color="DivA",x="8×11, 5%"),data=t[t$AlignmentName=="small-10-aa-RV100-BBA0039-DivA",])+
  stat_summary(aes(color="Divvier",x="8×11, 5%"),data=t[t$AlignmentName=="small-10-aa-RV100-BBA0039-Divvier",])+
ggsave("Figures/ErrParam_Figures/small-10-aa_ErrLenNumErrAlns_sizechange.pdf",width = 7.5,height =4)

qplot(AlignmentName,time/3600,data=t)+scale_y_log10()
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
d2=summ_roc(d[d$E %in% c( "small-10-aa_ErrLen","small-10-aa_NumErrAlns"),c(1,10:18,20)], ErrLenT+n~.)
d2$m="TAPER"
d3=summ_roc(t5, E+ErrLenT+n~.)
d3$m=factor(d3$E,labels = c("DivA","Divvier"))
d2=rbind(d2,d3[,2:12])
A = data.frame(x=d2$FP/(d2$FP+d2$TN),y=d2$TP/(d2$TP+d2$FN), n=d2$n, ErrLen=d2$ErrLenT,m =d2$m)
#B = data.frame(x=d2$FP0/(d2$FP0+d2$TN0),y=as.vector(matrix(1.1,nrow=nrow(d2))), n=d2$n, ErrLen=d2$ErrLenT)
ggplot(data=A, aes(x, y, color=as.factor(ErrLen),shape=m)) + 
  geom_point(alpha=.8)+
  theme_classic()+theme(legend.position = c(.73,.27),legend.box = "horizontal")+ #,legend.direction = "horizontal")+
  scale_shape(name=element_blank())+
  scale_color_brewer(name=element_blank(),palette = "Paired")+
  scale_x_continuous(name="FPR",labels=percent)+
  scale_y_continuous("Recall",labels=percent,breaks = c(0.2,0.4,0.6,0.8,1))+
  #geom_text(aes(label=paste(n,ErrLen,sep=", ")),nudge_y = 0.014,nudge_x = 0.003,size=2.6)+
  #coord_cartesian(xlim=c(0.0124,0.01425))
  #ggtitle("small-10-aa with Varying Error Lengths: ROC")
ggsave("Figures/ErrParam_Figures/small10aa_ErrLenNumErr_ROC.pdf", width=4, height=4)


ggplot(data=A, aes(x, y, shape=n, size=ErrLen,color=m)) + geom_point(alpha=.99)+
  theme_classic()+
  theme(legend.position = c(.7,0.25), legend.margin = margin(0,0,0,0,"pt"), legend.spacing = unit(5,"pt"), legend.box = "horizontal",
        axis.text = element_text(size=13),axis.title = element_text(size=14), 
        plot.tag.position = c(0.01, 0.975),plot.tag = element_text(size=15,face = "bold"))+labs(tag = "d)")+
  scale_shape_manual(name="Error Length",values=c(6,5,1,4,3,2,19,17,18,16,15,14,90))+
  scale_x_continuous(name="FPR",labels=percent)+
  scale_size_manual(name="Error Frequency",values=sqrt(c(2,4,12,30,60)))+
  scale_y_continuous("Recall",labels=percent)+
  scale_color_brewer(name="Method",palette = "Dark2")+
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




##########################################


hr = read.csv('CSV_Files/hacket_removed.csv',sep=" ",header=F)
head(hr)
hr$V14

htree = read.table("CSV_Files/Hacket-distances.csv", check.names=FALSE)
hm = merge(hr,data.frame(t(htree[1,]),taxon=rownames(t(htree[1,]))),by.x="V12",by.y="taxon")
head(hm)
nrow(hm)
nrow(hr)

source("ggplot_smooth_func.R")

ggplot(aes(x=Alligator, y=V14),data=hm[,])+
  theme_classic()+
  xlab("Distance to Alligator")+ylab("Mean number of nucleotides removed")+
  #stat_summary(fun.y=sum,color="darkgreen")+
  stat_summary(color="blue",fun.data = mean_se,size=.33,alpha=0.5)+
  facet_wrap((ifelse(V4==V6,1,V6))~V5,scales="free_y",ncol=2)+
  stat_smooth(se=F,method="lm",color="red")+
  #geom_jitter(size=0.1)+
  stat_smooth_func(geom="text",method="lm",hjust=0,parse=TRUE,color="red",size=3,xpos=0) +
  theme(axis.text.x=element_text(angle=90,size=7))+
  ggsave("Figures/ErrParam_Figures/Hacket-removed-len.pdf",width=8.8,height=11)


############################################
ggplot(read.csv('CSV_Files/gatsey-err-len.csv'),aes(x=L))+stat_ecdf()+
  geom_vline(xintercept = 20,color="red",linetype=2)+theme_bw()+
  xlab("Error length")+ylab("ECDF")+ggsave("Figures/ErrParam_Figures/bio-errlen.pdf",width=5,height = 5)


##########################################
require(reshape2); require(ggplot2); require(scales)

sl=(read.csv('CSV_Files/st-likelihood.txt',he=F,sep=' '))
#slold=(read.csv('CSV_Files/old-st-likelihood.txt',he=F,sep=' '))


rc=read.csv('CSV_Files/removed-count-default.txt',head=F,sep=' '); 
gtl=read.csv('CSV_Files/genetree-like.txt',sep=" ",header=F); 
m=merge(dcast(V1~V2,data=sl[,c(1,2,6)]),rc);
names(m)[5:6]=c("removed","all");
m=merge(m,dcast(V1~V2,data=gtl[,c(1,2,9)]),by="V1")


summary(with(m,removed/all))

ggplot(aes(x=removed/all,y=(before-after)/before),data=m)+
  geom_point(aes(shape="ST"),alpha=0.66)+
  geom_point(aes(y=(before-after.randomremoved)/before,shape="ST"),color="red",alpha=0.5)+
  scale_x_continuous(labels=percent,name="Portition of nucleotides removed")+
  scale_y_continuous(labels=percent,name="Increase in likelihood")+
  theme_bw()+
  geom_text(aes(label=V1),size=0,nudge_y=0.0033,data=m[m$removed/m$all>.002,])+
  geom_point(aes(y=(before.gt-after.gt)/before.gt,shape="GT"),alpha=0.66)+
  geom_point(aes(y=(before.gt-randomremoved.gt)/before.gt,shape="GT"),color="red",alpha=0.5)+
  scale_shape(name="")+theme(legend.position=c(.1,.8))+
  ggsave("Figures/ErrParam_Figures/st-likelihood.pdf",width=5,height=3.5);

ggplot(aes(x=reorder(V1,removed/all),y=(before-after)/before),data=m)+
  geom_point(aes(shape="ST"),alpha=0.66)+
  geom_point(aes(y=(before-after.randomremoved)/before),color="red",alpha=0.5)+
  scale_y_continuous(labels=percent,name="Increase in likelihood")+
  theme_bw()+
  geom_text(aes(label=V1),size=0,nudge_y=0.0033,data=m[m$removed/m$all>.002,])+
  geom_point(aes(y=(before.gt-after.gt)/before.gt,shape="GT"),alpha=0.66)+
  scale_shape(name="")+theme(legend.position=c(.1,.8),axis.text.x=element_text(angle=90))+
  scale_x_discrete(name="genes (ordered by portion removed)")+
  ggsave("Figures/ErrParam_Figures/likelihood-by-gene.pdf",width=6.5,height=3.5)

qplot(reorder(V2,`.`), `.`,data=recast(read.csv('CSV_Files/removed-length.txt',sep=' ',h=F),V2~.,fun.ag=sum,measure.var="V3"))+
  theme_bw()+xlab("species")+ylab("total nucleotides removed")+theme(axis.text.x=element_text(angle=90))+
  ggsave("Figures/ErrParam_Figures/removed-len.pdf",width=8,height=4)
