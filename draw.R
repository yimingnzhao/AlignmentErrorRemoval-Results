require(ggplot2); require(scales); require(reshape2)
d=(read.csv("./CSV_Files/res.csv", sep=",", header=F))
names(d) <- c("E", "DR", "X", "Diameter", "PD", "N", "ErrLen", "NumErrSeqDiv", "Rep", "FP0", "FN0", "TP0", "TN0", "FP", "FN", "TP", "TN")
d=d[d$E!="small-10-aa_NumErrAlns",]
nlabels = c("1","2%","5%","10%","20%")

d$n=with(d,as.factor(round(100/((NumErrSeqDiv==0&grepl("ErrLen$",E))*20+(NumErrSeqDiv!=N|grepl("ErrLen$",E))*NumErrSeqDiv+(NumErrSeqDiv==N&!grepl("ErrLen$",E))*100))))
levels(d$n) <- c(levels(d$n)[1],paste(levels(d$n)[0:-1],"%",sep=""),"~5%")
d[grepl("General$",d$E),"n"]="~5%"

d$ErrLen = (d$ErrLen==0)*8+d$ErrLen
d$ErrLenT = paste(d$ErrLen, intToUtf8(215), "11",sep="")
d[grepl("General$",d$E),"ErrLenT"]="~50"
d$ErrLenT = factor(d$ErrLenT,levels=c("2×11","3×11", "4×11", "8×11", "16×11", "32×11", "64×11","~50" ))


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






# ROC all 
options(digits = 2)
d$E2=sub("_.*","",d$E)
d2=summ_roc(d[grepl("Err",d$E) & (d$N > 19|d$N==0)  ,], E2+n+ErrLenT+cut(Diameter, breaks = c(0, 0.1, 0.2, 0.5, 0.8, 1), right = F)~.)
A = data.frame(x=d2$FP/(d2$FP+d2$TN),y=d2$TP/(d2$TP+d2$FN),E =d2$E2, ErrLenT=d2$ErrLenT,  n=d2$n, DR=d2$`cut(Diameter, breaks = c(0, 0.1, 0.2, 0.5, 0.8, 1), right = F)`)
B = data.frame(x=d2$FP0/(d2$FP0+d2$TN0),y=as.vector(matrix(1,nrow=nrow(d2))), E =d2$E2, n=d2$n, ErrLenT=d2$ErrLenT, DR=d2$`cut(Diameter, breaks = c(0, 0.1, 0.2, 0.5, 0.8, 1), right = F)`)
ggplot(data=A, aes(x, y, shape=interaction(n,ErrLenT,sep=", "),color=as.factor(DR))) + 
  geom_point(alpha=1)+facet_wrap(~sub("_.*","",E),ncol=1,scales="free_x")+
  geom_path(aes(group=interaction(DR,n)),linetype=2)+geom_path(aes(group=interaction(DR,ErrLenT)),linetype=1)+
  theme_bw()+theme(legend.position = "right",legend.text.align = 1)+
  scale_shape_manual(name="Err Freq%, Len",values=c(1,2,5,15,17,8,19,18,3,4,6,10))+
  scale_color_brewer(name="Diameter",palette = "Dark2")+
  scale_x_continuous(name="FPR",labels=percent)+
  scale_y_continuous("Recall",labels=percent)+
  geom_linerange(aes(x=x,ymin=0.995,ymax=1.005,color=as.factor(DR)),data=B,linetype=1,size=1)
ggsave("Figures/ErrParam_Figures/All_both_ROC.pdf", width=6, height=10)



# 16S.B - Varying both

# ROC for 16S.B with varying error lengths and fixed percentage of erroneous sequences
options(digits = 2)
d2=summ_roc(d[d$E %in% c( "16S.B_ErrLen","16S.B_NumErrAlns", "16S.B_General") & d$N > 19 ,], n+ErrLenT+cut(Diameter, breaks = c(0, 0.1, 0.2, 0.5, 0.8, 1), right = F)~.)
A = data.frame(x=d2$FP/(d2$FP+d2$TN),y=d2$TP/(d2$TP+d2$FN), ErrLenT=d2$ErrLenT,  n=d2$n, DR=d2$`cut(Diameter, breaks = c(0, 0.1, 0.2, 0.5, 0.8, 1), right = F)`)
B = data.frame(x=d2$FP0/(d2$FP0+d2$TN0),y=as.vector(matrix(1,nrow=nrow(d2))), n=d2$n, ErrLenT=d2$ErrLenT, DR=d2$`cut(Diameter, breaks = c(0, 0.1, 0.2, 0.5, 0.8, 1), right = F)`)
ggplot(data=A, aes(x, y, shape=interaction(ErrLenT,n,sep=", "),color=as.factor(DR))) + 
  geom_point(alpha=1)+
  geom_path(aes(group=interaction(DR,n)),data=A[A$n!="~5%",],linetype=2)+
  geom_path(aes(group=interaction(DR,ErrLenT)),data=A[A$n!="~5%",],linetype=1)+
  theme_bw()+theme(legend.position = "right",legend.text.align = 1)+
  scale_shape_manual(name="Length, Freq",values=c(15,17,1,2,5,8,9,7,6,19,18,3))+
  scale_color_brewer(name="Diameter",palette = "Dark2")+
  scale_x_continuous(name="FPR",labels=percent)+
  scale_y_continuous("Recall",labels=percent)+
  geom_linerange(aes(x=x,ymin=0.995,ymax=1.005,color=as.factor(DR)),data=B,linetype=1,size=1)
ggsave("Figures/ErrParam_Figures/16SB_ErrLenNumErr_ROC.pdf", width=8, height=6)



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



ggplot(data=A, aes(x, y, shape=interaction(ErrLenT,n,sep=", "),color=as.factor(DR))) + 
  geom_point(alpha=1)+
  geom_path(aes(group=DR),data=A[A$n!="~5%",],linetype=1)+
  theme_bw()+theme(legend.position = "bottom",legend.text.align = 1)+
  scale_shape_manual(name="",values=c(15,17,1,2,5,8,9,7,6,19,18,3))+
  scale_color_brewer(name="Diameter",palette = "Dark2")+
  scale_x_continuous(name="FPR",labels=percent)+
  scale_y_continuous("Recall",labels=percent)+facet_wrap(~ErrLenT=="8×11",labeller = function(x) list(c("Changing Error Length", "Changing Error Frequency")))+
  geom_linerange(aes(x=x,ymin=0.995,ymax=1.005,color=as.factor(DR)),data=B,linetype=1,size=1)
ggsave("Figures/ErrParam_Figures/16SB_ErrLenNumErr_ROC_faceted.pdf", width=10, height=6)

ggplot(aes(x=Diameter,y=TP/(TP+FN),color=interaction(ErrLenT,n,sep=", ")),data=d[d$E %in% c( "16S.B_ErrLen","16S.B_NumErrAlns") & d$N > 19,])+
  geom_point(alpha=0.4,size=.5)+
  theme_classic()+theme(legend.position = c(.75,.15),legend.direction = "horizontal", legend.text.align = 1)+
  geom_smooth()+scale_y_continuous("Recall")+
  scale_shape(name="")+facet_wrap(~E,labeller = function(x) list(E=c("Changing Error Length","Changing Error Frequency")))+
  scale_color_brewer(palette = "Paired",name="", labels = function(x) (paste(sub("1%","  1",x), intToUtf8(215), "11",sep="")))
ggsave("Figures/ErrParam_Figures/16S.B_ErrLenNumErr_Recall.pdf",width = 9,height = 4.5)

ggplot(aes(x=Diameter,y=FP/(TP+FP),color=interaction(ErrLen,n,sep=", ")),data=d[d$E %in% c( "16S.B_ErrLen","16S.B_NumErrAlns") & d$N > 19,])+
  geom_point(alpha=0.4,size=.5)+
  theme_classic()+theme(legend.position ="bottom",legend.direction = "horizontal", legend.text.align = 1)+
  geom_smooth()+scale_y_continuous("FDR")+
  scale_shape(name="")+facet_wrap(~E,labeller = function(x) list(E=c("Changing Error Length","Changing Error Frequency")))+
  scale_color_brewer(palette = "Paired",name="", labels = function(x) (paste(sub("1%","  1",x), intToUtf8(215), "11",sep="")))
ggsave("Figures/ErrParam_Figures/16S.B_ErrLenNumErr_FDR.pdf",width = 9,height = 4.5)


ggplot(aes(x=Diameter,y=FP/(FP+TN),color=interaction(ErrLen,n,sep=", ")),data=d[d$E %in% c( "16S.B_ErrLen","16S.B_NumErrAlns") & d$N > 19,])+
  geom_point(alpha=0.4,size=.5)+
  theme_classic()+theme(legend.position ="bottom",legend.direction = "horizontal", legend.text.align = 1)+
  geom_smooth()+scale_y_continuous("FPR",labels=percent)+
  scale_shape(name="")+facet_wrap(~E,labeller = function(x) list(E=c("Changing Error Length","Changing Error Frequency")))+
  scale_color_brewer(palette = "Paired",name="", labels = function(x) (paste(sub("1%","  1",x), intToUtf8(215), "11",sep="")))
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


ggplot(aes(x=n,y=TP/(TP+FN),color=ErrLenT),data=d[d$E %in% c( "Hackett_ErrLen","Hackett_NumErrAlns") ,])+
  geom_boxplot()+#geom_point(alpha=0.5,size=1)+
  theme_classic()+theme(legend.position = c(.5,.15),legend.direction = "horizontal", legend.text.align = 1)+
  geom_smooth(se=F,method="lm")+scale_y_continuous("Recall")+
  scale_shape(name="")+scale_x_discrete(name="Error frequency")+
  scale_color_brewer(palette = "Paired",name="Error Length")
ggsave("Figures/ErrParam_Figures/Hackett_ErrLenNumErr_Recall.pdf",width = 5,height = 4.5)


ggplot(aes(x=n,y=FP/(TP+FP),color=ErrLenT),data=d[d$E %in% c( "Hackett_ErrLen","Hackett_NumErrAlns") ,])+
  geom_boxplot()+#geom_point(alpha=0.5,size=1)+
  theme_classic()+theme(legend.position = c(.15,.25),legend.direction = "vertical", legend.text.align = 1)+
  geom_smooth(se=F,method="lm")+scale_y_continuous("FDR")+
  scale_shape(name="")+scale_x_discrete(name="Error frequency")+
  scale_color_brewer(palette = "Paired",name="Error Length")
ggsave("Figures/ErrParam_Figures/Hackett_ErrLenNumErr_FDR.pdf",width = 9,height = 4.5)


ggplot(aes(x=n,y=FP/(TN+FP),color=ErrLenT),data=d[d$E %in% c( "Hackett_ErrLen","Hackett_NumErrAlns") ,])+
  geom_boxplot()+#geom_point(alpha=0.5,size=1)+
  theme_classic()+theme(legend.position = c(.15,.25),legend.direction = "vertical", legend.text.align = 1)+
  geom_smooth(se=F,method="lm")+scale_y_continuous("FPR")+
  scale_shape(name="")+scale_x_discrete(name="Error frequency")+
  scale_color_brewer(palette = "Paired",name="Error Length")
ggsave("Figures/ErrParam_Figures/Hackett_ErrLenNumErr_FPR.pdf",width = 9,height = 4.5)


# ROC for Hackett with varying error lengths and fixed percentage of erroneous sequences
options(digits = 2)
d2=summ_roc(d[d$E %in% c( "Hackett_ErrLen","Hackett_NumErrAlns") ,], ErrLenT+n~.)
A = data.frame(x=d2$FP/(d2$FP+d2$TN),y=d2$TP/(d2$TP+d2$FN), n=d2$n, ErrLen=d2$ErrLenT)
B = data.frame(x=d2$FP0/(d2$FP0+d2$TN0),y=as.vector(matrix(1.1,nrow=nrow(d2))), n=d2$n, ErrLen=d2$ErrLenT)
ggplot(data=A, aes(x, y, color=as.factor(ErrLen), shape=n)) + geom_point(alpha=1)+
  theme_light()+theme(legend.position = "right")+
  scale_shape(name="Error Frequency")+scale_color_brewer(name="Error Length",palette = "Paired")+
  scale_x_continuous(name="FPR",labels=percent)+
  scale_y_continuous("Recall",labels=percent,breaks = c(0.2,0.4,0.6,0.8,1))+
  ggtitle("Hackett with Varying Error Lengths and Frequency: ROC")
ggsave("Figures/ErrParam_Figures/Hackett_ErrLenFreq_ROC.pdf", width=6, height=6)

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



ggplot(aes(x=Diameter,y=TP/(TP+FN),color=interaction(n,ErrLen,sep="%, ")),data=d[d$E %in% c( "small-10-aa_ErrLen","small-10-aa_NumErrAlns") & d$N > 19,])+
  geom_point(alpha=0.5,size=1)+
  theme_classic()+theme(legend.position = c(.75,.15),legend.direction = "horizontal", legend.text.align = 1)+
  geom_smooth(se=F,method="lm")+scale_y_continuous("Recall")+
  scale_shape(name="")+facet_wrap(~E,labeller = function(x) list(E=c("Changing Error Length","Changing Error Frequency")))+
  scale_color_brewer(palette = "Paired",name="", labels = function(x) (paste(sub("1%","  1",x), intToUtf8(215), "11",sep="")))
ggsave("Figures/ErrParam_Figures/small-10-aa_ErrLenNumErr_Recall.pdf",width = 9,height = 4.5)

ggplot(aes(x=Diameter,y=FP/(TP+FP),color=interaction(n,ErrLen,sep="%, ")),data=d[d$E %in% c( "small-10-aa_ErrLen","small-10-aa_NumErrAlns") & d$N > 19,])+
  geom_point(alpha=0.5,size=1)+
  theme_classic()+theme(legend.position = "bottom",legend.direction = "horizontal", legend.text.align = 1)+
  geom_smooth(se=F,method="lm")+scale_y_continuous("FDR")+
  scale_shape(name="")+facet_wrap(~E,labeller = function(x) list(E=c("Changing Error Length","Changing Error Frequency")))+
  scale_color_brewer(palette = "Paired",name="", labels = function(x) (paste(sub("1%","  1",x), intToUtf8(215), "11",sep="")))
ggsave("Figures/ErrParam_Figures/small-10-aa_ErrLenNumErr_FDR.pdf",width = 9,height = 4.5)

# small-10-aa - Varying error lengths and fixed percentage of erroneous sequences
ggplot(aes(x=Diameter,y=TP/(TP+FN),color=as.factor(ErrLen)),data=d[d$E=="small-10-aa_ErrLen" & d$N > 19,])+
  geom_point(alpha=0.5)+theme_classic()+geom_smooth()+scale_y_continuous("Recall")+
  scale_shape(name="")+scale_color_brewer(palette = "Paired",name="error len")+
  ggtitle("small-10-aa with Varying Error Lengths: Recall vs Diameter")
ggsave("Figures/ErrParam_Figures/small10aa_ErrLen_Recall.pdf",width = 6,height = 6)

# small-10-aa - Varying percentage of erroneous sequences and fixed error lengths
ggplot(aes(x=Diameter,y=TP/(TP+FN), group= as.factor(100/((NumErrSeqDiv!=N)*NumErrSeqDiv+(NumErrSeqDiv==N)*100)),
           color=as.factor(100/((NumErrSeqDiv!=N)*NumErrSeqDiv+(NumErrSeqDiv==N)*100)), shape=cut((FP/(FP+TN)),breaks=c(-1,0,0.001,0.1,1))),
       data=d[d$E=="small-10-aa_NumErrAlns" ,])+geom_point(alpha=0.5)+
  theme_classic()+geom_smooth(se=F)+scale_shape_manual(name="FPR",values=c(1,16,4,2))+
  scale_color_brewer(palette = "Paired",name="n",labels=nlabels)+
  scale_y_continuous(name="Recall")+coord_cartesian(ylim=c(0.35,1))+
  ggtitle("small-10-aa with Varying Number of Erroneous Sequences: Recall vs Diameter")
ggsave("Figures/ErrParam_Figures/small10aa_NumErrAlns_Recall.pdf",width = 6,height = 6)


# ROC for small-10-aa with varying error lengths and fixed percentage of erroneous sequences
options(digits = 2)
d2=summ_roc(d[d$E=="small-10-aa_ErrLen" & d$N > 19 & d$ErrLen!=0,], ErrLen+cut(Diameter, breaks = c(0, 0.1, 0.2, 0.5, 0.8, 1), right = F)~.)
A = data.frame(x=d2$FP/(d2$FP+d2$TN),y=d2$TP/(d2$TP+d2$FN), ErrLen=d2$ErrLen, DR=d2$`cut(Diameter, breaks = c(0, 0.1, 0.2, 0.5, 0.8, 1), right = F)`)
B = data.frame(x=d2$FP0/(d2$FP0+d2$TN0),y=as.vector(matrix(1.1,nrow=nrow(d2))), ErrLen=d2$ErrLen, DR=d2$`cut(Diameter, breaks = c(0, 0.1, 0.2, 0.5, 0.8, 1), right = F)`)
ggplot(data=A, aes(x, y, color=as.factor(ErrLen), shape=as.factor(DR))) + geom_point(alpha=1)+
  theme_light()+theme(legend.position = "right")+geom_point(data=B)+
  scale_shape(name="Diameter")+scale_color_brewer(name="Error Length",palette = "Paired",labels = function(x) (paste(x, intToUtf8(215), "11")))+
  scale_x_continuous(name="FPR",labels=percent)+
  scale_y_continuous("Recall",labels=percent,breaks = c(0.2,0.4,0.6,0.8,1))+
  ggtitle("small-10-aa with Varying Error Lengths: ROC")
ggsave("Figures/ErrParam_Figures/small10aa_ErrLen_ROC.pdf", width=6, height=6)

# ROC for small-10-aa with varying percentages of erroneous sequences and fixed error lengths
d2=d[d$E=="small-10-aa_NumErrAlns",]
d2$n=with(d2,as.factor(100/((NumErrSeqDiv!=N)*NumErrSeqDiv+(NumErrSeqDiv==N)*100)))
ggplot(aes(x=FP/(FP+TN),y=TP/(TP+FN), color=as.factor(n) ),data=summ_roc(d2,n~.))+
  geom_point(alpha=1)+
  theme_light()+theme(legend.position = c(.85,.25))+
  scale_color_brewer(name="n",labels=nlabels, palette="Paired")+
  scale_x_continuous(name="FPR",labels=percent)+scale_y_continuous("Recall")+
  ggtitle("small-10-aa with Varying Number of Erroneous Sequences: ROC")
ggsave("Figures/ErrParam_Figures/small10aa_NumErrAlns_ROC.pdf", width=6, height=6)



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






