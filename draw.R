require(ggplot2); require(scales); require(reshape2)
d=(read.csv("~/Linux Files/AlignmentErrorRemoval-Results/res.csv", sep=",", header=F))
names(d) <- c("E", "DR", "X", "Diameter", "PD", "N", "ErrLen", "NumErrSeqDiv", "Rep", "FP0", "FN0", "TP0", "TN0", "FP", "FN", "TP", "TN")

# General Results: Recall vs Diameter
ggplot(aes(x=Diameter, y=TP/(TP+FN)), data=d[d$E=="16S.B" & d$N > 10 & d$ErrLen==0,]) + theme_classic() + geom_point() + geom_smooth() + scale_y_continuous("Recall") + ggtitle("16S.B: Recall vs Diameter")
ggsave("16SB_General_Recall.pdf", width=6, height=6)
ggplot(aes(x=Diameter, y=TP/(TP+FN)), data=d[d$E=="Hackett" & d$N > 10 & d$ErrLen==0,]) + theme_classic() + geom_point() + geom_smooth() + scale_y_continuous("Recall") + ggtitle("Hackett: Recall vs Diameter")
ggsave("Hackett_General_Recall.pdf", width=6, height=6)
ggplot(aes(x=Diameter, y=TP/(TP+FN)), data=d[d$E=="small-10-aa-RV100-BBA0039" & d$N > 10 & d$ErrLen==0,]) + theme_classic() + geom_point() + geom_smooth() + scale_y_continuous("Recall") + ggtitle("small-10-aa: Recall vs Diameter")
ggsave("small10aa_General_Recall.pdf", width=6, height=6)

# Recall vs Diameter: 
ggplot(aes(x=Diameter,y=TP/(TP+FN),color=as.factor(ErrLen)),data=d[d$E=="16S.B" & d$N > 10 & d$ErrLen!=0,])+
  geom_point(alpha=0.1)+theme_classic()+geom_smooth(se=F)+scale_y_continuous("Recall")+
  scale_shape(name="")+scale_color_brewer(palette = "Paired",name="error len", labels = function(x) (paste(x, intToUtf8(215), "11")))+
  ggtitle("16S.B with Specific Error Lengths: Recall vs Diameter")
ggsave("16SB_ErrLen_Recall.pdf", width=6, height=6)

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

options(digits = 2)
d2=summ_roc(d[d$E=="16S.B" & d$N > 19 & d$ErrLen!=0,], ErrLen+cut(Diameter, breaks = c(0, 0.1, 0.2, 0.5, 0.8, 1), right = F)~.)
A = data.frame(x=d2$FP/(d2$FP+d2$TN),y=d2$TP/(d2$TP+d2$FN), ErrLen=d2$ErrLen, DR=d2$`cut(Diameter, breaks = c(0, 0.1, 0.2, 0.5, 0.8, 1), right = F)`)
B = data.frame(x=d2$FP0/(d2$FP0+d2$TN0),y=as.vector(matrix(1.1,nrow=nrow(d2))), ErrLen=d2$ErrLen, DR=d2$`cut(Diameter, breaks = c(0, 0.1, 0.2, 0.5, 0.8, 1), right = F)`)
ggplot(data=A, aes(x, y, color=as.factor(ErrLen), shape=as.factor(DR))) + geom_point(alpha=1)+
  theme_light()+theme(legend.position = "right")+geom_point(data=B)+
  scale_shape(name="Diameter")+scale_color_brewer(name="Error Length",palette = "Paired",labels = function(x) (paste(x, intToUtf8(215), "11")))+
  scale_x_continuous(name="FPR",labels=percent)+
  scale_y_continuous("Recall",labels=percent,breaks = c(0.2,0.4,0.6,0.8,1))+
  ggtitle("16S.B with Specific Error Lengths: ROC")
ggsave("16SB_ErrLen_ROC.pdf", width=6, height=6)

# ROC for General 16S.B
options(digits = 5)
d2=summ_roc(d[d$E=="16S.B" & d$N > 10 & d$ErrLen==0,], ErrLen+cut(Diameter, breaks = c(0, 0.1, 0.2, 0.5, 0.8, 1), right = F)~.)
A = data.frame(x=d2$FP/(d2$FP+d2$TN),y=d2$TP/(d2$TP+d2$FN), ErrLen=d2$ErrLen, DR=d2$`cut(Diameter, breaks = c(0, 0.1, 0.2, 0.5, 0.8, 1), right = F)`)
B = data.frame(x=d2$FP0/(d2$FP0+d2$TN0),y=as.vector(matrix(1.1,nrow=nrow(d2))), ErrLen=d2$ErrLen, DR=d2$`cut(Diameter, breaks = c(0, 0.1, 0.2, 0.5, 0.8, 1), right = F)`)
ggplot(data=A, aes(x, y, color=as.factor(DR))) + geom_point(alpha=1)+
  theme_light()+theme(legend.position = "right")+geom_point(data=B)+
  scale_color_brewer(name="Diameter",palette = "Paired",labels = function(x) (paste(x)))+
  scale_x_continuous(name="FPR",labels=percent)+
  scale_y_continuous("Recall",labels=percent,breaks = c(0.2,0.4,0.6,0.8,0.9,0.95,1,1.1))+
  ggtitle("16S.B: ROC")
ggsave("16SB_General_ROC.pdf", width=6, height=6)

# ROC for General Hackett
options(digits = 5)
d2=summ_roc(d[d$E=="Hackett" & d$N > 10,], ErrLen+cut(Diameter, breaks = c(0, 0.1, 0.2, 0.5, 0.8, 1), right = F)~.)
A = data.frame(x=d2$FP/(d2$FP+d2$TN),y=d2$TP/(d2$TP+d2$FN), ErrLen=d2$ErrLen, DR=d2$`cut(Diameter, breaks = c(0, 0.1, 0.2, 0.5, 0.8, 1), right = F)`)
B = data.frame(x=d2$FP0/(d2$FP0+d2$TN0),y=as.vector(matrix(1.1,nrow=nrow(d2))), ErrLen=d2$ErrLen, DR=d2$`cut(Diameter, breaks = c(0, 0.1, 0.2, 0.5, 0.8, 1), right = F)`)
ggplot(data=A, aes(x, y, color=as.factor(DR))) + geom_point(alpha=1)+
  theme_light()+theme(legend.position = "right")+geom_point(data=B)+
  scale_color_brewer(name="Diameter",palette = "Paired",labels = function(x) (paste(x)))+
  scale_x_continuous(name="FPR",labels=percent)+
  scale_y_continuous("Recall",labels=percent,breaks = c(0.2,0.4,0.6,0.8,0.9,0.95,1,1.1))+
  ggtitle("Hackett: ROC")
ggsave("Hackett_General_ROC.pdf", width=6, height=6)

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
ggsave("small10aa_General_ROC.pdf", width=6, height=6)






d=(read.csv("~/Linux Files/AlignmentErrorRemoval-Results/variedUnion.csv", sep=",", header=F))
names(d) <- c("E", "DR", "Diameter", "PD", "N", "ErrLen", "Rep", "FP0", "FN0", "TP0", "TN0", "FP", "FN", "TP", "TN")

ggplot(aes(x=Diameter,y=TP/(TP+FN),color=as.factor(ErrLen)),data=d[d$E=="16S.B_ErrLen" & d$N > 10 & d$ErrLen!=0,])+
  geom_point(alpha=0.1)+theme_classic()+geom_smooth(se=F)+scale_y_continuous("Recall")+
  scale_shape(name="")+scale_color_brewer(palette = "Paired",name="error len", labels = function(x) (paste(x, intToUtf8(215), "11")))+
  ggtitle("16S.B without Union: Recall vs Diameter")
ggsave("16SB_NoUnion_Recall.pdf", width=6, height=6)

ggplot(aes(x=Diameter,y=TP/(TP+FN),color=as.factor(ErrLen)),data=d[d$E=="16S.B_UnionErrLen2k" & d$N > 10 & d$ErrLen!=0,])+
  geom_point(alpha=0.1)+theme_classic()+geom_smooth(se=F)+scale_y_continuous("Recall")+
  scale_shape(name="")+scale_color_brewer(palette = "Paired",name="error len", labels = function(x) (paste(x, intToUtf8(215), "11")))+
  ggtitle("16S.B using Union of 2k: Recall vs Diameter")
ggsave("16SB_2k_Recall.pdf", width=6, height=6)

ggplot(aes(x=Diameter,y=TP/(TP+FN),color=as.factor(ErrLen)),data=d[d$E=="16S.B_UnionErrLen3k" & d$N > 10 & d$ErrLen!=0,])+
  geom_point(alpha=0.1)+theme_classic()+geom_smooth(se=F)+scale_y_continuous("Recall")+
  scale_shape(name="")+scale_color_brewer(palette = "Paired",name="error len", labels = function(x) (paste(x, intToUtf8(215), "11")))+
  ggtitle("16S.B using Union of 3k: Recall vs Diameter")
ggsave("16SB_3k_Recall.pdf", width=6, height=6)

ggplot(aes(x=Diameter,y=TP/(TP+FN),color=as.factor(ErrLen)),data=d[d$E=="16S.B_UnionErrLen4k" & d$N > 10 & d$ErrLen!=0,])+
  geom_point(alpha=0.1)+theme_classic()+geom_smooth(se=F)+scale_y_continuous("Recall")+
  scale_shape(name="")+scale_color_brewer(palette = "Paired",name="error len", labels = function(x) (paste(x, intToUtf8(215), "11")))+
  ggtitle("16S.B using Union of 4k: Recall vs Diameter")
ggsave("16SB_4k_Recall.pdf", width=6, height=6)

# False Discovery Rate vs Diameter
ggplot(aes(x=Diameter,y=FP/(TP+FP),color=as.factor(ErrLen)),data=d[d$E=="16S.B_ErrLen" & d$N > 10 & d$ErrLen!=0,])+
  geom_point(alpha=0.1)+theme_classic()+geom_smooth(se=F)+scale_y_continuous("FDR")+
  scale_shape(name="")+scale_color_brewer(palette = "Paired",name="error len", labels = function(x) (paste(x, intToUtf8(215), "11")))+
  ggtitle("16S.B without Union: FDR vs Diameter")
ggsave("16SB_NoUnion_FDR.pdf", width=6, height=6)

ggplot(aes(x=Diameter,y=FP/(TP+FP),color=as.factor(ErrLen)),data=d[d$E=="16S.B_UnionErrLen2k" & d$N > 10 & d$ErrLen!=0,])+
  geom_point(alpha=0.1)+theme_classic()+geom_smooth(se=F)+scale_y_continuous("FDR")+
  scale_shape(name="")+scale_color_brewer(palette = "Paired",name="error len", labels = function(x) (paste(x, intToUtf8(215), "11")))+
  ggtitle("16S.B using Union of 2k: FDR vs Diameter")
ggsave("16SB_2k_FDR.pdf", width=6, height=6)

ggplot(aes(x=Diameter,y=FP/(TP+FP),color=as.factor(ErrLen)),data=d[d$E=="16S.B_UnionErrLen3k" & d$N > 10 & d$ErrLen!=0,])+
  geom_point(alpha=0.1)+theme_classic()+geom_smooth(se=F)+scale_y_continuous("FDR")+
  scale_shape(name="")+scale_color_brewer(palette = "Paired",name="error len", labels = function(x) (paste(x, intToUtf8(215), "11")))+
  ggtitle("16S.B using Union of 3k: FDR vs Diameter")
ggsave("16SB_3k_FDR.pdf", width=6, height=6)

ggplot(aes(x=Diameter,y=FP/(TP+FP),color=as.factor(ErrLen)),data=d[d$E=="16S.B_UnionErrLen4k" & d$N > 10 & d$ErrLen!=0,])+
  geom_point(alpha=0.1)+theme_classic()+geom_smooth(se=F)+scale_y_continuous("FDR")+
  scale_shape(name="")+scale_color_brewer(palette = "Paired",name="error len", labels = function(x) (paste(x, intToUtf8(215), "11")))+
  ggtitle("16S.B using Union of 4k: FDR vs Diameter")
ggsave("16SB_4k_FDR.pdf", width=6, height=6)

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
ggsave("16SB_NoUnion_ROC.pdf", width=6, height=6)

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
ggsave("16SB_2k_ROC.pdf", width=6, height=6)

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
ggsave("16SB_3k_ROC.pdf", width=6, height=6)

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
ggsave("16SB_4k_ROC.pdf", width=6, height=6)






