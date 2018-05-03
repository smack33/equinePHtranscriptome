library(edgeR)
library(limma)
library(locfit)



setwd("/Users/savannahmack/Documents/chipper")

#read counts
files <- list.files(".","counts$")
counts <- readDGE(files)
counts$samples
head(counts$samples)
write.csv(counts$samples, "counts.csv")


## Example code
## apply the threathold on 50% of cases = ~12
keep<-rowSums(cpm(counts)>0.05) >= 12
counts <- counts[keep, , keep.lib.sizes=FALSE]
## normalization
counts <- calcNormFactors(counts)



#design matrix Diseased State (D=diseased, H=healthy) Season (W=winter, S=summer)
Condition2_<-c(rep("H", 23), rep("D",25))
Subject2_<- c("1","1","1","1","2","2","2","2","3","3","3","4","4","4","4","5","5","5","5","6","6","6","6","1","1","1","1","2","2","2","2","3","3","3","3","4","4","4","4","5","5","5","5","6","6","6","6","6")
Season2_<- c("W","W","S","S","W","W","S","S","W","W","S","S","S","W","W","S","S","W","W","S","S","W","W","W","W","S","S","W","W","S","S","W","W","S","S","W","W","S","S","W","W","S","S","W","W","S","S","S")

Targets2<-cbind(Subject2_, Condition2_, Season2_)
head(Targets2)
Targets2<-as.data.frame(Targets2)
#as.factor(Condition2_)
#relevel(as.factor(Condition2_), ref = "D")

#create factors
Condition2_ <- factor(Targets2$Condition2_, levels=c("H","D"))
Season2_ <- factor(Targets2$Season2_, levels=c("W","S"))
#Subject2_ <- factor(Targets2$Subject2_, levels=c("1","2","3","4","5","6","7","8","9","10","11","12"))


design2<- model.matrix(~Condition2_+Condition2_:Subject2_+Condition2_:Season2_)

design2
colnames(design2)

#estimate dispersion
y2 <- estimateGLMCommonDisp(counts,design2)
y2 <- estimateGLMTrendedDisp(y2,design2)
y2 <- estimateGLMTagwiseDisp(y2,design2)

#check to see if data frame matches the design!
data.frame(Condition2_, Subject2_, Season2_)

#fit to a linear model
fit2 <- glmFit(y2, design2)

#genes that respond differently in disease state (comparison1)
lrt <- glmLRT(fit2, coef="Condition2_D")
topTags(lrt)

#gene changes in disease state in summer vs winter (comparison 2)
lrt <- glmLRT(fit2, coef="Condition2_D:Season2_S")
topTags(lrt)

#genes chanegs in healthy state, summer vs winter (comparison 3)
lrt <- glmLRT(fit2, coef="Condition2_H:Season2_S")
topTags(lrt)

#genes that are affected by season depending on disease state (?) (comparison 4)
lrt <- glmLRT(fit2, contrast=c(0,0,0,0,0,0,0,0,0,0,0,0,-1,1))
topTags(lrt)

#genes that respond to summer in healthy and disease (comparison5 is just a combined table for 2 & 3)
lrt <- glmLRT(fit2, coef=13:14)
topTags(lrt)

#comparison 4 reversed (comparison 6)
lrt <- glmLRT(fit2, contrast=c(0,0,0,0,0,0,0,0,0,0,0,0,1,-1))
topTags(lrt)

#genes respond n a overall affect to summer (comparison 7)
lrt <- glmLRT(fit2, contrast=c(0,0,0,0,0,0,0,0,0,0,0,0,1/2,1/2))
topTags(lrt)





#make table 1
GLMresults1<-topTags(lrt, n=20000, p.value=1)
rel_tags<-GLMresults1$table
head(GLMresults1$table)

#combine annotation files
annot<-read.delim("RNAseqSupTrans.merge (1).reduced")
head(annot)
GLMresultsnew_names<-rownames(GLMresults1$table)
rel_tags<-cbind(GLMresultsnew_names, rel_tags)
anot_rel_all<-merge(rel_tags, annot, by.x = "GLMresultsnew_names", by.y="transcript.ID")
head(anot_rel_all)
write.csv(anot_rel_all, "comparison1.csv")


#make new table (2)
GLMresults2<-topTags(lrt, n=20000, p.value=1)
rel_tags<-GLMresults2$table
head(GLMresults2$table)
#combine w annotation file
GLMresultsnew_names<-rownames(GLMresults2$table)
rel_tags<-cbind(GLMresultsnew_names, rel_tags)
anot_rel_all<-merge(rel_tags, annot, by.x = "GLMresultsnew_names", by.y="transcript.ID")
head(anot_rel_all)
write.csv(anot_rel_all, "comparison2.csv")

#make new table (3)
GLMresults3<-topTags(lrt, n=20000, p.value=1)
rel_tags<-GLMresults3$table
head(GLMresults3$table)
#combine w annotation file
GLMresultsnew_names<-rownames(GLMresults3$table)
rel_tags<-cbind(GLMresultsnew_names, rel_tags)
anot_rel_all<-merge(rel_tags, annot, by.x = "GLMresultsnew_names", by.y="transcript.ID")
head(anot_rel_all)
write.csv(anot_rel_all, "comparison3.csv")

#make new table (4)
GLMresults4<-topTags(lrt, n=20000, p.value=1)
rel_tags<-GLMresults4$table
head(GLMresults4$table)
#combine w annotation file
GLMresultsnew_names<-rownames(GLMresults4$table)
rel_tags<-cbind(GLMresultsnew_names, rel_tags)
anot_rel_all<-merge(rel_tags, annot, by.x = "GLMresultsnew_names", by.y="transcript.ID")
head(anot_rel_all)
write.csv(anot_rel_all, "comparison4.csv")

#make new table (5)
GLMresults5<-topTags(lrt, n=20000, p.value=1)
rel_tags<-GLMresults5$table
head(GLMresults5$table)
#combine w annotation file
GLMresultsnew_names<-rownames(GLMresults5$table)
rel_tags<-cbind(GLMresultsnew_names, rel_tags)
anot_rel_all<-merge(rel_tags, annot, by.x = "GLMresultsnew_names", by.y="transcript.ID")
head(anot_rel_all)
write.csv(anot_rel_all, "comparison5.csv")

#make new table (6)
GLMresults6<-topTags(lrt, n=20000, p.value=1)
rel_tags<-GLMresults6$table
head(GLMresults6$table)
#combine w annotation file
GLMresultsnew_names<-rownames(GLMresults6$table)
rel_tags<-cbind(GLMresultsnew_names, rel_tags)
anot_rel_all<-merge(rel_tags, annot, by.x = "GLMresultsnew_names", by.y="transcript.ID")
head(anot_rel_all)
write.csv(anot_rel_all, "comparison6.csv")

#make new table (7)
GLMresults7<-topTags(lrt, n=20000, p.value=1)
rel_tags<-GLMresults7$table
head(GLMresults7$table)
#combine w annotation file
GLMresultsnew_names<-rownames(GLMresults7$table)
rel_tags<-cbind(GLMresultsnew_names, rel_tags)
anot_rel_all<-merge(rel_tags, annot, by.x = "GLMresultsnew_names", by.y="transcript.ID")
head(anot_rel_all)
write.csv(anot_rel_all, "comparison7.csv")



