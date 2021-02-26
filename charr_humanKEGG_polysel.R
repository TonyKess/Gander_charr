setwd("~/Desktop/Software/polysel/")
# set project variables
project.name<-"charr_mander"
data.path<-file.path("./data/charmander/")

code.path<-"./R"
empfdr.path<-"./empfdr"
results.path<-"./results"

source(file.path(code.path,'polysel.R'))

minsetsize<-5
result<-ReadSetObjTables(in.path=data.path,
                         set.info.file="SetInfo.txt",
                         set.obj.file="SetObj.txt",
                         obj.info.file="ObjInfo.txt",
                         minsetsize=minsetsize,
                         obj.in.set=F,
                         merge.similar.sets=T)


set.info<-result$set.info
obj.info<-result$obj.info
set.obj<-result$set.obj
set.info.lnk<-result$set.info.lnk

cat("Number of sets: ", nrow(set.info), "\n", sep="")
cat("Number of genes: ", nrow(obj.info), "\n", sep="")


print(set.info[1:5,],row.names=F, right=F)
print(obj.info[1:5,],row.names=F, right=F)
print(set.obj[1:5,],row.names=F, right=F)

for (sz in c(10,50,250)){
  CheckStatDistribution(obj.info,setsize=sz,
                        n.rand = 100000, 
                        xlab="sum(objStat)")
}

obj.stat<-obj.info[,c("objID", "objStat", "objBin")]
approx.null <- TRUE
use.bins <- FALSE
seq.rnd.sampling <- FALSE
nrand <- 1
test <- "highertail"
qvalue.method <- "smoother"

result<-EnrichmentAnalysis(set.info, set.obj, obj.stat,
                           nrand=nrand, approx.null=approx.null, 
                           seq.rnd.sampling=seq.rnd.sampling,
                           use.bins=use.bins, test=test,
                           do.pruning=FALSE, minsetsize=minsetsize,
                           project.txt=project.name, do.emp.fdr=FALSE,
                           qvalue.method=qvalue.method
)

set.scores.prepruning <- result$set.scores.prepruning
print(set.scores.prepruning[1:10,],row.names=F, right=F)

write.table(set.scores.prepruning, "~/Desktop/Charr_Reanalysis/Charr_polysel_human_KEGG.txt", col.names = T, row.names = F, sep = "\t", quote = F)
