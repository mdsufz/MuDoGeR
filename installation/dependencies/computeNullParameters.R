llikelihood = read.table("llikelihood.matrix")
means = apply(X=llikelihood,FUN=mean,MARGIN=1)
sigmas = apply(X=llikelihood,FUN=sd,MARGIN=1)

nullParam = cbind(means,sigmas)

write.table(nullParam,file="nullParameters.tsv",quote=F,sep = "\t",col.names=F)
