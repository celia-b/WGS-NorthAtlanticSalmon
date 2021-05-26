a=read.table("roh_500.hom.summary", header=T) 
plot((a$UNAFF/124*100), col=a$CHR, pch=20, xlab="chromosomes", ylab="% samples")