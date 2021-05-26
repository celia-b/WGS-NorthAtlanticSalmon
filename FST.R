a = read.table("salmon_fst_windows", header = FALSE)
plot(a$V6, col=a$V1, pch=20, xlab="chromosomes", ylab="mFST", ylim = c(0,0.2))
