y9_results<-read.delim("path to y9_80_20_repeat_level_results.txt")
y9_sig_counts<-tapply(y9_results$test_nominal_same_primary,y9_results$probe,sum)

y15_results<-read.delim("path to y15_80_20_repeat_level_results.txt")
y15_sig_counts<-tapply(y15_results$test_nominal_same_primary,y15_results$probe,sum)

y22_results<-read.delim("path to y22_80_20_repeat_level_results.txt")
y22_sig_counts<-tapply(y22_results$test_nominal_same_primary,y22_results$probe,sum)

pdf("Supplementary_Figure_6.pdf")
par(mfrow=c(1,3))
hist(y9_sig_counts,breaks=30,main = "",xlab = "Age 9")
hist(y15_sig_counts,breaks=30,main = "",xlab = "Age 15")
hist(y22_sig_counts,breaks=30,main = "",xlab = "Age 22")
dev.off()
par(mfrow=c(1,1))
