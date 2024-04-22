args <- commandArgs(trailingOnly = TRUE)

input_file<- args[1]
output_dir<- args[2]

df<-read.table(input_file)
files<-df$V5
genos<-lapply(files,read.table)
big_geno<-rbindlist(genos)
big_geno$V1<-paste0("chr",big_geno$V1)
colnames(big_geno)<-c("#CHROM","POS","Ch_v2_0_2","EuSprat")
output_file_name=paste0("../",output_dir,"/",unique(df$V1),".genes.geno")
write.table(big_geno,file=output_file_name,col.names = T,row.names=F,quote=F,sep="\t")