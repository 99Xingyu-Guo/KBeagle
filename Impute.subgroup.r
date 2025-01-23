#XYcluster is the file that needs to be converted and imputation. 
Impute.subgroup=function(XYcluster){
     print("subgroup will read file...")
     X=read.table(XYcluster,header = TRUE,na.strings="NA",sep="\t")
     print("subgroup has done in reading file.")
     X[X==0]="AA"
     X[X==1]="AG"
     X[X==2]="GG"
     X[is.na(X)]="NN"
	 print("subgroup is preparing to convert the file Format.")
     hmp=t(X)
     N=nrow(hmp)
     rs=seq(1,N)
     alleles=rep("A/G",N)
     chrom=rep(2,N)
     pos=seq(1,N)
     hmp5.11=matrix(NA,N,7)
     out.hmp=cbind(rs,alleles,chrom,pos,hmp5.11,hmp)
     colnames(out.hmp)[5:11]=c("strand","assembly","center","protLSID","assayLSID","panelLSID","QCcode")
     colnames(out.hmp)[12:ncol(out.hmp)]=paste0("V",1:(ncol(out.hmp)-11))
     write.table(out.hmp,"mid1.hmp",col.names=T,sep="\t",row.names=F,quote=F)
     rm(hmp,out.hmp,X)
	 print("subgroup has converted from num format to hapmap format.")
	 tassel_script = file.path(getwd(),"run_pipeline.pl")
     system(paste("perl",
         tassel_script,
         "-Xms512m -Xmx2g",
         "-fork1 -h",
         hapmap_file="mid1.hmp",
         "-export -exportType VCF -runfork1",
         sep=" "))
     print("subgroup has converted from hapmap format to VCF format.")
     system("java -Xmx50g -jar beagle.22Jul22.46e.jar gt=mid1.vcf out=mid3")
	 print("The file has been completed imputation.")
	 tryCatch({
         lj = gzfile("mid3.vcf.gz", "rt")
    }, error = function(e) {
         print("Error occurred while opening the file.")
         system("java -Xmx20g -jar beagle.22Jul22.46e.jar gt=mid1.vcf out=mid3")
         lj = gzfile("mid3.vcf.gz", "rt")
    })
     content = readLines(lj)  
     close(lj)
     output_file = "mid3.vcf"
     writeLines(content, output_file)
     vcf_path = file.path(getwd(),"mid3.vcf")
     snp_path = file.path(getwd(),"mid")
     command = paste("perl",
         "vcf2hpm2_wy_v1.pl", 
         shQuote(vcf_path), 
         shQuote(snp_path), 
         sep = " ")		 
     system(command)
	 print("subgroup has converted from VCF format to hapmap format.")
     X4=data.table::fread("mid2.txt",sep=" ")
     X5=as.matrix(X4[,-c(1:12)]) 
     X6=t(X5)
     index.AA=X6=="AA"
     index.AG=X6=="AG"
	 index.GA=X6=="GA"
     index.GG=X6=="GG"
     nr=nrow(X6)
     nc=ncol(X6)
     genotype.n=matrix(0,nr,nc)
     genotype.n[index.AA]=0
     genotype.n[index.AG]=1
	 genotype.n[index.GA]=1
     genotype.n[index.GG]=2
	 XYcluster=gsub(".txt",".kbeagle",XYcluster)
	 print(dim(genotype.n))
     write.table(genotype.n,file=XYcluster,quote=F,sep="\t",col.names=F,row.names=F)
     print("subgroup has finished.")
}
