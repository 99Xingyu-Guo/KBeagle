

#KBeagle_function
#We need to keep Tassel-related files, files to impute, etc. all in one folder.The tassel folder contains sTASSEL.jar, run_pipeline.pl, and lib.
#The first column in the imputed file is the individual ID name, and the remaining columns are the SNP.
Impute.subgroup <- function(XYcluster){
     if(!file.exists("beagle.22Jul22.46e.jar")){
      system("wget -O beagle.22Jul22.46e.jar https://github.com/99Xingyu-Guo/KBeagle/blob/main/beagle.22Jul22.46e.jar")
    }
     print("subgroup will read file...")
     X <- read.table(XYcluster, header <- TRUE, na.strings <- "NA", sep <- "\t")
     print("subgroup has done in reading file.")
     X[X==0] <- "AA"
     X[X==1] <- "AG"
     X[X==2] <- "GG"
     X[is.na(X)] <- "NN"
	 print("subgroup is preparing to convert the file Format.")
     hmp <- t(X)
     N <- nrow(hmp)
     rs <- seq(1, N)
     alleles <- rep("A/G", N)
     chrom <- rep(2, N)
     pos <- seq(1, N)
     hmp5.11 <- matrix(NA, N, 7)
     out.hmp <- cbind(rs, alleles, chrom, pos, hmp5.11, hmp)
     colnames(out.hmp)[5:11] <- c("strand", "assembly", "center", "protLSID", "assayLSID", "panelLSID", "QCcode")
     colnames(out.hmp)[12:ncol(out.hmp)] <- paste0("V", 1:(ncol(out.hmp)-11))
     write.table(out.hmp, "mid1.hmp", col.names <- T, sep <- "\t", row.names <- F, quote <- F)#输出填充前hmp格式文件
     rm(hmp, out.hmp, X)
	 print("subgroup has converted from num format to hapmap format.")
	 tassel_script <- file.path(getwd(), "run_pipeline.pl")
	 system(paste("perl", 
         tassel_script, 
         "-Xms512m -Xmx2g", 
         "-fork1 -h", 
         hapmap_file <- "mid1.hmp", 
         "-export -exportType VCF -runfork1", 
         sep <- " "))
     print("subgroup has converted from hapmap format to VCF format.")
     system("java -Xmx50g -jar beagle.22Jul22.46e.jar gt <- mid1.vcf out <- mid3")
	 print("The file has been completed imputation.")
	 tryCatch({
         lj <- gzfile("mid3.vcf.gz", "rt")
    }, error <- function(e) {
         print("Error occurred while opening the file.")
         system("java -Xmx10g -jar beagle.22Jul22.46e.jar gt <- mid1.vcf out <- mid3")
         lj <- gzfile("mid3.vcf.gz", "rt")
    })
     content <- readLines(lj)  
     close(lj)
     output_file <- "mid3.vcf"
     writeLines(content, output_file)
     vcf_path <- file.path(getwd(), "mid3.vcf")
     snp_path <- file.path(getwd(), "mid")
     command <- paste("perl", 
         "vcf2hpm2_wy_v1.pl", 
         shQuote(vcf_path), 
         shQuote(snp_path), 
         sep <- " ")		 
     system(command)
	 print("subgroup has converted from VCF format to hapmap format.")
     X4 <- data.table::fread("mid2.txt", sep <- " ")
     X5 <- as.matrix(X4[, -c(1:11)])
     X6 <- t(X5)
     index.AA <- X6=="AA"
     index.AG <- X6=="AG"
     index.GG <- X6=="GG"
     nr <- nrow(X6)
     nc <- ncol(X6)
     genotype.n <- matrix(0, nr, nc)
     genotype.n[index.AA] <- 0
     genotype.n[index.AG] <- 1
     genotype.n[index.GG] <- 2
	 XYcluster <- gsub(".txt", ".kbeagle", XYcluster)
	 print(dim(genotype.n))
     write.table(genotype.n, file <- XYcluster, quote <- F, sep <- "\t", col.names <- F, row.names <- F)
     print("subgroup has finished.")
}
KBeagle <- function(data_NA){
	 print("The packages about KBeagle are installing.")
	 library(data.table)
     library(factoextra)
     library(cluster)
     library(ggplot2)
     library(parallel)
	 library(foreach)
     library(doParallel)
	 set.seed(1)
	 if(!file.exists("beagle.22Jul22.46e.jar")){
         system("wget -O beagle.22Jul22.46e.jar https://github.com/99Xingyu-Guo/KBeagle/blob/main/beagle.22Jul22.46e.jar")
    }
	 if(!file.exists("vcf2hpm2_wy_v1.pl")){
         system("wget -O vcf2hpm2_wy_v1.pl https://github.com/99Xingyu-Guo/KBeagle/blob/main/vcf2hpm2_wy_v1.pl")
    }
     
	 print("The data_NA is converting to cluster file format.")
	 data1 <- read.table(data_NA, header = TRUE, sep = "\t")#读取缺失文件
	 data1[data1==0] <- "0"
	 data1[data1==1] <- "1"
	 data1[data1==2] <- "2"
	 data1[is.na(data1)] <- "1"
	 data <- data1[-1]
	 write.table(data, "data_KMeans", quote = FALSE, sep = "\t", col.name = TRUE, row.name = FALSE)
     data <- read.table("data_KMeans", sep = "\t", header = TRUE)
	 print("The data_NA has been converted to cluster file.")

     print("The number of clusters(KN) is determining.")
	 sil_width <- numeric(10)
     for(kn in 2:10){
      km <- kmeans(data,centers=kn,nstart=25)
      ss <- silhouette(km$cluster,dist(data))
      sil_width[kn] <- mean(ss[,3])
    }
     KN <- which.max(sil_width) 
 
	 b1 <- as.data.frame(data1[1])
	 cluster_data <- data.frame(IID = b1, data)
     data.kmeans <- kmeans(data, KN)
     t <- table(cluster_data$IID, data.kmeans$cluster)#展示为列表型，其中第一列为个体ID，之后的列表示个体是否属于该聚类，表示为0、1矩阵
     colnames(t) <- paste("type", 1:KN, sep = "")#将列名表示为“type 1”、“type 2”……“type KN”
     t1 <- cbind(data.frame(b1), apply(t, 2, as.numeric))#将ID列加入到聚类的t文件中，转换文件格式
     colnames(t1)[1] <- "IID"
	
	 print("The cluster sub files are Generating.")
     #生成聚类的文件
	 for(i in 1:KN){
         d <- t1[grep("1", t1[, i+1]), ]
         type <- merge(d, data1, all = FALSE, by = "IID")
         type_result <- type[, -c(2:(KN+1))]
         X <- type_result[, -1]
         fwrite(X, file <- paste("cluster_", i, ".txt", sep = ""), quote = FALSE, sep = "\t", col.name = TRUE, row.name = FALSE)
    }
	 
	 print("The related imputation files are transferring.")
	 for(n in 1:KN){
	     dir_name <- paste0("mid_", n)
		 if(!file.exists(dir_name)){
             dir.create(dir_name)
        }
		 raw_cluster_path = paste(getwd(), "/cluster_", n, ".txt", sep = "")
		 now_cluster_path = paste0(getwd(), "/mid_", n)
		 file.copy(from = raw_cluster_path, to = now_cluster_path)
		 raw_tassel_1 = paste0(getwd(), "/sTASSEL.jar")
		 file.copy(from = raw_tassel_1, to = now_cluster_path)
		 raw_tassel_2 = paste0(getwd(), "/run_pipeline.pl")
		 file.copy(from = raw_tassel_2, to = now_cluster_path)
		 raw_tassel_3 = paste0(getwd(), "/lib")
		 file.copy(from = raw_tassel_3, to = now_cluster_path, recursive = TRUE)
		 raw_perl = paste0(getwd(), "/vcf2hpm2_wy_v1.pl")
		 file.copy(from = raw_perl, to = now_cluster_path)
		 raw_beagle = paste0(getwd(), "/beagle.22Jul22.46e.jar")
		 file.copy(from = raw_beagle, to = now_cluster_path)
	}
     
	 print("The cluster sub files are entering multi-threaded imputation.")
     num_cores <- detectCores()
     cl <- makeCluster(num_cores)
	 registerDoParallel(cl)
     raw_path <- getwd()
     XYparallel <- foreach(j <- 1:KN) %dopar% {
       parallel_path <- paste0(raw_path, "/mid_", j)
       dir.create(parallel_path, showWarnings = FALSE)
       setwd(parallel_path)
       data_cluster <- file.path(getwd(), paste0("cluster_", j, ".txt"))
       result <- Impute.subgroup(XYcluster = data_cluster)
    }
	 stopCluster(cl)
	
	 setwd(raw_path)
     for(n in 1:KN){
	   raw_cluster_path <- paste0(getwd(), "/mid_", n, "/cluster_", n, ".kbeagle")
	   now_cluster_path <- getwd()
	   file.copy(from = raw_cluster_path, to = now_cluster_path)
    }
	 
	 print("The multi-threaded results are integrating.")
	 X.bgl2 <- NULL
     for(i in 1:KN){
        cluster_some <- data.table::fread(paste("cluster_", i, ".kbeagle", sep = ""), sep = "\t", head = TRUE)
        d_1 <- t1[grep("1", t1[, i+1]), ]
        X.bgl <- cbind(d_1[, 1], cluster_some)
        colnames(X.bgl) <- as.character(colnames(data1))
        X.bgl2 <- rbind(X.bgl2, X.bgl)
    }
	 colnames(X.bgl2) <- as.character(colnames(data1))
     X2 <- data.frame(X.bgl2)
	 loc <- match(b1$IID, X2$IID)
	 W <- X2[loc, ]
	 finished_file_path <- paste0(raw_path, "/KBeagle_wc")
     dir.create(finished_file_path, showWarnings = FALSE)
     setwd(finished_file_path)
	 fwrite(W, "KBeagle_finshed.txt", quote = F, sep = "\t", col.name = T, row.name = F)
     print("The multi-threaded results are integrating.")
	 setwd(raw_path)
	 for(n in 1:KN){
	     unlink(paste0(raw_path, "/mid_", n), recursive = TRUE)
		 file.remove(paste0(raw_path, "/cluster_", n, ".txt"))
		 file.remove(paste0(raw_path, "/cluster_", n, ".kbeagle"))
	}
     return("KBeagle_finshed.txt")
}
	
	




