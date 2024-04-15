library(GenomicRanges)
library(GenomicFeatures)
install.packages("foreach")
library(foreach)
install.packages("doMC")
library(doMC)

### reading in segmentation file as produced by crllm/cbs package and returns genomic ranges
readSegFile <- function(segFilePath){
  if(!file.exists(segFilePath)){stop("file does not exist!")}
  seg.data<-read.delim(segFilePath)
  paste("chr",seg.data$Chromosome, sep ="")->chr
  seg.data[,2]<-chr
  names(seg.data)[2]<-"chr"
  as.numeric(seg.data$Start.Position)->seg.data$Start.Position
  as.numeric(seg.data$End.Position)->seg.data$End.Position
  colnames(seg.data)[2:4]<-c("chromosome","start","end")
  makeGRangesFromDataFrame(seg.data, keep.extra.columns=TRUE,ignore.strand=TRUE)->grange
  return(grange)
}

##### reads a Gencode DB created using GenomicFeatures and returns a GrangesList where each slot represents a gencode gene
readGencodeDb <- function(gencode.db.path){
  gencode.db <- loadDb(file=gencode.db.path)
  GRList <- transcriptsBy(gencode.db, by = "gene")
  #GRList[match(rownames(assay(vsd)), names(GRList)),]->grange.list
  return(range(GRList))
}

######## convert seg.grange to matrix 
######## generates a list of matrices containing seg.means (average CN ration within that segment), 
######## CN (ratio adjusted to an actual copy number), CNcall (a categorical call AMP, DEL, etc), 
######## num.mark (number of markers i.e. SNP probes within the defining segment)
segMatrix <- function(exp.grange, seg.grange){
  #exp.grange = exp.grange
  #seg.grange = seg.grange
  #### set how many individuals
  levels(factor(seg.grange$Sample)) -> sample.ids
  #### gather the relevant data into lists
  #### have a way of filtering out small segments minimum number of num.markers
  seg.mean.list<-list()
  CN.list<-list()
  CNcall.list<-list()
  num.mark.list<-list()
  for(j in 1:length(sample.ids)){
    seg.mean<-rep(NA,length(exp.grange))
    #CN<-rep(NA,length(exp.grange))
    CNcall<-rep(NA,length(exp.grange))
    num.mark<-rep(NA,length(exp.grange))
    seg.grange[seg.grange$Sample==sample.ids[j],]->temp.seg.grange
    #### find appropriate overlaps
    #### need to deal with multiple overlaps
    findOverlaps(exp.grange, temp.seg.grange)->index
    as.data.frame(index)->index
    seg.mean[index[,1]]<-temp.seg.grange[index[,2]]$Seg.CN
    #CN[index[,1]]<-temp.seg.grange[index[,2]]$CN
    if(length(temp.seg.grange[index[,2]]$CNcall) == 0){
      next
    }
    CNcall[index[,1]]<-as.character(temp.seg.grange[index[,2]]$CNcall)
    num.mark[index[,1]]<-temp.seg.grange[index[,2]]$Num.markers
    #CN->CN.list[[j]]
    CNcall->CNcall.list[[j]]
    num.mark->num.mark.list[[j]]
    seg.mean->seg.mean.list[[j]]
  }
  do.call(cbind, seg.mean.list) -> seg.means.out
  colnames(seg.means.out) <- sample.ids
  rownames(seg.means.out) <- names(exp.grange)
  #do.call(cbind, CN.list) -> CN.out
  #colnames(CN.out) <- sample.ids
  #rownames(CN.out) <- names(exp.grange)
  do.call(cbind, CNcall.list) -> CNcall.out
  colnames(CNcall.out)<-sample.ids
  rownames(CNcall.out)<-names(exp.grange)
  do.call(cbind, num.mark.list)->num.mark.out
  colnames(num.mark.out)<-sample.ids
  rownames(num.mark.out)<-names(exp.grange)
  #output.list<-list(seg.means = seg.means.out, CN = CN.out, CNcall = CNcall.out, num.mark = num.mark.out)
  output.list<-list(seg.means = seg.means.out, CNcall = CNcall.out, num.mark = num.mark.out)
  return(output.list)
}



calculate.normal.vals <- function(matched.cn, matched.exp, cores = 20){
  registerDoMC(cores)
  normal.values <- foreach(i = 1:nrow(matched.cn), .combine = rbind)%dopar%{
    temp.mean <- mean(as.numeric(matched.exp[i,which(matched.cn[i,]=="Normal")]))
    temp.sd <- sd(as.numeric(matched.exp[i,which(matched.cn[i,]=="Normal")]))
    return(c(temp.mean, temp.sd))
  }
  rownames(normal.values) <- rownames(matched.exp)
  colnames(normal.values) <- c("mean", "sd")
  return(normal.values)
}


##### test individual copy number changes against normal equivalents
t.test.cna <- function(matched.cn, matched.exp, min.number.samples){
  
  
  apply(matched.cn,1,function(x){length(which(x=="Focal_Gain"))})>min.number.samples -> retain.focal.gain
  apply(matched.cn,1,function(x){length(which(x=="Focal_Amp"))})>min.number.samples -> retain.focal.amp
  apply(matched.cn,1,function(x){length(which(x=="Broad_Gain"))})>min.number.samples -> retain.broad.gain
  apply(matched.cn,1,function(x){length(which(x=="Focal_Del"))})>min.number.samples -> retain.focal.del
  apply(matched.cn,1,function(x){length(which(x=="Broad_Del"))})>min.number.samples -> retain.broad.del
  
  #focal.gain.values <- foreach(i = which(retain.focal.gain), .combine = rbind)%dopar%{
  #mean(as.numeric(matched.exp[i,which(matched.cn[i,]=="Focal_Gain")]))->temp.mean
  #sd(as.numeric(matched.exp[i,which(matched.cn[i,]=="Focal_Gain")]))->temp.sd
  #return(c(temp.mean,temp.sd))
  #}
  
  #####testing focal gains which have passed above min threshold of frequency
  if (length(which(retain.focal.gain))>1){
    focal.gain.test <- foreach(i = which(retain.focal.gain), .combine = rbind, .errorhandling=c('remove'))%dopar%{
      as.numeric(matched.exp[i,which(matched.cn[i,]=="Focal_Gain")]) -> matched.exp.altered
      as.numeric(matched.exp[i,which(matched.cn[i,]=="Normal")]) -> matched.exp.normal
      if(sd(matched.exp.normal)==0){return(c(NA,NA))}
      t.test(matched.exp.altered,matched.exp.normal) -> temp.test
      return(c(temp.test$statistic,temp.test$p.value))
    }
    rownames(focal.gain.test) <- rownames(matched.exp)[retain.focal.gain]
    colnames(focal.gain.test)<-c("t.stat","p.val")
  }else{focal.gain.test <- NA}
  
  
  if (length(which(retain.focal.amp))>1){
    #####testing amps which have passed above min threshold of frequency
    focal.amp.test <- foreach(i = which(retain.focal.amp), .combine = rbind, .errorhandling=c('remove'))%dopar%{
      as.numeric(matched.exp[i,which(matched.cn[i,]=="Focal_Amp")])->matched.exp.altered
      as.numeric(matched.exp[i,which(matched.cn[i,]=="Normal")])->matched.exp.normal
      if(sd(matched.exp.normal)==0){return(c(NA,NA))}
      t.test(matched.exp.altered,matched.exp.normal)->temp.test
      return(c(temp.test$statistic,temp.test$p.value))
    }
    rownames(focal.amp.test)<-rownames(matched.exp)[retain.focal.amp]
    colnames(focal.amp.test)<-c("t.stat","p.val")
  }else{focal.amp.test <- NA}
  
  
  if (length(which(retain.broad.gain))>1){
    #####testing broad gains which have passed above min threshold of frequency
    broad.gain.test <- foreach(i = which(retain.broad.gain), .combine = rbind, .errorhandling=c('remove'))%dopar%{
      as.numeric(matched.exp[i,which(matched.cn[i,]=="Broad_Gain")])->matched.exp.altered
      as.numeric(matched.exp[i,which(matched.cn[i,]=="Normal")])->matched.exp.normal
      if(sd(matched.exp.normal)==0){return(c(NA,NA))}
      t.test(matched.exp.altered,matched.exp.normal)->temp.test
      return(c(temp.test$statistic,temp.test$p.value))
    }
    rownames(broad.gain.test)<-rownames(matched.exp)[retain.broad.gain]
    colnames(broad.gain.test)<-c("t.stat","p.val")
  }else{broad.gain.test <- NA}
  
  
  if (length(which(retain.broad.del))>1){
    #####testing broad del which have passed above min threshold of frequency
    broad.del.test <- foreach(i = which(retain.broad.del), .combine = rbind, .errorhandling=c('remove'))%dopar%{
      as.numeric(matched.exp[i,which(matched.cn[i,]=="Broad_Del")])->matched.exp.altered
      as.numeric(matched.exp[i,which(matched.cn[i,]=="Normal")])->matched.exp.normal
      if(sd(matched.exp.normal)==0){return(c(NA,NA))}
      t.test(matched.exp.altered,matched.exp.normal)->temp.test
      return(c(temp.test$statistic,temp.test$p.value))
    }
    rownames(broad.del.test)<-rownames(matched.exp)[retain.broad.del]
    colnames(broad.del.test)<-c("t.stat","p.val")
  }else{broad.del.test <- NA}
  
  if (length(which(retain.focal.del))>1){
    #####testing focal.dels which have passed above min threshold of frequency
    focal.del.test <- foreach(i = which(retain.focal.del), .combine = rbind, .errorhandling=c('remove'))%dopar%{
      as.numeric(matched.exp[i,which(matched.cn[i,]=="Focal_Del")])->matched.exp.altered
      as.numeric(matched.exp[i,which(matched.cn[i,]=="Normal")])->matched.exp.normal
      if(sd(matched.exp.normal)==0){return(c(NA,NA))}
      t.test(matched.exp.altered,matched.exp.normal)->temp.test
      return(c(temp.test$statistic,temp.test$p.value))
    }
    rownames(focal.del.test)<-rownames(matched.exp)[retain.focal.del]
    colnames(focal.del.test)<-c("t.stat","p.val")
  }else{focal.del.test <- NA}
  
  list("focal.gain" = focal.gain.test, "focal.amp"= focal.amp.test, 
       "broad.gain" = broad.gain.test, "broad.del" = broad.del.test, 
       "focal.del" = focal.del.test) -> tests
  
  # add BH P.adjust
  tests <- foreach(i = 1:length(tests))%do%{
    if(is.na(tests)[i]){
      return(NA)
    }else{
      return(cbind(tests[[i]],padj = p.adjust(tests[[i]][,2], method = "BH")))
    }
  }
  
  return(tests)
}

chromosomal.plot<-function(normal.values, exp , exp.grange, sample.number, matched.seg.means, matched.cn, output.file = "chromosomes.pdf", span = 0.05){
  
  as.character(seqnames(exp.grange))->chromosome.locs
  as.numeric(start(exp.grange))->start.locs
  levels(as.factor(chromosome.locs))->all.chrom
  all.chrom <- all.chrom[all.chrom!="chrM"&all.chrom!="chrX"&all.chrom!="chrY"]
  pdf(file = output.file)
  for(i in c("chr1","chr2","chr3","chr4","chr5","chr6","chr7","chr8","chr9","chr10","chr11","chr12","chr13","chr14","chr15","chr16","chr17","chr18","chr19","chr20","chr21","chr22")){
    log(exp[,sample.number]/normal.values[,1],2) -> sample.exp.ratio
    sample.exp.ratio[chromosome.locs==i] -> sample.exp.ratio
    start.locs[chromosome.locs==i] -> start.locs.temp
    matched.seg.means[chromosome.locs==i,sample.number]->matched.seg.means.temp
    matched.cn[chromosome.locs==i,sample.number]->matched.cn.temp
    
    
    col.exp<-vector()
    for(j in 1:length(sample.exp.ratio)){
      col.exp[j] <- ifelse(sample.exp.ratio[j]>1|sample.exp.ratio[j]< -1,rgb(1,0,0,1), ifelse(sample.exp.ratio[j]>0.5|sample.exp.ratio[j]< -0.5,rgb(1,0,0,0.5), rgb(1,0,0,0.3)))
    }


    
    plot(sample.exp.ratio,start.locs.temp, pch = 20, cex = 0.2, col = col.exp, main = paste(colnames(exp)[sample.number],i), xlim = c(-5, 5), ylab = "bp", xlab = "ratio individual/normals")
    col.cn<-vector()
    for(j in 1:length(sample.exp.ratio)){
      col.cn[j] <- ifelse(matched.cn.temp[j]=="Focal_Del","#2739D7", ifelse(matched.cn.temp[j]=="Focal_Gain","#E41C1C", ifelse(matched.cn.temp[j]=="Broad_Gain","#E41C1C", ifelse(matched.cn.temp[j]=="Broad_Del","#2739D7", "#898888"))))
    }
    
    
    points(matched.seg.means.temp*10,start.locs.temp, pch = 20, cex = 0.2, col = col.cn)
    
    x.loess <- loess(x ~ y, span=span, data.frame(x=sample.exp.ratio[order(start.locs.temp)], y=start.locs.temp[order(start.locs.temp)]))
    x.predict <- predict(x.loess, data.frame(y=start.locs.temp[order(start.locs.temp)]))
    #lines(x.loess, col="red", cex = 0.5)
    lines(x.predict,y=start.locs.temp[order(start.locs.temp)], col = "purple")
    #lines(lowess(x=sample.exp.ratio[order(start.locs.temp)],y=start.locs.temp[order(start.locs.temp)], f= 2/3, delta = 10000)[[1]], start.locs.temp[order(start.locs.temp)], col="red")
    #lines(lowess(x=sample.exp.ratio[order(start.locs.temp)] f= 2/3, delta = 10000)[[1]], start.locs.temp[order(start.locs.temp)], col="red")
    abline(v = 0)
    
  }
  dev.off()
}


draw.cn.dep.plot <- function(exp.grange, all.tests, chr = "chr17", title = "Chromsome 17, Group 4 MB"){
  ### set gene characteristics, specifically chrosome and start locations and genes
  names(exp.grange) -> gene.names
  as.character(seqnames(exp.grange)) -> chromosome.locs
  as.numeric(start(exp.grange)) -> start.locs
  levels(as.factor(chromosome.locs)) -> all.chrom
  locs <- cbind(chromosome.locs,start.locs)
  rownames(locs) <- gene.names
  
  ### extract the relevant category of test and extract only the results which relate to the chromosome of interest
  all.tests[[1]] -> focal.gain.test
  cbind(locs[rownames(focal.gain.test),],focal.gain.test) -> focal.gain.test
  focal.gain.test[focal.gain.test[,1]==chr,] -> select.focal.gain.test
  
  
  all.tests[[3]] -> broad.gain.test
  cbind(locs[rownames(broad.gain.test),],broad.gain.test) -> broad.gain.test
  broad.gain.test[broad.gain.test[,1]==chr,] -> select.broad.gain.test
  
  all.tests[[2]] -> focal.amp.test
  cbind(locs[rownames(focal.amp.test),],focal.amp.test) -> focal.amp.test
  focal.amp.test[focal.amp.test[,1]==chr,] -> select.focal.amp.test
  
  all.tests[[4]] -> broad.del.test
  cbind(locs[rownames(broad.del.test),],broad.del.test) -> broad.del.test
  broad.del.test[broad.del.test[,1]==chr,] -> select.broad.del.test
  
  all.tests[[5]] -> focal.del.test
  cbind(locs[rownames(focal.del.test),],focal.del.test) -> focal.del.test
  focal.del.test[focal.del.test[,1]==chr,] -> select.focal.del.test
  

 ### compile all the results in list 
  select.test.list <- list(
    broad.del.test=select.broad.del.test,
    broad.gain.test=select.broad.gain.test,
    focal.amp.test=select.focal.amp.test,
    focal.del.test=select.focal.del.test, 
    focal.gain.test=select.focal.gain.test)

  ### check what is the maximum bp size needed on this chromsome
unlist(lapply(select.test.list,class)) -> idx
as.numeric(max(unlist(lapply(select.test.list[idx=="matrix"], function(x){max(x[,2],na.rm = T)})), na.rm = T)) -> x.max
  
  ## set colour according to significance threshold and plot
if(class(select.focal.gain.test)=="matrix"){
  color <- ifelse(as.numeric(select.focal.gain.test[,5])<0.05,"purple","dodgerblue")
  size <- ifelse(as.numeric(select.focal.gain.test[,5])<0.05,1,0.5)
  plot(c(1:10),c(1:10),type = "n", main = title, pch = 19, xlim = c(1,x.max),ylim=c(-20,21), xlab = "bp",ylab="CN dosage dependency")
  points(as.numeric(select.focal.gain.test[,2]),as.numeric(select.focal.gain.test[,3]), pch = 19, col = color, cex=size)
  abline(h=0)
  select.focal.gain.test[which(color=="purple"),]
}
if(class(select.broad.gain.test)=="matrix"){
  color <- ifelse(as.numeric(select.broad.gain.test[,5])<0.05,"orange","yellow")
  size <- ifelse(as.numeric(select.broad.gain.test[,5])<0.05,1,0.5)
  points(as.numeric(select.broad.gain.test[,2]),as.numeric(select.broad.gain.test[,3]), pch = 19, col = color, cex=size)
  abline(h=0)
  select.broad.gain.test[which(color=="orange"),]
}
if(class(select.focal.amp.test)=="matrix"){
  color <- ifelse(as.numeric(select.focal.amp.test[,5])<0.05,"darkred","red")
  size <- ifelse(as.numeric(select.focal.amp.test[,5])<0.05,1,0.5)
  points(as.numeric(select.focal.amp.test[,2]),as.numeric(select.focal.amp.test[,3]), pch = 19, col = color, cex=size)
  select.focal.amp.test[which(color=="darkred"),]
 } 
if(class(select.broad.del.test)=="matrix"){
  color <- ifelse(as.numeric(select.broad.del.test[,5])<0.05,"darkgreen","green")
  size <- ifelse(as.numeric(select.broad.del.test[,5])<0.05,1,0.5)
  points(as.numeric(select.broad.del.test[,2]),as.numeric(select.broad.del.test[,3]), pch = 19, col = color, cex=size)
  select.broad.del.test[which(color=="darkgreen"),]
  }
if(class(select.focal.del.test)=="matrix"){
  color <- ifelse(as.numeric(select.focal.del.test[,5])<0.05,"brown","beige")
  size <- ifelse(as.numeric(select.focal.del.test[,5])<0.05,1,0.5)
  points(as.numeric(select.focal.del.test[,2]),as.numeric(select.focal.del.test[,3]), pch = 19, col = color, cex=size)
  select.focal.del.test[which(color=="brown"),]
  }
  #### add legend
  legend("topright",c("focal amp - Significant","focal amp - Non-Significant","focal gain - Significant","focal gain - Non-Significant","broad gain - Significant","broad gain - Non-Significant","broad del - Significant","broad del - Non-Significant","focal del - Significant","focal del - Non-Significant"), pch = 19, col = c("darkred","red","purple","dodgerblue","orange","yellow","darkgreen","green","brown","beige"),bty="n")
return(select.test.list)
  }
