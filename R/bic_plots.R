#' Plot dispersion estimates
#'
#' Plot dispersion estimates for all conditions in experiment. Draw one
#' plot per PDF file.
#' 
#' @param cds          DESeq countsDataSet
#' @param out.dir      Directory in which plots should be saved; Default: $PWD
#' @param file.prefix  String to prepend to PDF file names (optional)
#' @export
bic.plot.dispersion.estimates <- function(cds,out.dir=getwd(),file.prefix=""){
  if(file.prefix != ""){
    file.prefix <- paste(file.prefix,"_",sep="")
  }
  for(cond in ls(cds@fitInfo)){
    pdf(file.path(out.dir,paste(file.prefix,"dispersion_estimates_",cond,".pdf",sep="")))
    plotDispEsts(cds,name=cond)
    dev.off()
  }  
}

#' Draw a heatmap of counts table
#' 
#' Draw heatmap of either raw counts or variance stabilization
#' transformed (?) data for all samples. Write plot to PDF file.
#' 
#' @param cds       DESeq countsDataSet object
#' @param file      PDF file to which heatmap should be saved (optional)
#' @param num.gns   Include this number of most highly expressed genes in 
#'                  the heatmap; Default: 100
#' @param transform logical indicating whether to transform counts; 
#'                  Default: FALSE
#' @export
bic.deseq.heatmap <- function(cds,file=NULL,num.gns=100,transform=FALSE){
  hmcol <- colorRampPalette(brewer.pal(9,"GnBu"))(100)
  dat <- NULL
  if(transform){
    cdsBlind <- estimateDispersions(cds, method="blind")
    vst <- varianceStabilizingTransformation(cdsBlind)
    dat <- exprs(vst)
    dat <- dat[-grep("alignment_not_unique|ambiguous|no_feature|not_aligned|too_low_aQual",
                           rownames(dat)),]
    select <- order(rowMeans(dat),decreasing=TRUE)[1:num.gns]
    dat <- dat[select,]
  } else {
    counts.cds <- counts(cds)
    counts.cds <- counts.cds[-grep("alignment_not_unique|ambiguous|no_feature|not_aligned|too_low_aQual",
                           rownames(counts.cds)),]
    select <- order(rowMeans(counts.cds), decreasing=TRUE)[1:num.gns]
    dat <- counts.cds[select,]
  }
  if(!is.null(file)){
    pdf(file)
  }
  heatmap.2(dat, col=hmcol, trace="none", margin=c(10,10),
            main=paste("Expression of the \n",num.gns," most highly expressed genes",sep=""),
            cexRow=0.6,cexCol=1.0,keysize=1.5,key.title=NA)
  par(cex.main=0.3)
  if(!is.null(file)){
    dev.off()
  }
}

#' Draw a heatmap showing sample-to-sample distances using variance 
#' stabilized data
#' 
#' Find potential sample mislabeling errors by viewing distances
#' between every pair of samples. Save plot to PDF file.
#' 
#' @param cds    DESeq countDataSet
#' @param conds  vector of sample conditions that was used to generate \code{cds}
#' @param file   PDF file to which heatmap whould be saved (optional)
#' @export
bic.sample.to.sample.distances <- function(cds,conds,file=NULL){
  #if(is.null(file)){
  #  file = "heatmap_sample_to_sample_distances.pdf"
  #}
  hmcol = colorRampPalette(brewer.pal(9, "GnBu"))(100)
  cdsBlind <- estimateDispersions(cds,method="blind")
  vst <- varianceStabilizingTransformation(cdsBlind)
  distances <- dist(t(exprs(vst)))
  dist.mat <- as.matrix(distances)
  rownames(dist.mat) = colnames(dist.mat) = with(pData(cdsBlind), paste(rownames(dist.mat), " : " ,conds,sep=""))
  if(!is.null(file)){
    pdf(file)
  }
  heatmap.2(dist.mat, trace="none", col=rev(hmcol), margin=c(13,13), key.title=NA)
  if(!is.null(file)){
    dev.off()
  }
}

#' Plot PCA
#' 
#' Run DESeq \code{plotPCA()} on variance stabilised tranformed data
#' 
#' @param cds  DESeq countDataSet
#' @param file PDF file to which plot should be saved
#' @export
bic.deseq.plot.pca <- function(cds,file=NULL){
  #if(is.null(file)){
  #  file="pca.pdf"
  #}
  cdsBlind <- estimateDispersions(cds,method="blind")
  vst <- varianceStabilizingTransformation(cdsBlind)
  if(!is.null(file)){
    pdf(file)
  }
  plotPCA(vst,intgroup="condition")
  if(!is.null(file)){
    dev.off()
  }
}

#' Plot MA
#'
#' Plot log2 fold changes versus mean of normalized counts for a DESeq
#' comparison
#' 
#' @param res   object returned from DESeq function nbinomTest
#' @param file  PDF file to which plot should be saved (optional)
#' @export
bic.plot.ma <- function(res, file=NULL){
  #if(is.null(file)){
  #  file = "ma_plot.pdf"
  #}
  if(!is.null(file)){
    pdf(file)
  }
  DESeq::plotMA(res)
  if(!is.null(file)){
    dev.off()
  }
}

#' Histogram of p-values
#' 
#' @param dat  data frame with at least one column with name 'pvals'
#' @param file PDF file to which plot should be saved (optional)
#' @export
bic.pval.histogram <- function(dat,file=NULL){
  #if(is.null(file)){
  #  file="pval_histogram.pdf"
  #}
  if(!is.null(file)){
    pdf(file)
  }
  hist(dat$pval, breaks=100, col="blue", border="slateblue", main="",xlab="p-value",ylab="Frequency")
  if(!is.null(file)){
    dev.off()
  }
}

#' Plot distibution of alignment across different genomic locations
#'
#' Bar chart showing for each sample the distribution of read
#' alignments across ribosomal, coding, UTR, intronic and 
#' intergenic bases
#' 
#' @param dat      data frame consisting of data from CollectRNASeqMetrics
#'                 output
#' @param col.pal  name of color palette; must be from list of RColorBrewer palettes
#'                 Default: "Set2"
#' @param file     PDF file to which plot should be saved (optional)
#' @param pct      plot percentages
#'
#' @export
bic.plot.alignment.distribution <- function(dat,pct=FALSE,col.pal="Set2",file=NULL){

  y <- data.frame(Sample = dat$SAMPLE, 
                 Ribosomal = dat$RIBOSOMAL_BASES, 
                 Coding = dat$CODING_BASES,
                 UTR = dat$UTR_BASES,
                 Intronic = dat$INTRONIC_BASES,
                 Intergenic = dat$INTERGENIC_BASES)

  y.m <- melt(y,id.var="Sample")

  position <- "stack"
  y.lbl <- "Total Bases (millions)"
  if(pct){ 
    position <- "fill" 
    y.lbl <- "Percent Bases"
  } 

  ### counts
  p <- ggplot(y.m, aes(x = Sample, y = value/1000000, fill = variable)) + 
    geom_bar(stat="identity", position=position, width=0.7, color="black")+
    theme(axis.text.x = element_text(angle=45,size=9,hjust=1,color="black"),
        legend.position="bottom",
        legend.title = element_blank()
     ) + 
    scale_fill_brewer(palette=col.pal) +      
    #scale_fill_manual(name="Type", values = c3, labels = c("Ribosomal", "Coding", "UTR", "Intronic", "Intergenic")) + 
    labs(title="Alignment Distribution") + 
    xlab("") +  
    ylab(y.lbl)
    if(pct){
      p <- p + scale_y_continuous(labels=percent)
    } 
  if(!is.null(file)){
    pdf(file)
  }
  print(p)
  if(!is.null(file)){
    dev.off()
  }
}

#' Plot coverage uniformity
#' 
#' Bar chart showing median CV coverage, median 5' and 3' bias,
#' and median 5' to 3' bias for each sample.
#'
#' @param dat      data frame consisting of data from CollectRNASeqMetrics
#'                 output
#' @param pct      logical indicating that plot should show percentages
#' @param stacked  logical indicating that bar chart should be stacked; Default: TRUE
#' @param col.pal  name of color palette; must be from list of RColorBrewer palettes
#'                 Default: "Set2"
#' @param file PDF file to which plot should be saved (optional)
#'
#' @export
bic.plot.5prime3prime.bias <- function(dat,pct=FALSE,stacked=TRUE,col.pal="Set2",file=NULL){

  y <- data.frame(Sample = dat$SAMPLE,
                  cvCoverage = dat$MEDIAN_CV_COVERAGE,
                  fivePrimeBias = dat$MEDIAN_5PRIME_BIAS,
                  threePrimeBias = dat$MEDIAN_3PRIME_BIAS,
                  fivetoThreePrimeBias = dat$MEDIAN_5PRIME_TO_3PRIME_BIAS)
  suppressMessages(y.m <- melt(y))

  position <- "stack"
  if(!stacked){
    position <- position_dodge(width=0.7)
  }
  if(pct){
    position <- "fill"
  }

  p <- ggplot(y.m, aes(x = Sample, y = value, fill = variable)) + 
    geom_bar(stat = "identity", position = position, width = 0.7, color = "black") + 
    theme(axis.text.x = element_text(angle=45,size=9,hjust=1,color="black"),
          legend.position="right",
          legend.title = element_blank()
       ) +    
    scale_fill_brewer(palette=col.pal, labels = c("Median CV of Coverage","Median 5\' Bias","Median 3\' Bias","Median 5\' to 3\' Bias")) +
    labs(title="Coverage Uniformity") + 
    xlab("") + 
    ylab("")
  if(pct){
    p <- p + scale_y_continuous(labels=percent)
  }

  if(!is.null(file)){
    pdf(file)
  }
  print(p)
  if(!is.null(file)){
    dev.off()
  }

}

#' Plot normalized coverage
#'
#' Line chart showing normalized coverage for each sample
#'
#' @param dat  data frame containing combined histograms from collectrnaseqmetrics
#' @param col.pal  name of color palette; must be from list of RColorBrewer palettes
#'                 Default: "Set2"
#' @param file PDF file to which plot should be saved (optional)
#' @export
bic.plot.coverage <- function(dat,col.pal="Set2",file=NULL){

  suppressMessages(x.m <- melt(dat, id.vars="position",col.pal=col.pal))
  x.m$position <- as.integer(as.character(x.m$position))
  p <-  ggplot(x.m, aes(x = position, y = value, color = variable)) +  
        geom_line(size=0.7) + 
        theme(legend.position="right",
            legend.text = element_text(size=9),
            legend.title = element_blank()
            ) +
        labs(title="Normalized Coverage") + 
         xlab("Read position") + 
         ylab("") +
         scale_fill_brewer(palette=col.pal) 
         #scale_color_hue(colnames(dat)[-1]) + 
         #guides(colour = guide_legend(override.aes = list(size=5)))
  if(!is.null(file)){
    pdf(file)
  }
  print(p)
  if(!is.null(file)){
    dev.off()
  }
}


#' Plot numbers of mapped and unmapped reads for each sample
#'
#' Stacked bar chart showing either absolute values or percentages of mapped 
#' and unmapped reads for R1 and R2 of each sample. Takes in data frame
#' containing PICARD AlignmentSummaryMetrics, or at minimum a data frame with
#' columns: CATEGORY, PF_READS, PF_READS_ALIGNED, SAMPLE, where CATEGORY column
#' contains at least categories FIRST_OF_PAIR and SECOND_OF_PAIR.
#'
#' @param dat  data frame of PICARD AlignmentSummaryMetrics, or at least 
#'             containing columns CATEGORY, PF_READS, PF_READS_ALIGNED and 
#'             SAMPLE, where CATEGORY contains at least values FIRST_OF_PAIR
#'             and SECOND_OF_PAIR
#' @param pct  logical indicating whether to show percentages rather than absolute
#'             read counts
#' @param col.pal  name of color palette; must be from list of RColorBrewer palettes
#'                 Default: "Set2"
#' @param file PDF file to which plot should be saved (optional)
#'
#' @export
bic.plot.alignment.summary <- function(dat,pct=FALSE,col.pal="Set2",file=NULL){
  dat <- dat[-which(dat$CATEGORY=="PAIR"),]
  dat$UNMAPPED <- dat$PF_READS-dat$PF_READS_ALIGNED
  dat <- dat[,c("CATEGORY","PF_READS","UNMAPPED","SAMPLE")]
  dat$CATEGORY <- revalue(dat$CATEGORY,c("FIRST_OF_PAIR"="R1","SECOND_OF_PAIR"="R2"))  
  colnames(dat) <- c("Category","MappedReads","UnmappedReads","Sample")
  suppressMessages(dat.m <- melt(dat))
  
  position <- "stack"
  cat.lbl.y <- -2
  y.lbl <- "Reads (xMillion)"
  if(pct){ 
    position <- "fill" 
    cat.lbl.y <- -0.04
    y.lbl <- ""
  }

  p <- ggplot(dat.m[which(dat.m$variable!="Total"),], aes(x=Category, y=value/1000000, fill=variable)) + 
   geom_bar(stat="identity", position=position) +
   facet_wrap( ~ Sample, nrow=1, strip.position="bottom") + 
   geom_text(data=dat.m[which(dat.m$variable=="MappedReads"),], mapping=aes(x=Category, y=cat.lbl.y, label=Category), vjust=0) +
   theme(axis.text.x = element_blank(),
      axis.ticks=element_blank(),
      legend.position="right",
      legend.title = element_blank(),
      strip.text.x = element_text(size = 9,color="black",angle=90, hjust=0.5, vjust=1)#,
      #plot.title = element_text(hjust = 0.5)
    ) +
  scale_fill_brewer(palette=col.pal) + 
  labs(title="Alignment Summary") + 
  xlab("") +
  ylab(y.lbl)
  if(pct){
    p = p + scale_y_continuous(labels=percent) +
    ylab(y.lbl)
  } 
  if(!is.null(file)){
    pdf(file)
  }
  print(p)
  if(!is.null(file)){
    dev.off()
  }
}

#' Write PDF file containing heierarchical clustering tree
#'
#' Quickly plot heirarchical clustering of data with several default
#' parameters like width, height, font size, etc.
#' 
#' @param dat            matrix containing data to plot
#' @param file.name      file name. Default: "tmp.pdf"
#' @param title          Title of plot
#' @param sample.labels  Vector of sample labels. Default: column names in matrix
#' @param width          plot width (see plot docs)
#' @param height         plot height (see plot docs)
#' @param lwd            line width
#' @param cex.main       size of title
#' @param cex.lab        size of labels
#' @param cex            font size
#' @param xlab           label of x axis
#' @param ylab           label of y axis
#' @export
bic.pdf.hclust<-function(dat,file.name="tmp.pdf",title="",width=26,height=16,lwd=3,
    cex.main=3,cex.lab=3,cex=3,xlab="",ylab="",sample.labels=NULL)
{
  if(is.null(sample.labels)){
    sample.labels=colnames(dat)
  }
  pdf(file.name,width=width,height=height)
  plot(hclust(dist(t(dat))), lwd=lwd, main=paste(title), cex.main=cex.main,
        xlab=xlab,ylab=ylab, cex.lab=cex.lab, cex=cex, labels=sample.labels)
  dev.off()
}


#' Plot heirachical clustering of all samples
#'
#' Plot clustering of samples using the heierarchical clustering method
#' 
#' @param norm.counts  data matrix containing normalized counts, where
#'                     a column represents a sample and a row represents
#'                     a gene. may or may not contain "GeneID" and "GeneSymbol"
#'                     columns.
#' @param file.name    plot will be saved in this file; Default: 
#'                     $PWD/counts_scaled_hclust.pdf
#' @export
bic.hclust.samples <- function(norm.counts,file.name=NULL){

  if("GeneID" %in% colnames(norm.counts) | 
     "GeneSymbol" %in% colnames(norm.counts)){
    norm.counts <- norm.counts[,-grep("GeneID|GeneSymbol",colnames(norm.counts))]
  }
  norm.counts <- bic.matrix2numeric(norm.counts)

  if(length(colnames(norm.counts)) < 3){
    cat("Less than three samples; can not run cluster analysis\n")
    return()
  }

  if(is.null(file.name)){
    file.name <- "counts_scaled_hclust.pdf"
  }
  
  pseudo <- min(norm.counts[norm.counts > 0])
  counts2hclust <- log2(norm.counts + pseudo)

  sink("/dev/null")
  tryCatch({
      bic.pdf.hclust(counts2hclust,file.name=file.name,title="All counts scaled using DESeq method")
     }, error = function() {
        cat("Error writing pdf file.")
     }, finally = {
       sink()
     }
  )

}

#' Plot MDS clustering of all samples
#'
#' Plot clustering of samples using the MDS method
#' 
#' @param norm.counts  data matrix containing normalized counts, where
#'                     a column represents a sample and a row represents
#'                     a gene
#' @param conds        vector of conditions to be appended to sample IDs
#'                     for labeling
#' @param file         PDF file to which plot should be saved (optional)
#' @param labels       logical indicating whether to include sample labels on plot; 
#'                     Default: TRUE
#' @export
bic.mds.clust.samples <- function(norm.counts,file=NULL,conds=NULL,labels=TRUE){

  if("GeneID" %in% colnames(norm.counts) |
     "GeneSymbol" %in% colnames(norm.counts)){
    norm.counts <- norm.counts[,-grep("GeneID|GeneSymbol",colnames(norm.counts))]
  }

  norm.counts <- bic.matrix2numeric(norm.counts)

  if(length(colnames(norm.counts)) < 3){
    cat("Less than three samples; can not run cluster analysis\n")
    return
  }

  if(is.null(conds)){
    conds=rep('s',length(colnames(norm.counts)))
  }

  pseudo <- min(norm.counts[norm.counts > 0])
  counts2hclust <- log2(norm.counts + pseudo)
  md <- cmdscale(dist(t(counts2hclust)),2)

  if(!is.null(file)){
    pdf(file,width=18,height=12)
  }
  if(labels){
    plot(md, col=as.factor(conds),main="",lwd=2.5,cex=1.5)
    legend("topleft", levels(as.factor(conds)),col=as.factor(levels(as.factor(conds))),pch=1,cex=1.2)
    text(md,colnames(counts2hclust),cex=0.9)
  } else {
    plot(md, col=as.factor(conds),pch=as.numeric(as.factor(conds)),main="",lwd=2.5,cex=1.5)
    legend("topleft", levels(as.factor(conds)),col=as.factor(levels(as.factor(conds))),pch=levels(as.numeric(as.factor(conds))),cex=1.2)
  }
  if(!is.null(file)){
    dev.off()
  }
}

#' Create PDF of generic red/green/black heatmap showing relative
#' expression values
#'
#' Given the results of DESeq comparison of two conditions and the
#' full normalized count matrix, generate a heatmap showing the relative 
#' expression values of "top 100" differentially expressed genes. If
#' there are fewer than 100 genes in the results file, all of them will
#' be included in the heatmap. Color scheme is red/black/green.
#'
#' @param norm.counts.matrix  matrix of normalized counts, where a row is a gene
#'                            and a column is a sample. May contain "GeneID" and/or
#'                            "GeneSymbol" column. If GeneSymbol column is present, 
#'                            it will be used for heatmap labeling. If not, GeneID 
#'                            column will be used.
#' @param condA               the first condition named in the bic.deseq.file name
#' @param condB               the second condition named in the bic.deseq.file name
#' @param genes               a vector of genes to be included in the heatmap;
#'                            if normalized counts matrix includes "GeneSymbol" column,
#'                            genes must be from that column. Otherwise, they must
#'                            be from the "GeneID" column.                      
#' @param file            name of PDF file to which heatmap should be written. (optional) 
#' @export
bic.standard.heatmap <- function(norm.counts.matrix,condA,condB,genes=NULL,file=NULL){

  dat <- norm.counts.matrix
  if ("GeneSymbol" %in% colnames(dat)){
    idHeader <- "GeneSymbol"
  } else {
    idHeader <- "GeneID"
  }
  
  ## extract data for the DE genes
  htmp.dat <- dat[genes,grep(paste(idHeader,condA,condB,sep="|"),colnames(dat))]

  ## remove duplicate genes
  if(length(htmp.dat[which(duplicated(htmp.dat[,idHeader])),idHeader])>0){
    htmp.dat <- htmp.dat[-which(duplicated(as.vector(htmp.dat[,idHeader]))),]
  }

  ## remove any rows that have NAs
  htmp.dat <- htmp.dat[complete.cases(htmp.dat),]
  ## as long as gene symbols are unique, assign them
  ## as rownames; if not, we'll have to average values
  ## for genes that occur multiple times
  rownames(htmp.dat) <- htmp.dat[,idHeader]
  htmp.dat <- htmp.dat[,-1]
  htmp.dat <- bic.matrix2numeric(htmp.dat)

  ## replace any zeros with ones before taking log2
  htmp.dat[htmp.dat==0] <- 1
  htmp.dat <- as.matrix(log2(htmp.dat))

  ## take only the top 100 genes in order to make heatmap legible
  if(dim(htmp.dat)[1]>100){
    htmp.dat <- htmp.dat[1:100,]
  }

  #if(is.null(out.file)){
  #  out.file <- paste("ResDESeq_",condA,"_vs_",condB,"_heatmap.pdf",sep="")
  #}
 
  if(dim(htmp.dat)[1]>1 && dim(htmp.dat)[2]>1){
    ## make heatmap pdf
    if(!is.null(file)){
      pdf(file,width=25,height=16)
    }
    heatmap.2(htmp.dat - apply(htmp.dat, 1, mean), 
              trace='none', 
              col=colorpanel(16,"green","black","red"),
              cexRow=0.9,
              cexCol=1.2, 
              dendrogram="both",
              main=paste("Top Differentially Expressed Genes ",condA," vs ",condB,sep=""), 
              symbreaks=TRUE, 
              keysize=0.5, 
              margin=c(20,10))
    if(!is.null(file)){
      dev.off()
    }
  }
 
}
