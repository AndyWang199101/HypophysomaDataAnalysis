#### infer CNV of single cells from RNA-seq
setwd( "~/Documents/projects/hypophysoma/" )

### read in genome annotation
options(stringsAsFactors = F)
gencode <- read.table( "gencode.v27.annotation.gtf",header = F,sep = '\t',skip = 5 )
#colnames(ref.flat) <- c( "geneName","name","chrom","strand","txStrat","txEnd","cdsStart","cdsEnd","exonCounts","exonStrats","exonEnds" )
colnames(gencode) <- c( "chrom","source","feature","start","end","score","strand","frame","attribute" )
idx <- which(gencode$feature == "gene")
gencode.genes <- gencode[idx,]

genes.names <- unlist( lapply(
    gencode.genes$attribute
    , function(i){
        tmp <- strsplit(i,split = ';')[[1]][3]
        strsplit(tmp,split = ' ')[[1]][3]
    }
))
gencode.genes <- cbind( geneName=genes.names,gencode.genes )
rm( gencode )

##### split genes by chromosomes
## assume this file has ordered genes iby their start position
chrs <- unlist(lapply( c(1:22,'X','Y'),function(i){paste0('chr',i)} ))
chrs_genes <- lapply(chrs, function(i){ 
        idx = which(gencode.genes$chrom==i)
        data.frame(
            chrom = array( i,dim = length(idx) ),
            geneName = gencode.genes$geneName[idx],
            Start = gencode.genes$start[idx],
            End = gencode.genes$end[idx]
        )
    })

##### read in gene expression data
## Methods from cell paper: Single-cell transcriptimoc analysis of primary and metastatic tumor ecosystems in head and neck cancer, 2017
##    Take log2(TPM_like +1 ) as E(xpression)
##    Initial CNVs (CNV0) were estimated by sorting the analyzed genes by their chromosomal location and applying a moving average to
#   the relative expression values, with a sliding window of 100 genes within each chromosome.
#   - relative expression: Er = E-E_a (scale by gene)  

## Analyze tumor data
load( "codes/tpm.rda" )
# TPM <- 2^(log_TPM)-1
# E_TPM <- apply( TPM, MARGIN = 1, mean)
# E_LOG <- log2(E_TPM+1)
Er.genes <- t(scale( t(log_TPM),center = T,scale = F ))

genes.names <- read.table( file="filtered_gene_bc_matrices/hg19/genes.tsv",head=F,sep = '\t')
emsg <- rownames(Er.genes)
gene.symbols <- c()
for(g in emsg){
    idx <- which(genes.names$V1 == g)
    gene.symbols <- c(gene.symbols,genes.names$V2[idx])
}
rownames( Er.genes )<- gene.symbols

# read in gene names

# ceiling expression between [-3,3]
Er.genes[Er.genes>3] <- 3
Er.genes[Er.genes < -3] <- -3

## moving agerage

# delete genes not measured
count.windows <- c()
k = 100
for( i in 1:length(chrs_genes) ){
    tmp <- chrs_genes[[i]]
    genes.in.chr <- tmp$geneName
    idx <- which( genes.in.chr %in% gene.symbols )
    tmp <- tmp[idx,]
    chrs_genes[[i]] <- tmp
    count.windows <- c(count.windows,length(idx)-k+1)
}

num.of.cells <- dim(Er.genes)[2]
#CNV0 <- matrix(0,nrow = dim(Er.genes)[2],ncol = sum(count.windows))

CNV0.chr <- list()
for( chr in 1:length(chrs_genes) ){
    genes.chr <- chrs_genes[[chr]]$geneName
    num.of.genes <- length(genes.chr)
    expr.chr <- Er.genes[genes.chr,]
    mov.avg <- apply(
        expr.chr,2,function(expr.vec){
            temp <- sum( expr.vec[1:k] )
            res <- c(temp/k)
            for( j in 2:(num.of.genes-k+1) ){
                temp <- temp - expr.vec[j-1] + expr.vec[j+k-1]
                res <- c(res,temp/k)
            }
            res
        }
    )
    CNV0.chr[[chr]] <- mov.avg
}
save(CNV0.chr,file = "CNV0.rda" )

#load( "CNV0.rda" )

CNV0 <- c()
for( chr in 1:(length(chrs_genes)-1) ){
    CNV0 <- rbind(CNV0,CNV0.chr[[chr]])
}

CNV0[which(is.na(CNV0))] <- 0
## check magnilant cells
cell.score <- apply( CNV0,2,mean )
avg.score <- rowMeans(CNV0)
cell.cor <- apply( CNV0,2,function(vec){cor(vec,avg.score)} )

cell.order <- sort( cell.score,decreasing = T,index.return = T )
cell.order.idx <- cell.order$ix

cell.CNV0.matrix <- CNV0[,cell.order.idx]

data.plot <- data.frame(
    cnv.score = cell.score,
    corr = cell.cor
)
library(ggplot2)
p <- ggplot( data.plot,aes(x=cnv.score,y=cell.cor) )
p +geom_point()+labs( x="CNV.score",y="Cor.with.Avg" )+theme_bw()+
    theme(plot.title=element_text(size=20),
          axis.title.y=element_text(size = 25, vjust=+0.2),
          axis.title.x=element_text(size = 25, vjust=-0.2),
          axis.text.y=element_text(size = 25),
          axis.text.x=element_text(size = 25),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank())
#dev.off()

