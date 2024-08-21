#!/usr/bin/env Rscript
###############################################################################
# script calculates PRCC coefficients only for the subset of boundaries,
# which control reactions with non-zero flux range, i.e. cardinality>1 and
# range > 1e-5.
#------------------------------------------------------------------------------
Sys.setenv('R_MAX_VSIZE'=64000000000)
suppressPackageStartupMessages(library(data.table))
source('prcc_columnwise.R')
args = commandArgs(trailingOnly=TRUE)
oop<-options(boot.parallel="multicore",boot.ncpus=8L)

infname <- args[1]
outfname <- args[2]
fitCol <- as.numeric(args[3])
chunkSize <- as.numeric(args[4])
chunk <- as.numeric(args[5])
fnum <- as.numeric(args[6])

cat(format(Sys.time(), "%b %d %X"),'\ninfname',infname,'\n')
cat('outfname',outfname,'\n')
cat('fitCol',fitCol,'chunkSize',chunkSize,'chunk',chunk,'N files:',fnum,'\n')
outfname<-paste0(outfname,'_nz_',chunk,'_',chunkSize,'_',fnum,'.rds')
if(file.exists(outfname)){
    cat('File exists',outfname,'\n')
    stop('File exists',outfname,'\n')
}else{
    l<-list()
    saveRDS(l,file=outfname)
    frname <- paste0(infname,'/flux_range_',fnum,'.csv.gz')
    if(!file.exists(frname)){
        cat('File does not exists',frname,'\n')
        stop('File does not exists',frname,'\n')
    }
    fr<-fread(frname)
    idx<-grep('^.+_fl$',fr$varname)
    zdx<-which((fr$card[idx]>1) & (fr$range[idx]>1e-5))
    ldx<-which(fr$varname %in% sub('_fl$','_l',fr$varname[idx[zdx]]))
    udx<-which(fr$varname %in% sub('_fl$','_u',fr$varname[idx[zdx]]))
    cnm<-c(fr$varname[ldx],fr$varname[udx])
    cat("Prepare: length(cnm)",length(cnm),'\n')
    if((chunkSize*(chunk-1)+1) > length(cnm)+1){
        cat("Start position:",(chunkSize*(chunk-1)+1),' > Number of parameters:',(length(cnm)+1),'\n')
        file.remove(outfname)
        stop("Start position:",(chunkSize*(chunk-1)+1),' > Number of parameters:',(length(cnm)+1),'\n')
    }
    #fl<-dir(infname,pattern='*.csv')
    #fl<-paste0('fba_sobol_',1:fnum,'_8192.csv')
    fl<-paste0('fba_sobol_solution_',1:fnum,'_8192.csv.gz')


    for(f in fl){
        dt<-fread(paste0(infname,'/',f))
        cat(f,"dim(dt)",dim(dt),'\n')
        l[[f]]<-dt
    }
    dt<-do.call(rbind,l)
    rm(l)
    df<-as.data.frame(dt)
    rm(dt)
    cat("Final: dim(df)",dim(df),'\n')
    Y<-df[,fitCol]
    X<-df[,names(df) %in% cnm]
    orignames<-names(X)
    rm(df)
    gc(verbose = TRUE, full = TRUE)
    #names(X)<-paste0('C',1:dim(X)[2])
    cat("Final: dim(X)",dim(X),'\n')

    N<-(dim(X)[2])
    L<-(dim(X)[1])
    Nchunks<-ceiling(N/chunkSize)
    p<- N-1
    js<-(chunkSize*(chunk-1)+1):min(chunkSize*chunk,N+1)
    cat(format(Sys.time(), "%b %d %X"),'js=[',range(js),']\n')

    getPRCC<-function(.x){
        cat(format(Sys.time(), "%b %d %X"),'Start:',.x,'\n')
        if(.x<=N){
            r1<-pccJ(X,Y,j=.x,rank=TRUE,nboot=0)$PRCC
            names(r1)[1]<-'original'
            s<-significancePVal(r1$original,L,p)
            r1$T<-s$T
            r1$pval<-s$pval
            r1$rownum<-L
            cat(format(Sys.time(), "%b %d %X"),.x,'r1:',as.character(r1),'\n')
        }else{
            r1<-data.frame(original=0, varname='', varnum=1, T=0.1, pval=0.1, rownum=1)[FALSE,]
        }
        return(r1)
    }
    out<-lapply(js,getPRCC)
    saveRDS(out,file=outfname)
    options(oop)
    cat(format(Sys.time(), "%b %d %X"),'Done.\n')
}
