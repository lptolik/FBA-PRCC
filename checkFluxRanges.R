#!/usr/bin/env Rscript
Sys.setenv('R_MAX_VSIZE'=64000000000)
suppressPackageStartupMessages(library(data.table))
args = commandArgs(trailingOnly=TRUE)

resfn <- 'flux_range.csv'
fnum <- 16 #64

mdir <- args[1] #Sobol results repository /flash/GoryaninU/anatoly/FBA_sensitivity/carveme_sobol
cat(format(Sys.time(), "%b %d %X"),'mdir=',mdir,'\n')
oodir <- args[2] #Sobol results repository /flash/GoryaninU/anatoly/FBA_sensitivity/carveme_sobol
cat(format(Sys.time(), "%b %d %X"),'oodir=',oodir,'\n')
dirN <- args[3] #Position of the folder to process
cat(format(Sys.time(), "%b %d %X"),'dirN=',dirN,'\n')
if(grepl('^ *[0-9]+ *',dirN)){
    dirN<-as.numeric(dirN)
}else{
    stop("Folder position should be integer.")
}
fnum <-args[4] #Number of results files to combine for analysis
cat(format(Sys.time(), "%b %d %X"),'fnum=',fnum,'\n')
if(grepl('^ *[0-9]+ *',fnum)){
    fnum<- as.numeric(fnum)
}else{
    stop('Aggregation size should be integer')
}
resfn <- paste0('flux_range_',fnum,'.csv.gz')

dl<-list.dirs(mdir,recursive=FALSE)
dl<-sort(dl)
cat(format(Sys.time(), "%b %d %X"),'len(dl)=',length(dl),'\n')
while(dirN<=length(dl)){
    cat(format(Sys.time(), "%b %d %X"),'dirN=',dirN,'\n')
    infname<-dl[dirN]
    cat(format(Sys.time(), "%b %d %X"),'infname=',infname,'\n')
    ddir<-file.path(infname)
    oddir<-sub(mdir,oodir,infname)
    if(!dir.exists(oddir)){
        dir.create(oddir,showWarnings = FALSE,recursive = TRUE)
    }
    odir<-file.path(oddir,resfn)
    cat(format(Sys.time(), "%b %d %X"),'ddir=',ddir,'\n')
    cat(format(Sys.time(), "%b %d %X"),'odir=',odir,'\n')
    dirN<-dirN+1
    if(file.exists(odir) ){
        cat(format(Sys.time(), "%b %d %X"),'Folder skipped ',infname,':',odir,' exists.',(dirN-1),'\n')
        next
    }else{
        fl<-paste0('fba_sobol_fluxes_',1:fnum,'_8192.csv.gz')
        if(file.exists(paste0(infname,'/',fl[fnum]))){
            system(paste('touch ',odir))
            l<-list()
            for(f in fl){
                dt<-fread(paste0(infname,'/',f))
                dt<-dt[is.finite(Solution)]
                cat(f,"dim(dt)",dim(dt),'\n')
                l[[f]]<-dt
            }
            dt<-do.call(rbind,l)
            cat(format(Sys.time(), "%b %d %X"),'dim(dt)=',dim(dt),'\n')
            res<-as.data.frame(t(sapply(dt,
                                        function(.x){c(min=min(.x),
                                                       q1=quantile(.x,0.25),
                                                       median=median(.x),
                                                       q3=quantile(.x,0.75),
                                                       max=max(.x),
                                                       sd=sd(.x),
                                                       range=diff(range(.x)),
                                                       card=length(table(.x)))})))
            res$varname<-rownames(res)
            cat(format(Sys.time(), "%b %d %X"),'dim(res)=',dim(res),'\n')
            fwrite(res,file=odir)
        }else{
            cat(format(Sys.time(), "%b %d %X"),"File ",
                paste0(infname,'/',fl[fnum]),' does not exists\n')
        }
    }
    resbl<-paste0('sbatch checkFluxRanges.sh ',mdir,' ',oodir,' ',(dirN),' ',fnum)
    cat(format(Sys.time(), "%b %d %X"),resbl,'\n')
    #system(resbl)
}
