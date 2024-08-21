#!/usr/bin/env Rscript
Sys.setenv('R_MAX_VSIZE'=64000000000)
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(IRanges))
tstmp<-format(Sys.time(), "%b %d %X")

args = commandArgs(trailingOnly=TRUE)

model <- args[1]
if(!file.exists(model)){
    stop("Model '",model,"' does not exists.\n")
}
outdir<- args[2]
fnum  <- args[3]
if(grepl('^ *[0-9]+ *',fnum)){
    fnum<- as.numeric(fnum)
}else{
    stop('Number of chunks should be integer')
}

chunkSize <- 10
cat(format(Sys.time(), "%b %d %X"),'model=',model,'\n','outdir=',outdir,'\n')
missingRanges<-function(fl.,chunk=fnum){
	ir<-IRanges(start=which(!1:chunk %in% unique(fl.)),width=1)
	rir<-reduce(ir)
	sir<-rir[width(rir)==1]
	lir<-rir[width(rir)>1]
	if(length(lir)>0){
		str<-c(paste0(start(sir)),paste0(start(lir),'-',end(lir)))
	}else{
		str<-paste0(start(sir))
	}
	res<-paste0(str,collapse=',')
	return(res)
}

mn<-sub('^.+/([^/]+).(xml|mat|yml).*$','\\1',model)
cat(format(Sys.time(), "%b %d %X"),'mn=',mn,'\n')

if(!dir.exists(file.path(outdir,mn))){
	dir.create(file.path(outdir,mn))
}
#fl<-dir('/flash/GoryaninU/anatoly/FBA_sensitivity/carveme_sobol',pattern='fba_sobol_fluxes.*.csv',recursive=TRUE)
fl<-dir(file.path(outdir,mn),pattern='fba_sobol_fluxes.*.csv',recursive=TRUE)

if(length(fl)==0){
  sbl<-paste0('sbatch --array=1-',fnum,' calc_FBA_ecoSobol_dump.sh ',model,' ',file.path(outdir,mn))
	cat(format(Sys.time(), "%b %d %X"),sbl,'\n')
	res<-system(sbl,intern=TRUE)
	if(grepl('Submitted batch job *',res)){
		resbl<-paste0('sbatch --depend=afterany:',sub('Submitted batch job *','',res),' makeRanges64.sh ',model,' ',outdir,' ',fnum)
		cat(format(Sys.time(), "%b %d %X"),resbl,'\n')
		system(resbl)
	}
  q(save='no')
}

#fl<-dir('/flash/GoryaninU/anatoly/FBA_sensitivity/carveme_sobol',pattern='fba_sobol_fluxes.*.csv',recursive=TRUE)
#fl<-dir(file.path(outdir,mn),pattern='fba_sobol_fluxes.*.csv',recursive=TRUE)
fldf<-data.frame(folder=fl,num=as.numeric(sub('fba_sobol_fluxes_','',sub('_8192.csv.*','',fl))))
fldf$Len<-NA
fldf$Len[grep('.*gz$',fl)]<-8192
cat(format(Sys.time(), "%b %d %X"),'len(fl)=',length(fl),'dim(fldf)=',dim(fldf),'\n')

for(i in grep('.*gz$',fl,invert=TRUE)){
cat(format(Sys.time(), "%b %d %X"),i,fl[i],'\n')
d<-fread(file.path(outdir,mn,fl[i]))
fldf$Len[i]<-dim(d)[1]
cat(format(Sys.time(), "%b %d %X"),i,dim(d)[1],'Done.\n')
}

df<- fldf[(fldf$Len==8192),]
#save.image(paste0(mn,'_fldf.Rdata'))

cat(format(Sys.time(), "%b %d %X"),'dim(df)=',dim(df),'unique(df$num)=',unique(df$num),'\n')
doneFlag<-TRUE
if(any(!1:fnum %in% unique(df$num))){
	doneFlag<-FALSE
	l<-missingRanges(df$num)
	sbl<-paste0('sbatch --array=',l,' calc_FBA_ecoSobol_dump.sh ',model,' ',file.path(outdir,mn))
	cat(format(Sys.time(), "%b %d %X"),sbl,'\n')
	res<-system(sbl,intern=TRUE)
	if(grepl('Submitted batch job *',res)){
		resbl<-paste0('sbatch --depend=afterany:',sub('Submitted batch job *','',res),' makeRanges64.sh ',model,' ',outdir,' ',fnum)
		cat(format(Sys.time(), "%b %d %X"),resbl,'\n')
		system(resbl)
	}
}

for(i in grep('.*gz$',df$folder,invert=TRUE)){
 comline<-paste0('gzip ',file.path(outdir,mn,df$folder[i]))
cat(format(Sys.time(), "%b %d %X"),i,comline,'\n')
 system(comline)
comline1<-sub('fba_sobol_fluxes','fba_sobol_solution',comline)
cat(format(Sys.time(), "%b %d %X"),i,comline1,'\n')
 system(comline1)
cat(format(Sys.time(), "%b %d %X"),i,'Done.\n')
}

if(doneFlag){
    ddir<-file.path(outdir,mn)
    odir<-file.path(ddir,'sobol_data')
    dt<-fread(paste0(ddir,'/','fba_sobol_solution_1_8192.csv.gz'))
    fcol <- grep('Solution',names(dt))
    fr<-fread(paste0(ddir,'/','flux_range_',fnum,'.csv.gz'))
    idx<-grep('^.+_fl$',fr$varname)
    zdx<-which((fr$card[idx]>1) & (fr$range[idx]>1e-5))
    ldx<-which(fr$varname %in% sub('_fl$','_l',fr$varname[idx[zdx]]))
    udx<-which(fr$varname %in% sub('_fl$','_u',fr$varname[idx[zdx]]))
    cnm<-c(fr$varname[ldx],fr$varname[udx])
    narr <- ceiling((1+length(cnm))/chunkSize)
    mems <-256
    if(fnum<=16){
        mems <- 32
    }else if(fnum<=50){
        mems <- 64
    }else if(fnum <= 90){
        mems <- 128
    }
    #resfl<-dir(ddir,pattern=paste0('sobol_data_.+',fnum))
    resfl<-dir(ddir,pattern=paste0('sobol_data_nz_.+',fnum))
    if(length(resfl) < narr){
        num<-as.numeric(sub(paste0('_10_',fnum,'.rds$'),'',sub('^sobol_data_','',resfl)))
        l<-missingRanges(num,narr)
        sbl<-paste0('sbatch --mem=',mems,'G --array=',l,
                    ' prcc_j_sobol_NZrange10.sh ',ddir,' ',odir,' ',fcol,' ',fnum)
#                    ' prcc_j_sobol_N10.sh ',ddir,' ',odir,' ',fcol,' ',fnum)
        cat(format(Sys.time(), "%b %d %X"),sbl,'\n')
        res<-system(sbl,intern=TRUE)
        cat(format(Sys.time(), "%b %d %X"),res,'\n')
    }
}
