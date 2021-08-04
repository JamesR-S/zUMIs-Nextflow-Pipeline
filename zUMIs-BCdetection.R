#!/usr/bin/env Rscript
library(methods)
library(data.table)
library(yaml)
library(ggplot2)

##########################
suppressMessages(require(R.utils))
opt<-commandArgs(trailingOnly = T,asValues=T)


setwd(opt$out_dir)
source(paste0(opt$zUMIs_directory,"/barcodeIDFUN.R"))
options(datatable.fread.input.cmd.message=FALSE)
data.table::setDTthreads(threads=opt$num_threads)
if(!is.null(opt$barcode_sharing)){
  if(opt$barcode_sharing == ""){
    opt$barcode_sharing <- NULL
  }
}

#######################################################################
#######################################################################
##### Barcode handling & chunking

#read file with barcodecounts
# bc is the vector of barcodes to keep
bccount<-cellBC(bcfile      = opt$barcode_file,
           bcnum       = opt$barcode_num,
           bcauto      = opt$automatic,
           bccount_file= paste0(opt$out_dir,"/", opt$project, ".BCstats.txt"),
           outfilename = paste0(opt$out_dir,"/zUMIs_output/stats/",opt$project,".detected_cells.pdf"))

fwrite(bccount,file=paste0(opt$out_dir,"/zUMIs_output/",opt$project,"kept_barcodes.txt"))

#check if binning of adjacent barcodes should be run
if(opt$BarcodeBinning > 0 | !is.null(opt$barcode_sharing)){
  binmap <- BCbin(bccount_file = paste0(opt$out_dir,"/", opt$project, ".BCstats.txt"),
                  bc_detected  = bccount)
  fwrite(binmap,file=paste0(opt$out_dir,"/zUMIs_output/",opt$project,".BCbinning.txt"))
  #update the number reads in BCcount table
  binmap_additional <- binmap[, .(addtl = sum(n)), by = trueBC]
  bccount[match(binmap_additional$trueBC,XC),n := n + binmap_additional$addtl]
  fwrite(bccount,file=paste0(opt$out_dir,"/zUMIs_output/",opt$project,"kept_barcodes_binned.txt"))
}

##############################################################
q()
