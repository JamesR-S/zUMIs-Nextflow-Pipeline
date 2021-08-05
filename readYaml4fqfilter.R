#!/usr/bin/env Rscript
options(warn = -1)

suppressMessages(require(R.utils))
opt<-commandArgs(trailingOnly = T,asValues=T)

input<-read.csv(opt$input)

fq_pair <- as.list(readLines(file(opt$pairs_file)))

input <- input[grep((paste(fq_pair,collapse="|")),input$name),]

sequence_files <- as.list(NULL)
for (i in 1:nrow(input)){
  name <- input$name[i]
  cDNA <- if (is.na(input$cDNA[i])) NULL else paste0("cDNA(",input$cDNA[i],")")
  BC <- if (is.na(input$BC[i])) NULL else paste0("BC(",input$BC[i],")")
  UMI <- if (is.na(input$UMI[i])) NULL else paste0("UMI(",input$UMI[i],")")
  if (cDNA == "cDNA()") cDNA <- NULL
  if (BC == "BC()") BC <- NULL
  if (UMI == "UMI()") UMI <- NULL
  base_definition <- c(cDNA,BC,UMI)
  if (is.null(base_definition)) base_definition <- NULL
  filter_cutoffs <- input$filter_cutoffs[i]
  correct_frameshift <- input$correct_frameshift[i]
  if (is.na(filter_cutoffs)) filter_cutoffs <- NULL
  if (is.na(correct_frameshift)) correct_frameshift <- NULL
  sequence_files[[paste0("file",i)]] <- list(name,base_definition,filter_cutoffs,correct_frameshift)
  names(sequence_files[[paste0("file",i)]]) <- (c("name","base_definition","filter_cutoffs","correct_frameshift"))
}

BC_filter <- NULL
BC_filter$num_bases <- opt$BCfbases
BC_filter$phred <- opt$BCfphred

UMI_filter <- NULL
UMI_filter$num_bases <- opt$UMIfbases
UMI_filter$phred <- opt$UMIfphred

print( paste(sapply(sequence_files,function(x) { gsub("[[:space:]]", "", x$name) }),collapse=" ") )
print( paste(sapply(sequence_files,function(x) { paste(x$base_definition, collapse=";")}),collapse=" "))
print( opt$project)
print( opt$num_threads)
print( paste(BC_filter))
print( paste(UMI_filter))
print( paste(sapply(sequence_files,function(x) { paste(x$find_pattern)}),collapse=" "))
print( paste(sapply(sequence_files,function(x) { paste(x$correct_frameshift)}),collapse=" "))

q()