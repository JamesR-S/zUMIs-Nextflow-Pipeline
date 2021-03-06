#!/usr/bin/env nextflow

process copyResourceFiles {
    // Copy resource files into working directory -UNIX
    cpus = '4'
    memory = '20G'
    time = '00:05:00'

    output:
      file "${params.fastaName}" into fasta_ch
      file "${params.gtfName}" into gtf_ch1
      file "${params.gtfName}" into gtf_ch2

    """
    cp ${params.genome} .
    cp ${params.gtf} .

    gunzip *.gz

    """

}

process generateIndex {
    // Generate index with STAR -UNIX
    cpus = '6'
    memory = "${params.mem_limit}G"
    time = '02:00:00'

    input:
        val fasta from fasta_ch
        val gtf from gtf_ch1

    output:
        path "${params.genomeName}" into genome_index_ch

    """
    date
    Overhang=`expr ${params.read_length} - 1`

    if [ ! -d "${params.genomeName}" ]; then
    mkdir ${params.genomeName}
    fi

    STAR --runThreadN 6 \\
    --runMode genomeGenerate \\
    --genomeDir ${params.genomeName} \\
    --genomeFastaFiles ${fasta} \\
    --sjdbGTFfile ${gtf} \\
    --sjdbOverhang \${Overhang}
    date
    """
}


process detect_read_pairs {
    // Rscript created to detect read pairs so they can be split accross multiple processes
    cpus = '4'
    memory = '10G'
    time = '00:10:00'

    output:
    file "*_pair.txt" into read_pairs_ch1
    file "*_pair.txt" into read_pairs_ch2
    file "*_pair.txt" into read_pairs_ch3

    """
    #!/usr/bin/env Rscript
    Sys.time()
    input <- read.csv("${params.input_csv}") # read in iput csv
    Basesubset <-subset(input,BC != "" | UMI != "")\$name
    cDNAsubset <-subset(input,cDNA != "" )\$name

    compare_reads <- function(x,n,len){ # function to iterate over the header lines of the fq files looking for match
    line2 <- system(paste0("zcat -f < ", cDNAsubset2[n]," | head -n 1"), intern = TRUE)
    line2 <- as.list(strsplit(line2,split="\\\\s+")[[1]])
    if (identical(line1[[1]],line2[[1]])){
        read_pair <- c(x,cDNAsubset2[n])
        lapply(read_pair, write, paste0(basename(x),"_pair.txt"), append=TRUE, ncolumns=2)
    }else if (n<len){ 
        n <-n+1
        compare_reads(x,n,len)
    }else {
        if(x %in% Basesubset & x %in% cDNAsubset) { # if the file has both cdna and BC and a match cannot be identified assume SE read
            SE_read <- x
            lapply(SE_read, write, paste0(basename(x),"_pair.txt"), append=TRUE, ncolumns=2)
        }else { # else if match not identified and missing cDNA or barcode error out
            print("Read pair could not be found for fastq file and it did not have complete base definition. Please check samples")
            quit(status=1)
        }
    }
    }

    for (x in Basesubset){
    cDNAsubset2 <- cDNAsubset[cDNAsubset != x]
    len <- length(cDNAsubset2)
    n <- 1
    line1 <- system(paste0("zcat -f < ", x," | head -n 1"), intern = TRUE)
    line1 <- as.list(strsplit(line1,split="\\\\s+")[[1]])
    compare_reads(x,n,len)
    }
    Sys.time()
    q()
    """
}

process splitfq {
    // split fastq process - UNIX
    cpus = "${params.num_threads}"
    memory = '20G'
    time = '00:15:00'

    input:
    each file from read_pairs_ch3

    output:
    file '*_tempfiles.txt' into tempfiles_ch1
    file '*_tempfiles.txt' into tempfiles_ch2
    // uses pigz to read gzipped files and splits files into chunks to match number of available cores
    """
    date
    IFS=\$'\\r\\n' GLOBIGNORE='*' command eval  'reads=(\$(cat ${file}))'

    fullsize=\$(stat -L --printf="%s" \${reads[1]})

    basefile=\$(basename ${file} | cut -f 1 -d '.')

    if [[ \${reads[1]} =~ \\.gz\$ ]]; then
      pigz -dc \${reads[1]} | head -n 4000000 | pigz > ${params.projectDir}${params.outputDir}.\${basefile}.1mio.check.fq.gz
      smallsize=\$(stat --printf="%s" ${params.projectDir}${params.outputDir}.\${basefile}.1mio.check.fq.gz)
      rm ${params.projectDir}${params.outputDir}.\${basefile}.1mio.check.fq.gz
      nreads=\$(expr \${fullsize} \\* 1000000 / \${smallsize})
      n=`expr \$nreads / ${params.num_threads}`
      n=`expr \$n + 1`
      nl=`expr \$n \\* 4`
      for i in \${reads[@]} ; do
        pref=\$(basename \${i} | cut -f 1 -d '.')
        if [ ! -d "${params.projectDir}${params.outputDir}.\${pref}_tmpMerge" ]; then
            mkdir ${params.projectDir}${params.outputDir}.\${pref}_tmpMerge
        fi
        pigz -dc -p ${params.num_threads} \${i} | split --lines=\$nl --filter=''pigz' -p '${params.num_threads}' > \$FILE.gz' - ${params.projectDir}${params.outputDir}.\${pref}_tmpMerge/\${pref}${params.projectName}
        ls ${params.projectDir}${params.outputDir}.\${pref}_tmpMerge/\${pref}* | sed "s|${params.projectDir}${params.outputDir}.\${pref}_tmpMerge/\${pref}||"  > \${pref}_tempfiles.txt
        done
    else
      cat \${reads[1]} | head -n 4000000 > ${params.projectDir}${params.outputDir}.\${basefile}.1mio.check.fq
      smallsize=\$(stat --printf="%s" ${params.projectDir}${params.outputDir}.\${basefile}.1mio.check.fq)
      rm ${params.projectDir}${params.outputDir}.\${basefile}.1mio.check.fq
      nreads=\$(expr \${fullsize} \\* 1000000 / \${smallsize})
      n=`expr \$nreads / ${params.num_threads}`
      n=`expr \$n + 1`
      nl=`expr \$n \\* 4`
      for i in \${reads[@]} ; do
        pref=\$(basename \${i} | cut -f 1 -d '.')
        if [ ! -d "${params.projectDir}${params.outputDir}.\${pref}_tmpMerge" ]; then
            mkdir ${params.projectDir}${params.outputDir}.\${pref}_tmpMerge
        fi
        split --lines=\$nl --filter=''pigz' -p '${params.num_threads}' > \$FILE.gz' \${i} ${params.projectDir}${params.outputDir}.\${pref}_tmpMerge/\${pref}${params.projectName}
        ls ${params.projectDir}${params.outputDir}.\${pref}_tmpMerge/\${pref}* | sed "s|${params.projectDir}${params.outputDir}.\${pref}_tmpMerge/\${pref}||"  > \${pref}_tempfiles.txt
      done
    fi
    date
    """
}

process fqFilter {
    // filter fq file - Base script in unix but calls perl script
    cpus = "${params.num_threads}"
    memory = '20G'
    time { 2.hour * task.attempt }

    input:
    val tfiles from tempfiles_ch2.collect()
    each file from read_pairs_ch2

    output:
    file "*.bamlist.txt" into bamlists_ch
    env read_layout into read_layout_ch1
    env read_layout into read_layout_ch2
    env read_layout into read_layout_ch3
    file "*.BCstats.txt" into bcstats_ch
    // The script iterates over each temp file and filters the barcodes and umis based on parameters in config.
    // the perl script it calls was slightly modified from original provided with zUMIs pipeline to accomodate yaml removal
    """
    date
    filename=\$(basename ${file})

    parentfq=\${filename%"_pair.txt"}

    basefile=\$(echo \${parentfq} | cut -f 1 -d '.')

    tmpfiles=(${tfiles.join(" ")})

    tlookup=\$(printf -- '%s\n' "\${tmpfiles[@]}" | grep \$basefile)

    IFS=\$'\\r\\n' GLOBIGNORE='*' command eval  'temp_files=(\$(cat \${tlookup}))'

    for x in "\${temp_files[@]}" 
    do 
        perl ${params.projectDir}${params.binDir}fqfilter_v2.pl ${params.input_csv} samtools Rscript pigz ${params.projectDir}${params.binDir} \${x} ${params.UMI_phred} \\
        ${params.projectName} ${params.num_threads} \${parentfq} \${basefile} ${params.BC_num_bases} ${params.BC_phred} ${params.UMI_num_bases} ${params.projectDir}${params.outputDir} ${file}
    done
    wait
    ls ${params.projectDir}${params.outputDir}.\${basefile}_tmpMerge/${params.projectName}.*.filtered.tagged.bam > \${basefile}_${params.projectName}.bamlist.txt
    
    f=`head -n1 \${basefile}_${params.projectName}.bamlist.txt`
    flag=`samtools view \$f | head -n1 | cut -f2`

    if [[ \$flag == 4 ]]; then
    read_layout="SE"
    else
    read_layout="PE"
    fi

    date
    """
}

process collectBCstats {
    // simple unix process to collate the barcode counts produced by each of the read pair filtering processes.
    cpus = "${params.num_threads}"
    memory = '20G'
    time = '00:15:00'

    input:
    val files from bcstats_ch.collect()

    output:
    file "${params.projectName}.BCstats.txt" into cat_bcstats_ch
    file "${params.projectName}.BCstats.txt" into cat_bcstats_ch2

    """
    bcfiles=(${files.join(" ")})

    for f in "\${bcfiles[@]}"
    do
        cat \${f} >> ${params.projectName}.BCstats.txt
    done
    """

}

process barcode_detect {
    // Barcode detection Rscript - uses inflection package to determine optimal number of barcode to retain
    cpus = "${params.num_threads}"
    memory = "${params.mem_limit}G"
    time { 1.hour * task.attempt }

    input:
    val file from cat_bcstats_ch

    output:
    file "*.detected_cells.pdf" into detected_cells_ch
    file "*_kept_barcodes.txt" into kept_barcodes_ch1
    file "*_kept_barcodes.txt" into kept_barcodes_ch2
    file "*_kept_barcodes.txt" into kept_barcodes_ch3 
    file("*.BCbinning.txt") optional true into bc_binning_ch
    file("*kept_barcodes_binned.txt") optional true into kept_bc_binned_ch1
    file("*kept_barcodes_binned.txt") optional true into kept_bc_binned_ch2


    """
    #!/usr/bin/env Rscript
    Sys.time()
    .libPaths(c("/usr/local/lib/R/site-library","/usr/lib/R/site-library","/usr/lib/R/library"))

    library(methods)
    library(data.table)
    library(ggplot2)
    library(inflection)
    library(ggrastr)
    library(cowplot)
    library(stringdist)
    library(mclust)

    ##########################
    # Set parameters from config file
    barcode_file <- "${params.barcode_file}"
    if (barcode_file == "") barcode_file <- NULL
    barcode_num <- "${params.barcode_num}"
    if (barcode_num == "") barcode_num <- NULL
    bcauto <- ${params.barcode_automatic}
    out_dir <- "${params.projectDir}${params.outputDir}"
    project <- "${params.projectName}"
    zUMIs_directory <- "${params.projectDir}${params.binDir}"
    n_reads_per_cell <- ${params.n_reads_per_cell}
    barcode_binning <- "${params.barcode_binning}"
    num_threads <- ${params.num_threads}
    barcode_sharing <- "${params.barcode_sharing}"
    if (barcode_sharing == "") barcode_sharing <- NULL
    # read in csv file
    input <- read.csv("${params.input_csv}")
    # convert cs file input into original yaml format
    sequence_files <- as.list(NULL)
    for (i in 1:nrow(input)){
        name <- input\$name[i]
        base_definition <- c(if (is.na(input\$cDNA[i])) NULL else paste0("cDNA(",input\$cDNA[i],")"),if (is.na(input\$BC[i])) NULL else paste0("BC(",input\$BC[i],")"),if (is.na(input\$UMI[i])) NULL else paste0("UMI(",input\$UMI[i],")"))
        complete<- NULL
        complete\$name <- name
        complete\$base_definition <- base_definition
        complete\$filter_cutoffs <- input\$filter_cutoffs[i]
        complete\$correct_frameshift <- input\$correct_frameshift[i]
        sequence_files[[paste0("file",i)]] <- complete
    }


    source(paste0(zUMIs_directory,"/barcodeIDFUN.R"))
    options(datatable.fread.input.cmd.message=FALSE)
    data.table::setDTthreads(threads=num_threads)
    if(!is.null(barcode_sharing)){
    if(barcode_sharing == ""){
        barcode_sharing <- NULL
    }
    }

    #######################################################################
    # Barcode detection

    #read BC count file
    bccount<-cellBC(bcfile      = barcode_file,
            bcnum       = barcode_num,
            bcauto      = bcauto,
            bccount_file= "${file}",
            outfilename = paste0(project,".detected_cells.pdf"))

    fwrite(bccount,file=paste0(project,"_kept_barcodes.txt"))

    #check if binning of adjacent barcodes should be run
    if(barcode_binning > 0 | !is.null(barcode_sharing)){
    binmap <- BCbin(bccount_file = "${file}",
                    bc_detected  = bccount)
    fwrite(binmap,file=paste0(project,".BCbinning.txt"))
    #update the number reads in BCcount table
    binmap_additional <- binmap[, .(addtl = sum(n)), by = trueBC]
    bccount[match(binmap_additional\$trueBC,XC),n := n + binmap_additional\$addtl]
    fwrite(bccount,file=paste0(project,"kept_barcodes_binned.txt"))
    }
    Sys.time() 
    ##############################################################
    q()    

    """


}
// output can vary depending on options - optional output channels filled with dummy value if output not produced.
bc_binning_ch = bc_binning_ch.ifEmpty('EMPTY')
kept_bc_binned_ch1 = kept_bc_binned_ch1.ifEmpty('EMPTY')

process BC_bin_bam_collate {
    // If barcodes were binned this process amends the binned barcodes within the bam files to reflect new sequences
    // otherwise this process just collates the bam files into single dir ready for mapping.
    cpus = "${params.num_threads}"
    memory = "${params.mem_limit}G"
    time { 1.hour * task.attempt }

    input:
    val BCbinning from bc_binning_ch
    val tfiles from tempfiles_ch1.collect()
    each file from read_pairs_ch1


    output:
    env tmpFILT into collated_bam_ch

    """
    date
    filename=\$(basename ${file})

    parentfq=\${filename%"_pair.txt"}

    basefile=\$(echo \${parentfq} | cut -f 1 -d '.')

    tmpfiles=(${tfiles.join(" ")})

    tlookup=\$(printf -- '%s\n' "\${tmpfiles[@]}" | grep \$basefile)

    IFS=\$'\\r\\n' GLOBIGNORE='*' command eval  'temp_files=(\$(cat \${tlookup}))'

    [ -d ${params.projectDir}${params.outputDir}.tmpFiltered ] || mkdir ${params.projectDir}${params.outputDir}.tmpFiltered
    # if statement checks if value in bc binning channel is file - if bc binning was run then it will just be a dummy value and return False
    if [[ -f "${BCbinning}" ]] ; then
        for x in "\${temp_files[@]}" ; do
            rawbam="${params.projectDir}${params.outputDir}.\${basefile}_tmpMerge/${params.projectName}.\${x}.raw.tagged.bam"
            fixedbam="${params.projectDir}${params.outputDir}.\${basefile}_tmpMerge/${params.projectName}.\${x}.filtered.tagged.bam"
            mv \${fixedbam} \${rawbam}
            perl ${params.projectDir}${params.binDir}correct_BCtag.pl \${rawbam} \${fixedbam} ${BCbinning} samtools &
        done
        wait
    fi

    for x in "\${temp_files[@]}" ; do
        mv ${params.projectDir}${params.outputDir}.\${basefile}_tmpMerge/${params.projectName}.\${x}.filtered.tagged.bam ${params.projectDir}${params.outputDir}.tmpFiltered/${params.projectName}.\${basefile}.\${x}.filtered.tagged.bam
    done

    while read p; do
        fname=\$(basename \${p} | cut -f 1 -d '.')
        rm -r ${params.projectDir}${params.outputDir}.\${fname}_tmpMerge
    done <${file}

    tmpFILT="${params.projectDir}${params.outputDir}.tmpFiltered"
    date
    """
}

process Mapping {
    // mapping process (rscript) maps bam files using STAR mapper. if resources available this process unes on multiple threads
    cpus = "${params.num_threads}"
    memory = "${params.mem_limit}G"
    time { 1.hour * task.attempt }

    input:
    val tmpFILT from collated_bam_ch.last()
    path index from genome_index_ch
    val layout from read_layout_ch1.last()
    file gtf from gtf_ch2

    publishDir "${params.projectDir}${params.outputDir}", mode: 'copy', pattern: "*.filtered.tagged.Log.final.out"

    output:
    file "*.filtered.tagged.Log.final.out" into mapping_log_ch
    file "*.filtered.tagged.Aligned.out.bam" into aligned_bam_ch
    file "*.filtered.tagged.Aligned.toTranscriptome.out.bam" into txbam_ch


    """
    #!/usr/bin/env Rscript --vanilla

    .libPaths(c("/usr/local/lib/R/site-library","/usr/lib/R/site-library","/usr/lib/R/library"))

    suppressMessages(require(data.table))
    suppressMessages(require(R.utils))
    options(datatable.fread.input.cmd.message=FALSE)
    inp <- commandArgs(trailingOnly = T, asValues=T)
    # executables 
    samtools <- "samtools"
    STAR_exec <- "STAR"

    # Read in parameters from config file
    mem_limit <- ${params.mem_limit}
    out_dir <- "${params.projectDir}${params.outputDir}"
    project <- "${params.projectName}"
    STAR_index <- "${index}"
    num_threads <- ${params.num_threads}
    GTF_file <- "${gtf}"
    read_layout <- "${layout}"
    twoPass <- ${params.two_pass}
    additional_STAR_params <- "${params.additional_STAR_params}"
    if (identical(additional_STAR_params,"NONE")){
    additional_STAR_params <- NULL
    }
    additional_fq <- "${params.additional_fq}"
    if (identical(additional_fq,"NONE")){
    additional_fq <- NULL
    }else{
        additional_fq <- strsplit(additional_fq,",")
    }
    if(is.null(mem_limit)){
    mem_limit <- 100
    }else if(mem_limit == 0){
    mem_limit <- 100
    }

    # merge chunked bam files
    tmpfolder <- "${tmpFILT}"
    filtered_bams <- list.files(path = tmpfolder, pattern=paste(project,".*.filtered.tagged.bam",sep=""),full.names=T)
    sammerge_command <- paste(samtools,"cat -o",paste0(out_dir,project,".filtered.tagged.unmapped.bam"),paste0(filtered_bams,collapse=" "))
    system(sammerge_command)


    # check multithreading resources

    genome_size <- system(command = paste("du -sh",STAR_index,"| cut -f1"), intern = TRUE)
    genome_size <- as.numeric(gsub(pattern = "G",replacement = "", x = genome_size))
    if(is.na(genome_size)){
    genome_size <- 25 #set average genome size if there was a problem detecting
    }
    num_star_instances <- floor(mem_limit/genome_size)
    if(num_star_instances < 1){
    num_star_instances = 1 #set the number of STAR instances to 1 if it is 0
    }
    if(num_star_instances > num_threads){
    num_star_instances = num_threads
    } 

    # GTF file setup 

    if ( is.null(additional_fq[1]) | length(additional_fq)==0) {
    gtf_to_use <- GTF_file
    param_additional_fa <- NULL
    system(paste0("cp ",gtf_to_use," ",out_dir,project,".final_annot.gtf"))
    }else{
    for (i in additional_fq) {
        system(paste(samtools,"faidx",i))
        assign(paste("fai",i,sep="_"),data.table::fread(input = paste("cut -f1,2 ",i,".fai",sep=""),stringsAsFactors = F,data.table = F))
    }

    ref_df <- do.call("rbind", mget(ls(pattern = "fai_")))

    user_gtf <- data.frame(
        V1 = ref_df\$V1,
        V2 = "User",
        V3 = "exon",
        V4 = 1,
        V5 = ref_df\$V2,
        V6 = ".",
        V7 = "+",
        V8 = ".",
        V9 = paste('gene_id "',ref_df\$V1,'"; transcript_id "',ref_df\$V1,'"; exon_number "1"; gene_name "',ref_df\$V1,'"; gene_biotype "User"; transcript_name "',ref_df\$V1,'"; exon_id "',ref_df\$V1,'"',sep = ""),
        stringsAsFactors = F
    )

    write.table(user_gtf,file = paste(out_dir,"additional_sequence_annot.gtf",sep = ""),sep = "\t",quote = F,row.names = F,col.names = F)

    system(command = paste("cat ",GTF_file," ",paste(out_dir,"additional_sequence_annot.gtf",sep = "")," > ",out_dir,project,".final_annot.gtf",sep=""))

    gtf_to_use <- paste(out_dir,project,".final_annot.gtf",sep="")
    param_additional_fa <- paste("--genomeFastaFiles",paste(additional_fq,collapse = " "))
    }

    #GTF_file_final <- gtf_to_use

    # Detect read length from from top 1000 reads - used for sjdb ooverhang
    filtered_bam <- paste(out_dir,project,".filtered.tagged.unmapped.bam",sep="")

    cDNA_peek <- data.table::fread(cmd = paste(samtools,"view",filtered_bam,"| cut -f10 | uniq | head -n 1000"),stringsAsFactors = F,data.table = T, header = F)

    getmode <- function(v) {
    uniqv <- unique(v)
    uniqv[which.max(tabulate(match(v, uniqv)))]
    }

    cDNA_read_length <- getmode(nchar(cDNA_peek\$V1))


    # Setup STAR commands
    samtools_load_cores <- ifelse(num_threads>8,2,1)
    avail_cores <- num_threads - samtools_load_cores #reserve threads for samtools file opening
    
    avail_cores <- floor(avail_cores / num_star_instances)


    if(avail_cores < 2){
    avail_cores = 1
    }


    param_defaults <- paste("--readFilesCommand ",samtools," view -@",samtools_load_cores," --outSAMmultNmax 1 --outFilterMultimapNmax 50 --outSAMunmapped Within --outSAMtype BAM Unsorted --quantMode TranscriptomeSAM")
    param_misc <- paste("--genomeDir",STAR_index,
                        "--sjdbGTFfile",gtf_to_use,
                        "--runThreadN",avail_cores,
                        "--sjdbOverhang", cDNA_read_length-1,
                        "--readFilesType SAM",read_layout)

    STAR_command <- paste(STAR_exec,param_defaults,param_misc,additional_STAR_params,param_additional_fa)
    if(twoPass==TRUE){
    STAR_command <- paste(STAR_command,"--twopassMode Basic")
    }

    #Execute STAR mapping process 
    if(num_star_instances>1){
    map_tmp_dir <- "tmpMap"
    dir.create(path = map_tmp_dir,showWarnings = FALSE)
    input_split <- split(filtered_bams, ceiling(seq_along(filtered_bams) / ceiling(length(filtered_bams) / num_star_instances)))
    input_split <- sapply(input_split, paste0, collapse = ",")
    STAR_preset <- STAR_command
    STAR_command <- lapply(seq(num_star_instances), function(x){
        paste(STAR_preset,
        "--readFilesIn",input_split[x],
        "--outFileNamePrefix",paste0(map_tmp_dir,"/tmp.",project,".",x,"."))
    })
    STAR_command <- paste(unlist(STAR_command), collapse = " & ")
    system(paste(STAR_command,"& wait"), wait = TRUE)
    
    # if multithreading was used files need to be collated
    out_logs <- list.files(map_tmp_dir, pattern = paste0("tmp.",project,".*.Log.final.out"), full = TRUE)
    merge_logs <- paste("cat",paste(out_logs, collapse = " "),">",paste0(project,".filtered.tagged.Log.final.out"))
    out_bams <- list.files(map_tmp_dir, pattern = paste0("tmp.",project,".*.Aligned.out.bam"), full = TRUE)
    merge_bams <- paste(samtools,"cat -o",paste0(project,".filtered.tagged.Aligned.out.bam"),paste(out_bams, collapse = " "))
    out_txbams <- list.files(map_tmp_dir, pattern = paste0("tmp.",project,".*.Aligned.toTranscriptome.out.bam"), full = TRUE)
    merge_txbams <- paste(samtools,"cat -o",paste0(project,".filtered.tagged.Aligned.toTranscriptome.out.bam"),paste(out_txbams, collapse = " "))
    system(paste(merge_logs,";",merge_bams,";",merge_txbams, "; wait"), wait = TRUE)
    system(paste0("echo 'Cleaning temp files' ; rm -r ", map_tmp_dir, " ; rm -r ",out_dir, ".tmpFiltered ; echo '...'" ), wait = TRUE)
    while( dir.exists(map_tmp_dir)) {
        system("sleep 1")
        print("...")
    } 
    while( dir.exists(paste0(out_dir,".tmpFiltered"))) {
        system("sleep 1")
        print("...")
    }
    while( !dir.exists(paste0(out_dir,".tmpFiltered"))) {
    while( !dir.exists(map_tmp_dir)) {
    print(paste0("Cleanup complete, ending process at ",Sys.time()))
    quit(save = "no")
    }}
    }
    if(num_star_instances==1){
    STAR_command <- paste(STAR_command,
        "--readFilesIn",paste0(filtered_bams,collapse=","),
        "--outFileNamePrefix",paste(project,".filtered.tagged.",sep="")
    )
    system(paste(STAR_command,"& wait"))
    system(paste0("echo 'Cleaning temp files' ; rm -r ",out_dir, ".tmpFiltered ; wait'" ), wait = TRUE)
    while( !dir.exists(paste0(out_dir,".tmpFiltered"))) {
        print(paste0("Cleanup complete, ending process at ",Sys.time()))
        quit(save = "no")
    }}

    """


}


process Counting {
    // Count reads mapped to gtf file features
    cpus = "${params.num_threads}"
    memory = "${params.mem_limit}G"

    time { 3.hour * task.attempt }

    errorStrategy { task.exitStatus!=0 ? 'retry' : 'terminate' }
    maxRetries 3
    
    input:
    val barcodes from kept_barcodes_ch1
    val keptBC from kept_bc_binned_ch1
    val readLO from read_layout_ch2.last()
    file alignedBAM from aligned_bam_ch

    publishDir "${params.projectDir}${params.outputDir}zUMIs_output", mode: 'copy'

    output:
    file "*.filtered.Aligned.GeneTagged.UBcorrected.sorted.bam" into UB_corrected_bam_ch1
    file "*.filtered.Aligned.GeneTagged.UBcorrected.sorted.bam" into UB_corrected_bam_ch2
    file "*.intronProbability.rds" optional true into intronProbability_ch
    file "*.dgecounts.rds" into dgecounts_ch1
    file "*.dgecounts.rds" into dgecounts_ch2

    """
    #!/usr/bin/env Rscript
    library(methods)
    library(data.table)
    library(yaml)
    library(ggplot2)
    suppressMessages(require(R.utils))
    suppressPackageStartupMessages(library(Rsamtools))
    Sys.time()
    # Read in parameters from config file
    opt <- NULL
    opt\$zUMIs_directory <- "${params.projectDir}${params.binDir}"
    opt\$counting_opts\$Ham_Dist <- ${params.ham_dist}
    opt\$barcodes\$BarcodeBinning <- ${params.barcode_binning}
    opt\$reference\$exon_extension <- ${params.exon_extension}
    opt\$reference\$extension_length <- ${params.extension_length}
    opt\$reference\$scaffold_length_min <- ${params.scaffold_length_min}
    opt\$counting_opts\$multi_overlap <- ${params.multi_overlap}
    opt\$counting_opts\$strand <- ${params.strand}
    opt\$counting_opts\$primaryHit <- ${params.primary_hit}
    opt\$counting_opts\$introns <- ${params.include_introns}
    opt\$counting_opts\$intronProb <- ${params.intron_probability}
    opt\$counting_opts\$downsampling <- ${params.downsampling}
    opt\$barcodes\$demultiplex <- ${params.demultiplex}
    opt\$out_dir <- "${params.projectDir}${params.outputDir}"
    opt\$project <- "${params.projectName}"
    opt\$num_threads <- ${params.num_threads}
    opt\$mem_limit <- ${params.mem_limit}
    opt\$read_layout <- "${readLO}"
    opt\$keptbc <- "${barcodes}"

    # read in input csv
    input <- read.csv("${params.input_csv}")
    # convert to original yaml format
    sequence_files <- as.list(NULL)
    for (i in 1:nrow(input)){
        name <- input\$name[i]
        base_definition <- c(if (is.na(input\$cDNA[i])) NULL else paste0("cDNA(",input\$cDNA[i],")"),if (is.na(input\$BC[i])) NULL else paste0("BC(",input\$BC[i],")"),if (is.na(input\$UMI[i])) NULL else paste0("UMI(",input\$UMI[i],")"))
        complete<- NULL
        complete\$name <- name
        complete\$base_definition <- base_definition
        complete\$filter_cutoffs <- input\$filter_cutoffs[i]
        complete\$correct_frameshift <- input\$correct_frameshift[i]
        sequence_files[[paste0("file",i)]] <- complete
    }

    # load in required functions from other script files
    source(paste0(opt\$zUMIs_directory,"/runfeatureCountFUN.R"))
    source(paste0(opt\$zUMIs_directory,"/misc/featureCounts.R"))
    source(paste0(opt\$zUMIs_directory,"/UMIstuffFUN.R"))
    source(paste0(opt\$zUMIs_directory,"/barcodeIDFUN.R"))
    options(datatable.fread.input.cmd.message=FALSE)
    print(Sys.time())
    samtoolsexc <- "samtools"
    data.table::setDTthreads(threads=1)

    #Check Rsubread version
    fcounts_clib <- paste0(opt\$zUMIs_directory,"/misc/fcountsLib2")

    opt <- fixMissingOptions(opt)

    # check for non-UMI
    UMIcheck <- check_nonUMIcollapse(sequence_files)
    if(UMIcheck == "nonUMI"){
    opt\$counting_opts\$Ham_Dist <- 0
    }
    #is the data Smart-seq3?
    smart3_flag <- ifelse(any(grepl(pattern = "ATTGCGCAATG",x = unlist(sequence_files))), TRUE, FALSE)


    #read in barcode count file

    #check binning options from config
    if(opt\$barcodes\$BarcodeBinning > 0){
    bccount <- fread("${keptBC}")
    }else{
    bccount <- fread("${barcodes}")
    }
    bccount<-splitRG(bccount=bccount, mem= opt\$mem_limit, hamdist = opt\$counting_opts\$Ham_Dist)

    #Run feature counts

    abamfile<-"${alignedBAM}"
    outbamfile <-paste0(opt\$project,".filtered.Aligned.GeneTagged.bam")

    # gene annotation
    saf<-.makeSAF(gtf = paste0(opt\$out_dir,"/",opt\$project,".final_annot.gtf"),
                extension_var = opt\$reference\$exon_extension,
                exon_extension = opt\$reference\$extension_length,
                buffer_length = (opt\$reference\$extension_length / 2),
                scaff_length = opt\$reference\$scaffold_length_min,
                multi_overlap_var = opt\$multi_overlap,
                samtoolsexc = samtoolsexc)
    try(gene_name_mapping <- .get_gene_names(gtf = paste0(opt\$out_dir,"/",opt\$project,".final_annot.gtf"), threads = opt\$num_threads), silent = TRUE)
    try(data.table::fwrite(gene_name_mapping, file = paste0(opt\$out_dir,"zUMIs_output/expression/",opt\$project,".gene_names.txt"), sep ="\t", quote = FALSE), silent = TRUE)
    ##

    if(smart3_flag & opt\$counting_opts\$strand == 1){
    #split bam file
    print("Preparing Smart-seq3 data for stranded gene assignment...")
    print(Sys.time())
    tmp_bams <- split_bam(bam = abamfile, cpu = opt\$num_threads, samtoolsexc=samtoolsexc)

    #check stranding
    fnex_int<-.runFeatureCount(tmp_bams[1], saf=saf\$exons, strand=0, type="ex", primaryOnly = opt\$counting_opts\$primaryHit, cpu = opt\$num_threads, mem = opt\$mem_limit, fcounts_clib = fcounts_clib, multi_overlap_var = opt\$counting_opts\$multi_overlap)
    fnex_umi<-.runFeatureCount(tmp_bams[2], saf=saf\$exons, strand=1, type="ex", primaryOnly = opt\$counting_opts\$primaryHit, cpu = opt\$num_threads, mem = opt\$mem_limit, fcounts_clib = fcounts_clib, multi_overlap_var = opt\$counting_opts\$multi_overlap)
    ffiles_int <- paste0(fnex_int,".tmp")
    ffiles_umi <- paste0(fnex_umi,".tmp")
    # intron mapped read counting
    if(opt\$counting_opts\$introns){
        fnin_int<-.runFeatureCount(ffiles_int, saf=saf\$introns, strand=0, type="in", primaryOnly = opt\$counting_opts\$primaryHit, cpu = opt\$num_threads, mem = opt\$mem_limit, fcounts_clib = fcounts_clib, multi_overlap_var = opt\$counting_opts\$multi_overlap)
        fnin_umi<-.runFeatureCount(ffiles_umi, saf=saf\$introns, strand=1, type="in", primaryOnly = opt\$counting_opts\$primaryHit, cpu = opt\$num_threads, mem = opt\$mem_limit, fcounts_clib = fcounts_clib, multi_overlap_var = opt\$counting_opts\$multi_overlap)
        ffiles_int <- paste0(fnin_int,".tmp")
        ffiles_umi <- paste0(fnin_umi,".tmp")
    }
    join_bam_cmd <- paste(samtoolsexc, "cat -o", outbamfile, ffiles_int, ffiles_umi)
    system(join_bam_cmd)
    system(paste0("rm ",tmp_bams[1],"* ",tmp_bams[2],"*"))
    }else{
    fnex<-.runFeatureCount(abamfile,
                            saf=saf\$exons,
                            strand=opt\$counting_opts\$strand,
                            type="ex",
                            primaryOnly = opt\$counting_opts\$primaryHit,
                            cpu = opt\$num_threads,
                            mem = opt\$mem_limit,
                            fcounts_clib = fcounts_clib,
                            multi_overlap_var = opt\$counting_opts\$multi_overlap)
    ffiles<-paste0(fnex,".tmp")

    if(opt\$counting_opts\$introns){
        fnin  <-.runFeatureCount(ffiles,
                                saf=saf\$introns,
                                strand=opt\$counting_opts\$strand,
                                type="in",
                                primaryOnly = opt\$counting_opts\$primaryHit,
                                cpu = opt\$num_threads,
                                mem = opt\$mem_limit,
                                fcounts_clib = fcounts_clib,
                                multi_overlap_var = opt\$counting_opts\$multi_overlap)
        system(paste0("rm ",fnex,".tmp"))
        ffiles<-paste0(fnin,".tmp")
    }

    system(paste("mv",ffiles,outbamfile))
    }

    if(is.null(opt\$mem_limit)){
    mempercpu <- max(round(100/opt\$num_threads,0),1)
    }else{
    mempercpu <- max(round(opt\$mem_limit/opt\$num_threads,0),1)
    }

    # binning similar UMIs
    if(opt\$counting_opts\$Ham_Dist == 0){
    sortbamfile <-paste0(opt\$project,".filtered.Aligned.GeneTagged.sorted.bam")
    print(Sys.time())
    print("Coordinate sorting final bam file...")
    sort_cmd <- paste0(samtoolsexc," sort -O 'BAM' -@ ",opt\$num_threads," -m ",mempercpu,"G -o ",sortbamfile," ",outbamfile)
    system(sort_cmd)
    system(paste0("rm ",outbamfile))
    }else{
    #run hamming distance collapsing here and write output into bam file
    if(!dir.exists( paste0(opt\$out_dir,"zUMIs_output/molecule_mapping/") )){
        dir.create( paste0(opt\$out_dir,"zUMIs_output/molecule_mapping/") )
    }

    tmpbamfile <- outbamfile
    outbamfile <- paste0(opt\$project,".filtered.Aligned.GeneTagged.sorted.bam")
    print(Sys.time())
    print("Coordinate sorting intermediate bam file...")
    sort_cmd <- paste0(samtoolsexc," sort -O 'BAM' -@ ",opt\$num_threads," -m ",mempercpu,"G -o ",outbamfile," ",tmpbamfile)
    system(sort_cmd)
    index_cmd <- paste(samtoolsexc,"index -@",opt\$num_threads,outbamfile)
    system(index_cmd)
    system(paste0("rm ",tmpbamfile))
    print(Sys.time())
    
    #check read layout
    if(is.null(opt\$read_layout)){
        opt\$read_layout <- check_read_layout(outbamfile)
    }
    
    for(i in unique(bccount\$chunkID)){
        print( paste( "Hamming distance collapse in barcode chunk", i, "out of",length(unique(bccount\$chunkID)) ))
        reads <- reads2genes_new(featfile = outbamfile,
                                bccount  = bccount,
                                inex     = opt\$counting_opts\$introns,
                                chunk    = i,
                                cores    = opt\$num_threads)
        reads <- reads[!UB==""] #make sure only UMI-containing reads go further
        u <- umiCollapseHam(reads,bccount, HamDist=opt\$counting_opts\$Ham_Dist)
    }
    print("Correcting UMI barcode tags...")
    sortbamfile <- correct_UB_tags_new(outbamfile, opt\$project)
    file.remove(outbamfile)
    bccount<-splitRG(bccount=bccount, mem= opt\$mem_limit, hamdist = 0) # allow more reads to be in RAM fur subsequent steps
    }
    index_cmd <- paste(samtoolsexc,"index -@",opt\$num_threads,sortbamfile)
    system(index_cmd)
    print(Sys.time())

    #check if read layout
    if(is.null(opt\$read_layout)){
    opt\$read_layout <- check_read_layout(sortbamfile)
    }

    #Adaptive Dowwnsampling

    data.table::setDTthreads(threads=opt\$num_threads)

    subS<-setDownSamplingOption( opt\$counting_opts\$downsampling,
                                bccount= bccount,
                                filename=paste(opt\$out_dir,"zUMIs_output/stats/",opt\$project,
                                                ".downsampling_thresholds.pdf",sep=""))
    print("Here are the detected subsampling options:")
    if(is.null(row.names(subS))){
    print("Automatic downsampling")
    }else{
    print(row.names(subS))
    }
    if( opt\$counting_opts\$introns ){
    mapList<-list("exon"="exon",
                    "inex"=c("intron","exon"),
                    "intron"="intron")
    }else{
    mapList<-list("exon"="exon")
    }


    #assign reads to UMI and Genes

    for(i in unique(bccount\$chunkID)){
        print( paste( "Working on barcode chunk", i, "out of",length(unique(bccount\$chunkID)) ))
        print( paste( "Processing",length(bccount[chunkID==i]\$XC), "barcodes in this chunk..." ))
        reads <- reads2genes_new(featfile = sortbamfile,
                                bccount  = bccount,
                                inex     = opt\$counting_opts\$introns,
                                chunk    = i,
                                cores    = opt\$num_threads)

        tmp<-collectCounts(  reads =reads,
                            bccount=bccount[chunkID==i],
                            subsample.splits=subS[which(max(bccount[chunkID==i]\$n) >= subS[,1]), , drop = FALSE],
                            mapList=mapList
                            )

        if(i==1){
        allC<-tmp
        }else{
        allC<-bindList(alldt=allC,newdt=tmp)
        }
    }

    if( UMIcheck == "UMI"  ){
    if(smart3_flag){
        final<-list( umicount  = convert2countM(alldt=allC,what="umicount"),
                    readcount = convert2countM(allC,"readcount"),
                    readcount_internal = convert2countM(allC,"readcount_internal"))
    }else{
        final<-list( umicount  = convert2countM(alldt=allC,what="umicount"),
                    readcount = convert2countM(allC,"readcount"))
    }
    }else{
    final<-list(readcount = convert2countM(allC,"readcount"))
    }

    if(UMIcheck == "nonUMI" | smart3_flag == TRUE ){
    tx_len <- .get_tx_lengths( paste0(opt\$out_dir,"/",opt\$project,".final_annot.gtf") )
    rpkms <- if(smart3_flag){
        RPKM.calc(final\$readcount_internal\$exon\$all, tx_len)
    }else{
        RPKM.calc(final\$readcount\$exon\$all, tx_len)
    }

    final\$rpkm <- list(exon = list(all = rpkms))
    }

    saveRDS(final,file=paste(opt\$project,".dgecounts.rds",sep=""))


    #Intron Probability calculation
    if(opt\$counting_opts\$intronProb == TRUE){
    if(opt\$counting_opts\$introns == FALSE){
        print("Intron information is needed to calculate Intron-Probability! Please change yaml settings counting_opts->introns to yes.\\n Skipping probability calculation...")
    }else{
        library(extraDistr)
        print("Starting Intron-Probability Calculations...")

        # Extractingreads from sam files again, this could be sped up by integrating code further upstream of zUMIs-dge2
        genesWithIntronProb<-.intronProbability(bccount=bccount,
                                                featfile=sortbamfile,
                                                inex=opt\$counting_opts\$introns,
                                                cores=opt\$num_threads,
                                                samtoolsexc=samtoolsexc,
                                                allC=allC,
                                                saf=saf)
        saveRDS(genesWithIntronProb, file = paste0(opt\$project,".intronProbability.rds"))
    }
    }

    #Demultiplex bam files
    if(opt\$barcodes\$demultiplex){
    print("Demultiplexing output bam file by cell barcode...")
    demultiplex_bam(opt, sortbamfile, nBCs = length(unique(bccount\$XC)), bccount = bccount, samtoolsexc = samtoolsexc, BCkeepfile= "${barcodes}")
    }

    #################

    print(Sys.time())
    print(paste("I am done!! Look what I produced...",opt\$out_dir,"zUMIs_output/",sep=""))
    print(gc())
    Sys.time()
    q(save = "no")

    """
}
// Check if loompy output should be produced
if(params.loompy) {
    process rds2loom {
        // Convert RDS count file to loompy format
        cpus = "${params.num_threads}"
        memory = "${params.mem_limit}G"
        time { 2.hour * task.attempt }
        
        input:
        file counts from dgecounts_ch1

        publishDir "${params.projectDir}${params.outputDir}zUMIs_output", mode: 'copy'

        output:
        file "*.loom" into loom_files_ch

        """
        #!/usr/bin/env Rscript
        require(yaml)
        require(Matrix)
        suppressMessages(require(R.utils))
        Sys.time()

        opt <- NULL

        opt\$out_dir <- "${params.projectDir}${params.outputDir}"
        opt\$project <- "${params.projectName}"
        #  Read in dgecounts rds file and convert to loompy file - relies on loomR package
        rds_to_loom <- function(zUMIsRDS){
        rds <- readRDS(zUMIsRDS)
        for(type in names(rds)){
            for(quant in names(rds[[type]])){
            outfile <- paste0(opt\$project,".",type,".",quant,".all.loom")
            loomR::create(filename = outfile, data = as.matrix(rds[[type]][[quant]][["all"]]), do.transpose = TRUE,overwrite = TRUE)
            
            for(d in names(rds[[type]][[quant]][["downsampling"]])){
                outfile <- paste0(opt\$project,".",type,".",quant,".",d,".loom")
                loomR::create(filename = outfile, data = as.matrix(rds[[type]][[quant]][["downsampling"]][[d]]), do.transpose = TRUE,overwrite = TRUE)
                
            }
            }
        }
        }

        suppressMessages(require("loomR"))


        ##########################

        rds_loc <- "${counts}"
        rds_to_loom(rds_loc)

        Sys.time()
        q(save = "no")

        """
    }
}
// Run RNA velocity analysis?
if(params.velocyto) {

    process velocyto {
        // velocyto R package used to run RNA velocity analysis       
        cpus = "${params.num_threads}"
        memory = "${params.mem_limit}G"
        
        time { 3.hour * task.attempt }

        errorStrategy { task.exitStatus!=0 ? 'retry' : 'terminate' }
        maxRetries 3

        input:
        file bam from UB_corrected_bam_ch1
        file barcodes from kept_barcodes_ch3

        """
        #!/usr/bin/env Rscript
        suppressMessages(require(R.utils))
        Sys.time()

        ##########################

        samtoolsexc <- "samtools"

        opt <- NULL

        opt\$out_dir <- "${params.projectDir}${params.outputDir}"
        opt\$project <- "${params.projectName}"
        opt\$num_threads <- ${params.num_threads}
        opt\$mem_limit <- ${params.mem_limit}
        opt\$counting_opts\$primaryHit <- ${params.primary_hit}
        #Read in input csv file
        input<-read.csv("${params.input_csv}")
        #Convert into original yaml format
        sequence_files <- as.list(NULL)
        for (i in 1:nrow(input)){
        name <- input\$name[i]
        cDNA <- if (is.na(input\$cDNA[i])) NULL else paste0("cDNA(",input\$cDNA[i],")")
        BC <- if (is.na(input\$BC[i])) NULL else paste0("BC(",input\$BC[i],")")
        UMI <- if (is.na(input\$UMI[i])) NULL else paste0("UMI(",input\$UMI[i],")")
        if (cDNA == "cDNA()") cDNA <- NULL
        if (BC == "BC()") BC <- NULL
        if (UMI == "UMI()") UMI <- NULL
        base_definition <- c(cDNA,BC,UMI)
        if (is.null(base_definition)) base_definition <- NULL
        filter_cutoffs <- input\$filter_cutoffs[i]
        correct_frameshift <- input\$correct_frameshift[i]
        if (is.na(filter_cutoffs)) filter_cutoffs <- NULL
        if (is.na(correct_frameshift)) correct_frameshift <- NULL
        sequence_files[[paste0("file",i)]] <- list(name,base_definition,filter_cutoffs,correct_frameshift)
        names(sequence_files[[paste0("file",i)]]) <- (c("name","base_definition","filter_cutoffs","correct_frameshift"))
        }

        if(is.null(opt\$mem_limit)){
          mempercpu <- round(100/opt\$num_threads,0)
        }else{
          mempercpu <- round(opt\$mem_limit/opt\$num_threads,0)
          if(mempercpu==0){
                mempercpu <- 1
          }
        }

        featfile <- "${bam}"


        print(Sys.time())
        #load in bam file and reformat for velocyto
        print("Preparing bam file for velocyto...")
        retag_cmd <- paste0(samtoolsexc," view -@ 2 -h ",featfile," | sed 's/BC:Z:/CB:Z:/' | sed 's/GE:Z:/GX:Z:/'")
        velobam <- paste0(opt\$out_dir,"/",opt\$project,".tagged.forVelocyto.bam")
        out_cmd <- paste0(samtoolsexc," view -b -@ ",opt\$num_threads," -o ",velobam," - " )
        system(paste(retag_cmd,out_cmd,sep=" | "))

        print(Sys.time())

        #create BC whitelist from complete barcodes
        bc<-data.table::fread("${barcodes}",select = 1, header = T)
        bcpath<-paste0(opt\$out_dir,"zUMIs_output/",opt\$project,".BCwhitelist.txt")
        data.table::fwrite(bc,file = bcpath,col.names = F,row.names = F)

        #prepare annotation
        gtf <- paste0(opt\$out_dir,"/",opt\$project,".final_annot.gtf")

        #check if umis present in data
        UMI_check <- lapply(sequence_files, 
              function(x) {
                if(!is.null(x\$base_definition)) {
                  if(any(grepl("^UMI",x\$base_definition))) return("UMI method detected.")
                }
              })

        umi_decision <- ifelse(length(unlist(UMI_check))>0,"","--without-umi")
        mm_decision <- ifelse(opt\$counting_opts\$primaryHit,"--multimap","")

        #run velocyto

        print("Attempting to run RNA velocity...")
        velo_check <- suppressWarnings(system("which velocyto",intern =T))
        if(length(velo_check) == 0){
          print("No velocyto installation found in path. Please install it via pip.")
        }else{
          velo_cmd <- paste(velo_check[1],"run -vv --umi-extension Gene -b",bcpath,"-o",paste0(opt\$out_dir,"zUMIs_output/velocity/"),umi_decision,mm_decision, "-e",opt\$project,"--samtools-threads",opt\$num_threads,"--samtools-memory",mempercpu*1000,velobam,gtf,sep=" ")
          
          try(system(paste0(velo_cmd," > ",opt\$out_dir,opt\$project,".velocityo.log.out 2>&1")))
          print("RNA velocity done!")
          
        }

        print(Sys.time())

        q(save = "no")

        """

    }
}

process summary_stats {
    // Create summary statistics on read composition, barcode counts, etc.
    cpus = "${params.num_threads}"
    memory = "${params.mem_limit}G"
    time { 1.hour * task.attempt }

    input:
    file bam from UB_corrected_bam_ch2
    val read_layout from read_layout_ch3.first()
    file barcodes from kept_barcodes_ch2
    file binned_BC from kept_bc_binned_ch2
    file BCstats from cat_bcstats_ch2
    file counts from dgecounts_ch1

    output:

    """
    #!/usr/bin/env Rscript
    print("I am loading useful packages for plotting...")
    print(Sys.time())

    library(methods)
    library(data.table)
    library(yaml)
    library(ggplot2)
    library(Matrix)
    suppressMessages(library(dplyr))
    suppressMessages(library(Rsamtools))
    suppressMessages(library(cowplot))
    suppressMessages(require(R.utils))
    options(datatable.fread.input.cmd.message=FALSE)
 
    # Read in parameters from config file
    opt <- NULL
    opt\$zUMIs_directory <- "${params.projectDir}${params.binDir}"
    opt\$barcodes\$BarcodeBinning <- ${params.barcode_binning}
    opt\$barcodes\$nReadsperCell <- ${params.n_reads_per_cell}
    opt\$counting_opts\$introns <- ${params.include_introns}
    opt\$out_dir <- "${params.projectDir}${params.outputDir}"
    opt\$project <- "${params.projectName}"
    opt\$num_threads <- ${params.num_threads}
    opt\$mem_limit <- ${params.mem_limit}
    opt\$read_layout <- "${read_layout}"
    opt\$BCcount <- "${BCstats}"
    # read in input csv
    input<-read.csv("${params.input_csv}")

    #convert csv to yaml format
    sequence_files <- as.list(NULL)
    for (i in 1:nrow(input)){
    name <- input\$name[i]
    cDNA <- if (is.na(input\$cDNA[i])) NULL else paste0("cDNA(",input\$cDNA[i],")")
    BC <- if (is.na(input\$BC[i])) NULL else paste0("BC(",input\$BC[i],")")
    UMI <- if (is.na(input\$UMI[i])) NULL else paste0("UMI(",input\$UMI[i],")")
    if (cDNA == "cDNA()") cDNA <- NULL
    if (BC == "BC()") BC <- NULL
    if (UMI == "UMI()") UMI <- NULL
    base_definition <- c(cDNA,BC,UMI)
    if (is.null(base_definition)) base_definition <- NULL
    filter_cutoffs <- input\$filter_cutoffs[i]
    correct_frameshift <- input\$correct_frameshift[i]
    if (is.na(filter_cutoffs)) filter_cutoffs <- NULL
    if (is.na(correct_frameshift)) correct_frameshift <- NULL
    sequence_files[[paste0("file",i)]] <- list(name,base_definition,filter_cutoffs,correct_frameshift)
    names(sequence_files[[paste0("file",i)]]) <- (c("name","base_definition","filter_cutoffs","correct_frameshift"))
    }

    samtoolsexc <- "samtools"
    # set up graphic colour schemes
    featColors<-c("#1A5084", "#914614" ,"#118730","grey33","tan1","#631879FF","gold1","grey73","firebrick3")
    names(featColors)<-c("Exon","Intron+Exon","Intron","Unmapped","Ambiguity","MultiMapping","Intergenic","Unused BC","User")


    source(paste0(opt\$zUMIs_directory,"/statsFUN.R"))

    data.table::setDTthreads(threads=opt\$num_threads)

    user_seq<- getUserSeq(paste0(opt\$out_dir,"/",opt\$project,".final_annot.gtf")) 
    bc<-data.table::fread("${barcodes}",select = 1, header = T)
    AllCounts<-readRDS("${counts}")


    featfile <- "${bam}"

    #check if PE / SE flag is set correctly
    if(is.null(opt\$read_layout)){
      opt\$read_layout <- check_read_layout(featfile)
    }
    # UMI fragment counts
    if(any(grepl(pattern = "ATTGCGCAATG",x = unlist(sequence_files)))){
      print("Counting UMI fragments...")
      script_filepath <- paste0(opt\$zUMIs_directory,"/misc/countUMIfrags.py")
      bam_filepath <- featfile
      if(opt\$barcodes\$BarcodeBinning > 0){
        bc_filepath <- "${binned_BC}"
      }else{
        bc_filepath <- "${barcodes}"
      }
      system(paste(script_filepath,'--bam',bam_filepath,'--bcs',bc_filepath,'--p',opt\$num_threads,"&"))
    }

    # Gene counts

    genecounts <- suppressWarnings(dplyr::bind_rows( lapply(names(AllCounts\$readcount), function(i){
                          countGenes(AllCounts\$readcount[[i]][["all"]], user_seq=user_seq) %>%
                          mutate(type=case_when( i == "exon" ~ "Exon",
                                                  i == "inex" ~ "Intron+Exon",
                                                  i == "intron" ~ "Intron")) })))

    umicounts <- suppressWarnings(dplyr::bind_rows( lapply(names(AllCounts\$umicount), function(i){
                            countUMIs(AllCounts\$umicount[[i]][["all"]], user_seq=user_seq) %>%
                            mutate(type=case_when( i == "exon" ~ "Exon",
                                                    i == "inex" ~ "Intron+Exon",
                                                    i == "intron" ~ "Intron"))})))

    med<-genecounts %>% dplyr::group_by(type) %>% dplyr::summarise(n=round(median(Count)))
    write.table(genecounts,file = paste0(opt\$out_dir,"zUMIs_output/stats/",opt\$project,".genecounts.txt"),sep="\t",row.names = F,col.names = T)

    ag <- countBoxplot(cnt = genecounts,
                      ylab= "Number of Genes",
                      fillcol=featColors[unique(genecounts\$type)],
                      lab = med)

    if(length(umicounts) > 0){
      medUMI<-try(umicounts %>% dplyr::group_by(type) %>% dplyr::summarise(n=round(median(Count))))
      try(write.table(umicounts,file = paste0(opt\$out_dir,"zUMIs_output/stats/",opt\$project,".UMIcounts.txt"),sep="\t",row.names = F,col.names = T))
      bg <- try(countBoxplot(cnt = umicounts,
                            ylab= "Number of UMIs",
                            fillcol=featColors[unique(umicounts\$type)],
                            lab = medUMI))
      cp<-try(cowplot::plot_grid(ag,bg,ncol = 2))
    }else{
      cp <- ag
    }

    try(ggsave(cp,filename = paste(opt\$out_dir,"zUMIs_output/stats/",opt\$project,".geneUMIcounts.pdf",sep=""),width = 10,height = 5))

    ## Calculate total read count per gene
    reads_per_gene <- sumGene(counts = AllCounts)
    data.table::fwrite(reads_per_gene, file = paste0(opt\$out_dir,"zUMIs_output/stats/",opt\$project,".reads_per_gene.txt"), quote = F, sep = "\t")

    ## Calculate read count per cell

    typeCount <- sumstatBAM( featfile = featfile,
                            cores = opt\$num_threads,
                            outdir= opt\$out_dir,
                            user_seq = user_seq,
                            bc = bc,
                            inex = opt\$counting_opts\$introns,
                            outfile = paste0(opt\$out_dir,"zUMIs_output/stats/",opt\$project,".bc.READcounts.rds"),
                            samtoolsexc=samtoolsexc)

    #only print per BC mapping stats if there are fewer than 200 BCs
    tc<-data.frame(typeCount)
    tc\$type<-factor(tc\$type, levels=rev(c("Exon","Intron","Intergenic","Ambiguity","Unmapped","User")))
    write.table(tc,file = paste(opt\$out_dir,"zUMIs_output/stats/",opt\$project,".readspercell.txt",sep=""),sep="\t",row.names = F,col.names = T)


    bar<-totReadCountBarplot(typeCount = typeCount,
                            fillcol= featColors)


    box<-totReadBoxplot(typeCount = typeCount,
                        fillcol= featColors)

    d<-plot_grid(cp,bar,box,ncol = 1,rel_heights  = c(0.3,0.2,0.5))

    ggsave(d,filename = paste(opt\$out_dir,"zUMIs_output/stats/",opt\$project,".features.pdf",sep=""),width = 12,height = 9)


    ###############
    system("wait")
    gc()
    Sys.time()
    q(save = "no")

    """
}