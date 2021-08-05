#!/usr/bin/env nextflow

process copyResourceFiles {
 
    cpus = '4'
    mem = '20G'
    time = '00:05:00'

    output:
      file "${params.fastaName}" into fasta_ch
      file "${params.gtfName}" into gtf_ch

    """
    cp ${params.genome} .
    cp ${params.gtf} .

    gunzip *.gz

    """

}

process generateIndex {

    cpus = '6'
    mem = '20G'
    time = '02:00:00'

    input:
        val fasta from fasta_ch
        val gtf from gtf_ch

    output:
        path "${params.genomeName}" into genome_index_ch

    """

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

    """
}

process readFastqFiles {

    cpus = '4'
    mem = '20G'
    time = '00:05:00'

    output:
        file '*' into fastqs_ch

    """

    nrow=\$(wc -l ${params.input_csv} | cut -d' ' -f1)
    length=`expr \$nrow + 1`
    for i in \$( seq 2 \$length )
    do
      awk -v rownum="\${i}" -F "\\"*,\\"*" 'FNR==rownum {print \$1}' ${params.input_csv} | xargs -I{} cp {} .
    done

    """
}

process splitfq {

    cpus = '32'
    mem = '20G'
    time = '00:15:00'

    input:
    each file from fastqs_ch

    output:
    file '*_tempfiles.txt' into tempfiles_ch1
    file '*_tempfiles.txt' into tempfiles_ch2


    """

    fullsize=\$(stat -L --printf="%s" ${file})

    basefile=\$(basename ${file} | cut -f 1 -d '.')

    if [ ! -d "${params.projectDir}${params.outputDir}.\${basefile}_tmpMerge" ]; then
        mkdir ${params.projectDir}${params.outputDir}.\${basefile}_tmpMerge
    fi

    if [[ ${file} =~ \\.gz\$ ]]; then
      pigz -dc ${file} | head -n 4000000 | pigz > ${params.projectDir}${params.outputDir}.\${basefile}_tmpMerge/\${basefile}.1mio.check.fq.gz
      smallsize=\$(stat --printf="%s" ${params.projectDir}${params.outputDir}.\${basefile}_tmpMerge/\${basefile}.1mio.check.fq.gz)
      rm ${params.projectDir}${params.outputDir}.\${basefile}_tmpMerge/\${basefile}.1mio.check.fq.gz
      nreads=\$(expr \${fullsize} \\* 1000000 / \${smallsize})
      n=`expr \$nreads / ${params.num_threads}`
      n=`expr \$n + 1`
      nl=`expr \$n \\* 4`
      pigz -dc -p ${params.num_threads} ${file} | split --lines=\$nl --filter=''pigz' -p '${params.num_threads}' > \$FILE.gz' - ${params.projectDir}${params.outputDir}.\${basefile}_tmpMerge/\${basefile}${params.projectName}
    else
      cat ${file} | head -n 4000000 > ${params.projectDir}${params.outputDir}.\${basefile}_tmpMerge/\${basefile}.1mio.check.fq
      smallsize=\$(stat --printf="%s" ${params.projectDir}${params.outputDir}.\${basefile}_tmpMerge/\${basefile}.1mio.check.fq)
      rm ${params.projectDir}${params.outputDir}.\${basefile}_tmpMerge/\${basefile}.1mio.check.fq
      nreads=\$(expr \${fullsize} \\* 1000000 / \${smallsize})
      n=`expr \$nreads / ${params.num_threads}`
      n=`expr \$n + 1`
      nl=`expr \$n \\* 4`
      split --lines=\$nl --filter=''pigz' -p '${params.num_threads}' > \$FILE.gz' ${file} ${params.projectDir}${params.outputDir}.\${basefile}_tmpMerge/\${basefile}${params.projectName}
    fi

    ls ${params.projectDir}${params.outputDir}.\${basefile}_tmpMerge/\${basefile}* | sed "s|${params.projectDir}${params.outputDir}.\${basefile}_tmpMerge/\${basefile}||"  > \${basefile}_tempfiles.txt
    """
}

process fqFilter {

    cpus = '32'
    mem = '20G'
    time = '00:30:00'

    input:
    each file from tempfiles_ch1

    output:
    file "*.bamlist.txt" into bamlists_ch
    env read_layout into read_layout_ch
    file "*.BCstats.txt" into bcstats_ch


    """

    IFS=\$'\\r\\n' GLOBIGNORE='*' command eval  'temp_files=(\$(cat ${file}))'

    filename=\$(basename ${file})

    parentfq=\${filename%"_tempfiles.txt"}

    basefile=\$(echo \${parentfq} | cut -f 1 -d '.')

    for x in "\${temp_files[@]}" 
    do 
        perl ${params.projectDir}${params.binDir}fqfilter_v2.pl ${params.input_csv} samtools Rscript pigz ${params.projectDir}${params.binDir} \${x} ${params.UMI_phred} \\
        ${params.projectName} ${params.num_threads} \${parentfq} \${basefile} ${params.BC_num_bases} ${params.BC_phred} ${params.UMI_num_bases} ${params.projectDir}${params.outputDir}
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


    """
}

process collectBCstats {

    cpus = '32'
    mem = '20G'
    time = '00:05:00'

    input:
    val files from bcstats_ch.collect()

    output:
    file "${params.projectName}.BCstats.txt" into cat_bcstats_ch

    """
    bcfiles=(${files.join(" ")})

    for f in "\${bcfiles[@]}"
    do
        cat \${f} >> ${params.projectName}.BCstats.txt
    done
    """

}

process barcode_detect {

    cpus = '32'
    mem = '40G'
    time = '04:30:00'

    input:
    file file from cat_bcstats_ch

    output:
    file "*.detected_cells.pdf" into detected_cells_ch
    file "*_kept_barcodes.txt" into kept_barcodes_ch
    file("*.BCbinning.txt") optional true into bc_binning_ch
    file("*kept_barcodes_binned.txt") optional true into kept_bc_binned_ch


    """
    #!/usr/bin/env Rscript

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

    
    barcode_sharing <- "${params.barcode_sharing}"
    barcode_file <- "${params.barcode_file}"
    if (barcode_file == "") barcode_file <- NULL
    barcode_num <- "${params.barcode_num}"
    if (barcode_num == "") barcode_num <- NULL
    bcauto <- as.logical("${params.barcode_automatic}")
    out_dir <- "${params.projectDir}${params.outputDir}"
    project <- "${params.projectName}"
    zUMIs_directory <- "${params.projectDir}${params.binDir}"
    n_reads_per_cell <- "${params.n_reads_per_cell}"
    barcode_binning <- "${params.barcode_binning}"
    num_threads <- "${params.num_threads}"
    barcode_sharing <- "${params.barcode_sharing}"
    if (barcode_sharing == "") barcode_sharing <- NULL

    input <- read.csv("${params.input_csv}")

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
    #######################################################################
    ##### Barcode handling & chunking

    #read file with barcodecounts
    # bc is the vector of barcodes to keep
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

    ##############################################################
    q()    

    """


}

bc_binning_ch = bc_binning_ch.ifEmpty('EMPTY')
kept_bc_binned_ch = kept_bc_binned_ch.ifEmpty('EMPTY')

process BC_bin_bam_collate {

    cpus = '32'
    mem = '40G'
    time = '01:30:00'

    input:
    val BCbinning from bc_binning_ch
    each file from tempfiles_ch2

    output:
    env tmpFILT into collated_bam_ch

    """
    IFS=\$'\\r\\n' GLOBIGNORE='*' command eval  'temp_files=(\$(cat ${file}))'

    filename=\$(basename ${file})

    parentfq=\${filename%"_tempfiles.txt"}

    basefile=\$(echo \${parentfq} | cut -f 1 -d '.')

    [ -d ${params.projectDir}${params.outputDir}.tmpFiltered ] || mkdir ${params.projectDir}${params.outputDir}.tmpFiltered

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
    tmpFILT="${params.projectDir}${params.outputDir}.tmpFiltered"
    """
}

process Mapping {

    cpus = '32'
    mem = '40G'
    time = '04:30:00'

    input:
    val tmpFILT from collated_bam_ch.last()
    val barcodes from kept_barcodes_ch
    val barcodes_binned from kept_bc_binned_ch
    path index from genome_index_ch
    val layout from read_layout_ch.last()

    output:
    file "*.filtered.tagged.Log.final.out" into mapping_log_ch
    file "*.filtered.tagged.Aligned.out.bam" into aligned_bam_ch
    file "*.filtered.tagged.Aligned.toTranscriptome.out.bam" into txbam_ch


    """
    #!/usr/bin/env Rscript

    .libPaths(c("/usr/local/lib/R/site-library","/usr/lib/R/site-library","/usr/lib/R/library"))

    suppressMessages(require(data.table))
    suppressMessages(require(R.utils))
    options(datatable.fread.input.cmd.message=FALSE)
    Sys.time()
    inp <- commandArgs(trailingOnly = T, asValues=T)

    samtools <- "samtools"
    STAR_exec <- "STAR"

    mem_limit = ${params.mem_limit}
    out_dir = "${params.projectDir}${params.outputDir}"
    project = "${params.projectName}"
    STAR_index = "${index}"
    num_threads = ${params.num_threads}
    GTF_file = "${params.gtf}"
    read_layout = "${layout}"
    twoPass = ${params.two_pass}
    additional_STAR_params = "${params.additional_STAR_params}"
    if (identical(additional_STAR_params,"NONE")){
    additional_STAR_params <- NULL
    }
    additional_fq = "${params.additional_fq}"
    if (identical(additional_fq,"NONE")){
    additional_fq <- NULL
    }else{
        additional_fq = strsplit(additional_fq,",")
    }
    if(is.null(mem_limit)){
    mem_limit <- 100
    }else if(mem_limit == 0){
    mem_limit <- 100
    }

    # collect filtered bam files ----------------------------------------------
    tmpfolder <- "${tmpFILT}"
    filtered_bams <- list.files(path = tmpfolder, pattern=paste(project,".*.filtered.tagged.bam",sep=""),full.names=T)
    #also merge the unmapped bam files:
    sammerge_command <- paste(samtools,"cat -o",paste0(out_dir,project,".filtered.tagged.unmapped.bam"),paste0(filtered_bams,collapse=" "))
    system(sammerge_command)


    # check if multiple STAR instances can be run -----------------------------

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

    # GTF file setup ----------------------------------------------------------
    #in case of additional sequences, we need to create a custom GTF

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
    #yaml::write_yaml(inp,file = paste(out_dir,project,".postmap.yaml",sep=""))

    # Detect read length ------------------------------------------------------
    #check the first 100 reads to detect the read length of the cDNA read
    filtered_bam <- paste(out_dir,project,".filtered.tagged.unmapped.bam",sep="")

    cDNA_peek <- data.table::fread(cmd = paste(samtools,"view",filtered_bam,"| cut -f10 | uniq | head -n 1000"),stringsAsFactors = F,data.table = T, header = F)

    getmode <- function(v) {
    uniqv <- unique(v)
    uniqv[which.max(tabulate(match(v, uniqv)))]
    }

    cDNA_read_length <- getmode(nchar(cDNA_peek\$V1))


    # Setup STAR mapping ------------------------------------------------------
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

    #finally, run STAR
    if(num_star_instances>1){
    map_tmp_dir <- paste0(out_dir,"zUMIs_output/.tmpMap/")
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
    system(paste(STAR_command,"& wait"))
    
    #after parallel instance STAR, collect output data in the usual file places
    out_logs <- list.files(map_tmp_dir, pattern = paste0("tmp.",project,".*.Log.final.out"), full = TRUE)
    merge_logs <- paste("cat",paste(out_logs, collapse = " "),">",paste0(project,".filtered.tagged.Log.final.out"))
    out_bams <- list.files(map_tmp_dir, pattern = paste0("tmp.",project,".*.Aligned.out.bam"), full = TRUE)
    merge_bams <- paste(samtools,"cat -o",paste0(project,".filtered.tagged.Aligned.out.bam"),paste(out_bams, collapse = " "))
    out_txbams <- list.files(map_tmp_dir, pattern = paste0("tmp.",project,".*.Aligned.toTranscriptome.out.bam"), full = TRUE)
    merge_txbams <- paste(samtools,"cat -o",paste0(project,".filtered.tagged.Aligned.toTranscriptome.out.bam"),paste(out_txbams, collapse = " "))
    system(paste(merge_logs,"&",merge_bams,"&",merge_txbams,"& wait"))
    system(paste0("rm -r ", map_tmp_dir, "tmp.", project, ".*"))
    }else{
    STAR_command <- paste(STAR_command,
        "--readFilesIn",paste0(filtered_bams,collapse=","),
        "--outFileNamePrefix",paste(project,".filtered.tagged.",sep="")
    )

        system(paste(STAR_command,"& wait"))

    }

    """


}


