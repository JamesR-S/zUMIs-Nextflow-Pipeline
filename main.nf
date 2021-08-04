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
    #!/usr/bin/env Rscript --no-environ

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

bc_binning_ch.ifEmpty('EMPTY')

process BC_bin/bam_collate {

    cpus = '32'
    mem = '40G'
    time = '04:30:00'

    input:
    val BCbinning from bc_binning_ch
    each file from tempfiles_ch1

    output:
    file "*.detected_cells.pdf" into detected_cells_ch
    file "*_kept_barcodes.txt" into kept_barcodes_ch
    file("*.BCbinning.txt") optional true into bc_binning_ch
    file("*kept_barcodes_binned.txt") optional true into kept_bc_binned_ch


    """
    IFS=\$'\\r\\n' GLOBIGNORE='*' command eval  'temp_files=(\$(cat ${file}))'

    filename=\$(basename ${file})

    parentfq=\${filename%"_tempfiles.txt"}

    basefile=\$(echo \${parentfq} | cut -f 1 -d '.')

    if [[ -f "${BCbinning}" ]] ; then
        for x in "\${temp_files[@]}" ; do
            rawbam="${params.projectDir}${params.outputDir}.\${basefile}_tmpMerge/${params.projectName}.${x}.raw.tagged.bam"
            fixedbam="${params.projectDir}${params.outputDir}.\${basefile}_tmpMerge/${params.projectName}.${x}.filtered.tagged.bam"
            mv ${fixedbam} ${rawbam}
            perl ${params.projectDir}${params.binDir}correct_BCtag.pl \${rawbam} \${fixedbam} ${BCbinning} samtools &
        done
        wait
    fi

    for x in "\${temp_files[@]}" ; do
        mv ${params.projectDir}${params.outputDir}.\${basefile}_tmpMerge/${params.projectName}.${x}.filtered.tagged.bam ${params.projectDir}${params.outputDir}.tmpFiltered/${params.projectName}.\${basefile}.${x}.filtered.tagged.bam
    done
    """
    
