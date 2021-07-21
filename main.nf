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

    gunzip ${params.projectDir}${params.resourcesDir}/*.gz

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

    Overhang=$((${params.read_length} - 1))

    if [ ! -d "${params.genomeName}" ]; then
    mkdir ${params.genomeName}
    fi

    STAR --runThreadN 6 \
    --runMode genomeGenerate \
    --genomeDir ${params.genomeName} \
    --genomeFastaFiles ${fasta} \
    --sjdbGTFfile ${gtf} \
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

    nrow=\$(wc -l ${input_csv})
    for i in \$( seq 2 \$nrow )
    do
      awk -v rownum=\${i} -F "\"*,\"*" 'FNR==rownum {print \$1}' ${input_csv} | xargs -I{} ln {} .
    done

    """
}

process splitfq {

    cpus = '${params.num_threads}'
    mem = '20G'
    time = '00:05:00'

    input:
    each file from fastqs_ch

    output:
    file '*_tempfiles.txt' into tempfiles_ch


    """

    fullsize=\$(stat -L --printf="%s" ${file})

    basefile=\$(basename ${file} | cut -f 1 -d '.')

    if [! -d "${params.projectDir}${params.outputDir}.\${basefile}_tmpMerge"] ; then
        mkdir ${params.projectDir}${params.outputDir}.\${basefile}_tmpMerge
    fi

    n=`expr $nreads / ${params.num_threads}`
    n=`expr $n + 1`
    nl=`expr $n \* 4`

    if [[ ${file} =~ \.gz$ ]] ; then
      pigz -dc ${file} | head -n 4000000 | pigz > ${params.projectDir}${params.outputDir}.\${basefile}_tmpMerge/\${basefile}.1mio.check.fq.gz
      smallsize=$(stat --printf="%s" ${params.projectDir}${params.outputDir}.\${basefile}_tmpMerge/\${basefile}.1mio.check.fq.gz)
      rm ${params.projectDir}${params.outputDir}.\${basefile}_tmpMerge/\${basefile}.1mio.check.fq.gz
      nreads=$(expr \${fullsize} \* 1000000 / \${smallsize})
      n=`expr $nreads / ${params.num_threads}`
      n=`expr $n + 1`
      nl=`expr $n \* 4`
      pigz -dc -p ${params.num_threads} ${file} | split --lines=\$nl --filter=''pigz' -p '${params.num_threads}' > \$FILE.gz' - ${params.projectDir}${params.outputDir}.\${basefile}_tmpMerge/\${basefile}${params.projectName}
    else
      cat ${file} | head -n 4000000 > ${params.projectDir}${params.outputDir}.\${basefile}_tmpMerge/\${basefile}.1mio.check.fq
      smallsize=$(stat --printf="%s" ${params.projectDir}${params.outputDir}.\${basefile}_tmpMerge/\${basefile}.1mio.check.fq)
      rm ${params.projectDir}${params.outputDir}.\${basefile}_tmpMerge/\${basefile}.1mio.check.fq
      nreads=$(expr \${fullsize} \* 1000000 / \${smallsize})
      n=`expr $nreads / ${params.num_threads}`
      n=`expr $n + 1`
      nl=`expr $n \* 4`
      split --lines=\$nl --filter=''pigz' -p '${params.num_threads}' > \$FILE.gz' ${file} ${params.projectDir}${params.outputDir}.\${basefile}_tmpMerge/\${basefile}${params.projectName}
    fi

    ls ${params.projectDir}${params.outputDir}.\${basefile}_tmpMerge/\${basefile}* | sed "s|${params.projectDir}${params.outputDir}.\${basefile}_tmpMerge/\${basefile}*||"  > \${basefile}_tempfiles.txt
    """
}

process fqFilter {

    cpus = '${params.num_threads}'
    mem = '20G'
    time = '00:05:00'

    input:
    each file from tempfiles_ch

    output:
    file "*.bamlist.txt" into bamlists_ch
    env read_layout into read_layout_ch
    file "*.bcstats.txt" into bcstats_ch


    """
    tIFS=\$IFS
    tGLOBIGNORE=\$GLOBIGNORE
    IFS=\$'\r\n' GLOBIGNORE='*' command eval  'temp_files=(\$(cat ${file}))'
    IFS=\$tIFS
    GLOBIGNORE=\$tGLOBIGNORE

    filename=\$(basename ${file})

    parentfq=\${filename%"_tempfiles.txt"}

    for x in \${temp_files} ; do perl ${zumisdir}/fqfilter_v2.pl ${params.input_csv} samtools rscript pigz ${params.projectDir}${params.binDir} ${x} ${params.projectDir}${params.outputDir} ${params.projectName} ${params.num_threads} \${parentfq} & done
    wait
    basefile=\$(\${filename} | cut -f 1 -d '.')
    ls ${params.projectDir}${params.outputDir}.\${basefile}_tmpMerge/${params.projectName}.*.filtered.tagged.bam > \${basefile}_${params.projectName}.bamlist.txt
    
    f=`head -n1 \${basefile}_${params.projectName}.bamlist.txt`
    flag=`$samtoolsexc view $f | head -n1 | cut -f2`

    if [[ $flag == 4 ]]; then
    read_layout="SE"
    else
    read_layout="PE"
    fi


    """
}

process collectBCstats {

    cpus = '${params.num_threads}'
    mem = '20G'
    time = '00:05:00'

    input:
    val files from bcstats_ch.collect()

    output:
    file "${params.projectName}.BCstats.txt" into cat_bcstats_ch

    """
    bcfiles=(${files.join(" ")})

    for f in "\${bcfiles}"
    do
        cat ${f} >> ${params.projectName}.BCstats.txt
    """

}

process barcode_detect {

    cpus = '${params.num_threads}'
    mem = '20G'
    time = '00:05:00'

    input:
    file file from cat_bcstats_ch

    output:
    file "*.detected_cells.pdf" into detected_cells_ch
    file "*.kept_barcodes.txt" into kept_barcodes_ch
    file("*.BCbinning.txt") optional true into bc_binning_ch
    file("*kept_barcodes_binned.txt") optional true into kept_bc_binned_ch


    """
    #!/usr/bin/env rscript

    library(methods)
    library(data.table)
    library(yaml)
    library(ggplot2)

    ##########################

    opt\$barcode_sharing <- "${params.barcode_sharing}"
    opt\$barcode_file <- "${params.barcode_file}"
    opt\$barcode_num <- "${params.barcode_num}"
    opt\$automatic <- ${params.barcode_automatic}
    opt\$out_dir <- "${params.projectDir}${params.outputDir}"
    opt\$project <- "${params.projectName}"


    source(paste0(opt\$zUMIs_directory,"/barcodeIDFUN.R"))
    options(datatable.fread.input.cmd.message=FALSE)
    data.table::setDTthreads(threads=opt\$num_threads)
    if(!is.null(opt\$barcode_sharing)){
    if(opt\$barcode_sharing == ""){
        opt\$barcode_sharing <- NULL
    }
    }

    #######################################################################
    #######################################################################
    ##### Barcode handling & chunking

    #read file with barcodecounts
    # bc is the vector of barcodes to keep
    bccount<-cellBC(bcfile      = opt\$barcode_file,
            bcnum       = opt\$barcode_num,
            bcauto      = opt\$automatic,
            bccount_file= "${file}",
            outfilename = paste0(opt\$project,".detected_cells.pdf"))

    fwrite(bccount,file=paste0(opt\$project,"kept_barcodes.txt"))

    #check if binning of adjacent barcodes should be run
    if(opt\$BarcodeBinning > 0 | !is.null(opt\$barcode_sharing)){
    binmap <- BCbin(bccount_file = "${file}",
                    bc_detected  = bccount)
    fwrite(binmap,file=paste0(opt\$project,".BCbinning.txt"))
    #update the number reads in BCcount table
    binmap_additional <- binmap[, .(addtl = sum(n)), by = trueBC]
    bccount[match(binmap_additional\$trueBC,XC),n := n + binmap_additional\$addtl]
    fwrite(bccount,file=paste0(opt\$project,"kept_barcodes_binned.txt"))
    }

    ##############################################################
    q()    

    """


}