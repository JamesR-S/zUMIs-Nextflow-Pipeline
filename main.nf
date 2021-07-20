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
      awk -v rownum=\${i} -F "\\"*,\\"*" 'FNR==rownum {print \$1}' ${input_csv} | xargs -I{} ln {} .
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

    if [! -d "${params.projectDir}${params.outputDir}.\$(echo "${file}" | cut -f 1 -d '.')_tmpMerge"] ; then
        mkdir ${params.projectDir}${params.outputDir}.\$(echo "${file}" | cut -f 1 -d '.')_tmpMerge
    fi

    n=`expr $nreads / ${params.num_threads}`
    n=`expr $n + 1`
    nl=`expr $n \* 4`

    if [[ ${file} =~ \.gz$ ]] ; then
      pigz -dc ${file} | head -n 4000000 | pigz > ${params.projectDir}${params.outputDir}.\$(echo "${file}" | cut -f 1 -d '.')_tmpMerge/${file}.1mio.check.fq.gz
      smallsize=$(stat --printf="%s" ${params.projectDir}${params.outputDir}.\$(echo "${file}" | cut -f 1 -d '.')_tmpMerge/${file}.1mio.check.fq.gz)
      rm ${params.projectDir}${params.outputDir}.\$(echo "${file}" | cut -f 1 -d '.')_tmpMerge/${file}.1mio.check.fq.gz
      nreads=$(expr \${fullsize} \* 1000000 / \${smallsize})
      n=`expr $nreads / ${params.num_threads}`
      n=`expr $n + 1`
      nl=`expr $n \* 4`
      pigz -dc -p ${params.num_threads} ${file} | split --lines=\$nl --filter=''pigz' -p '${params.num_threads}' > \$FILE.gz' - ${params.projectDir}${params.outputDir}.\$(echo "${file}" | cut -f 1 -d '.')_tmpMerge/\$(basename ${file} .gz)${params.projectName}
    else
      cat ${file} | head -n 4000000 > ${params.projectDir}${params.outputDir}.\$(echo "${file}" | cut -f 1 -d '.')_tmpMerge/${file}.1mio.check.fq
      smallsize=$(stat --printf="%s" ${params.projectDir}${params.outputDir}.\$(echo "${file}" | cut -f 1 -d '.')_tmpMerge/${file}.1mio.check.fq)
      rm ${params.projectDir}${params.outputDir}.\$(echo "${file}" | cut -f 1 -d '.')_tmpMerge/${file}.1mio.check.fq
      nreads=$(expr \${fullsize} \* 1000000 / \${smallsize})
      n=`expr $nreads / ${params.num_threads}`
      n=`expr $n + 1`
      nl=`expr $n \* 4`
      split --lines=\$nl --filter=''pigz' -p '${params.num_threads}' > \$FILE.gz' ${file} ${params.projectDir}${params.outputDir}.\$(echo "${file}" | cut -f 1 -d '.')_tmpMerge/${file}${params.projectName}
    fi

    ls ${params.projectDir}${params.outputDir}.\$(echo "${file}" | cut -f 1 -d '.')_tmpMerge > \$(echo "${file}" | cut -f 1 -d '.')_tempfiles.txt
    """



}