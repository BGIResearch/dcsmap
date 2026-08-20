process BWA_MEM2_MAP_AND_EXTRACT {

    tag "$sample_name"
    cpus params.threads
    publishDir { "${params.outdir}/logs/${task.process}" }, mode: params.publish_mode, overwrite: true, pattern: '.command.*'

    input:
    val tools_root
    val fq1
    val fq2
    val ref_fasta
    val sample_name
    val platform
    val extract_model
    val threads

    output:
    path "${sample_name}.bwa-mem2.extract_1.fastq.gz", emit: extract_fq1
    path "${sample_name}.bwa-mem2.extract_2.fastq.gz", emit: extract_fq2
    path "${sample_name}.bwa-mem2.not-extract.bam",   emit: out_bam
    path "${sample_name}.bwa-mem2.not-extract.bam.bai", emit: out_bam_bai
    path '.command.log', emit: log
    path '.command.sh',  emit: script

    script:
    """
    set -euo pipefail
    export PATH=${tools_root}/libexec:\$PATH

    in_ref=${ref_fasta}
    sample=${sample_name}
    read_group="@RG\\tID:\${sample}\\tLB:\${sample}\\tSM:\${sample}\\tPL:${platform}"
    extract_prefix=\${sample}.bwa-mem2.extract
    out_bam=\${sample}.bwa-mem2.not-extract.bam

    echo "[INFO] `date "+%Y-%m-%d %H:%M:%S"` Start bwa-mem2 mapping and extract ..."
    /usr/bin/time -f "[bwa-mem2] real %e\\tuser %U\\tsys %S\\tcpu_percentage %P\\tmaxrss %M" \\
        bwa-mem2 mem -S -r 3 -y 0 -Y -t ${threads} -R \${read_group} \${in_ref} ${fq1} ${fq2} |\\
        mawk '/^@/ {print;next} {print \$0"\\tMX:Z:BWA-MEM2"}' |\\
        extract-bam -b 4000000 -t 12 -p \${extract_prefix} -M ${extract_model} |\\
        samtools sort -l 1 -OBAM --threads 4 -m1g -o \${out_bam}##idx##\${out_bam}.bai --write-index
    echo "[INFO] `date "+%Y-%m-%d %H:%M:%S"` End bwa-mem2 mapping and extract."
    """
}
