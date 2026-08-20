// ============================================================================
// BWA-MEM2 alignment module
//   bwa-mem2 mem -> samtools sort : produces a sorted, indexed BAM.
// ============================================================================

process BWA_MEM2 {
    tag "${sample_id}"
    cpus params.threads
    // publishDir "${params.outdir}/mapping", mode: 'copy'
    publishDir { "${params.outdir}/logs/${task.process}" }, mode: 'copy', overwrite: true, pattern: '.command.*'

    input:
    tuple val(sample_id), path(fq1), path(fq2)
    path ref
    path ref_idx

    output:
    tuple val(sample_id), path("${sample_id}.sorted.bam"), path("${sample_id}.sorted.bam.bai"), emit: bam
    path '.command.log', emit: log
    path '.command.sh',  emit: script

    script:
    def rg = "@RG\\tID:${sample_id}\\tSM:${sample_id}\\tPL:ILLUMINA\\tLB:lib1"
    """
    set -euo pipefail

    echo "[INFO] `date "+%Y-%m-%d %H:%M:%S"` Start bwa-mem2 mapping and samtools sort ..."
    /usr/bin/time -f "[bwa-mem2+samtools-sort] real %e\\tuser %U\\tsys %S\\tcpu_percentage %P\\tmaxrss %M" \\
        bwa-mem2 mem -t ${task.cpus} -R '${rg}' ${ref} ${fq1} ${fq2} | \\
        samtools sort -@ ${task.cpus} -m1g -OBAM -o ${sample_id}.sorted.bam##idx##${sample_id}.sorted.bam.bai --write-index
    echo "[INFO] `date "+%Y-%m-%d %H:%M:%S"` End bwa-mem2 mapping and samtools sort."
    """
}
