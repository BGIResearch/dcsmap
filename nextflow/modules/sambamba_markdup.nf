// ============================================================================
// Sambamba markdup module
//   Marks PCR/optical duplicates in the sorted BAM.
// ============================================================================

process SAMBAMBA_MARKDUP {
    tag "${sample_id}"
    cpus params.threads
    publishDir "${params.outdir}/markdup", mode: 'copy'
    publishDir { "${params.outdir}/logs/${task.process}" }, mode: 'copy', overwrite: true, pattern: '.command.*'

    input:
    tuple val(sample_id), path(bam), path(bai)

    output:
    tuple val(sample_id), path("${sample_id}.markdup.bam"), path("${sample_id}.markdup.bam.bai"), emit: bam
    path '.command.log', emit: log
    path '.command.sh',  emit: script

    script:
    """
    set -euo pipefail

    mkdir -p tmp_markdup

    echo "[INFO] `date "+%Y-%m-%d %H:%M:%S"` Start sambamba markdup ..."
    /usr/bin/time -f "[sambamba-markdup] real %e\\tuser %U\\tsys %S\\tcpu_percentage %P\\tmaxrss %M" \\
        sambamba markdup -t ${task.cpus} --tmpdir=tmp_markdup ${bam} ${sample_id}.markdup.bam
    echo "[INFO] `date "+%Y-%m-%d %H:%M:%S"` End sambamba markdup."
    """
}
