// ============================================================================
// GATK Base Quality Score Recalibration (BQSR) module
//   BaseRecalibrator -> ApplyBQSR -> sambamba index
//   Produces an analysis-ready recalibrated BAM (sample_id.bqsr.bam + bai).
// ============================================================================

process GATK_BQSR {
    tag "${sample_id}"
    cpus params.threads
    // publishDir "${params.outdir}/gatk/bqsr", mode: 'copy'
    publishDir { "${params.outdir}/logs/${task.process}" }, mode: 'copy', overwrite: true, pattern: '.command.*'

    input:
    tuple val(sample_id), path(bam), path(bai)
    path ref
    path ref_idx
    path ref_dict
    path known_sites_files
    path known_sites_indices

    output:
    tuple val(sample_id), path("${sample_id}.bqsr.bam"), path("${sample_id}.bqsr.bam.bai"), emit: bam
    path '.command.log', emit: log
    path '.command.sh',  emit: script

    script:
    """
    set -euo pipefail

    ks_args=""
    for ks in ${known_sites_files}; do
        ks_args="\${ks_args} --known-sites \${ks}"
    done

    echo "[INFO] `date "+%Y-%m-%d %H:%M:%S"` Start GATK BaseRecalibrator ..."
    /usr/bin/time -f "[gatk-BaseRecalibrator] real %e\\tuser %U\\tsys %S\\tcpu_percentage %P\\tmaxrss %M" \\
        ${params.gatk} BaseRecalibrator \\
            -R ${ref} \\
            -I ${bam} \\
            \${ks_args} \\
            -O ${sample_id}.recal_data.table
    echo "[INFO] `date "+%Y-%m-%d %H:%M:%S"` End GATK BaseRecalibrator."

    echo "[INFO] `date "+%Y-%m-%d %H:%M:%S"` Start GATK ApplyBQSR ..."
    /usr/bin/time -f "[gatk-ApplyBQSR] real %e\\tuser %U\\tsys %S\\tcpu_percentage %P\\tmaxrss %M" \\
        ${params.gatk} ApplyBQSR \\
            -R ${ref} \\
            -I ${bam} \\
            --bqsr-recal-file ${sample_id}.recal_data.table \\
            -O ${sample_id}.bqsr.bam
    echo "[INFO] `date "+%Y-%m-%d %H:%M:%S"` End GATK ApplyBQSR."

    echo "[INFO] `date "+%Y-%m-%d %H:%M:%S"` Start sambamba index ..."
    /usr/bin/time -f "[sambamba-index] real %e\\tuser %U\\tsys %S\\tcpu_percentage %P\\tmaxrss %M" \\
        sambamba index -t ${task.cpus} ${sample_id}.bqsr.bam
    echo "[INFO] `date "+%Y-%m-%d %H:%M:%S"` End sambamba index."
    """
}
