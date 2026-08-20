// ============================================================================
// GATK HaplotypeCaller module
//   1. GATK_HAPLOTYPECALLER : per-chromosome GVCF generation (-ERC GVCF)
//   2. MERGE_GVCFS          : GatherVcfs + IndexFeatureFile to merge per-chrom
//                             GVCFs into a single sample GVCF (ordered by ref dict)
// ============================================================================

process GATK_HAPLOTYPECALLER {
    tag "${sample_id}:${interval}"
    cpus params.threads
    publishDir "${params.outdir}/gatk/haplotypecaller/per_chrom", mode: 'copy'
    publishDir { "${params.outdir}/logs/${task.process}/${interval}" }, mode: 'copy', overwrite: true, pattern: '.command.*'

    input:
    tuple val(sample_id), path(bam), path(bai), val(interval)
    path ref
    path ref_idx
    path ref_dict

    output:
    tuple val(sample_id), path("${sample_id}.${interval}.g.vcf.gz"), path("${sample_id}.${interval}.g.vcf.gz.tbi"), emit: gvcf
    path '.command.log', emit: log
    path '.command.sh',  emit: script

    script:
    """
    set -euo pipefail

    echo "[INFO] `date "+%Y-%m-%d %H:%M:%S"` Start GATK HaplotypeCaller (${interval}) ..."
    /usr/bin/time -f "[gatk-HaplotypeCaller] real %e\\tuser %U\\tsys %S\\tcpu_percentage %P\\tmaxrss %M" \\
        ${params.gatk} HaplotypeCaller \\
            -R ${ref} \\
            -I ${bam} \\
            -L ${interval} \\
            -O ${sample_id}.${interval}.g.vcf.gz \\
            -ERC GVCF \\
            --native-pair-hmm-threads ${task.cpus}
    echo "[INFO] `date "+%Y-%m-%d %H:%M:%S"` End GATK HaplotypeCaller (${interval})."
    """
}

process MERGE_GVCFS {
    tag "${sample_id}"
    cpus params.threads
    publishDir "${params.outdir}/gatk/haplotypecaller", mode: 'copy'
    publishDir { "${params.outdir}/logs/${task.process}" }, mode: 'copy', overwrite: true, pattern: '.command.*'

    input:
    tuple val(sample_id), path(gvcfs), path(tbis)

    output:
    tuple val(sample_id), path("${sample_id}.g.vcf.gz"), path("${sample_id}.g.vcf.gz.tbi"), emit: gvcf
    path '.command.log', emit: log
    path '.command.sh',  emit: script

    script:
    """
    set -euo pipefail

    input_args=""
    for g in ${gvcfs}; do
        input_args="\${input_args} -I \${g}"
    done

    echo "[INFO] `date "+%Y-%m-%d %H:%M:%S"` Start GATK GatherVcfs ..."
    /usr/bin/time -f "[gatk-GatherVcfs] real %e\\tuser %U\\tsys %S\\tcpu_percentage %P\\tmaxrss %M" \\
        ${params.gatk} GatherVcfs \\
            \${input_args} \\
            -O ${sample_id}.g.vcf.gz
    echo "[INFO] `date "+%Y-%m-%d %H:%M:%S"` End GATK GatherVcfs."

    echo "[INFO] `date "+%Y-%m-%d %H:%M:%S"` Start GATK IndexFeatureFile ..."
    /usr/bin/time -f "[gatk-IndexFeatureFile] real %e\\tuser %U\\tsys %S\\tcpu_percentage %P\\tmaxrss %M" \\
        ${params.gatk} IndexFeatureFile -I ${sample_id}.g.vcf.gz
    echo "[INFO] `date "+%Y-%m-%d %H:%M:%S"` End GATK IndexFeatureFile."
    """
}
