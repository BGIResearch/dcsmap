// ============================================================================
// GATK GenotypeGVCFs module
//   GenotypeGVCFs -> IndexFeatureFile
//   Produces the final VCF (sample_id.vcf.gz), used directly by downstream
//   benchmarking (RTG/hap.py) against GIAB truth sets.
// ============================================================================

process GATK_GENOTYPE_GVCFS {
    tag "${sample_id}"
    cpus params.threads
    publishDir "${params.outdir}/gatk/genotype", mode: 'copy'
    publishDir { "${params.outdir}/logs/${task.process}" }, mode: 'copy', overwrite: true, pattern: '.command.*'

    input:
    tuple val(sample_id), path(gvcf), path(tbi)
    path ref
    path ref_idx
    path ref_dict

    output:
    tuple val(sample_id), path("${sample_id}.vcf.gz"), path("${sample_id}.vcf.gz.tbi"), emit: vcf
    path '.command.log', emit: log
    path '.command.sh',  emit: script

    script:
    """
    set -euo pipefail

    echo "[INFO] `date "+%Y-%m-%d %H:%M:%S"` Start GATK GenotypeGVCFs ..."
    /usr/bin/time -f "[gatk-GenotypeGVCFs] real %e\\tuser %U\\tsys %S\\tcpu_percentage %P\\tmaxrss %M" \\
        ${params.gatk} GenotypeGVCFs \\
            -R ${ref} \\
            -V ${gvcf} \\
            -O ${sample_id}.vcf.gz
    echo "[INFO] `date "+%Y-%m-%d %H:%M:%S"` End GATK GenotypeGVCFs."

    echo "[INFO] `date "+%Y-%m-%d %H:%M:%S"` Start GATK IndexFeatureFile (genotyped vcf) ..."
    /usr/bin/time -f "[gatk-IndexFeatureFile] real %e\\tuser %U\\tsys %S\\tcpu_percentage %P\\tmaxrss %M" \\
        ${params.gatk} IndexFeatureFile -I ${sample_id}.vcf.gz
    echo "[INFO] `date "+%Y-%m-%d %H:%M:%S"` End GATK IndexFeatureFile (genotyped vcf)."
    """
}
