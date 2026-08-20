process RTG_VCF_EVAL {

    publishDir { "${params.outdir}/${out_subdir}" }, mode: params.publish_mode, overwrite: true, pattern: 'rtg_benchmark_*'
    publishDir { "${params.outdir}/logs/${task.process}/${caller}" }, mode: params.publish_mode, overwrite: true, pattern: '.command.*'
    tag "$sample"

    input:
    tuple path(vcf), val(tag), val(rtg_outname), val(caller), val(out_subdir)
    val sample
    val giab_version
    val giab_stratify_version
    val rtg_sif
    val bind_path

    output:
    path "${rtg_outname}", emit: result
    path '.command.log', emit: log
    path '.command.sh',  emit: script

    script:
    """
    set -euo pipefail
    echo "[INFO] `date "+%Y-%m-%d %H:%M:%S"` Start RTG vcfeval benchmark ..."
    singularity exec -B ${bind_path} -B \$PWD --workdir \$PWD ${rtg_sif} \\
        rtg_benchmark.sh \\
            -b ${vcf} \\
            -e "${tag}" \\
            -s ${sample} \\
            -o "${rtg_outname}" \\
            --giab_version ${giab_version} \\
            --giab_stratify_version ${giab_stratify_version}
    echo "[INFO] `date "+%Y-%m-%d %H:%M:%S"` End RTG vcfeval benchmark."
    """
}
