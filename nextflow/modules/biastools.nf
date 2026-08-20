process BIASTOOLS {

    publishDir params.outdir, mode: params.publish_mode, overwrite: true, pattern: 'biastools'
    publishDir { "${params.outdir}/logs/${task.process}" }, mode: params.publish_mode, overwrite: true, pattern: '.command.*'
    tag "$sample"

    input:
    tuple path(bam), path(bai), val(tag)
    val ref_fasta
    val benchmark_vcf
    val sample
    val depth
    val threads
    val biastools_bin

    output:
    path "biastools", emit: result
    path '.command.log', emit: log
    path '.command.sh',  emit: script

    script:
    """
    set -euo pipefail
    mkdir -p biastools
    echo "[INFO] `date "+%Y-%m-%d %H:%M:%S"` Start biastools map-bias analysis ..."
    /usr/bin/time -f "[biastools] real %e\\tuser %U\\tsys %S\\tcpu_percentage %P\\tmaxrss %M" \\
    ${biastools_bin} --analyze --real -t ${threads} \\
        -o biastools \\
        -g ${ref_fasta} \\
        -v ${benchmark_vcf} \\
        -s ${sample} \\
        -r "${tag}" \\
        --bam ${bam} \\
        -d ${depth}
    echo "[INFO] `date "+%Y-%m-%d %H:%M:%S"` End biastools map-bias analysis."
    """
}
