process MERGE_BAMS {

    publishDir params.outdir, mode: params.publish_mode, overwrite: true,
        pattern: '*.merged.sort.bam*', enabled: params.publish_merged_bam ? true : false
    publishDir { "${params.outdir}/logs/${task.process}" }, mode: params.publish_mode, overwrite: true, pattern: '.command.*'
    tag "$sample_name"
    cpus params.threads

    input:
    val tools_root
    val sample_name
    path bam1
    path bam2
    val threads

    output:
    path "${sample_name}.merged.sort.bam",      emit: merged_bam
    path "${sample_name}.merged.sort.bam.bai",  emit: merged_bam_bai
    path '.command.log', emit: log
    path '.command.sh',  emit: script

    script:
    """
    set -euo pipefail
    export PATH=${tools_root}/libexec:\$PATH

    echo "[INFO] `date "+%Y-%m-%d %H:%M:%S"` Start merge bams ..."
    /usr/bin/time -f "[samtools-merge] real %e\\tuser %U\\tsys %S\\tcpu_percentage %P\\tmaxrss %M" \\
    samtools merge \\
        -f -p -c -@ ${threads} --write-index -OBAM \\
        -o ${sample_name}.merged.sort.bam##idx##${sample_name}.merged.sort.bam.bai \\
        ${bam1} ${bam2}
    echo "[INFO] `date "+%Y-%m-%d %H:%M:%S"` merge bams finished."
    """
}
