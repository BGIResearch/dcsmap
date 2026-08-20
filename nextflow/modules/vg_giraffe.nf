process VG_GIRAFFE {

    tag "$sample_name"
    cpus params.threads
    publishDir { "${params.outdir}/logs/${task.process}" }, mode: params.publish_mode, overwrite: true, pattern: '.command.*'

    input:
    val tools_root
    path fq1
    path fq2
    val sample_name
    val platform
    path gbz
    path dist
    path min
    path zipcodes
    val ref_contigs
    val threads
    val ref_path_prefix

    output:
    path "giraffe/${sample_name}.project_to_linear.bam",      emit: surject_bam
    path "giraffe/${sample_name}.project_to_linear.bam.bai",  emit: surject_bai
    path '.command.log', emit: log
    path '.command.sh',  emit: script

    script:
    """
    set -euo pipefail
    export PATH=${tools_root}/libexec:\$PATH

    input_fq1=${fq1}
    input_fq2=${fq2}
    in_sample_name=${sample_name}
    in_gbz_file=${gbz}
    in_dist_file=${dist}
    in_min_file=${min}
    in_zipcodes_file=${zipcodes}
    in_ref_contigs=${ref_contigs}

    giraffe_outdir=giraffe
    mkdir -p \${giraffe_outdir}

    printf -v read_group 'ID:%s\\tLB:%s\\tSM:%s\\tPL:${platform}' "\${in_sample_name}" "\${in_sample_name}" "\${in_sample_name}"
    threads=${threads}
    in_preset="default"

    output_surject_bam=\${giraffe_outdir}/\${in_sample_name}.project_to_linear.bam
    echo "[INFO] `date "+%Y-%m-%d %H:%M:%S"` Start vg giraffe mapping ..."
    /usr/bin/time -f "[vg-giraffe] real %e\\tuser %U\\tsys %S\\tcpu_percentage %P\\tmaxrss %M" \\
    vg giraffe \\
        --progress \\
        --read-group "\${read_group}" \\
        --sample "\${in_sample_name}" \\
        -b \${in_preset} \\
        --output-format bam \\
        --left-align \\
        -f \${input_fq1} -f \${input_fq2} \\
        -Z \${in_gbz_file} \\
        -d \${in_dist_file} \\
        -m \${in_min_file} \\
        -z \${in_zipcodes_file} \\
        --ref-paths \${in_ref_contigs} \\
        -t \${threads} |\\
        samtools addreplacerg -r "\${read_group}" --input-fmt BAM - |\\
        sed "s/${ref_path_prefix}//g" |\\
        mawk '/^@/ {print;next} {print \$0"\\tMX:Z:Giraffe"}' |\\
        samtools sort --threads 4 -m1g --write-index -OBAM -o \${output_surject_bam}##idx##\${output_surject_bam}.bai
    echo "[INFO] `date "+%Y-%m-%d %H:%M:%S"` End vg giraffe mapping."
    """
}
