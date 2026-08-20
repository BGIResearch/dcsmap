process DEEPVARIANT {

    publishDir params.outdir, mode: params.publish_mode, overwrite: true, pattern: 'deepvariant/*'
    publishDir { "${params.outdir}/logs/${task.process}" }, mode: params.publish_mode, overwrite: true, pattern: '.command.*'
    tag "$sample_name"
    cpus params.threads

    input:
    path bam
    path bai
    val ref_fasta
    val sample_name
    val threads
    val disable_small_model_arg
    val make_examples_extra_arg
    val dv_sif
    val bind_path

    output:
    path "deepvariant/${sample_name}.deepvariant.g.vcf.gz",       emit: gvcf
    path "deepvariant/${sample_name}.deepvariant.g.vcf.gz.tbi",   emit: gvcf_tbi
    path "deepvariant/${sample_name}.deepvariant.vcf.gz",         emit: vcf
    path "deepvariant/${sample_name}.deepvariant.vcf.gz.tbi",     emit: vcf_tbi
    path "deepvariant/${sample_name}.deepvariant.pass.vcf.gz",   emit: pass_vcf
    path "deepvariant/${sample_name}.deepvariant.pass.vcf.gz.tbi", emit: pass_vcf_tbi
    path '.command.log', emit: log
    path '.command.sh',  emit: script

    script:
    """
    set -euo pipefail

    ln -sf ${bam} input.bam
    ln -sf ${bai} input.bam.bai

    threads=${threads}
    in_bam=input.bam
    in_ref=${ref_fasta}
    sample_name=${sample_name}
    outdir="deepvariant"
    mkdir -p \${outdir}

    echo "[INFO] `date "+%Y-%m-%d %H:%M:%S"` Start DeepVariant run_deepvariant ..."
    tm_start=\$(date +%s)
    /usr/bin/time -f "[deepvariant] real %e\\tuser %U\\tsys %S\\tcpu_percentage %P\\tmaxrss %M" \\
    singularity exec -B ${bind_path} -B \$PWD --workdir \$PWD \\
        ${dv_sif} /opt/deepvariant/bin/run_deepvariant \\
        --model_type=WGS \\
        --ref \${in_ref} \\
        --reads \${in_bam} \\
        --intermediate_results_dir \${outdir}/intermediate_results \\
        --output_vcf \${outdir}/\${sample_name}.deepvariant.vcf.gz \\
        --output_gvcf \${outdir}/\${sample_name}.deepvariant.g.vcf.gz \\
        --num_shards=\${threads} ${disable_small_model_arg} \\
        --logging_dir=\${outdir}/logs \\
        ${make_examples_extra_arg}
    tm_end=\$(date +%s)
    elapsed=\$((\${tm_end} - \${tm_start}))
    echo "[INFO] `date "+%Y-%m-%d %H:%M:%S"` End DeepVariant run_deepvariant. Elapsed time: \${elapsed} seconds."

    bcftools view -f "PASS," \${outdir}/\${sample_name}.deepvariant.vcf.gz -Oz -o \${outdir}/\${sample_name}.deepvariant.pass.vcf.gz --write-index=tbi
    """
}
