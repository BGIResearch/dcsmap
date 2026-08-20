process VG_HAPLOTYPE_SAMPLING {

    tag "haplotype_sampling"
    cpus params.threads
    publishDir { "${params.outdir}/logs/${task.process}" }, mode: params.publish_mode, overwrite: true, pattern: '.command.*'

    input:
    val tools_root
    val fq1
    val fq2
    val gbz
    val hapl
    val threads
    val set_reference

    output:
    path "haplotype_sampling/haplotype_sampling/haplotype_sampling.gbz",      emit: haplotype_gbz
    path "haplotype_sampling/haplotype_sampling/haplotype_sampling.dist",     emit: haplotype_dist
    path "haplotype_sampling/haplotype_sampling/haplotype_sampling.min",      emit: haplotype_min
    path "haplotype_sampling/haplotype_sampling/haplotype_sampling.zipcodes", emit: haplotype_zipcodes
    path '.command.log', emit: log
    path '.command.sh',  emit: script

    script:
    """
    set -euo pipefail
    export PATH=${tools_root}/libexec:\$PATH

    threads=${threads}
    input_fq1=${fq1}
    input_fq2=${fq2}
    in_gbz_file=${gbz}
    in_hap_index=${hapl}
    outdir="haplotype_sampling"

    kmc_output_prefix=kmc
    kmer_len=29
    kmc_max_ram=64
    kmc_outdir=\${outdir}/kmc
    mkdir -p \${kmc_outdir}/kmc_tmp

    echo \${input_fq1} > \${kmc_outdir}/kmc_input_list.txt
    echo \${input_fq2} >> \${kmc_outdir}/kmc_input_list.txt

    echo "[INFO] `date "+%Y-%m-%d %H:%M:%S"` Start kmer counting ..."
    /usr/bin/time -f "[kmc] real %e\\tuser %U\\tsys %S\\tcpu_percentage %P\\tmaxrss %M" \\
    kmc -k\${kmer_len} -t\${threads} -m\${kmc_max_ram} -okff @\${kmc_outdir}/kmc_input_list.txt \${kmc_outdir}/\${kmc_output_prefix} \${kmc_outdir}/kmc_tmp
    rm -rf \${kmc_outdir}/kmc_tmp \${kmc_outdir}/kmc_input_list.txt
    echo "[INFO] `date "+%Y-%m-%d %H:%M:%S"` End kmer counting."

    in_kmer_info=\${kmc_outdir}/\${kmc_output_prefix}.kff
    haplotype_number=32
    present_discount=0.9
    het_adjust=0.05
    absent_score=0.8
    haplotype_sampling_outdir=\${outdir}/haplotype_sampling
    mkdir -p \${haplotype_sampling_outdir}

    unset OMP_NUM_THREADS
    echo "[INFO] `date "+%Y-%m-%d %H:%M:%S"` Start haplotype sampling ..."
    /usr/bin/time -f "[vg-haplotypes] real %e\\tuser %U\\tsys %S\\tcpu_percentage %P\\tmaxrss %M" \\
    vg haplotypes -v 2 -t \${threads} \\
            --num-haplotypes \${haplotype_number} \\
            --present-discount \${present_discount} \\
            --het-adjustment \${het_adjust} \\
            --absent-score \${absent_score} \\
            --include-reference \\
            --set-reference ${set_reference} \\
            --diploid-sampling \\
            -i \${in_hap_index} \\
            -k \${in_kmer_info} \\
            -g \${haplotype_sampling_outdir}/haplotype_sampling.gbz \${in_gbz_file}
    echo "[INFO] `date "+%Y-%m-%d %H:%M:%S"` End haplotype sampling."

    echo "[INFO] `date "+%Y-%m-%d %H:%M:%S"` Start creating indexes ..."
    /usr/bin/time -f "[vg-index] real %e\\tuser %U\\tsys %S\\tcpu_percentage %P\\tmaxrss %M" \\
    vg index -t \${threads} -j \${haplotype_sampling_outdir}/haplotype_sampling.dist \${haplotype_sampling_outdir}/haplotype_sampling.gbz
    in_minimizer_k=29
    in_minimizer_w=11
    /usr/bin/time -f "[vg-minimizer] real %e\\tuser %U\\tsys %S\\tcpu_percentage %P\\tmaxrss %M" \\
    vg minimizer -p -t \${threads} -k \${in_minimizer_k} -w \${in_minimizer_w} --weighted --save-memory \\
        -o \${haplotype_sampling_outdir}/haplotype_sampling.min -z \${haplotype_sampling_outdir}/haplotype_sampling.zipcodes \\
        -d \${haplotype_sampling_outdir}/haplotype_sampling.dist \${haplotype_sampling_outdir}/haplotype_sampling.gbz
    echo "[INFO] `date "+%Y-%m-%d %H:%M:%S"` End creating indexes."
    echo "[INFO] `date "+%Y-%m-%d %H:%M:%S"` Haplotype sampling workflow is done."
    """
}
