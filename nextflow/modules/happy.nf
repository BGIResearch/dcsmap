process HAPPY {

    publishDir { "${params.outdir}/${out_subdir}" }, mode: params.publish_mode, overwrite: true, pattern: 'happy_out'
    publishDir { "${params.outdir}/logs/${task.process}/${caller}" }, mode: params.publish_mode, overwrite: true, pattern: '.command.*'
    tag "$sample"

    input:
    tuple path(vcf), val(caller), val(out_subdir)
    val sample
    val happy_ref
    val ref_assembly
    val giab_version
    val happy_home
    val happy_sif
    val bind_path

    output:
    path "happy_out", emit: result
    path '.command.log', emit: log
    path '.command.sh',  emit: script

    script:
    // -r is optional: when happy_ref is empty, run_hap.py.sh defaults the
    // reference FASTA by assembly (GRCh38 / CHM13). -a selects the truth
    // set assembly and -v selects the GIAB benchmark version (v4.2.1 or
    // v5.0q; CHM13 requires v5.0q + HG002).
    def ref_arg = happy_ref ? "-r ${happy_ref}" : ''
    """
    set -euo pipefail
    mkdir -p happy_out
    echo "[INFO] `date "+%Y-%m-%d %H:%M:%S"` Start hap.py evaluation ..."
    bash ${happy_home}/run_hap.py.sh \\
        -c ${vcf} \\
        -s ${sample} \\
        ${ref_arg} \\
        -a ${ref_assembly} \\
        -v ${giab_version} \\
        -o happy_out \\
        -b ${bind_path} \\
        -i ${happy_sif}
    echo "[INFO] `date "+%Y-%m-%d %H:%M:%S"` End hap.py evaluation."
    """
}
