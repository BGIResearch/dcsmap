process REALIGN {

    tag "$sample_name"
    cpus params.threads
    publishDir { "${params.outdir}/logs/${task.process}" }, mode: params.publish_mode, overwrite: true, pattern: '.command.*'
    publishDir params.outdir, mode: params.publish_mode, overwrite: true,
        pattern: "*/*.abra2.bam*", enabled: params.publish_abra2_bam ?: false

    input:
    val tools_root
    val java_home
    path bam
    path bai
    val sample_name
    val ref_fasta
    val ref_fasta_fai
    val ref_fasta_dict
    val threads
    val out_subdir

    output:
    path "${out_subdir}/${sample_name}.abra2.bam",     emit: abra2_bam
    path "${out_subdir}/${sample_name}.abra2.bam.bai", emit: abra2_bai
    path '.command.log', emit: log
    path '.command.sh',  emit: script

    script:
    def gkl_flag = (params.use_gkl == false) ? '' : '--gkl'
    """
    set -euo pipefail
    export PATH=${tools_root}/libexec:\$PATH

    java_home=${java_home}
    if [ -n "\${java_home}" ]; then
        export JAVA_HOME="\${java_home}"
        export PATH="\${java_home}/bin:\$PATH"
    fi
    gatk_jar=${tools_root}/jar/GenomeAnalysisTK.jar
    abra2_jar=${tools_root}/jar/abra2-2.24-jar-with-dependencies.jar

    ln -sf ${ref_fasta} ref.fasta
    ln -sf ${ref_fasta_fai} ref.fasta.fai
    ln -sf ${ref_fasta_dict} ref.fasta.dict
    ln -sf ${ref_fasta_dict} ref.dict

    input_bam=${bam}
    in_sample_name=${sample_name}
    in_ref=ref.fasta
    threads=${threads}

    realign_outdir=${out_subdir}
    mkdir -p \${realign_outdir}

    echo "[INFO] `date "+%Y-%m-%d %H:%M:%S"` Start prepare realign targets ..."
    realign_targets_bed=\${realign_outdir}/\${in_sample_name}.realign_targets.bed
    /usr/bin/time -f "[gatk-RealignerTargetCreator] real %e\\tuser %U\\tsys %S\\tcpu_percentage %P\\tmaxrss %M" \\
    java -Xms24g -Xmx64g -XX:+UseG1GC -jar \${gatk_jar} \\
            -T RealignerTargetCreator \\
            -drf DuplicateRead \\
            --disable_bam_indexing \\
            -nt \${threads} \\
            -R \${in_ref} \\
            -I \${input_bam} \\
            --out \${realign_outdir}/\${in_sample_name}.forIndelRealigner.intervals
    awk -F '[:-]' 'BEGIN { OFS = "\\t" } { if( \$3 == "") { print \$1, \$2-1, \$2 } else { print \$1, \$2-1, \$3}}' \${realign_outdir}/\${in_sample_name}.forIndelRealigner.intervals > \${realign_targets_bed}
    echo "[INFO] `date "+%Y-%m-%d %H:%M:%S"` End prepare realign targets."

    realigned_bam=\${realign_outdir}/\${in_sample_name}.abra2.bam
    echo "[INFO] `date "+%Y-%m-%d %H:%M:%S"` Start abra2 indel realignment ..."
    /usr/bin/time -f "[abra2] real %e\\tuser %U\\tsys %S\\tcpu_percentage %P\\tmaxrss %M" \\
    java -Xms24g -Xmx100g -XX:+UseG1GC -jar \${abra2_jar} \\
            --targets \${realign_targets_bed} \\
            --in \${input_bam} \\
            --out \${realigned_bam} \\
            --ref \${in_ref} \\
            --index \\
            --threads \${threads} \\
            ${gkl_flag}

    # abra2 index naming differs by build: the original abra2 writes
    # <sample>.abra2.bai, while the optimized abra2 writes <sample>.abra2.bam.bai.
    # Normalize to <sample>.abra2.bam.bai so downstream steps see a consistent name.
    if [ ! -e "\${realigned_bam}.bai" ]; then
        legacy_bai="\${realigned_bam%.bam}.bai"
        if [ -e "\${legacy_bai}" ]; then
            mv "\${legacy_bai}" "\${realigned_bam}.bai"
        else
            echo "[ERROR] abra2 index not found: \${realigned_bam}.bai or \${legacy_bai}" >&2
            exit 1
        fi
    fi
    echo "[INFO] `date "+%Y-%m-%d %H:%M:%S"` End abra2 indel realignment."
    """
}
