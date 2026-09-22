#include "Scripts.h"

#include <string>

namespace {

std::string joinPath(const std::string& a, const std::string& b) {
    if (a.empty()) return b;
    if (a.back() == '/') return a + b;
    return a + "/" + b;
}

void replaceAll(std::string& s, const std::string& ph, const std::string& val) {
    size_t pos = 0;
    while ((pos = s.find(ph, pos)) != std::string::npos) {
        s.replace(pos, ph.size(), val);
        pos += val.size();
    }
}

std::string jobDir(const Config& c, const std::string& name) {
    return joinPath(c.work_root, "task-" + name);
}

}  // namespace

// T1: linear reference align (BWA-MEM2) + extract
std::string buildLinearAlignExtract(const Config& c, std::vector<std::string>& outputs) {
    const std::string jd = jobDir(c, "linear_align_extract");
    outputs = {
        joinPath(jd, c.sample_name + ".linear.extract_1.fastq.gz"),
        joinPath(jd, c.sample_name + ".linear.extract_2.fastq.gz"),
        joinPath(jd, c.sample_name + ".linear.not-extract.bam"),
        joinPath(jd, c.sample_name + ".linear.not-extract.bam.bai"),
    };
    std::string s = R"BASH(set -euo pipefail
export PATH="{{TOOLS_ROOT}}/libexec:$PATH"

in_ref="{{REF_FASTA}}"
sample="{{SAMPLE}}"
printf -v read_group '@RG\tID:%s\tLB:%s\tSM:%s\tPL:{{PLATFORM}}' "$sample" "$sample" "$sample"
extract_prefix="${sample}.linear.extract"
out_bam="${sample}.linear.not-extract.bam"

echo "[INFO] $(date '+%Y-%m-%d %H:%M:%S') Start linear align (bwa-mem2) mapping and extract ..."
timed "linear-align/bwa-mem2" \
    bwa-mem2 mem -S -r 3 -y 0 -Y -t {{THREADS}} -R "$read_group" "$in_ref" "{{FQ1}}" "{{FQ2}}" | \
    "${AWK}" '/^@/ {print;next} {print $0"\tMX:Z:BWA-MEM2"}' | \
    extract-bam -b 4000000 -t 12 -p "$extract_prefix" -M "{{EXTRACT_MODEL}}" | \
    samtools sort -l 1 -OBAM --threads {{SAMTOOLS_SORT_THREADS}} -m1g -o "$out_bam"##idx##"$out_bam".bai --write-index
echo "[INFO] $(date '+%Y-%m-%d %H:%M:%S') End linear align (bwa-mem2) mapping and extract."
)BASH";
    replaceAll(s, "{{TOOLS_ROOT}}", c.tools_root);
    replaceAll(s, "{{REF_FASTA}}", c.ref_fasta);
    replaceAll(s, "{{SAMPLE}}", c.sample_name);
    replaceAll(s, "{{PLATFORM}}", c.platform);
    replaceAll(s, "{{THREADS}}", std::to_string(c.threads));
    replaceAll(s, "{{SAMTOOLS_SORT_THREADS}}", std::to_string(c.samtools_sort_threads));
    replaceAll(s, "{{FQ1}}", c.fq1);
    replaceAll(s, "{{FQ2}}", c.fq2);
    replaceAll(s, "{{EXTRACT_MODEL}}", c.extract_model);
    return s;
}

// T2: VG haplotype sampling
std::string buildVgHaplotype(const Config& c, std::vector<std::string>& outputs) {
    const std::string jd = jobDir(c, "vg_haplotype");
    const std::string hs = joinPath(joinPath(jd, "haplotype_sampling"), "haplotype_sampling");
    outputs = {
        joinPath(hs, "haplotype_sampling.gbz"),
        joinPath(hs, "haplotype_sampling.dist"),
        joinPath(hs, "haplotype_sampling.min"),
        joinPath(hs, "haplotype_sampling.zipcodes"),
    };
    std::string s = R"BASH(set -euo pipefail
export PATH="{{TOOLS_ROOT}}/libexec:$PATH"

threads={{VG_HAPLOTYPE_THREADS}}
input_fq1="{{FQ1}}"
input_fq2="{{FQ2}}"
in_gbz_file="{{GBZ}}"
in_hap_index="{{HAPL}}"
outdir="haplotype_sampling"
kmc_output_prefix=kmc
kmer_len=29
kmc_max_ram=64
kmc_outdir="${outdir}/kmc"
mkdir -p "${kmc_outdir}/kmc_tmp"
echo "${input_fq1}" > "${kmc_outdir}/kmc_input_list.txt"
echo "${input_fq2}" >> "${kmc_outdir}/kmc_input_list.txt"

echo "[INFO] $(date '+%Y-%m-%d %H:%M:%S') Start kmer counting ..."
timed "kmc" \
    kmc -k${kmer_len} -t${threads} -m${kmc_max_ram} -okff @"${kmc_outdir}/kmc_input_list.txt" "${kmc_outdir}/${kmc_output_prefix}" "${kmc_outdir}/kmc_tmp"
rm -rf "${kmc_outdir}/kmc_tmp" "${kmc_outdir}/kmc_input_list.txt"
echo "[INFO] $(date '+%Y-%m-%d %H:%M:%S') End kmer counting."

in_kmer_info="${kmc_outdir}/${kmc_output_prefix}.kff"
haplotype_number=32
present_discount=0.9
het_adjust=0.05
absent_score=0.8
haplotype_sampling_outdir="${outdir}/haplotype_sampling"
mkdir -p "${haplotype_sampling_outdir}"
unset OMP_NUM_THREADS

echo "[INFO] $(date '+%Y-%m-%d %H:%M:%S') Start haplotype sampling ..."
timed "vg-haplotypes" \
    vg haplotypes -v 2 -t ${threads} \
        --num-haplotypes ${haplotype_number} \
        --present-discount ${present_discount} \
        --het-adjustment ${het_adjust} \
        --absent-score ${absent_score} \
        --include-reference \
        --set-reference {{SET_REFERENCE}} \
        --diploid-sampling \
        -i ${in_hap_index} \
        -k ${in_kmer_info} \
        -g ${haplotype_sampling_outdir}/haplotype_sampling.gbz ${in_gbz_file}
echo "[INFO] $(date '+%Y-%m-%d %H:%M:%S') End haplotype sampling."

echo "[INFO] $(date '+%Y-%m-%d %H:%M:%S') Start creating indexes ..."
timed "vg-index" \
    vg index -t ${threads} -j ${haplotype_sampling_outdir}/haplotype_sampling.dist ${haplotype_sampling_outdir}/haplotype_sampling.gbz
in_minimizer_k=29
in_minimizer_w=11
timed "vg-minimizer" \
    vg minimizer -p -t ${threads} -k ${in_minimizer_k} -w ${in_minimizer_w} --weighted --save-memory \
        -o ${haplotype_sampling_outdir}/haplotype_sampling.min -z ${haplotype_sampling_outdir}/haplotype_sampling.zipcodes \
        -d ${haplotype_sampling_outdir}/haplotype_sampling.dist ${haplotype_sampling_outdir}/haplotype_sampling.gbz
echo "[INFO] $(date '+%Y-%m-%d %H:%M:%S') End creating indexes."
echo "[INFO] $(date '+%Y-%m-%d %H:%M:%S') Haplotype sampling workflow is done."
)BASH";
    replaceAll(s, "{{TOOLS_ROOT}}", c.tools_root);
    replaceAll(s, "{{VG_HAPLOTYPE_THREADS}}", std::to_string(c.vg_haplotype_threads));
    replaceAll(s, "{{FQ1}}", c.fq1);
    replaceAll(s, "{{FQ2}}", c.fq2);
    replaceAll(s, "{{GBZ}}", c.gbz);
    replaceAll(s, "{{HAPL}}", c.hapl);
    replaceAll(s, "{{SET_REFERENCE}}", c.set_reference);
    return s;
}

// T3: VG giraffe
std::string buildVgGiraffe(const Config& c,
                           const std::string& extract_fq1,
                           const std::string& extract_fq2,
                           const std::string& hap_gbz,
                           const std::string& hap_dist,
                           const std::string& hap_min,
                           const std::string& hap_zipcodes,
                           std::vector<std::string>& outputs) {
    const std::string jd = jobDir(c, "vg_giraffe");
    outputs = {
        joinPath(joinPath(jd, "giraffe"), c.sample_name + ".project_to_linear.bam"),
        joinPath(joinPath(jd, "giraffe"), c.sample_name + ".project_to_linear.bam.bai"),
    };
    std::string s = R"BASH(set -euo pipefail
export PATH="{{TOOLS_ROOT}}/libexec:$PATH"

input_fq1="{{EXTRACT_FQ1}}"
input_fq2="{{EXTRACT_FQ2}}"
in_sample_name="{{SAMPLE}}"
in_gbz_file="{{HAP_GBZ}}"
in_dist_file="{{HAP_DIST}}"
in_min_file="{{HAP_MIN}}"
in_zipcodes_file="{{HAP_ZIPCODES}}"
in_ref_contigs="{{GRAPH_REF_CONTIGS}}"
giraffe_outdir=giraffe
mkdir -p "${giraffe_outdir}"
printf -v read_group 'ID:%s\tLB:%s\tSM:%s\tPL:{{PLATFORM}}' "${in_sample_name}" "${in_sample_name}" "${in_sample_name}"
threads={{THREADS}}
in_preset="default"
output_surject_bam="${giraffe_outdir}/${in_sample_name}.project_to_linear.bam"

echo "[INFO] $(date '+%Y-%m-%d %H:%M:%S') Start vg giraffe mapping ..."
timed "vg-giraffe" \
    vg giraffe \
        --progress \
        --read-group "${read_group}" \
        --sample "${in_sample_name}" \
        -b ${in_preset} \
        --output-format bam \
        --left-align \
        -f ${input_fq1} -f ${input_fq2} \
        -Z ${in_gbz_file} \
        -d ${in_dist_file} \
        -m ${in_min_file} \
        -z ${in_zipcodes_file} \
        --ref-paths ${in_ref_contigs} \
        -t ${threads} | \
    samtools addreplacerg -r "${read_group}" --input-fmt BAM - | \
    sed "s/{{REF_PATH_PREFIX}}//g" | \
    "${AWK}" '/^@/ {print;next} {print $0"\tMX:Z:Giraffe"}' | \
    samtools sort --threads {{THREADS}} -m1g --write-index -OBAM -o ${output_surject_bam}##idx##${output_surject_bam}.bai
echo "[INFO] $(date '+%Y-%m-%d %H:%M:%S') End vg giraffe mapping."
)BASH";
    replaceAll(s, "{{TOOLS_ROOT}}", c.tools_root);
    replaceAll(s, "{{EXTRACT_FQ1}}", extract_fq1);
    replaceAll(s, "{{EXTRACT_FQ2}}", extract_fq2);
    replaceAll(s, "{{SAMPLE}}", c.sample_name);
    replaceAll(s, "{{HAP_GBZ}}", hap_gbz);
    replaceAll(s, "{{HAP_DIST}}", hap_dist);
    replaceAll(s, "{{HAP_MIN}}", hap_min);
    replaceAll(s, "{{HAP_ZIPCODES}}", hap_zipcodes);
    replaceAll(s, "{{GRAPH_REF_CONTIGS}}", c.graph_ref_contigs);
    replaceAll(s, "{{PLATFORM}}", c.platform);
    replaceAll(s, "{{THREADS}}", std::to_string(c.threads));
    replaceAll(s, "{{REF_PATH_PREFIX}}", c.ref_path_prefix);
    return s;
}

// MERGE_BAMS
std::string buildMergeBams(const Config& c,
                           const std::string& bam1,
                           const std::string& bam2,
                           const std::string& final_out_bam,
                           std::vector<std::string>& outputs) {
    const std::string jd = jobDir(c, "merge_bams");
    const std::string out_bam = !final_out_bam.empty()
        ? final_out_bam
        : joinPath(jd, c.sample_name + ".merged.sort.bam");
    outputs = { out_bam, out_bam + ".bai" };
    std::string s = R"BASH(set -euo pipefail
export PATH="{{TOOLS_ROOT}}/libexec:$PATH"

out_bam="{{OUT_BAM}}"
out_bai="${out_bam}.bai"

echo "[INFO] $(date '+%Y-%m-%d %H:%M:%S') Start merge bams ..."
timed "samtools-merge" \
    samtools merge \
        -f -p -c -@ {{THREADS}} --write-index -OBAM \
        -o "${out_bam}"##idx##"${out_bai}" \
        "{{BAM1}}" "{{BAM2}}"
echo "[INFO] $(date '+%Y-%m-%d %H:%M:%S') merge bams finished."
)BASH";
    replaceAll(s, "{{TOOLS_ROOT}}", c.tools_root);
    replaceAll(s, "{{OUT_BAM}}", out_bam);
    replaceAll(s, "{{THREADS}}", std::to_string(c.threads));
    replaceAll(s, "{{BAM1}}", bam1);
    replaceAll(s, "{{BAM2}}", bam2);
    return s;
}

// REALIGN (GATK RealignerTargetCreator + abra2, always --gkl)
std::string buildRealign(const Config& c,
                         const std::string& in_bam,
                         const std::string& in_bai,
                         const std::string& out_subdir,
                         const std::string& final_out_bam,
                         std::vector<std::string>& outputs) {
    (void)in_bai;
    const std::string jd = jobDir(c, "realign");

    std::string realigned_bam;
    std::string mkdir_block;
    std::string targets_bed;
    std::string intervals_file;

    // Targets/intervals always live in <jobdir>/<out_subdir>; only the bam may be
    // redirected to the caller-provided output path.
    const std::string work_subdir = joinPath(jd, out_subdir);
    targets_bed = joinPath(work_subdir, c.sample_name + ".realign_targets.bed");
    intervals_file = joinPath(work_subdir, c.sample_name + ".forIndelRealigner.intervals");

    if (!final_out_bam.empty()) {
        // final task: write abra2 bam directly to outdir path
        realigned_bam = final_out_bam;
        std::string parent = realigned_bam;
        size_t slash = parent.find_last_of('/');
        parent = (slash == std::string::npos) ? "." : parent.substr(0, slash);
        mkdir_block = "realign_outdir=\"" + parent + "\"\nmkdir -p \"${realign_outdir}\" \"" +
                      work_subdir + "\"";
    } else {
        // intermediate: <jobdir>/<out_subdir>/<sample>.abra2.bam
        realigned_bam = joinPath(work_subdir, c.sample_name + ".abra2.bam");
        mkdir_block = "realign_outdir=\"" + out_subdir + "\"\nmkdir -p \"${realign_outdir}\"";
    }
    outputs = { realigned_bam, realigned_bam + ".bai" };

    std::string s = R"BASH(set -euo pipefail
export PATH="{{TOOLS_ROOT}}/libexec:$PATH"

java_home="{{JAVA_HOME}}"
if [ -n "${java_home}" ]; then
    export JAVA_HOME="${java_home}"
    export PATH="${java_home}/bin:$PATH"
fi
gatk_jar="{{TOOLS_ROOT}}/jar/GenomeAnalysisTK.jar"
abra2_jar="{{TOOLS_ROOT}}/jar/abra2-2.24-jar-with-dependencies.jar"

ln -sf "{{REF_FASTA}}" ref.fasta
ln -sf "{{REF_FASTA_FAI}}" ref.fasta.fai
ln -sf "{{REF_FASTA_DICT}}" ref.fasta.dict
ln -sf "{{REF_FASTA_DICT}}" ref.dict

input_bam="{{IN_BAM}}"
in_sample_name="{{SAMPLE}}"
in_ref=ref.fasta
threads={{THREADS}}
{{MKDIR_BLOCK}}
realign_targets_bed="{{TARGETS_BED}}"

echo "[INFO] $(date '+%Y-%m-%d %H:%M:%S') Start prepare realign targets ..."
timed "gatk-RealignerTargetCreator" \
    java -Xms24g -Xmx64g -XX:+UseG1GC -jar ${gatk_jar} \
        -T RealignerTargetCreator \
        -drf DuplicateRead \
        --disable_bam_indexing \
        -nt ${threads} \
        -R ${in_ref} \
        -I ${input_bam} \
        --out "{{INTERVALS_FILE}}"
"${AWK}" -F '[:-]' 'BEGIN { OFS = "\t" } { if( $3 == "") { print $1, $2-1, $2 } else { print $1, $2-1, $3}}' "{{INTERVALS_FILE}}" > ${realign_targets_bed}
echo "[INFO] $(date '+%Y-%m-%d %H:%M:%S') End prepare realign targets."

realigned_bam="{{REALIGNED_BAM}}"
echo "[INFO] $(date '+%Y-%m-%d %H:%M:%S') Start abra2 indel realignment ..."
timed "abra2" \
    java -Xms24g -Xmx80g -XX:+UseG1GC -jar ${abra2_jar} \
        --targets ${realign_targets_bed} \
        --in ${input_bam} \
        --out ${realigned_bam} \
        --ref ${in_ref} \
        --index \
        --threads ${threads} \
        --gkl
if [ ! -e "${realigned_bam}.bai" ]; then
    legacy_bai="${realigned_bam%.bam}.bai"
    if [ -e "${legacy_bai}" ]; then
        mv "${legacy_bai}" "${realigned_bam}.bai"
    else
        echo "[ERROR] abra2 index not found: ${realigned_bam}.bai or ${legacy_bai}" >&2
        exit 1
    fi
fi
echo "[INFO] $(date '+%Y-%m-%d %H:%M:%S') End abra2 indel realignment."
)BASH";

    replaceAll(s, "{{TOOLS_ROOT}}", c.tools_root);
    replaceAll(s, "{{JAVA_HOME}}", c.java_home);
    replaceAll(s, "{{REF_FASTA}}", c.ref_fasta);
    replaceAll(s, "{{REF_FASTA_FAI}}", c.ref_fasta_fai);
    replaceAll(s, "{{REF_FASTA_DICT}}", c.ref_fasta_dict);
    replaceAll(s, "{{IN_BAM}}", in_bam);
    replaceAll(s, "{{SAMPLE}}", c.sample_name);
    replaceAll(s, "{{THREADS}}", std::to_string(c.threads));
    replaceAll(s, "{{MKDIR_BLOCK}}", mkdir_block);
    replaceAll(s, "{{TARGETS_BED}}", targets_bed);
    replaceAll(s, "{{INTERVALS_FILE}}", intervals_file);
    replaceAll(s, "{{REALIGNED_BAM}}", realigned_bam);
    return s;
}
