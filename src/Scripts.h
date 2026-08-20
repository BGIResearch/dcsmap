#pragma once

#include <string>
#include <vector>

// Resolved configuration shared by all task builders. All paths are absolute.
// `work_root` is the mktemp work dir (work.XXXXXX). Each task writes into
// `<work_root>/task-<task_name>` (the Task base class handles that).
struct Config {
    // inputs
    std::string fq1;
    std::string fq2;
    std::string ref_fasta;          // required (with .fai + .dict)
    std::string ref_fasta_fai;       // = ref_fasta + ".fai"
    std::string ref_fasta_dict;      // = ref_fasta + ".dict"
    std::string gbz;
    std::string hapl;
    std::string graph_ref_contigs;   // vg ref-paths file (was --ref-contigs)
    std::string extract_model;

    // labels
    std::string sample_name = "SAMPLE";
    std::string platform = "Illumina";

    // resources
    int threads = 32;
    int samtools_sort_threads = 4;   // threads for samtools sort (4 when --parallel, else = threads)
    int vg_haplotype_threads = 16;   // threads for vg_haplotype (threads/2 when --parallel, else = threads)
    std::string tools_root;          // absolute (libexec + jar root)
    std::string java_home;           // absolute JDK root (verified Java 8)

    // derived from --graph-ref-contigs
    std::string set_reference;       // <ref>, used for vg haplotypes --set-reference
    std::string ref_path_prefix;     // <ref>#0#, stripped from vg giraffe surject BAM contigs

    // runtime
    std::string work_root;           // work.XXXXXX (absolute)
    std::string out_bam;             // final bam path (absolute); work dir is its dirname
};

// Each builder returns the bash script body and fills `outputs` with the
// absolute paths the task is expected to produce (used for verification).

// T1: BWA-MEM2 map + extract. Produces extract fastq pair + not-extract bam/bai.
std::string buildBwaMem2Extract(const Config& c, std::vector<std::string>& outputs);

// T2: VG haplotype sampling (kmc + vg haplotypes + vg index + vg minimizer).
std::string buildVgHaplotype(const Config& c, std::vector<std::string>& outputs);

// T3: VG giraffe mapping on the haplotype-sampled subgraph.
// `extract_fq1/fq2` are outputs of T1; `hap_gbz/dist/min/zipcodes` are outputs of T2.
std::string buildVgGiraffe(const Config& c,
                           const std::string& extract_fq1,
                           const std::string& extract_fq2,
                           const std::string& hap_gbz,
                           const std::string& hap_dist,
                           const std::string& hap_min,
                           const std::string& hap_zipcodes,
                           std::vector<std::string>& outputs);

// MERGE_BAMS: samtools merge --write-index.
// `final_out_bam` is empty for the intermediate case (bam written in job dir);
// when non-empty the merged bam is written directly to that absolute path.
std::string buildMergeBams(const Config& c,
                           const std::string& bam1,
                           const std::string& bam2,
                           const std::string& final_out_bam,
                           std::vector<std::string>& outputs);

// REALIGN: GATK RealignerTargetCreator + abra2 (always --gkl).
// `out_subdir` is the relative subdir under the job dir where abra2 writes
// (e.g. "giraffe" for dcsmap intermediate, "merged_realign" for dcsmap-m1 final).
// `final_out_bam` is empty for the intermediate case; when non-empty the abra2
// bam is written directly to that absolute path (and out_subdir is ignored for
// placement).
std::string buildRealign(const Config& c,
                         const std::string& in_bam,
                         const std::string& in_bai,
                         const std::string& out_subdir,
                         const std::string& final_out_bam,
                         std::vector<std::string>& outputs);
