#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

// ============================================================================
// BWA-MEM2 Germline Variant Calling & Benchmarking Pipeline (entry workflow)
//
//   Always on (shared preprocessing):
//     modules/bwa_mem2.nf          : BWA_MEM2
//     modules/sambamba_markdup.nf  : SAMBAMBA_MARKDUP
//
//   GATK path      (params.use_gatk      default: true):
//     modules/gatk_bqsr.nf         : GATK_BQSR
//     modules/haplotypecaller.nf   : GATK_HAPLOTYPECALLER (per-chrom) -> MERGE_GVCFS
//     modules/genotype_gvcfs.nf    : GATK_GENOTYPE_GVCFS
//
//   DeepVariant    (params.use_deepvariant default: false):
//     modules/deepvariant.nf       : DEEPVARIANT
//
//   Benchmarking   (params.use_rtg / params.use_happy / params.use_biastools default: false):
//     modules/rtg.nf               : RTG_VCF_EVAL   (per enabled caller, VCF-based)
//     modules/happy.nf             : HAPPY          (per enabled caller, VCF-based)
//     modules/biastools.nf         : BIASTOOLS      (single run on markdup BAM, BAM-based)
//     VCF benchmarks (RTG/happy) of each enabled caller (GATK and/or DeepVariant)
//     are separated under <outdir>/benchmark/<caller>/{rtg,happy}/.
//     biastools runs once on the markduplicate BAM -> <outdir>/biastools/.
//
//   Logs / scripts (.command.sh, .command.log) are unified under:
//     <outdir>/logs/<process>/                          (single-run tasks)
//     <outdir>/logs/<process>/<caller>/                 (RTG_VCF_EVAL, HAPPY)
//     <outdir>/logs/<process>/<chrom>/                  (GATK_HAPLOTYPECALLER)
//   so that processes invoked multiple times do not overwrite each other.
// ============================================================================

include { BWA_MEM2 } from './modules/bwa_mem2'
include { SAMBAMBA_MARKDUP } from './modules/sambamba_markdup'
include { GATK_BQSR } from './modules/gatk_bqsr'
include { GATK_HAPLOTYPECALLER; MERGE_GVCFS } from './modules/haplotypecaller'
include { GATK_GENOTYPE_GVCFS } from './modules/genotype_gvcfs'
include { DEEPVARIANT } from './modules/deepvariant'
include { RTG_VCF_EVAL } from './modules/rtg'
include { HAPPY } from './modules/happy'
include { BIASTOOLS } from './modules/biastools'

// --- Core params -----------------------------------------------------------
params.threads      = 16
params.gatk         = 'gatk'
params.outdir       = './results'
params.publish_mode = 'copy'

// --- Inputs ----------------------------------------------------------------
params.sample_id = null
params.ref       = null
params.fq1       = null
params.fq2       = null

// --- GATK params -----------------------------------------------------------
params.use_gatk    = true
params.use_bqsr    = true   // BaseRecalibrator+ApplyBQSR; set false to call directly on the markdup BAM (e.g. CHM13, no known-sites)
params.known_sites = null   // required only when use_gatk=true AND use_bqsr=true

// --- DeepVariant params ----------------------------------------------------
params.use_deepvariant     = false
params.deepvariant_sif     = null
params.disable_small_model = false
params.use_make_examples_extra_args = false
params.bind_path           = ''

// --- RTG params ------------------------------------------------------------
params.use_rtg               = false
params.rtg_sif               = null
params.rtg_tag               = 'BWA-MEM2'
params.giab_version          = 'v4.2.1'
params.giab_stratify_version = 'v2.0'

// --- hap.py params ---------------------------------------------------------
params.use_happy      = false
params.happy_home     = null
params.happy_sif      = null
params.happy_ref      = null
params.happy_assembly = 'GRCh38'   // truth-set assembly: 'GRCh38' or 'CHM13' (CHM13 forces v5.0q)

// --- biastools params ------------------------------------------------------
params.use_biastools  = false
params.biastools_bin  = 'biastools'
params.biastools_tag  = 'BWA-MEM2'
params.biastools_depth = 30

// ============================================================================
// Helper: parse .dict file to extract chromosome/contig names
// ============================================================================

def getChromosomes(dict_file) {
    def chroms = []
    new File(dict_file).eachLine { line ->
        if (line.startsWith('@SQ')) {
            def matcher = line =~ /SN:(\S+)/
            if (matcher.find()) {
                chroms << matcher.group(1)
            }
        }
    }
    return chroms
}

// ============================================================================
// Workflow
// ============================================================================

workflow {
    // --- Input validation --------------------------------------------------
    if (!params.fq1 || !params.fq2) {
        exit 1, "params.fq1 and params.fq2 are required"
    }
    if (!params.ref) {
        exit 1, "params.ref is required"
    }
    if (!params.sample_id) {
        exit 1, "params.sample_id is required"
    }
    if (params.use_gatk && params.use_bqsr && !params.known_sites) {
        exit 1, "params.known_sites is required when use_gatk=true and use_bqsr=true (comma-separated VCFs for BQSR)"
    }
    if (params.use_deepvariant && !params.deepvariant_sif) {
        exit 1, "params.deepvariant_sif is required when use_deepvariant=true"
    }
    if ((params.use_deepvariant || params.use_rtg || params.use_happy) && !params.bind_path) {
        exit 1, "params.bind_path is required when use_deepvariant/use_rtg/use_happy=true (data root bound into singularity containers)"
    }
    if (params.use_rtg) {
        if (!params.use_gatk && !params.use_deepvariant) {
            exit 1, "params.use_gatk or params.use_deepvariant must be true when use_rtg=true (a caller VCF is needed to benchmark)"
        }
        if (!params.rtg_sif) {
            exit 1, "params.rtg_sif is required when use_rtg=true"
        }
    }
    if (params.use_happy) {
        if (!params.use_gatk && !params.use_deepvariant) {
            exit 1, "params.use_gatk or params.use_deepvariant must be true when use_happy=true (a caller VCF is needed to benchmark)"
        }
        if (!params.happy_home || !params.happy_sif || !params.happy_ref) {
            exit 1, "params.happy_home, params.happy_sif and params.happy_ref are required when use_happy=true"
        }
    }
    if (params.use_biastools) {
        if (!params.happy_home) {
            exit 1, "params.happy_home is required when use_biastools=true (GIAB benchmark VCF source)"
        }
        if (!params.biastools_bin) {
            exit 1, "params.biastools_bin is required when use_biastools=true"
        }
    }

    // --- Shared reference / index / reads ----------------------------------
    ref           = file(params.ref)
    ref_idx_files = channel.fromPath("${params.ref}.*").collect()
    reads_ch      = channel.of(tuple(params.sample_id, file(params.fq1), file(params.fq2)))

    // --- Always on: alignment + mark duplicate -----------------------------
    BWA_MEM2(reads_ch, ref, ref_idx_files)
    SAMBAMBA_MARKDUP(BWA_MEM2.out.bam)

    // --- Collect PASS VCFs from each enabled caller for benchmarking ----------
    // benchmark_ch emits tuple(caller_name, vcf) per enabled caller.
    // (biastools is BAM-based and caller-independent -> runs once on the
    //  markduplicate BAM, not per caller.)
    benchmark_ch = channel.empty()

    // GIAB benchmark truth VCF (used by biastools; happy uses its own ref path)
    def benchmark_vcf = params.happy_home \
        ? "${params.happy_home}/data/giab_small_variant_benchmark/giab_${params.giab_version}/${params.sample_id}_${params.giab_version}_benchmark.vcf.gz" \
        : null

    // --- GATK germline variant calling -------------------------------------
    if (params.use_gatk) {
        ref_dict = file(params.ref.replaceAll(/\.fasta$|\.fa$/, '.dict'))

        // BQSR is optional: when use_bqsr=false (e.g. CHM13 has no known-sites),
        // HaplotypeCaller runs directly on the markdup BAM.
        def hc_bam_ch
        if (params.use_bqsr) {
            known_sites_list = params.known_sites.tokenize(',')
            known_sites_files   = known_sites_list.collect { ks -> file(ks) }
            known_sites_indices = known_sites_list.collect { ks -> file("${ks}.tbi") }
            GATK_BQSR(SAMBAMBA_MARKDUP.out.bam, ref, ref_idx_files, ref_dict, known_sites_files, known_sites_indices)
            hc_bam_ch = GATK_BQSR.out.bam
        } else {
            hc_bam_ch = SAMBAMBA_MARKDUP.out.bam
        }

        // HaplotypeCaller per chromosome
        chromosomes_ch = channel.fromList(getChromosomes(ref_dict.toString()))
        hc_input_ch = hc_bam_ch
            .combine(chromosomes_ch)
            .map { sample_id, bam, bai, chrom -> tuple(sample_id, bam, bai, chrom) }
        GATK_HAPLOTYPECALLER(hc_input_ch, ref, ref_idx_files, ref_dict)

        // Merge per-chromosome GVCFs (order must match ref dict)
        def chrom_order = getChromosomes(ref_dict.toString())
        merge_input_ch = GATK_HAPLOTYPECALLER.out.gvcf
            .groupTuple(by: 0)
            .map { sample_id, gvcfs, tbis ->
                def order_map = chrom_order.withIndex().collectEntries { chrom, idx -> [(chrom): idx] }
                def sorted_gvcfs = gvcfs.sort { g -> order_map[g.name.replaceAll(/.*\.(\S+)\.g\.vcf\.gz/, '$1')] }
                def sorted_tbis  = tbis.sort  { t -> order_map[t.name.replaceAll(/.*\.(\S+)\.g\.vcf\.gz\.tbi/, '$1')] }
                tuple(sample_id, sorted_gvcfs, sorted_tbis)
            }
        MERGE_GVCFS(merge_input_ch)

        // Genotype merged GVCFs -> final VCF (used directly by downstream benchmarking)
        GATK_GENOTYPE_GVCFS(MERGE_GVCFS.out.gvcf, ref, ref_idx_files, ref_dict)

        benchmark_ch = benchmark_ch.concat(
            GATK_GENOTYPE_GVCFS.out.vcf.map { _sid, vcf, _tbi -> tuple('gatk', vcf) }
        )
    }

    // --- DeepVariant variant calling ---------------------------------------
    if (params.use_deepvariant) {
        dv_bam_ch = SAMBAMBA_MARKDUP.out.bam.map { _sample_id, bam, _bai -> bam }
        dv_bai_ch = SAMBAMBA_MARKDUP.out.bam.map { _sample_id, _bam, bai -> bai }
        def dv_small_model_arg = params.disable_small_model ? '--disable_small_model=true' : ''
        def dv_make_examples_extra_arg = params.use_make_examples_extra_args ? '--make_examples_extra_args="min_mapping_quality=0,keep_legacy_allele_counter_behavior=true,normalize_reads=true"' : ''

        DEEPVARIANT(
            dv_bam_ch, dv_bai_ch,
            params.ref,
            params.sample_id, params.threads,
            dv_small_model_arg, dv_make_examples_extra_arg, params.deepvariant_sif, params.bind_path
        )

        benchmark_ch = benchmark_ch.concat(
            DEEPVARIANT.out.pass_vcf.map { vcf -> tuple('deepvariant', vcf) }
        )
    }

    // --- Benchmarking (against GIAB truth set), per enabled caller ----------
    // Each caller's PASS VCF is benchmarked independently. Results are
    // separated under <outdir>/benchmark/<caller>/{rtg,happy}/.
    if (params.use_rtg) {
        rtg_in_ch = benchmark_ch.map { caller, vcf ->
            tuple(
                vcf,
                "${params.rtg_tag}_${caller}",
                "rtg_benchmark_${caller}_${params.giab_version}_${params.giab_stratify_version}",
                caller,
                "benchmark/${caller}/rtg"
            )
        }
        RTG_VCF_EVAL(
            rtg_in_ch,
            params.sample_id, params.giab_version, params.giab_stratify_version,
            params.rtg_sif, params.bind_path
        )
    }

    if (params.use_happy) {
        happy_in_ch = benchmark_ch.map { caller, vcf ->
            tuple(vcf, caller, "benchmark/${caller}/happy")
        }
        HAPPY(
            happy_in_ch,
            params.sample_id, params.happy_ref, params.happy_assembly, params.giab_version,
            params.happy_home, params.happy_sif, params.bind_path
        )
    }

    // --- Mapping-bias analysis (biastools), single run on the markdup BAM ---
    // biastools measures alignment mapping bias, which is a property of the
    // alignment (caller-independent), so it runs once on the markduplicate BAM.
    // Results are published under <outdir>/biastools/.
    if (params.use_biastools) {
        biastools_in_ch = SAMBAMBA_MARKDUP.out.bam
            .map { _sample_id, bam, bai -> tuple(bam, bai, params.biastools_tag) }
        BIASTOOLS(
            biastools_in_ch,
            params.ref, benchmark_vcf,
            params.sample_id, params.biastools_depth, params.threads, params.biastools_bin
        )
    }
}
