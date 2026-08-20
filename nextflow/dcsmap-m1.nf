#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

// ============================================================================
// DCSMAP-M1: merge-first variant of DCSMAP.
//
// Differs from dcsmap.nf only in step ordering: dcsmap.nf realigns the VG
// Giraffe surject BAM with abra2 first and then merges it with the BWA-MEM2
// BAM, whereas DCSMAP-M1 merges the BWA-MEM2 BAM and the VG Giraffe surject
// BAM first and runs a single abra2 realignment pass on the merged BAM.
// ============================================================================

include { BWA_MEM2_MAP_AND_EXTRACT } from './modules/bwa-mem2-extract'
include { VG_HAPLOTYPE_SAMPLING    } from './modules/vg_haplotype'
include { VG_GIRAFFE               } from './modules/vg_giraffe'
include { MERGE_BAMS               } from './modules/merge_bams'
include { REALIGN                  } from './modules/realign'
include { DEEPVARIANT              } from './modules/deepvariant'
include { RTG_VCF_EVAL             } from './modules/rtg'
include { HAPPY                    } from './modules/happy'
include { BIASTOOLS                } from './modules/biastools'

params.fq1                   = null
params.fq2                   = null
params.sample_name           = 'SAMPLE'
params.platform              = 'Illumina'
params.threads               = 32

params.tools_root            = './tools'
params.java_home             = ''
params.extract_model         = null

params.ref_fasta             = null
params.gbz                   = null
params.hapl                  = null
params.ref_contigs           = null

// Reference genome type selector. When set to 'GRCh38' or 'CHM13', the
// ref_fasta / gbz / hapl / ref_contigs paths are auto-derived from the
// corresponding HPRC data directory below, unless they are provided
// explicitly (explicit values take precedence over ref_type-derived ones).
// ref_type also drives: (1) the --set-reference value passed to vg
// haplotypes and (2) the vg giraffe surject BAM contig prefix to strip,
// both of which must match the reference path-name prefix embedded in the
// graph ('GRCh38' or 'CHM13'). Defaults to 'GRCh38' when unset.
params.ref_type              = null

params.deepvariant_sif       = null
params.disable_small_model   = false
params.use_make_examples_extra_args = true

params.use_rtg               = false
params.use_happy             = false
params.use_biastools         = false

params.rtg_sif               = null
params.rtg_tag               = 'DCSMAP-M1-DV'
params.giab_version          = 'v4.2.1'
params.giab_stratify_version = 'v2.0'

params.happy_home            = null
params.happy_sif             = null
params.happy_ref             = null
params.happy_giab_version    = null   // optional override for hap.py GIAB version (v4.2.1 / v5.0q); CHM13 forces v5.0q

params.biastools_bin         = 'biastools'
params.biastools_tag         = 'DCSMAP-M1'
params.biastools_depth       = 30

params.bind_path             = ''
params.outdir                = './results'
params.publish_mode          = 'link'
params.publish_abra2_bam      = true   // M1: abra2 BAM is the final output (realigned merged BAM)
params.publish_merged_bam    = false   // M1: merged BAM is only an intermediate input to REALIGN
params.use_gkl               = true


workflow {

    // Effective reference resources. Start from the explicit params; when
    // params.ref_type is set, auto-derive ref_fasta / gbz / hapl / ref_contigs
    // from the corresponding HPRC data directory (explicit values win) and
    // use ref_type as the vg haplotypes --set-reference value.
    def ref_fasta     = params.ref_fasta
    def gbz           = params.gbz
    def hapl          = params.hapl
    def ref_contigs   = params.ref_contigs
    def set_reference = params.ref_type ?: 'GRCh38'

    if (params.ref_type) {
        def ref_type_dirs = [
            'GRCh38': '/zdswhst2/ST_BIOINTEL/P24Z32400N0004/yangqi4/data/hprc',
            'CHM13' : '/zdswhst2/ST_BIOINTEL/P24Z32400N0004/yangqi4/data/hprc-v2-chm13',
        ]
        def ref_type_prefix = [
            'GRCh38': 'hprc-v2.0-mc-grch38',
            'CHM13' : 'hprc-v2.0-mc-chm13',
        ]
        if (!(params.ref_type in ref_type_dirs)) {
            exit 1, "params.ref_type must be 'GRCh38' or 'CHM13' (got: ${params.ref_type})"
        }
        def data_root = ref_type_dirs[params.ref_type]
        def prefix    = ref_type_prefix[params.ref_type]
        if (!ref_fasta)   ref_fasta   = "${data_root}/${prefix}.ref.fasta"
        if (!gbz)         gbz         = "${data_root}/${prefix}.gbz"
        if (!hapl)        hapl        = "${data_root}/${prefix}.hapl"
        if (!ref_contigs) ref_contigs = "${data_root}/${prefix}.ref.pathnames"
    }

    // hap.py assembly / GIAB version. The assembly follows the resolved
    // reference (set_reference). For CHM13 the v5.0q truth set is required
    // (only available for HG002); for GRCh38 the configured giab_version
    // (default v4.2.1) is used unless params.happy_giab_version overrides it.
    def happy_assembly    = set_reference
    def happy_giab_version = params.happy_giab_version
        ?: (params.ref_type == 'CHM13' ? 'v5.0q' : params.giab_version)

    // Strip the vg giraffe surject BAM contig prefix ("<set_reference>#0#")
    // so it matches the plain contig names of ref_fasta.
    def ref_path_prefix = "${set_reference}#0#"

    if (!params.fq1 || !params.fq2) {
        exit 1, "params.fq1 and params.fq2 are required"
    }
    if (!ref_fasta) {
        exit 1, "params.ref_fasta is required (or set params.ref_type to 'GRCh38'/'CHM13')"
    }
    if (!gbz || !hapl || !ref_contigs) {
        exit 1, "params.gbz, params.hapl and params.ref_contigs are required (or set params.ref_type to 'GRCh38'/'CHM13')"
    }
    if (!params.extract_model) {
        exit 1, "params.extract_model is required"
    }
    if (!params.deepvariant_sif) {
        exit 1, "params.deepvariant_sif is required (DeepVariant is a core step)"
    }
    if (!params.bind_path) {
        exit 1, "params.bind_path is required (data root bound into singularity containers)"
    }
    if (params.use_rtg && !params.rtg_sif) {
        exit 1, "params.rtg_sif is required when use_rtg=true"
    }
    if (params.use_happy && (!params.happy_home || !params.happy_sif || !params.happy_ref)) {
        exit 1, "params.happy_home, params.happy_sif and params.happy_ref are required when use_happy=true"
    }
    if (params.use_biastools && !params.happy_home) {
        exit 1, "params.happy_home is required when use_biastools=true (benchmark VCF source)"
    }

    def ref_fasta_fai  = "${ref_fasta}.fai"
    def ref_fasta_dict = "${ref_fasta}.dict"

    def dv_small_model_arg = params.disable_small_model ? '--disable_small_model=true' : ''
    def dv_make_examples_extra_arg = params.use_make_examples_extra_args ? '--make_examples_extra_args="min_mapping_quality=0,keep_legacy_allele_counter_behavior=true,normalize_reads=true"' : ''
    def rtg_outname        = "rtg_benchmark_${params.giab_version}_${params.giab_stratify_version}"
    def benchmark_vcf      = params.happy_home \
        ? "${params.happy_home}/data/giab_small_variant_benchmark/giab_${params.giab_version}/${params.sample_name}_${params.giab_version}_benchmark.vcf.gz" \
        : null

    BWA_MEM2_MAP_AND_EXTRACT(
        params.tools_root,
        params.fq1, params.fq2, ref_fasta,
        params.sample_name, params.platform, params.extract_model, params.threads
    )

    VG_HAPLOTYPE_SAMPLING(
        params.tools_root,
        params.fq1, params.fq2, gbz, hapl, params.threads,
        set_reference
    )

    VG_GIRAFFE(
        params.tools_root,
        BWA_MEM2_MAP_AND_EXTRACT.out.extract_fq1, BWA_MEM2_MAP_AND_EXTRACT.out.extract_fq2,
        params.sample_name, params.platform,
        VG_HAPLOTYPE_SAMPLING.out.haplotype_gbz, VG_HAPLOTYPE_SAMPLING.out.haplotype_dist,
        VG_HAPLOTYPE_SAMPLING.out.haplotype_min, VG_HAPLOTYPE_SAMPLING.out.haplotype_zipcodes,
        ref_contigs, params.threads, ref_path_prefix
    )

    // Merge first: BWA-MEM2 BAM + VG Giraffe surject BAM (pre-realign).
    MERGE_BAMS(
        params.tools_root,
        params.sample_name,
        BWA_MEM2_MAP_AND_EXTRACT.out.out_bam, VG_GIRAFFE.out.surject_bam, params.threads
    )

    // Then realign once on the merged BAM.
    REALIGN(
        params.tools_root, params.java_home,
        MERGE_BAMS.out.merged_bam, MERGE_BAMS.out.merged_bam_bai,
        params.sample_name,
        ref_fasta, ref_fasta_fai, ref_fasta_dict, params.threads,
        'merged_realign'
    )

    DEEPVARIANT(
        REALIGN.out.abra2_bam, REALIGN.out.abra2_bai,
        ref_fasta,
        params.sample_name, params.threads,
        dv_small_model_arg, dv_make_examples_extra_arg, params.deepvariant_sif, params.bind_path
    )

    if (params.use_rtg) {
        rtg_input_ch = DEEPVARIANT.out.pass_vcf.map { vcf -> tuple(vcf, params.rtg_tag, rtg_outname, 'deepvariant', 'benchmark/deepvariant/rtg') }
        RTG_VCF_EVAL(
            rtg_input_ch,
            params.sample_name,
            params.giab_version, params.giab_stratify_version,
            params.rtg_sif, params.bind_path
        )
    }

    if (params.use_happy) {
        happy_input_ch = DEEPVARIANT.out.pass_vcf.map { vcf -> tuple(vcf, 'deepvariant', 'benchmark/deepvariant/happy') }
        HAPPY(
            happy_input_ch,
            params.sample_name, params.happy_ref, happy_assembly, happy_giab_version,
            params.happy_home, params.happy_sif, params.bind_path
        )
    }

    if (params.use_biastools) {
        biastools_in_ch = REALIGN.out.abra2_bam
            .combine(REALIGN.out.abra2_bai)
            .map { bam, bai -> tuple(bam, bai, params.biastools_tag) }
        BIASTOOLS(
            biastools_in_ch,
            ref_fasta, benchmark_vcf,
            params.sample_name, params.biastools_depth, params.threads, params.biastools_bin
        )
    }
}
