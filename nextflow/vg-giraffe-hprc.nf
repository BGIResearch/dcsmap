#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { VG_GIRAFFE   } from './modules/vg_giraffe'
include { REALIGN      } from './modules/realign'
include { DEEPVARIANT  } from './modules/deepvariant'
include { RTG_VCF_EVAL } from './modules/rtg'
include { HAPPY        } from './modules/happy'
include { BIASTOOLS    } from './modules/biastools'

params.fq1                   = null
params.fq2                   = null
params.sample_name           = 'SAMPLE'
params.platform              = 'Illumina'
params.threads               = 32

params.tools_root            = './tools'
params.java_home             = ''

params.ref_fasta             = null
params.gbz                   = null
params.dist                  = null
params.min                   = null
params.zipcodes              = null
params.ref_contigs           = null
// Reference genome type: 'GRCh38' (default) or 'CHM13'. Used to infer the
// vg giraffe surject BAM contig prefix to strip ("GRCh38#0#" / "CHM13#0#")
// so the BAM matches the plain contig names of ref_fasta.
params.ref_type              = null

params.deepvariant_sif       = null
params.disable_small_model   = false
params.use_make_examples_extra_args = true

params.use_rtg               = false
params.use_happy             = false
params.use_biastools         = false

params.rtg_sif               = null
params.rtg_tag               = 'VG-GIRAFFE-DV-hprc'
params.giab_version          = 'v4.2.1'
params.giab_stratify_version = 'v2.0'

params.happy_home            = null
params.happy_sif             = null
params.happy_ref             = null
params.happy_assembly        = 'GRCh38'   // truth-set assembly: 'GRCh38' or 'CHM13' (CHM13 forces v5.0q)

params.biastools_bin         = 'biastools'
params.biastools_tag         = 'VG-GIRAFFE-hprc'
params.biastools_depth       = 30

params.bind_path             = ''
params.outdir                = './results'
params.publish_mode          = 'link'
params.publish_abra2_bam     = true
params.use_gkl               = true


workflow {

    if (!params.fq1 || !params.fq2) {
        exit 1, "params.fq1 and params.fq2 are required"
    }
    if (!params.ref_fasta) {
        exit 1, "params.ref_fasta is required"
    }
    if (!params.gbz || !params.dist || !params.min || !params.zipcodes || !params.ref_contigs) {
        exit 1, "params.gbz, params.dist, params.min, params.zipcodes and params.ref_contigs are required"
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

    def ref_fasta      = params.ref_fasta
    def ref_fasta_fai  = "${params.ref_fasta}.fai"
    def ref_fasta_dict = "${params.ref_fasta}.dict"

    // Infer the vg giraffe surject BAM contig prefix from ref_type
    // ("GRCh38#0#" / "CHM13#0#") so the BAM matches ref_fasta contig names.
    def ref_path_prefix = params.ref_type ? "${params.ref_type}#0#" : 'GRCh38#0#'

    def fq1      = file(params.fq1)
    def fq2      = file(params.fq2)
    def gbz      = file(params.gbz)
    def dist     = file(params.dist)
    def min      = file(params.min)
    def zipcodes = file(params.zipcodes)

    def dv_small_model_arg = params.disable_small_model ? '--disable_small_model=true' : ''
    def dv_make_examples_extra_arg = params.use_make_examples_extra_args ? '--make_examples_extra_args="min_mapping_quality=0,keep_legacy_allele_counter_behavior=true,normalize_reads=true"' : ''
    def rtg_outname        = "rtg_benchmark_${params.giab_version}_${params.giab_stratify_version}"
    def benchmark_vcf      = params.happy_home \
        ? "${params.happy_home}/data/giab_small_variant_benchmark/giab_${params.giab_version}/${params.sample_name}_${params.giab_version}_benchmark.vcf.gz" \
        : null

    // No VG_HAPLOTYPE_SAMPLING: map straight onto the full input pangenome graph.
    VG_GIRAFFE(
        params.tools_root,
        fq1, fq2,
        params.sample_name, params.platform,
        gbz, dist, min, zipcodes,
        params.ref_contigs, params.threads, ref_path_prefix
    )

    REALIGN(
        params.tools_root, params.java_home,
        VG_GIRAFFE.out.surject_bam, VG_GIRAFFE.out.surject_bai,
        params.sample_name,
        ref_fasta, ref_fasta_fai, ref_fasta_dict, params.threads,
        'giraffe'
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
            params.sample_name, params.happy_ref, params.happy_assembly, params.giab_version,
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
