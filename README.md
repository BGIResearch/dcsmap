# dcsmap

[中文文档](doc/README_zh.md)

`dcsmap` is a fast and accurate pangenome-aware read aligner developed as part of DCSTools. It employs a hybrid alignment strategy that combines linear-reference alignment with pangenome graph alignment. Compared with the vg toolkit, `dcsmap` achieves a 2.5-fold improvement in alignment speed(from ~5.5h to ~2.17h on a 32c128g server) while reducing SNP/Indel errors by 7–12% when evaluated using the same downstream DeepVariant variant-calling pipeline. A machine-learning model is used to identify low-confidence alignments from the linear alignment results, which are then selectively realigned using vg giraffe, improving alignment accuracy while avoiding the computational cost of graph-based alignment for all reads. `dcsmap` incorporates optimized vg toolkit and ABRA2 for efficient haplotype sampling, graph alignment and indel realignment.

## Build

```bash
mkdir build && cd build
cmake ..
make -j
```

The binary `dcsmap` is output to `build/`. It is **purely statically linked**
(`-static`), so it has no runtime shared-library dependencies and can be copied
to any Linux machine of the same architecture.

### Example

```bash
dcsmap \
    --fq1  sample.R1.fastq.gz \
    --fq2  sample.R2.fastq.gz \
    --out-bam  /out/sample.bam \
    --ref-fasta  /data/hprc-v2.0-mc-grch38.ref.fasta \
    --gbz  /data/hprc-v2.0-mc-grch38.gbz \
    --hapl  /data/hprc-v2.0-mc-grch38.hapl \
    --graph-ref-contigs  /data/hprc-v2.0-mc-grch38.ref.pathnames \
    --extract-model  /data/extract.model \
    --sample-name HG002 \
    --threads 32
```

`dcsmap` relies on the DCSTools environment: it invokes `vg`, `samtools`,
`bwa-mem2`, `abra2.jar`, etc. from a tools root directory (containing `libexec/`
and `jar/` subdirectories) and requires Java 8 for ABRA2.

- `--tools-root`: defaults to the parent of the `dcsmap` executable's directory
  (e.g. if `dcsmap` is at `<DCS_HOME>/libexec/dcsmap`, the default is
  `<DCS_HOME>`). If `dcsmap` is placed elsewhere, set `DCS_HOME` env
  var or pass `--tools-root` explicitly. `--tools-root` takes precedence over
  `DCS_HOME`.
- `--java-home`: defaults to the `JAVA_HOME` env var. Must point to a Java 8
  installation (required by ABRA2). Pass `--java-home` explicitly to override on
  a per-invocation basis.

```bash
# Optional: override via env vars instead of CLI flags
export DCS_HOME=/path/to/dcstools
export JAVA_HOME=/path/to/jdk8
export PATH="${DCS_HOME}/libexec:${PATH}"
```

**Graph and linear reference inputs**

- `--gbz`: download an HPRC minigraph-cactus gbz graph (see [HPRC Graphs](#hprc-graphs)).
- `--hapl`, `--graph-ref-contigs` (`<prefix>.ref.pathnames`), and
  `--ref-fasta` (`<prefix>.ref.fasta` + `.fai` + `.dict`): all three are produced
  by `scripts/build_hapl_fasta.sh` from the gbz (see [Index preparation](#index-preparation), step 1).
- The bwa-mem2 index files (`.0123`, `.amb`, `.ann`, `.bwt.2bit.64`, `.pac`)
  used by `linear_align_extract` are looked up next to `--ref-fasta`; build them
  with `bwa-mem2 index` (see [Index preparation](#index-preparation), step 2).
- `--extract-model`: the machine-learning model used by `extract-bam` to flag
  low-confidence linear alignments for graph realignment. Obtain it from the
  DCSTools distribution (`$DCS_HOME/share/dcsmap/extract.model` or similar).

## Usage

```
Program: dcsmap
version: 1.0.1

Usage: dcsmap [-option]

Required:
  --fq1 <file>                          input read1 fastq path
  --fq2 <file>                          input read2 fastq path
  --out-bam <file>                      output bam path (work dir defaults to its dirname)
  --ref-fasta <file>                    reference fasta path (with .fai + .dict + aligner index)
  --gbz <file>                          vg gbz graph path
  --hapl <file>                         vg hapl index path
  --graph-ref-contigs <file>            vg ref-paths file
  --extract-model <file>                model file for extract-bam

Options:
  --sample-name <str>                   sample name (default: SAMPLE)
  --platform <str>                      sequencing platform (default: DNBSEQ)
  --threads <int>                       number of threads to use (default: 32)
  --mode <str>                          workflow mode (default: dcsmap)
                                        available options: {dcsmap, dcsmap-m1}
  --tools-root <dir>                    tools root dir (libexec + jar)
                                        default: DCS_HOME env or parent-of-exe-dir
  --java-home <dir>                     JAVA home (must be Java 8)
                                        default: JAVA_HOME env
  --work-dir <dir>                      work dir (default: out-bam-dir/work.XXXXXX)
  --parallel <bool>                     run linear_align_extract and vg_haplotype in parallel (default: true)
                                        set to false to run vg_haplotype first, then linear_align_extract
  --clean <bool>                        clean work dir on success, keeping command.sh/logs/rc (default: true)
  -h, --help                            display help message
  --version                             display version message
```

## Workflow

![dcsmap alignment workflow](./doc/image/dcsmap_alignment_workflow.png)

### dcsmap mode

- `linear_align_extract` and `vg_haplotype` run **in parallel** by default
  (`--parallel true`). With `--parallel false`, `vg_haplotype` runs first,
  then `linear_align_extract`.
- `realign` runs ABRA2 on the giraffe BAM.
- `merge_bams` merges the ABRA2 BAM and the linear not-extract BAM into the
  final `--out-bam`.

### dcsmap-m1 mode

- `merge_bams` merges the linear not-extract BAM and the giraffe BAM.
- `realign` runs ABRA2 on the merged BAM to produce the final `--out-bam`.

### Tasks

| Task                  | Description                                        |
|-----------------------|----------------------------------------------------|
| `linear_align_extract`| BWA-MEM2 alignment + extract-bam (extract fastq + not-extract bam) |
| `vg_haplotype`        | vg haplotype sampling (produces personalized-gbz, dist, min, zipcodes) |
| `vg_giraffe`          | vg giraffe mapping on extract fastq                |
| `merge_bams`          | Merge BAMs (samtools merge)                        |
| `realign`             | ABRA2 realignment                                  |

## Work directory

- Default: `<dirname(--out-bam)>/work.XXXXXX` (created with `mktemp`).
- Override with `--work-dir`.
- Each task runs in its own subdirectory: `<work_root>/task-<name>`.
- Each task writes `command.sh`, `command.sh.log.o`, `command.sh.log.e`, and
  `command.sh.rc`.
- On success, if `--clean true` (default), intermediate files are removed but
  `command.sh`, logs, and rc are preserved for traceability.

## HPRC Graphs

dcsmap uses minigraph-cactus (mc) pangenome graphs from the
[Human Pangenome Reference Consortium (HPRC)](https://humanpangenome.org/).
Download the gbz graph for your reference coordinate system from the S3 URLs
below. See `data/hprc_graph_urls.tsv` for the full list.

| Coordinate | Version | Size | URL |
|------------|---------|------|-----|
| GRCh38 | v2 | 5.4 GB | https://s3-us-west-2.amazonaws.com/human-pangenomics/pangenomes/freeze/release2/minigraph-cactus/v2.0/hprc-v2.0-mc-grch38/hprc-v2.0-mc-grch38.gbz |
| CHM13 | v2 | 5.7 GB | https://s3-us-west-2.amazonaws.com/human-pangenomics/pangenomes/freeze/release2/minigraph-cactus/v2.0/hprc-v2.0-mc-chm13/hprc-v2.0-mc-chm13.gbz |

> **Note:** See `data/hprc_graph_urls.tsv` for additional versions (v1.1, v2-eval) and filtered (d46) graphs.

## Index preparation

dcsmap requires two sets of indices: the vg graph-side indices (hapl, ref fasta,
pathnames) and the bwa-mem2 linear-side index. Prepare them in two steps:

### 1. Build vg graph indices (scripts/build_hapl_fasta.sh)

Generates the vg hapl index, reference fasta (+ .fai / .dict), and the ref
pathnames file from a vg gbz graph. Reads `DCS_HOME` from the environment
and expects `vg`, `samtools` under `$DCS_HOME/libexec`.

```
Usage: build_hapl_fasta.sh <gbz> <ref_path> [ref_dict] [out_prefix]

  <gbz>        vg gbz graph, e.g. hprc-v2.0-mc-grch38.gbz
  <ref_path>   reference path prefix in the graph, e.g. GRCh38 or CHM13
  [ref_dict]   SAM dict of the linear reference (optional)
               If given, contigs are ordered by the dict's @SQ entries;
               otherwise natural chromosome order: chr1..chr22, chrX, chrY, chrM, then the rest.
               For GRCh38, use data/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.dict
               bundled in this project.
  [out_prefix] output path prefix (default: <gbz_dir>/<gbz_basename_without_.gbz>,
               i.e. outputs are written next to the input gbz)

Outputs:
  <prefix>.ref.fasta / .ref.fasta.fai / .ref.fasta.dict
  <prefix>.ref.pathnames
  <prefix>.hapl  (+ .snarls / .xg / .ri / .dist)
```

Set `SKIP_INDEX=1` to skip the snarls/xg/ri/dist/hapl build (useful for testing
the fasta + pathnames generation only).

1) GRCh38 linear reference coordinate system
```bash
export DCS_HOME=/path/to/dcstools
scripts/build_hapl_fasta.sh /data/hprc-v2.0-mc-grch38.gbz GRCh38 \
    data/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.dict
# Outputs written to /data/hprc-v2.0-mc-grch38.* (next to the input gbz)
```

2) CHM13 T2T linear reference coordinate system
```bash
export DCS_HOME=/path/to/dcstools
scripts/build_hapl_fasta.sh /data/hprc-v2.0-mc-chm13.gbz CHM13
# Outputs written to /data/hprc-v2.0-mc-chm13.* (next to the input gbz)
# (no ref_dict needed — CHM13 has only chr1-22/X/Y/M, natural order is used)
```

**build vg hapl index resource consumption**

Measured on the HPRC v2 GRCh38 Default graph (`hprc-v2.0-mc-grch38.gbz`,
5.4 GB). 

| Step | Command | Wall time | CPU % | Peak RSS | Output size |
|------|---------|-----------|-------|----------|-------------|
| fasta + pathnames | `vg paths` / `samtools faidx` / `dict` | ~2m 20s | — | ~3 GB | 3.0 GB |
| snarls | `vg snarls` | 29m 51s | 7098% | 67.8 GB | 215 MB |
| xg | `vg convert -x --drop-haplotypes` | 9m 58s | 222% | 36.8 GB | 8.4 GB |
| ri | `vg gbwt -Z -r` | 3m 28s | 6085% | 46.4 GB | 9.5 GB |
| dist | `vg index -j` | 1h 00m | 96% | **262 GB** | **105.5 GB** |
| hapl | `vg haplotypes -H` | 1h 53m | 1078% | 63.9 GB | 19.9 GB |
| **Total** | | **~3h 40m** | — | — | **~148 GB** |

> **Note:** `vg index -j` (dist) is the bottleneck — single-threaded, ~262 GB
> peak memory, and ~105 GB output (71% of the total).

Recommended: 32+ cores, 300 GB memory, 300 GB SSD (the 300 GB memory floor is
driven by `vg index -j` dist, which peaks at ~262 GB RSS).

### 2. Build bwa-mem2 index

`dcsmap`'s linear aligner (bwa-mem2) looks up its index files (`.0123`, `.amb`,
`.ann`, `.bwt.2bit.64`, `.pac`) next to the `--ref-fasta` path. Build them from
the `<prefix>.ref.fasta` produced in step 1:

```bash
bwa-mem2 index <prefix>.ref.fasta
```

This writes the index files next to the fasta. The `.fai` and `.dict` already
exist from step 1.

**bwa-mem2 index resource consumption**

Measured on the same GRCh38 reference fasta (`hprc-v2.0-mc-grch38.ref.fasta`,
3.0 GB):

| Step | Wall time | CPU % | Peak RSS | Output size |
|------|-----------|-------|----------|-------------|
| `bwa-mem2 index` | 16m 51s | 99% | 69.3 GB | 16.9 GB |

Output breakdown: `.bwt.2bit.64` (9.4 GB), `.0123` (5.8 GB), `.pac` (738 MB),
`.amb` / `.ann` (< 30 KB).

Recommended: 2+ cores, 128 GB memory, 50 GB SSD (the ~70 GB peak RSS of
`bwa-mem2 index` dominates the memory requirement).

## License

The `dcsmap` workflow (this project) is released under the **BSD 2-Clause
License**. See the `LICENSE` file for the full text.

`dcsmap` is part of **DCSTools**, which is **commercially licensed software**.
Use, redistribution, or modification of the bundled DCSTools components
(including `extract-bam`) requires a valid commercial license from the DCSTools
copyright holder. No open-source license is granted for the DCSTools components
themselves. Contact the DCSTools maintainers for licensing terms and
acquisition.

The `dcsmap` workflow additionally bundles and invokes several third-party
open-source tools, each distributed under its own license. Users must comply
with the terms of every license listed below when redistributing or running
`dcsmap`.

| Tool | Used by | License |
|------|---------|---------|
| extract-bam | `linear_align_extract` | **Commercial (DCSTools)** |
| vg | `vg_haplotype`, `vg_giraffe` | MIT |
| samtools / htslib | `linear_align_extract`, `merge_bams` | MIT/BSD |
| bwa-mem2 | `linear_align_extract` | MIT |
| ABRA2 | `realign` | MIT |
| KMC | `extract-bam` k-mer counting | GPL-3.0 |
