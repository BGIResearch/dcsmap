#!/bin/bash
# build_hapl_fasta.sh — generate the vg hapl index and the reference fasta
# (+ .fai / .dict) and the ref pathnames file from a vg gbz graph.
#
# These are the graph-side inputs dcsmap needs (besides the gbz itself):
#   --hapl                <prefix>.hapl
#   --ref-fasta           <prefix>.ref.fasta        (+ .fai / .dict)
#   --graph-ref-contigs   <prefix>.ref.pathnames
#
# Usage: build_hapl_fasta.sh <gbz> <ref_path> [ref_dict] [out_prefix]
#   <gbz>        vg gbz graph, e.g. hprc-v2.0-mc-grch38.gbz
#   <ref_path>   reference path prefix in the graph, e.g. GRCh38 or CHM13.
#                Paths matching "<ref_path>#0#<contig>" are extracted.
#   [ref_dict]   SAM dict of the linear reference, e.g.
#                GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.dict
#                Its @SQ SN: entries define the contig set and order to keep.
#                If omitted, contigs are derived from the graph in natural
#                chromosome order: chr1..chr22, chrX, chrY, chrM, then the rest.
#   [out_prefix] output path prefix (default: <gbz_dir>/<gbz_basename_without_.gbz>,
#                i.e. outputs are written next to the input gbz)
#
# Outputs:
#   <prefix>.ref.fasta / .ref.fasta.fai / .ref.fasta.dict
#   <prefix>.ref.pathnames
#   <prefix>.hapl  (+ .snarls / .xg / .ri / .dist)
set -euo pipefail

: "${DCS_HOME:=/path/to/dcstools}"
export DCS_HOME
export PATH="${DCS_HOME}/libexec:${PATH}"

if [ $# -lt 2 ]; then
    echo "Usage: $0 <gbz> <ref_path> [ref_dict] [out_prefix]" >&2
    exit 1
fi

gbz=$1
ref_path=$2
ref_dict=${3:-}
# Default output prefix: same directory as the input gbz, with the gbz basename.
if [ $# -ge 4 ]; then
    prefix=$4
else
    gbz_dir=$(cd "$(dirname "${gbz}")" && pwd)
    prefix="${gbz_dir}/$(basename "${gbz}" .gbz)"
fi

echo "[INFO] $(date '+%F %T') gbz=${gbz} ref_path=${ref_path} ref_dict=${ref_dict:-<none>} prefix=${prefix}"

# --- full path listing from the graph (written to a file to avoid SIGPIPE
#     from `head`/`awk ... exit` closing the pipe under `set -o pipefail`) ---
vg paths -L -R -x "${gbz}" > "${prefix}.pathnames.raw"

# --- filter pathnames.raw to the requested ref_path, into pathnames.tmp.
#     Both ${prefix}.contigs and ${prefix}.ref.pathnames are derived from
#     this file so they always agree on the contig set. ---
grep "${ref_path}#0#" "${prefix}.pathnames.raw" > "${prefix}.pathnames.tmp"

# --- contigs to keep — ALWAYS the complete set of plain contigs from
#     pathnames.tmp. Ordering:
#     - if ref_dict is given: dict @SQ order first (intersected with the
#       graph), then any graph-only contigs in natural chromosome order;
#     - otherwise: natural chromosome order chr1..chr22, chrX, chrY, chrM,
#       then the rest. ---
awk -v rp="${ref_path}#0#" -v dict_file="${ref_dict:-}" '
    BEGIN {
        if (dict_file != "") {
            while ((getline line < dict_file) > 0) {
                if (line ~ /^@SQ/) {
                    nf = split(line, f, "\t")
                    for (i = 1; i <= nf; i++) {
                        if (f[i] ~ /^SN:/) {
                            sub(/^SN:/, "", f[i])
                            dict[f[i]] = 1
                            dict_order[++d] = f[i]
                        }
                    }
                }
            }
            close(dict_file)
        }
    }
    {
        c = $0; sub(rp, "", c)
        if (index(c, "[") > 0) next
        graph[c] = 1
    }
    END {
        if (dict_file != "") {
            for (i = 1; i <= d; i++) {
                c = dict_order[i]
                if (c in graph) { print c; emitted[c] = 1 }
            }
        }
        m = 0
        for (k in graph) {
            if (!(k in emitted)) remaining[++m] = k
        }
        for (i = 1; i <= m; i++) {
            key = remaining[i]; sub(/^chr/, "", key)
            if (key ~ /^[0-9]+$/)      sortkey[i] = sprintf("%03d", key + 0)
            else if (key == "X")       sortkey[i] = "023"
            else if (key == "Y")       sortkey[i] = "024"
            else if (key == "M")       sortkey[i] = "025"
            else                       sortkey[i] = "999" remaining[i]
        }
        for (i = 1; i <= m; i++)
            for (j = i + 1; j <= m; j++)
                if (sortkey[i] > sortkey[j]) {
                    t = sortkey[i]; sortkey[i] = sortkey[j]; sortkey[j] = t
                    t = remaining[i]; remaining[i] = remaining[j]; remaining[j] = t
                }
        for (i = 1; i <= m; i++) print remaining[i]
    }
' "${prefix}.pathnames.tmp" > "${prefix}.contigs"

# --- reference path names from the graph, one-to-one with the contigs list,
#     emitted in contigs-list order (dict @SQ order or natural chromosome order). ---
awk -v rp="${ref_path}#0#" '
    NR==FNR { dict[$0] = 1; order[++n] = $0; next }
    {
        c = $0; sub(rp, "", c)
        if (c in dict) paths[c] = $0
    }
    END {
        for (i = 1; i <= n; i++) {
            c = order[i]
            if (c in paths) print paths[c]
        }
    }
' "${prefix}.contigs" "${prefix}.pathnames.tmp" > "${prefix}.ref.pathnames"
rm -f "${prefix}.pathnames.tmp"

# --- extract reference fasta from the graph, subset/reorder to contigs order ---
vg paths -x "${gbz}" -F -Q "${ref_path}" "${gbz}" \
    | sed "s/${ref_path}#0#//g" > "${prefix}.ref.tmp.fasta"
samtools faidx "${prefix}.ref.tmp.fasta"
samtools faidx "${prefix}.ref.tmp.fasta" $(paste -sd' ' "${prefix}.contigs") > "${prefix}.ref.fasta"
samtools faidx "${prefix}.ref.fasta"
samtools dict "${prefix}.ref.fasta" > "${prefix}.ref.fasta.dict"
rm -f "${prefix}.ref.tmp.fasta" "${prefix}.ref.tmp.fasta.fai" "${prefix}.pathnames.raw"

# --- build the vg hapl index (depends on snarls / xg / ri / dist) ---
if [ "${SKIP_INDEX:-0}" != "1" ]; then
    echo "[INFO] $(date '+%F %T') building vg graph indexes (snarls/xg/ri/dist)..."
    /usr/bin/time -v vg snarls "${gbz}" > "${prefix}.snarls"
    /usr/bin/time -v vg convert -x --drop-haplotypes "${gbz}" > "${prefix}.xg"
    /usr/bin/time -v vg gbwt -Z "${gbz}" -r "${prefix}.ri"
    /usr/bin/time -v vg index -j "${prefix}.dist" "${gbz}"
    /usr/bin/time -v vg haplotypes -H "${prefix}.hapl" "${gbz}"
else
    echo "[INFO] $(date '+%F %T') SKIP_INDEX=1 — skipping snarls/xg/ri/dist/hapl"
fi

echo "[INFO] $(date '+%F %T') done."
echo "  ${prefix}.ref.fasta (+ .fai / .dict)"
echo "  ${prefix}.ref.pathnames"
if [ "${SKIP_INDEX:-0}" != "1" ]; then
    echo "  ${prefix}.hapl  (+ .snarls / .xg / .ri / .dist)"
fi
