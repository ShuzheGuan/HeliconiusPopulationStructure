#!/usr/bin/env python3

# Purpose: split non-CDS VCFs into BPP locus VCFs.
# Input: per-contig non-CDS VCF files and an exon BED file.
# Output: one directory of locus-level VCFs per contig.
# Software: Python, cyvcf2, bgzip, bcftools

from pathlib import Path
import subprocess
from cyvcf2 import VCF


INPUT_VCF_DIR = Path("<path/to/per_contig_ncds_vcf_dir>")
OUTPUT_DIR = Path("<path/to/ncds_locus_vcf_output_dir>")
EXON_BED = Path("<path/to/exon_coordinates.bed>")

POSITION_SPACING_BETWEEN_LOCI = 2000
POSITION_SPAN_PER_LOCUS_LIMIT = 1000
MIN_RECORDS_PER_LOCUS = 100
EXCLUDE_EXON_OVERLAPS = True
CONTIG_FIELD_IN_FILENAME = 2


def contig_from_filename(path):
    parts = path.name.removesuffix(".vcf.gz").split("_")
    return parts[CONTIG_FIELD_IN_FILENAME - 1]


def load_exons(bed_path):
    exons = {}
    with bed_path.open() as handle:
        for line in handle:
            if not line.strip() or line.startswith("#"):
                continue
            chrom, start, end = line.split()[:3]
            exons.setdefault(chrom, []).append((int(start), int(end)))

    for chrom in exons:
        exons[chrom].sort()
    return exons


def overlaps_exon(contig, start, end, exon_map):
    return any(start < exon_end and exon_start < end for exon_start, exon_end in exon_map.get(contig, []))


def compress_and_index(vcf_path):
    subprocess.run(["bgzip", "-f", str(vcf_path)], check=True)
    subprocess.run(["bcftools", "index", "-f", "-t", str(vcf_path) + ".gz"], check=True)


def write_locus(records, header, out_dir, contig, locus_id, exon_map):
    if len(records) < MIN_RECORDS_PER_LOCUS:
        return False

    start = records[0].POS - 1
    end = records[-1].POS
    if EXCLUDE_EXON_OVERLAPS and overlaps_exon(contig, start, end, exon_map):
        return True

    span = end - start
    out_vcf = out_dir / f"{contig}_locus{locus_id}_{span}_{len(records)}.vcf"
    with out_vcf.open("w") as out:
        out.write(header)
        for record in records:
            out.write(str(record))
    compress_and_index(out_vcf)
    return True


def process_vcf(vcf_path, exon_map):
    contig = contig_from_filename(vcf_path)
    out_dir = OUTPUT_DIR / contig
    out_dir.mkdir(parents=True, exist_ok=True)

    vcf = VCF(str(vcf_path))
    header = vcf.raw_header
    locus_id = 1
    locus = []
    locus_start = None
    last_locus_end = -10**18

    for record in vcf:
        pos = record.POS
        if locus_start is None:
            if pos >= last_locus_end + POSITION_SPACING_BETWEEN_LOCI:
                locus = [record]
                locus_start = pos
            continue

        if pos - locus_start + 1 <= POSITION_SPAN_PER_LOCUS_LIMIT:
            locus.append(record)
            continue

        if write_locus(locus, header, out_dir, contig, locus_id, exon_map):
            last_locus_end = locus[-1].POS
            locus_id += 1

        if pos >= last_locus_end + POSITION_SPACING_BETWEEN_LOCI:
            locus = [record]
            locus_start = pos
        else:
            locus = []
            locus_start = None

    if locus:
        write_locus(locus, header, out_dir, contig, locus_id, exon_map)


OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
exon_map = load_exons(EXON_BED)

for vcf_path in sorted(INPUT_VCF_DIR.glob("*_ncds.vcf.gz")):
    process_vcf(vcf_path, exon_map)
