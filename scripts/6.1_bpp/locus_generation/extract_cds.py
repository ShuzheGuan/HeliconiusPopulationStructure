#!/usr/bin/env python3

# Purpose: split CDS VCFs into exon-based BPP locus VCFs.
# Input: per-contig CDS VCF files and an exon BED file.
# Output: one directory of locus-level VCFs per contig.
# Software: Python, cyvcf2, bgzip, bcftools

from pathlib import Path
import subprocess
from cyvcf2 import VCF


INPUT_VCF_DIR = Path("<path/to/per_contig_cds_vcf_dir>")
OUTPUT_DIR = Path("<path/to/cds_locus_vcf_output_dir>")
EXON_BED = Path("<path/to/exon_coordinates.bed>")

POSITION_SPACING_BETWEEN_LOCI = 2000
MIN_RECORDS_PER_LOCUS = 100
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


def compress_and_index(vcf_path):
    subprocess.run(["bgzip", "-f", str(vcf_path)], check=True)
    subprocess.run(["bcftools", "index", "-f", "-t", str(vcf_path) + ".gz"], check=True)


def write_locus(records, header, out_dir, contig, locus_id, span):
    out_vcf = out_dir / f"{contig}_locus{locus_id}_{span}_{len(records)}.vcf"
    with out_vcf.open("w") as out:
        out.write(header)
        for record in records:
            out.write(str(record))
    compress_and_index(out_vcf)


def process_vcf(vcf_path, exon_map):
    contig = contig_from_filename(vcf_path)
    out_dir = OUTPUT_DIR / contig
    out_dir.mkdir(parents=True, exist_ok=True)

    vcf = VCF(str(vcf_path))
    header = vcf.raw_header
    locus_id = 1
    last_end_kept = -10**18

    for start, end in exon_map.get(contig, []):
        if start < last_end_kept + POSITION_SPACING_BETWEEN_LOCI:
            continue

        records = [record for record in vcf(f"{contig}:{start + 1}-{end}")]
        if len(records) < MIN_RECORDS_PER_LOCUS:
            continue

        write_locus(records, header, out_dir, contig, locus_id, end - start)
        last_end_kept = end
        locus_id += 1


OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
exon_map = load_exons(EXON_BED)

for vcf_path in sorted(INPUT_VCF_DIR.glob("*_cds.vcf.gz")):
    process_vcf(vcf_path, exon_map)
