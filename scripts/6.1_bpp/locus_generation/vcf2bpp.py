#!/usr/bin/env python3

"""Convert locus-level VCF files to BPP sequence format."""

from collections import defaultdict
from pathlib import Path
from cyvcf2 import VCF


INPUT_VCF_ROOT = Path("<path/to/locus_vcf_root_dir>")
OUTPUT_BPP_ROOT = Path("<path/to/bpp_locus_output_dir>")


IUPAC = {
    frozenset({"A"}): "A",
    frozenset({"C"}): "C",
    frozenset({"G"}): "G",
    frozenset({"T"}): "T",
    frozenset({"A", "G"}): "R",
    frozenset({"C", "T"}): "Y",
    frozenset({"G", "C"}): "S",
    frozenset({"A", "T"}): "W",
    frozenset({"G", "T"}): "K",
    frozenset({"A", "C"}): "M",
}


def allele_base(index, ref, alts):
    if index == 0:
        return ref
    if 1 <= index <= len(alts):
        return alts[index - 1]
    return None


def genotype_base(genotype, ref, alts):
    alleles = list(genotype)
    if alleles and isinstance(alleles[-1], bool):
        alleles = alleles[:-1]

    called = [int(a) for a in alleles[:2] if a is not None and a >= 0]
    if not called:
        return "N" * len(ref) if len(ref) > 1 and not alts else "N"

    bases = [allele_base(a, ref, alts) for a in called]
    bases = [b for b in bases if b is not None]
    if not bases:
        return "N"
    if len(set(bases)) == 1:
        return bases[0]
    if all(len(b) == 1 and b in "ACGT" for b in bases):
        return IUPAC.get(frozenset(bases), "N")
    return "N"


def convert_vcf(vcf_path, out_path):
    vcf = VCF(str(vcf_path))
    samples = vcf.samples
    seqs = defaultdict(list)

    for variant in vcf:
        alts = [a for a in (variant.ALT or []) if a and a != "."]
        for sample, genotype in zip(samples, variant.genotypes):
            seqs[sample].append(genotype_base(genotype, variant.REF, alts))

    if not samples or not seqs:
        return

    length = len("".join(seqs[samples[0]]))
    out_path.parent.mkdir(parents=True, exist_ok=True)

    with out_path.open("w") as out:
        out.write(f"{len(samples)} {length}\n")
        for sample in samples:
            out.write(f"^{sample}  {''.join(seqs[sample])}\n")


for vcf_path in sorted(INPUT_VCF_ROOT.rglob("*.vcf.gz")):
    rel = vcf_path.relative_to(INPUT_VCF_ROOT)
    out_path = OUTPUT_BPP_ROOT / rel.with_suffix("").with_suffix(".bpp.txt")
    convert_vcf(vcf_path, out_path)
