"""
vcf_parser.py
-------------
Parses VCF v4.x files and returns structured variant records.
Handles multi-allelic INFO fields, missing values, and malformed lines gracefully.
"""

import re
from pathlib import Path
from dataclasses import dataclass, field
from typing import Optional


@dataclass
class Variant:
    chrom: str
    pos: int
    rsid: str
    ref: str
    alt: str
    qual: Optional[float]
    filter_status: str
    gene: str
    allele_freq: Optional[float]
    depth: Optional[int]
    clinical_significance: str
    disease_name: str
    genotype: str
    zygosity: str          # homozygous_ref | heterozygous | homozygous_alt
    conditions: list       # parsed list from CLNDN (split on '/')

    @property
    def is_variant(self) -> bool:
        """True if the sample actually carries the ALT allele."""
        return self.genotype not in ("0/0", "0|0", "./.")

    @property
    def risk_tier(self) -> str:
        sig = self.clinical_significance.lower().replace(" ", "_")
        if "likely_pathogenic" in sig:
            return "Moderate"
        if "pathogenic" in sig:
            return "High"
        if "risk_factor" in sig:
            return "Low"
        return "Benign / VUS"


def _parse_info(info_str: str) -> dict:
    result = {}
    for token in info_str.split(";"):
        if "=" in token:
            k, v = token.split("=", 1)
            result[k.strip()] = v.strip()
        else:
            result[token.strip()] = True
    return result


def _parse_genotype(fmt_str: str, sample_str: str) -> str:
    keys = fmt_str.split(":")
    vals = sample_str.split(":")
    fmt_map = dict(zip(keys, vals))
    return fmt_map.get("GT", "./.")


def _zygosity(gt: str) -> str:
    alleles = re.split(r"[/|]", gt)
    if all(a == "0" for a in alleles):
        return "homozygous_ref"
    if len(set(alleles)) == 1:
        return "homozygous_alt"
    return "heterozygous"


def parse_vcf(filepath: str | Path) -> list[Variant]:
    """
    Parse a VCF file and return a list of Variant objects.

    Args:
        filepath: Path to the .vcf file.

    Returns:
        List of Variant dataclass instances.

    Raises:
        FileNotFoundError: If the VCF file does not exist.
        ValueError: If the file has no valid #CHROM header line.
    """
    path = Path(filepath)
    if not path.exists():
        raise FileNotFoundError(f"VCF file not found: {path}")

    variants = []
    header_found = False

    with open(path, "r") as fh:
        for raw_line in fh:
            line = raw_line.strip()
            if not line:
                continue
            if line.startswith("##"):
                continue
            if line.startswith("#CHROM"):
                header_found = True
                continue
            if not header_found:
                continue

            cols = line.split("\t")
            if len(cols) < 8:
                continue  # skip malformed lines

            chrom, pos_str, rsid, ref, alt, qual_str, filt = cols[:7]
            info_str = cols[7]
            fmt_str  = cols[8] if len(cols) > 8 else "GT"
            samp_str = cols[9] if len(cols) > 9 else "./."

            info = _parse_info(info_str)

            try:
                pos = int(pos_str)
            except ValueError:
                continue

            try:
                qual = float(qual_str)
            except ValueError:
                qual = None

            gene    = info.get("GENE", "UNKNOWN")
            af_raw  = info.get("AF")
            dp_raw  = info.get("DP")
            clnsig  = info.get("CLNSIG", "Unknown").replace("_", " ")
            clndn   = info.get("CLNDN", "Unknown").replace("_", " ")

            try:
                af = float(af_raw) if af_raw else None
            except ValueError:
                af = None

            try:
                dp = int(dp_raw) if dp_raw else None
            except ValueError:
                dp = None

            gt = _parse_genotype(fmt_str, samp_str)
            conditions = [c.strip() for c in clndn.split("/")]

            variants.append(Variant(
                chrom=chrom,
                pos=pos,
                rsid=rsid if rsid != "." else f"chr{chrom}:{pos}",
                ref=ref,
                alt=alt,
                qual=qual,
                filter_status=filt,
                gene=gene,
                allele_freq=af,
                depth=dp,
                clinical_significance=clnsig,
                disease_name=clndn,
                genotype=gt,
                zygosity=_zygosity(gt),
                conditions=conditions,
            ))

    if not header_found:
        raise ValueError(f"No #CHROM header found in {path}. Is this a valid VCF?")

    return variants
