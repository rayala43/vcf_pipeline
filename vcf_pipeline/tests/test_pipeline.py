"""
tests/test_pipeline.py
-----------------------
Unit tests for the VCF Variant Annotation Pipeline.
Run with: python -m pytest tests/ -v
"""

import sys
import json
import tempfile
import textwrap
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent.parent))

from src.vcf_parser       import parse_vcf, Variant, _zygosity
from src.condition_filter import filter_by_condition, group_by_condition, resolve_aliases, list_available_conditions
from src.report_generator import generate_report


# ── Helpers ──────────────────────────────────────────────────────────────────

MINIMAL_VCF = textwrap.dedent("""\
    ##fileformat=VCFv4.2
    ##INFO=<ID=GENE,Number=1,Type=String,Description="Gene">
    ##INFO=<ID=AF,Number=A,Type=Float,Description="AF">
    ##INFO=<ID=DP,Number=1,Type=Integer,Description="Depth">
    ##INFO=<ID=CLNSIG,Number=.,Type=String,Description="ClinSig">
    ##INFO=<ID=CLNDN,Number=.,Type=String,Description="Disease">
    ##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
    #CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE
    2\t227099598\trs1801282\tC\tG\t120\tPASS\tGENE=PPARG;AF=0.12;DP=45;CLNSIG=Likely_pathogenic;CLNDN=Type_2_Diabetes\tGT\t0/1
    10\t114758349\trs10740055\tA\tG\t110\tPASS\tGENE=TCF7L2;AF=0.22;DP=52;CLNSIG=Pathogenic;CLNDN=Type_2_Diabetes\tGT\t1/1
    1\t55505647\trs562556\tA\tG\t140\tPASS\tGENE=PCSK9;AF=0.23;DP=60;CLNSIG=Pathogenic;CLNDN=Coronary_artery_disease\tGT\t0/1
    9\t107545839\trs3814960\tG\tA\t85\tPASS\tGENE=ABO;AF=0.21;DP=40;CLNSIG=Benign;CLNDN=Blood_type\tGT\t0/1
    6\t20679709\trs7754840\tC\tG\t95\tPASS\tGENE=CDKAL1;AF=0.31;DP=38;CLNSIG=risk_factor;CLNDN=Type_2_Diabetes\tGT\t0/0
""")


def _write_vcf(content: str) -> Path:
    """Write content to a temp VCF file and return its path."""
    tmp = tempfile.NamedTemporaryFile(suffix=".vcf", mode="w", delete=False)
    tmp.write(content)
    tmp.close()
    return Path(tmp.name)


# ── Parser tests ─────────────────────────────────────────────────────────────

class TestVcfParser:

    def test_parses_correct_count(self):
        path = _write_vcf(MINIMAL_VCF)
        variants = parse_vcf(path)
        assert len(variants) == 5

    def test_rsid_and_gene_extracted(self):
        path = _write_vcf(MINIMAL_VCF)
        v = parse_vcf(path)[0]
        assert v.rsid == "rs1801282"
        assert v.gene == "PPARG"

    def test_position_is_int(self):
        path = _write_vcf(MINIMAL_VCF)
        v = parse_vcf(path)[0]
        assert isinstance(v.pos, int)
        assert v.pos == 227099598

    def test_allele_freq_parsed(self):
        path = _write_vcf(MINIMAL_VCF)
        v = parse_vcf(path)[0]
        assert abs(v.allele_freq - 0.12) < 1e-6

    def test_zygosity_heterozygous(self):
        assert _zygosity("0/1") == "heterozygous"

    def test_zygosity_hom_alt(self):
        assert _zygosity("1/1") == "homozygous_alt"

    def test_zygosity_hom_ref(self):
        assert _zygosity("0/0") == "homozygous_ref"

    def test_is_variant_true_for_het(self):
        path = _write_vcf(MINIMAL_VCF)
        variants = {v.rsid: v for v in parse_vcf(path)}
        assert variants["rs1801282"].is_variant is True

    def test_is_variant_false_for_hom_ref(self):
        path = _write_vcf(MINIMAL_VCF)
        variants = {v.rsid: v for v in parse_vcf(path)}
        assert variants["rs7754840"].is_variant is False

    def test_risk_tier_pathogenic(self):
        path = _write_vcf(MINIMAL_VCF)
        variants = {v.rsid: v for v in parse_vcf(path)}
        assert variants["rs10740055"].risk_tier == "High"

    def test_risk_tier_likely_pathogenic(self):
        path = _write_vcf(MINIMAL_VCF)
        variants = {v.rsid: v for v in parse_vcf(path)}
        assert variants["rs1801282"].risk_tier == "Moderate"

    def test_risk_tier_benign(self):
        path = _write_vcf(MINIMAL_VCF)
        variants = {v.rsid: v for v in parse_vcf(path)}
        assert variants["rs3814960"].risk_tier == "Benign / VUS"

    def test_conditions_split(self):
        """CLNDN with slash splits into multiple conditions."""
        vcf = MINIMAL_VCF.replace(
            "CLNDN=Coronary_artery_disease",
            "CLNDN=Coronary_artery_disease/Stroke"
        )
        path = _write_vcf(vcf)
        variants = {v.rsid: v for v in parse_vcf(path)}
        assert len(variants["rs562556"].conditions) == 2

    def test_missing_vcf_raises(self):
        import pytest
        with pytest.raises(FileNotFoundError):
            parse_vcf("/nonexistent/path/file.vcf")


# ── Filter tests ─────────────────────────────────────────────────────────────

class TestConditionFilter:

    def _variants(self):
        return parse_vcf(_write_vcf(MINIMAL_VCF))

    def test_filter_diabetes_returns_correct_count(self):
        result = filter_by_condition(self._variants(), ["diabetes"])
        assert len(result) == 2   # rs1801282 (het) + rs10740055 (hom-alt), hom-ref excluded

    def test_filter_cad(self):
        result = filter_by_condition(self._variants(), ["cad"])
        genes = {v.gene for v in result}
        assert "PCSK9" in genes

    def test_only_carriers_default_true(self):
        """Hom-ref variant (rs7754840) must be excluded by default."""
        result = filter_by_condition(self._variants(), ["diabetes"])
        rsids = {v.rsid for v in result}
        assert "rs7754840" not in rsids

    def test_include_non_carriers(self):
        result = filter_by_condition(
            self._variants(), ["diabetes"], only_carriers=False
        )
        rsids = {v.rsid for v in result}
        assert "rs7754840" in rsids

    def test_benign_excluded_by_default(self):
        result = filter_by_condition(self._variants(), ["blood_type"])
        assert len(result) == 0

    def test_include_benign(self):
        result = filter_by_condition(
            self._variants(), ["blood type"], include_benign=True, only_carriers=True
        )
        assert any(v.gene == "ABO" for v in result)

    def test_alias_resolution_diabetes(self):
        tokens = resolve_aliases("diabetes")
        assert "diabetes" in tokens

    def test_alias_resolution_unknown_falls_back(self):
        tokens = resolve_aliases("my_custom_condition")
        assert tokens == ["my_custom_condition"]

    def test_list_conditions_not_empty(self):
        assert len(list_available_conditions()) > 0

    def test_group_by_condition_structure(self):
        variants = filter_by_condition(self._variants(), ["diabetes", "cad"])
        groups = group_by_condition(variants, ["diabetes", "cad"])
        assert "Diabetes" in groups
        assert "Cad" in groups

    def test_empty_conditions_returns_empty(self):
        result = filter_by_condition(self._variants(), [])
        assert result == []


# ── Report generator tests ────────────────────────────────────────────────────

class TestReportGenerator:

    def test_report_creates_file(self):
        path = _write_vcf(MINIMAL_VCF)
        variants = filter_by_condition(parse_vcf(path), ["diabetes"])
        groups = group_by_condition(variants, ["diabetes"])
        with tempfile.TemporaryDirectory() as td:
            out = Path(td) / "report.html"
            result = generate_report(groups, "sample.vcf", out)
            assert result.exists()
            assert result.stat().st_size > 500

    def test_report_contains_patient_id(self):
        path = _write_vcf(MINIMAL_VCF)
        variants = filter_by_condition(parse_vcf(path), ["diabetes"])
        groups = group_by_condition(variants, ["diabetes"])
        with tempfile.TemporaryDirectory() as td:
            out = Path(td) / "report.html"
            generate_report(groups, "sample.vcf", out, patient_id="TEST_P999")
            html = out.read_text()
            assert "TEST_P999" in html

    def test_report_contains_gene(self):
        path = _write_vcf(MINIMAL_VCF)
        variants = filter_by_condition(parse_vcf(path), ["diabetes"])
        groups = group_by_condition(variants, ["diabetes"])
        with tempfile.TemporaryDirectory() as td:
            out = Path(td) / "report.html"
            generate_report(groups, "sample.vcf", out)
            html = out.read_text()
            assert "TCF7L2" in html
