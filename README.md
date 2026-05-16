# 🧬 VCF Variant Annotation Pipeline

A production-style genomics pipeline that parses VCF (Variant Call Format) files, filters variants by user-defined medical conditions at runtime, and generates **interactive HTML clinical reports** — no external bioinformatics tools required.

Built from real-world experience processing clinical genomics data at GenepoweRx (K&H Personalized Clinical Research).

---

## Demo

```
python main.py \
  --vcf data/sample_vcfs/sample_mixed.vcf \
  --conditions "diabetes" "cad" "breast cancer" \
  --output reports/report.html \
  --patient-id PATIENT_001
```

Open `reports/report.html` in any browser — no server needed.

---

## Features

- **Zero-dependency VCF parser** — handles VCF v4.x, multi-allelic INFO fields, missing values, malformed lines
- **Runtime condition filtering** — specify any condition at the command line; built-in aliases map common names to search tokens
- **Zygosity classification** — homozygous alt / heterozygous / homozygous ref
- **Risk tiering** — Pathogenic → High, Likely Pathogenic → Moderate, Risk Factor → Low, Benign/VUS
- **Interactive HTML report** with:
  - Summary stat cards (total variants, risk breakdown, gene count)
  - Risk distribution donut chart
  - Variants-per-gene bar chart
  - Collapsible condition sections with full variant tables
  - Sidebar navigation
- **28 unit tests** — parser, filter, and report generator all covered

---

## Project Structure

```
vcf_pipeline/
├── main.py                        # CLI entry point
├── src/
│   ├── vcf_parser.py              # VCF parser → Variant dataclass
│   ├── condition_filter.py        # Runtime condition filtering & grouping
│   └── report_generator.py        # Self-contained HTML report builder
├── data/
│   └── sample_vcfs/
│       ├── sample_diabetes.vcf    # Diabetes-related variants (PPARG, TCF7L2, KCNJ11...)
│       ├── sample_cad.vcf         # CAD variants (PCSK9, APOE, LPA, APOB...)
│       └── sample_mixed.vcf       # Multi-condition (BRCA1/2, TCF7L2, APOE, HBB...)
├── reports/                       # Generated HTML reports land here
├── tests/
│   └── test_pipeline.py           # 28 unit tests (pytest)
└── README.md
```

---

## Installation

No pip dependencies for the core pipeline. Only `pytest` for running tests.

```bash
git clone https://github.com/YOUR_USERNAME/vcf-variant-annotation-pipeline
cd vcf-variant-annotation-pipeline
pip install pytest          # for tests only
```

Python 3.10+ required (uses `match` type annotations).

---

## Usage

### Basic

```bash
python main.py --vcf data/sample_vcfs/sample_diabetes.vcf --conditions "diabetes"
```

### Multiple conditions

```bash
python main.py \
  --vcf data/sample_vcfs/sample_mixed.vcf \
  --conditions "diabetes" "cad" "breast cancer" \
  --patient-id PATIENT_042 \
  --output reports/patient_042.html
```

### Include benign and non-carrier records

```bash
python main.py \
  --vcf data/sample_vcfs/sample_mixed.vcf \
  --conditions "diabetes" \
  --include-benign \
  --include-non-carriers
```

### List all built-in condition aliases

```bash
python main.py --list-conditions
```

**Built-in aliases:**
| Input | Matches |
|---|---|
| `diabetes` | Type 2 Diabetes, T2D, Maturity Onset Diabetes |
| `cad` | Coronary Artery Disease, Familial Hypercholesterolemia |
| `breast cancer` | Breast Cancer, BRCA |
| `alzheimer` | Alzheimer Disease |
| `sickle cell` | Sickle Cell Disease |
| `cystic fibrosis` | Cystic Fibrosis, CFTR |

Any free-text condition string also works — the pipeline searches for partial matches against `CLNDN` fields.

---

## Input VCF Format

Standard VCF v4.x. The pipeline uses these INFO fields:

| Field | Description |
|---|---|
| `GENE` | Gene symbol |
| `AF` | Allele frequency |
| `DP` | Read depth |
| `CLNSIG` | ClinVar clinical significance |
| `CLNDN` | ClinVar disease name (slash-separated for multiple) |

If `CLNDN` is absent, the variant is retained but marked Unknown condition.

---

## Running Tests

```bash
python -m pytest tests/ -v
```

All 28 tests cover:
- VCF parsing (correct counts, rsID/gene extraction, zygosity, risk tiers)
- Condition filtering (carrier/non-carrier, benign inclusion, alias resolution)
- Report generation (file creation, patient ID embedding, gene presence)

---

## Technical Stack

| Layer | Technology |
|---|---|
| Language | Python 3.10+ |
| Core parsing | stdlib only (`re`, `pathlib`, `dataclasses`) |
| Charts | Chart.js 4.x (CDN, in report) |
| Tests | pytest |
| Output | Self-contained HTML (no server required) |

---

## Roadmap / Planned Features

- [ ] ClinVar API integration — live lookup for rsID significance
- [ ] Multi-sample VCF support (family trios)
- [ ] HGVS notation parsing
- [ ] CSV/JSON export alongside HTML
- [ ] Streamlit web UI for drag-and-drop VCF upload
- [ ] gnomAD allele frequency annotation

---

## Background

Developed based on real clinical genomics workflows — extracting disease-specific gene/rsID combinations from patient VCF files, segregating by condition, and building automated outputs for clinical reporting. The sample VCFs use publicly known variant rsIDs from ClinVar and dbSNP.

---
