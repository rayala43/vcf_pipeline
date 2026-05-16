"""
report_generator.py
--------------------
Generates a self-contained, interactive HTML clinical report from filtered variants.
No external dependencies — pure Python string templating.
"""

from __future__ import annotations
import json
from datetime import datetime
from pathlib import Path
from src.vcf_parser import Variant


RISK_COLOR = {
    "High":         ("#dc2626", "#fef2f2"),   # red  text / bg
    "Moderate":     ("#d97706", "#fffbeb"),   # amber
    "Low":          ("#2563eb", "#eff6ff"),   # blue
    "Benign / VUS": ("#6b7280", "#f9fafb"),   # gray
}

ZYGOSITY_LABEL = {
    "homozygous_alt": "Hom Alt",
    "heterozygous":   "Het",
    "homozygous_ref": "Hom Ref",
}


def _badge(text: str, color: str, bg: str) -> str:
    return (
        f'<span style="display:inline-block;padding:2px 8px;border-radius:99px;'
        f'font-size:11px;font-weight:600;color:{color};background:{bg}">{text}</span>'
    )


def _risk_badge(tier: str) -> str:
    c, bg = RISK_COLOR.get(tier, ("#6b7280", "#f9fafb"))
    return _badge(tier, c, bg)


def _variant_row(v: Variant) -> str:
    zyg_label = ZYGOSITY_LABEL.get(v.zygosity, v.zygosity)
    af_str    = f"{v.allele_freq:.3f}" if v.allele_freq is not None else "—"
    dp_str    = str(v.depth) if v.depth is not None else "—"
    return f"""
        <tr>
          <td><code style="font-size:12px">{v.rsid}</code></td>
          <td><strong>{v.gene}</strong></td>
          <td>chr{v.chrom}:{v.pos:,}</td>
          <td>{v.ref} → {v.alt}</td>
          <td>{zyg_label}</td>
          <td>{af_str}</td>
          <td>{dp_str}</td>
          <td>{_risk_badge(v.risk_tier)}</td>
          <td style="max-width:200px;font-size:12px">{v.clinical_significance}</td>
        </tr>"""


def _section(condition: str, variants: list[Variant], idx: int) -> str:
    if not variants:
        return ""
    rows = "".join(_variant_row(v) for v in variants)
    counts = {"High": 0, "Moderate": 0, "Low": 0}
    for v in variants:
        if v.risk_tier in counts:
            counts[v.risk_tier] += 1

    summary_pills = (
        (_badge(f"{counts['High']} High", *RISK_COLOR["High"]) if counts["High"] else "") + " " +
        (_badge(f"{counts['Moderate']} Moderate", *RISK_COLOR["Moderate"]) if counts["Moderate"] else "") + " " +
        (_badge(f"{counts['Low']} Low risk", *RISK_COLOR["Low"]) if counts["Low"] else "")
    )

    return f"""
    <div class="section" id="sec-{idx}">
      <div class="section-header" onclick="toggleSection({idx})">
        <span class="section-title">🧬 {condition}</span>
        <span class="section-meta">{len(variants)} variant(s) &nbsp;{summary_pills}</span>
        <span class="chevron" id="chev-{idx}">▾</span>
      </div>
      <div class="section-body" id="body-{idx}">
        <div style="overflow-x:auto">
          <table>
            <thead>
              <tr>
                <th>rsID</th><th>Gene</th><th>Location</th><th>Allele</th>
                <th>Zygosity</th><th>AF</th><th>Depth</th><th>Risk</th><th>Clinical Significance</th>
              </tr>
            </thead>
            <tbody>{rows}</tbody>
          </table>
        </div>
      </div>
    </div>"""


def generate_report(
    groups: dict[str, list[Variant]],
    vcf_filename: str,
    output_path: str | Path,
    patient_id: str = "PATIENT_001",
) -> Path:
    """
    Build a self-contained HTML clinical report.

    Args:
        groups:        Output of condition_filter.group_by_condition()
        vcf_filename:  Original VCF file name (for display)
        output_path:   Where to write the .html file
        patient_id:    Patient identifier string

    Returns:
        Path to the generated HTML file.
    """
    timestamp   = datetime.now().strftime("%Y-%m-%d %H:%M")
    total_vars  = sum(len(v) for v in groups.values())
    conditions  = list(groups.keys())

    sections_html = "\n".join(
        _section(cond, variants, i)
        for i, (cond, variants) in enumerate(groups.items())
    )

    nav_items = "\n".join(
        f'<li><a href="#sec-{i}" onclick="openSection({i})">{cond}</a></li>'
        for i, cond in enumerate(conditions)
    )

    # ── Chart data for the summary donut ─────────────────────────────────────
    risk_counts = {"High": 0, "Moderate": 0, "Low": 0}
    gene_counts: dict[str, int] = {}
    for variants in groups.values():
        for v in variants:
            if v.risk_tier in risk_counts:
                risk_counts[v.risk_tier] += 1
            gene_counts[v.gene] = gene_counts.get(v.gene, 0) + 1

    top_genes    = sorted(gene_counts.items(), key=lambda x: -x[1])[:8]
    gene_labels  = json.dumps([g for g, _ in top_genes])
    gene_data    = json.dumps([c for _, c in top_genes])
    donut_data   = json.dumps(list(risk_counts.values()))

    html = f"""<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="UTF-8">
<meta name="viewport" content="width=device-width, initial-scale=1.0">
<title>VCF Variant Annotation Report — {patient_id}</title>
<script src="https://cdnjs.cloudflare.com/ajax/libs/Chart.js/4.4.1/chart.umd.min.js"></script>
<style>
  *{{box-sizing:border-box;margin:0;padding:0}}
  body{{font-family:'Segoe UI',system-ui,sans-serif;background:#f1f5f9;color:#1e293b;line-height:1.6}}
  a{{color:#2563eb;text-decoration:none}}

  /* ── Header ── */
  .header{{background:linear-gradient(135deg,#1e3a5f 0%,#1a5276 100%);color:#fff;padding:28px 40px}}
  .header h1{{font-size:22px;font-weight:700;letter-spacing:.3px}}
  .header .sub{{font-size:13px;opacity:.8;margin-top:4px}}
  .meta-chips{{display:flex;gap:12px;margin-top:16px;flex-wrap:wrap}}
  .chip{{background:rgba(255,255,255,.15);border-radius:99px;padding:4px 14px;font-size:12px}}

  /* ── Layout ── */
  .container{{max-width:1200px;margin:0 auto;padding:24px 20px;display:grid;grid-template-columns:220px 1fr;gap:24px}}
  @media(max-width:768px){{.container{{grid-template-columns:1fr}}}}

  /* ── Sidebar nav ── */
  .sidebar{{position:sticky;top:20px;height:fit-content}}
  .nav-card{{background:#fff;border-radius:12px;padding:16px;box-shadow:0 1px 4px rgba(0,0,0,.07)}}
  .nav-card h3{{font-size:12px;text-transform:uppercase;letter-spacing:.8px;color:#94a3b8;margin-bottom:10px}}
  .nav-card ul{{list-style:none}}
  .nav-card li a{{display:block;padding:6px 8px;border-radius:6px;font-size:13px;color:#475569;transition:background .15s}}
  .nav-card li a:hover{{background:#eff6ff;color:#2563eb}}

  /* ── Summary cards ── */
  .summary-grid{{display:grid;grid-template-columns:repeat(auto-fit,minmax(130px,1fr));gap:12px;margin-bottom:20px}}
  .stat-card{{background:#fff;border-radius:10px;padding:14px 16px;box-shadow:0 1px 4px rgba(0,0,0,.07)}}
  .stat-card .num{{font-size:28px;font-weight:700}}
  .stat-card .lbl{{font-size:11px;color:#94a3b8;text-transform:uppercase;letter-spacing:.5px}}
  .stat-card.red .num{{color:#dc2626}}
  .stat-card.amber .num{{color:#d97706}}
  .stat-card.blue .num{{color:#2563eb}}

  /* ── Charts row ── */
  .charts-row{{display:grid;grid-template-columns:1fr 1fr;gap:16px;margin-bottom:20px}}
  @media(max-width:640px){{.charts-row{{grid-template-columns:1fr}}}}
  .chart-card{{background:#fff;border-radius:10px;padding:16px;box-shadow:0 1px 4px rgba(0,0,0,.07)}}
  .chart-card h4{{font-size:12px;color:#94a3b8;text-transform:uppercase;letter-spacing:.5px;margin-bottom:12px}}

  /* ── Variant sections ── */
  .section{{background:#fff;border-radius:10px;margin-bottom:14px;box-shadow:0 1px 4px rgba(0,0,0,.07);overflow:hidden}}
  .section-header{{padding:14px 18px;cursor:pointer;display:flex;align-items:center;gap:10px;user-select:none;transition:background .15s}}
  .section-header:hover{{background:#f8fafc}}
  .section-title{{font-weight:600;font-size:15px;flex:1}}
  .section-meta{{font-size:12px;color:#94a3b8;display:flex;gap:6px;align-items:center;flex-wrap:wrap}}
  .chevron{{font-size:18px;color:#94a3b8;transition:transform .2s}}
  .chevron.open{{transform:rotate(180deg)}}
  .section-body{{display:none;padding:0 0 16px}}
  .section-body.open{{display:block}}

  /* ── Table ── */
  table{{width:100%;border-collapse:collapse;font-size:13px}}
  thead tr{{background:#f8fafc}}
  th{{padding:8px 12px;text-align:left;font-size:11px;text-transform:uppercase;letter-spacing:.5px;color:#64748b;border-bottom:1px solid #e2e8f0}}
  td{{padding:8px 12px;border-bottom:1px solid #f1f5f9;vertical-align:middle}}
  tr:last-child td{{border-bottom:none}}
  tr:hover td{{background:#fafbfc}}
  code{{background:#f1f5f9;padding:2px 5px;border-radius:4px;font-size:11px}}

  /* ── Footer ── */
  .footer{{text-align:center;font-size:12px;color:#94a3b8;padding:20px;margin-top:10px}}
</style>
</head>
<body>

<div class="header">
  <h1>🧬 VCF Variant Annotation Report</h1>
  <div class="sub">Genomic variant analysis with condition-specific filtering</div>
  <div class="meta-chips">
    <span class="chip">👤 {patient_id}</span>
    <span class="chip">📄 {vcf_filename}</span>
    <span class="chip">📅 {timestamp}</span>
    <span class="chip">🔬 {total_vars} variant(s) found</span>
  </div>
</div>

<div class="container">

  <!-- Sidebar -->
  <aside class="sidebar">
    <div class="nav-card">
      <h3>Conditions</h3>
      <ul>{nav_items}</ul>
    </div>
  </aside>

  <!-- Main content -->
  <main>

    <!-- Summary stats -->
    <div class="summary-grid">
      <div class="stat-card">
        <div class="num">{total_vars}</div>
        <div class="lbl">Total Variants</div>
      </div>
      <div class="stat-card red">
        <div class="num">{risk_counts["High"]}</div>
        <div class="lbl">High Risk</div>
      </div>
      <div class="stat-card amber">
        <div class="num">{risk_counts["Moderate"]}</div>
        <div class="lbl">Moderate Risk</div>
      </div>
      <div class="stat-card blue">
        <div class="num">{risk_counts["Low"]}</div>
        <div class="lbl">Low Risk</div>
      </div>
      <div class="stat-card">
        <div class="num">{len(conditions)}</div>
        <div class="lbl">Conditions</div>
      </div>
      <div class="stat-card">
        <div class="num">{len(gene_counts)}</div>
        <div class="lbl">Genes Hit</div>
      </div>
    </div>

    <!-- Charts -->
    <div class="charts-row">
      <div class="chart-card">
        <h4>Risk distribution</h4>
        <canvas id="donutChart" height="180"></canvas>
      </div>
      <div class="chart-card">
        <h4>Variants per gene</h4>
        <canvas id="geneChart" height="180"></canvas>
      </div>
    </div>

    <!-- Variant sections -->
    {sections_html}

  </main>
</div>

<div class="footer">
  Generated by VCF Variant Annotation Pipeline &nbsp;|&nbsp; For research use only &nbsp;|&nbsp; {timestamp}
</div>

<script>
// ── Charts ────────────────────────────────────────────────────────────────
const donutCtx = document.getElementById('donutChart').getContext('2d');
new Chart(donutCtx, {{
  type: 'doughnut',
  data: {{
    labels: ['High', 'Moderate', 'Low'],
    datasets: [{{
      data: {donut_data},
      backgroundColor: ['#fca5a5', '#fcd34d', '#93c5fd'],
      borderColor: ['#dc2626', '#d97706', '#2563eb'],
      borderWidth: 2
    }}]
  }},
  options: {{
    responsive: true,
    plugins: {{ legend: {{ position: 'bottom', labels: {{ font: {{ size: 11 }} }} }} }},
    cutout: '65%'
  }}
}});

const geneCtx = document.getElementById('geneChart').getContext('2d');
new Chart(geneCtx, {{
  type: 'bar',
  data: {{
    labels: {gene_labels},
    datasets: [{{
      label: 'Variants',
      data: {gene_data},
      backgroundColor: '#93c5fd',
      borderColor: '#2563eb',
      borderWidth: 1,
      borderRadius: 4
    }}]
  }},
  options: {{
    responsive: true,
    plugins: {{ legend: {{ display: false }} }},
    scales: {{
      y: {{ beginAtZero: true, ticks: {{ stepSize: 1, font: {{ size: 11 }} }} }},
      x: {{ ticks: {{ font: {{ size: 11 }} }} }}
    }}
  }}
}});

// ── Accordion ────────────────────────────────────────────────────────────
function toggleSection(i) {{
  const body = document.getElementById('body-' + i);
  const chev = document.getElementById('chev-' + i);
  body.classList.toggle('open');
  chev.classList.toggle('open');
}}
function openSection(i) {{
  const body = document.getElementById('body-' + i);
  const chev = document.getElementById('chev-' + i);
  body.classList.add('open');
  chev.classList.add('open');
}}
// Open first section by default
document.addEventListener('DOMContentLoaded', () => {{ openSection(0); }});
</script>
</body>
</html>"""

    out = Path(output_path)
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(html, encoding="utf-8")
    return out
