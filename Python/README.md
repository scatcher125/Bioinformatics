# QC Interpreter Agent

A LangGraph-based agentic workflow that parses sequencing QC metrics and uses an LLM to generate plain-English summaries for non-computational research collaborators. This tool is built to support Illumina 30x coverage Whole Exome Sequencing (WES) data alignment assesed using Alfred - an efficient and versatile BAM alignment QC tool. As input, it takes a `.json` file ().


## Motivation

In genomics research, QC reports are full of technical metrics (mapping rates, duplication rates, coverage depth) that are essential for bioinformaticians but opaque to bench scientists and clinicians. This agent bridges that gap by automatically interpreting QC outputs and producing clear, actionable summaries — keeping scientific discussions focused on biology rather than the pipeline.

## Architecture

Three LangGraph nodes in sequence:

```
[parse_qc] → [llm_summary] → [format_report]
```

1. **parse_qc_metrics** — Loads QC metrics from JSON, evaluates against established thresholds, and flags any issues
2. **generate_llm_summary** — Sends structured metrics to GPT-4o with a bioinformatics-aware prompt to generate plain-English interpretation.
3. **format_report** — Assembles a clean, human-readable report combining raw metrics and the LLM summary.

## Setup

```bash
pip install -r requirements.txt
export OPENAI_API_KEY="your-api-key-here"
```

If user does not already have an OPENAI API key, they can register for one by logging into `platform.openai.com` and generating a new secret API key.

## Usage

```bash
# Basic usage
python qc_agent.py --input INPUT_Alfred_JSON_files/case1_clean_pass.json

# Save report to file
python qc_agent.py --input INPUT_Alfred_JSON_files/case1_clean_pass.json --output case1_report.txt
```

## Input Format

The agent expects a JSON file with the following fields (all optional except `sample_id`):

```json
{
    "sample_id": "SAMPLE_003_low_ontarget",
    "Mapped": 82000000,
    "MappedFraction": 0.955,
    "DuplicateMarked": 9840000,
    "DuplicateFraction": 0.120,
    "MedianMAPQ": 60,
    "FractionInBed": 0.621,
    "MedianCoverage": 44.7,
    "SDCoverage": 41.2,
    "EnrichmentOverBed": 1.4,
    "MedianInsertSize": 178,
    "SDInsertSize": 55.1,
    "GCContent": 0.409
}
```

## QC Thresholds

| Metric | WARN | FAIL |
|--------|------|------|
| Mapped Fraction | <75% | <50% |
| Duplicate Fraction | >35% | >50% |
| Median MAPQ | <30 | <20 |
| Fraction in BED | >70% | >35% |
| Enrichment over BED | >2 | >0.9 |
| Median Coverage | 15 | 5 |
| GC Content | <35% or >60% | <20% or >80% |

## Example Output

```
================================================================================
  QC REPORT — SAMPLE_003_low_ontarget
  Overall Status: ⚠️ WARN
================================================================================

"sample_id": "SAMPLE_003_low_ontarget",
    "Mapped": 82000000,
    "MappedFraction": 0.955,
    "DuplicateMarked": 9840000,
    "DuplicateFraction": 0.120,
    "MedianMAPQ": 56,
    "FractionInBed": 0.621,
    "MedianCoverage": 44.7,
    "SDCoverage": 41.2,
    "EnrichmentOverBed": 1.4,
    "MedianInsertSize": 178,
    "SDInsertSize": 55.1,
    "GCContent": 0.409

METRICS SUMMARY
---------------
  Total Reads:         82,000,000
  Mapped Fraction:     95.5%
  Duplicate Marked:    9,840,000
  Duplicate Fraction:  12%
  Median MAPQ:         60
  Fraction In Bed:     62.1%
  Median Coverage:     44.7x
  SD Coverage:         41.2
  Enrichment Over Bed: 1.4
  Median Insert Size:  178
  SD Insert Size:      55.1
  GC Content:          40.9%

QUALITY FLAGS
-------------
  • WARN: Fraction In Bed 62.1% is below recommended (threshold: >70%)
  • WARN: EnrichmentOverBed 1.4 is below recommended (threshold: >2)

BIOLOGICAL INTERPRETATION
-----------------------------------
The reads in this sample possibly do not optimally capture exon regions.
================================================================================
```

## Extensions (future work)

- Support Whole genome sequencing (WGS)
- Batch processing across multiple samples
- Slack/email notification integration
- Interactive HTML report output
- RAG over internal QC history to contextualize current sample against cohort
