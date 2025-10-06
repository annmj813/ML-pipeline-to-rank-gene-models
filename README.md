# ML-pipeline-to-rank-gene-models
Rank gene models vs official annotation with a logistic regression pipeline (R): overlap filtering, z-score scaling, OOF tuning, per-locus ranking &amp; metrics.
---

## Features
- Z-score scaling (fit on train only; no imputation)
- 10× stratified 70/30 splits with OOF predictions
- Threshold tuning by **MCC** (and optional AUC)
- **Per-locus ranking** metrics (TP/FP/FN/TN)
- Simple config file; runs with one command

---

## Quick start

### 1) Install packages
```bash
Rscript install.R
