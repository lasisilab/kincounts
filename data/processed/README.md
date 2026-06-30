# Processed Data

This folder is for aggregate summaries derived from IPUMS USA fertility data. It should contain only shareable processed outputs, not row-level microdata.

Expected regenerated outputs include:

- cohort sample-size summaries;
- aggregate counts of women by cohort/source and children-ever-born count;
- weighted and unweighted fertility PMFs;
- direct-1950 versus proxy-1950 validation summaries;
- simulator-ready `fertility_pmf_YYYY.csv` files, if still needed by the React app.

To rebuild from a private IPUMS extract, place the extract and any companion metadata/syntax files in `data/raw/` and run:

```bash
Rscript code/preprocess_ipums.R
```
