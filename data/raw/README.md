# Raw Data

Place private IPUMS USA extract files here, for example:

```text
data/raw/usa_00008.dat.gz
```

If IPUMS provides a companion metadata, DDI/XML, or syntax file, keep it here too. The raw files are ignored by Git.

The preprocessing script expects variables used by the paper and validation checks:

- `YEAR`
- `SAMPLE`
- `AGE`
- `MARST`
- `BIRTHYR`
- `CHBORN`
- `PERWT`
- `SEX` if selected

Design/default variables such as `SERIAL`, `HHWT`, `CLUSTER`, `STRATA`, `GQ`, and `PERNUM` may also be present.

Raw IPUMS extracts are intentionally ignored by Git. To regenerate the shareable aggregate files after placing the extract here, run:

```bash
Rscript code/preprocess_ipums.R
```
