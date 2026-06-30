# Raw Data

Place private IPUMS USA extract files here, for example:

```text
data/raw/usa_00009.dat.gz
data/raw/usa_00009.xml
```

Keep the DDI/XML codebook alongside the `.dat.gz` extract — `analysis/preprocess.qmd` reads the extract via the codebook (it matches whichever extract is present, so the exact number does not matter). The raw files are ignored by Git.

The preprocessing expects variables used by the paper and validation checks:

- `YEAR`
- `SAMPLE`
- `AGE`
- `MARST`
- `BIRTHYR`
- `CHBORN`
- `PERWT`
- `SEX` if selected

Design/default variables such as `SERIAL`, `HHWT`, `CLUSTER`, `STRATA`, `GQ`, and `PERNUM` may also be present.

Raw IPUMS extracts are intentionally ignored by Git. To regenerate the shareable aggregate files after placing the extract here, render the analysis site:

```bash
quarto render
```
