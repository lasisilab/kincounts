# Raw Data

Place private IPUMS USA extract files here, for example:

```text
data/raw/usa_00009.dat.gz
data/raw/usa_00009.xml
```

Keep the DDI/XML codebook alongside the `.dat.gz` extract — `analysis/preprocess.qmd` reads the extract via the codebook (it matches whichever extract is present, so the exact number does not matter). The raw files are ignored by Git.

The current preprocessing page (`analysis/preprocess.qmd`) requires:

- `YEAR`
- `PERWT`
- `SEX`
- `AGE`
- `BIRTHYR`
- `MARST`
- `GQ`
- `CHBORN`

Design/default variables such as `SAMPLE`, `SERIAL`, `HHWT`, `CLUSTER`, `STRATA`, and `PERNUM` may also be present, but they are not required by the current preprocessing code.

Raw IPUMS extracts are intentionally ignored by Git. To regenerate the shareable aggregate files after placing the extract here, render the analysis site:

```bash
quarto render
```
