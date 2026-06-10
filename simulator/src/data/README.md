# Empirical fertility PMF data

Drop the five CSV files produced by the IPUMS processing script here:

  fertility_pmf_1950.csv
  fertility_pmf_1960.csv
  fertility_pmf_1970.csv
  fertility_pmf_1980.csv
  fertility_pmf_1990.csv

Expected format (two columns, produced by pandas value_counts):

  children,probability
  0,0.123
  1,0.234
  ...
  15,0.001

Once the files are present, rebuild the app (`npm run dev` / `npm run build`)
and the butterfly chart will switch from ZINB theoretical to empirical PMFs.
