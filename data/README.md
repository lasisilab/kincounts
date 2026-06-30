# Data

This folder separates restricted IPUMS source files from shareable derived files.

- `raw/`: local-only IPUMS extracts. Do not commit raw IPUMS microdata.
- `processed/`: aggregate fertility summaries used by the analyses and simulator. This folder is regenerated from the private raw extract and should contain no row-level IPUMS records.

IPUMS USA's terms of use restrict redistribution of IPUMS-USA data without permission. The raw extract needed to rebuild the processed summaries must therefore be downloaded by a registered IPUMS user and placed in `data/raw/`. Only aggregate summaries should be committed.
