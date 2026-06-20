# cran-comments for GLOWpipeline 0.1.0

## Test environment
- R 4.3.3 on x86_64-conda-linux-gnu
- check command: `R CMD check --as-cran --no-manual`
- NOTE: no LaTeX (pdflatex) was available in this environment, so the PDF manual was not built (`--no-manual`). The Rd files themselves pass all Rd checks; the PDF manual builds on a LaTeX-equipped machine (e.g. CRAN).

## R CMD check results
- 0 ERROR | 0 WARNING | 2 NOTE

- No ERRORs.

- No WARNINGs.

### NOTEs (verbatim from the check log)

```
* checking CRAN incoming feasibility ... NOTE
Maintainer: ‘Zheyang Wu <zheyangwu@wpi.edu>’

Unknown, possibly misspelled, fields in DESCRIPTION:
  ‘Remotes’
```

```
* checking for future file timestamps ... NOTE
unable to verify current time
```
