# Merging TRGT VCFs

The `trgt merge` command can be used to merge TRGT VCFs into a single multi-sample
VCF.

## Input requirements

Input TRGT VCFs must be sorted, BGZF-compressed, and indexed (`.tbi` or `.csi`).
All VCFs must be genotyped with the same repeat catalog. Currently, inputs need
to be accessible on the local filesystem. If any input VCF predates TRGT 1.0 the
merge will require a reference genome file to be added with `--genome` such that
the missing padding base can be reintroduced. Sample IDs must be unique across
files, when duplicates are intentional, add the flag `--force-samples` to keep
them.

## Output formatting

Without `--output` set, results will stream to standard output in uncompressed
VCF format. Set `--output` to write to a file. The file extension (`.vcf`,
`.vcf.gz`, `.bcf`, `.bcf.gz`) will automatically determine both format and
compression level unless you override it using `--output-type`. Use
`--print-header` to inspect the merged header without merging variants, this is
useful to verify sample order or merged metadata.

## How records are merged

Records are traversed contig by contig in reference order, and you can constrain
work to specific contigs with `--contig chr20,chrX`. For each site the reference
alleles must match across all inputs; otherwise the site fails with a clear
error (unless `--quit-on-errors` is omitted, in which case the record is skipped
and a warning is logged). FORMAT fields (`AL`, `ALLR`, `SD`, `MC`, `MS`, `AP`,
`AM`) are merged sample-by-sample, with missing samples filled using the correct
sentinel values so downstream tooling sees absence explicitly. Headers are
merged once up-front; the command rewrites historic field definitions when
needed (for example, collapsing `ALCI` into the current `ALLR` representation)
and injects version metadata unless `--no-version` is set.

## Example usage

Merge three TRGT VCFs into one compressed VCF:

```bash
trgt merge \
  --vcf results/sample1.trgt.vcf.gz results/sample2.trgt.vcf.gz results/sample3.trgt.vcf.gz \
  --output cohort.trgt.vcf.gz
```

If the results folders contains these three VCF files you can also use auto-expansions:

```bash
trgt merge \
  --vcf results/*.vcf.gz \
  --output cohort.trgt.vcf.gz
```

Running the merge from a manifest file, while also limiting processing to chromosome 12:

```bash
trgt merge \
  --vcf-list cohort_vcfs.txt \
  --contig chr12 \
  --output cohort.chr12.trgt.vcf.gz
```

Here `cohort_vcfs.txt` is a plain text file with one VCF path per line. In this
file leading and trailing whitespace is ignored, blank lines are allowed, and
lines starting with `#` are treated as comments. Both relative and absolute
paths are accepted as long as they exist when the merge runs.

## Troubleshooting

If an errors is reported with "Reference alleles do not match", it indicates that the input VCFs
used different TR catalogs for genotyping.
