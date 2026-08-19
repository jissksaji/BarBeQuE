# Primer Input

There are three ways to tell BarBeQuE which primers to benchmark. `--input` and `--primer_set` are
mutually exclusive; one of them is required unless you are running `--build_references`.

| Route | Parameter | Amplicon bounds come from |
| --- | --- | --- |
| Named primer sets | `--primer_set COI,Fish16S` | the FooDMe2 catalog, per set |
| Samplesheet | `--input primers.tsv` | the `min`/`max` columns, per row |
| Primer FASTA | `--input primers.fasta` or `--input primer_dir/` | `--primer_min` / `--primer_max`, globally |

All three end up as the same `primer / fwd / rev / min / max` rows, so nothing downstream of
`INPUT_CHECK` knows or cares which route you used.

## How `--input` is recognised

A directory is always treated as primer FASTAs. A file is read as FASTA or as a samplesheet based on
**its first non-blank line**, not its extension: a line starting with `>` means FASTA, anything else
means samplesheet. This is `WorkflowPipeline.isFastaInput` in `lib/WorkflowPipeline.groovy`, and it
is the same check used at startup to decide whether `--primer_min`/`--primer_max` are required.

Naming a samplesheet `primers.fasta`, or a FASTA `primers.txt`, therefore still does the right
thing.

## What happens to a FASTA

```text
--input primers.fasta      \
--input primer_dir/         >--->  PARSE_PRIMERS  --->  primers.tsv  --->  INPUT_CHECK  --->  benchmarking
--primer_set COI,Fish16S   /      (parse_primers.py)
```

`--input` is passed to the `PARSE_PRIMERS` process as a single staged path - one file, or the whole
directory. `bin/parse_primers.py` then:

1. Expands the input into a list of FASTAs. For a directory, that is every `.fa`, `.fasta`, and
   `.fna` file in it, sorted; everything else in the directory is ignored. A directory with none of
   those is an error.
2. Parses and validates every file, collecting all problems rather than stopping at the first.
3. Writes `primers.tsv` and `<id>.primer_warnings.txt`, both published to `primers/` in your
   results.

`--primer_set` takes the same route: each set is downloaded from the catalog and handed to the same
`PARSE_PRIMERS` process, with that set's own `min`/`max` from its catalog config.

## Grouping: prefixes decide what a primer set is

Records are grouped by the **prefix** of their header - the text in front of an `fwd`, `forward`,
`rev`, or `reverse` token (case-insensitive). Each prefix becomes its own primer set.

Anything *after* the token is a variant label and is ignored, so `MA_FWD_1` and `MA_FWD_2` are both
prefix `MA`. The token must be delimited by a separator or the start/end of the header, so a prefix
that merely contains the letters (`TREV_fwd`) is not mistaken for the token itself.

A file with no direction tokens anywhere and exactly two records is accepted as forward then
reverse, in that order.

Grouping by prefix rather than by length is what keeps two markers in one file apart:

```text
>MA_FWD   AAAACCCCGGGGTTTTAC        ->  markers_MA
>MA_REV   TTTTGGGGCCCCAAAAGT
>POL_FWD  AAAACCCCGGGGTTTTAG        ->  markers_POL
>POL_REV  TTTTGGGGCCCCAAAAGA
```

Those four primers are all the same length, but MA and POL stay separate because their prefixes
differ.

## Collapsing and splitting

`obipcr` accepts one forward and one reverse string per run, so each prefix has to resolve to a
single pair.

**Same length: collapsed.** Variants of equal length within a prefix become one IUPAC-degenerate
primer.

```text
>MA_fwd_1  ACGT
>MA_fwd_2  ATGT      ->  fwd AYGT
>MA_rev    TGCA          rev TGCA
```

**Different lengths: split.** They cannot be collapsed, so every forward is paired with every
reverse and each combination is benchmarked as its own primer set. This is a warning, not an error.

```text
>ITS2_fwd_1  AAAA
>ITS2_fwd_2  AAAAAA    ->  ITS2_1  fwd AAAA    rev TTTT
>ITS2_rev    TTTT          ITS2_2  fwd AAAAAA  rev TTTT
```

A single forward and a single reverse - the common case - passes through untouched. Nothing is
degenerated, and any IUPAC codes already in your sequences are preserved exactly.

## Set naming

Names come from the filename, gaining a suffix only where the file is ambiguous:

| File contains | Set names |
| --- | --- |
| one prefix, one length per direction | `ITS2` |
| several prefixes | `markers_MA`, `markers_POL` |
| one prefix that splits on length | `ITS2_1`, `ITS2_2` |
| several prefixes, one splitting | `markers_MA_1`, `markers_MA_2`, `markers_POL` |

Prefer one primer pair per file. Two pairs in one file are kept apart correctly, but separate files
give shorter, clearer names in the results.

## When a file is rejected

Every FASTA is validated before anything else runs, and every problem across every file is reported
at once. Nothing is benchmarked while any file is unusable, so a typo can never silently drop a
primer from the comparison.

A file is rejected when:

- it holds no FASTA records (empty, or no `>` lines)
- a record has an empty sequence
- a sequence contains a non-nucleotide character - only `ACGT`/`U` and IUPAC ambiguity codes are allowed
- a prefix has a forward but no reverse primer, or the other way round
- a record has no direction token, in a file that is not a plain two-record pair

## Choosing between a FASTA folder and a samplesheet

A FASTA folder has one global `--primer_min`/`--primer_max` for every file. Markers with very
different amplicon lengths - ITS2, rbcL, trnL - do not share a sensible window, so a single range
will either cut off real amplicons or admit junk for some of them.

Use a folder when your markers share one length window and you want the convenience. Use a
samplesheet when they do not, since it carries per-primer bounds:

```text
primer  fwd                         rev                     min  max
ITS2    CGAGTYTTTGAAYGCAAGTTG       YCCCGYCTGAYCTGRGGT      200  500
rbcL    ATGTCACCACAAACAGAGACTAAAGC  GTAAAATCAAGTCCACCRCG    500  800
trnL    CGAAATCGGTAGACGCTACG        GGGGATAGAGGGACTTGAAC    100  300
```

The samplesheet takes literal primer strings, so it does no collapsing. If a marker genuinely needs
several variants merged, run it through the FASTA route first and copy the consensus out of
`primers/primers.tsv`.

## Checking what actually ran

`primers/primers.tsv` in your results is the exact set of primers that was benchmarked, and
`primers/*.primer_warnings.txt` records any prefix that had to be split. Read these first when
results do not look like what you expected.
