# FMD runs: the RLBWT form shared by fmd2rlbwt, Orbit and TeraLCP

`fmd2rlbwt IN.fmd OUT` writes the BWT held in a ropebwt3 FMD as
`OUT.bwt.heads` and `OUT.bwt.len`, the form Orbit's `ms build-rlbwt` and
`ms tms-build` read. This note fixes the conventions that any other producer
of per-run data for the same BWT (TeraLCP's `-ominima`, its Phi export) must
follow so that its records line up with these runs.

## The BWT

An FMD is the multi-string BWT of a string collection over ropebwt3's
alphabet `$ A C G T N` (codes 0 to 5), with one `$` per string (each contig
and, by default, its reverse complement). ropebwt3 sorts the symbols in code
order, so N sorts after T. Every `$` ends its own string: no two sentinels
are equal, so a common prefix of two suffixes stops at the first `$`.

## Run numbering

Runs are the **maximal runs of the BWT in row order**: run 0 holds row 0, and
two adjacent runs never have the same symbol. The rld encoding inside an FMD
may hold adjacent entries with the same symbol; fmd2rlbwt merges them. Runs of
`$` are merged as well: several adjacent sentinels form one run, with one
head byte and its length, although the sentinels are distinct.

The sum of the lengths is n, the FMD's symbol count (`ropebwt3 stat`'s
"symbols"), and fmd2rlbwt checks it, and each symbol's total, against the
FMD header.

## Files

* `OUT.bwt.heads`: one byte per run.
* `OUT.bwt.len`: one 5-byte little-endian unsigned length per run, so a run
  is at most 2^40 - 1 long (fmd2rlbwt stops with an error otherwise; n itself
  is not limited).

Head bytes:

| symbol | code | head byte |
|--------|------|-----------|
| `$`    | 0    | `$` (0x24) |
| A      | 1    | `A` (0x41) |
| C      | 2    | `C` (0x43) |
| G      | 3    | `G` (0x47) |
| T      | 4    | `T` (0x54) |
| N      | 5    | 0x7f       |

The bytes increase with the symbol codes, so a reader that orders characters
by byte value (as Orbit does when it builds LF) gets ropebwt3's order. N is
0x7f rather than `N` (0x4e) because `N` sorts before `T`. Orbit's loader maps
`$` to its terminator code and 0x7f to `ms_alpha::TEXT_N`, a byte that no
read character equals, so N in the text matches nothing; it refuses the byte
`N`, whose sort position it cannot know.

## Per-run data from TeraLCP

A TLMSM minima file for these runs (format: TeraTools
`src/TeraLCP/MS_MINIMA_FORMAT.md`) must have one record per maximal run, in
the order above, with header `r` equal to the number of runs, `n` equal to the
BWT length, flag bit 2 clear, and `input_runs = r`. A run of several `$` is
one record: its top is the LCP at its first row and its interior values are
the LCPs at its other rows. LCP values follow the FMD's semantics: extension
stops at any `$`. (Whether N matches N inside an LCP does not change any
matching statistic, since no read character matches N; the numbering rules
above are what must agree.)

TeraLCP at TeraTools bce9ea7 numbers `-f fmd` records by rld entry, with
every run of `l` sentinels split into `l` runs of length 1. Until TeraLCP
writes the maximal numbering itself, fmd2rlbwt converts such a file while it
dumps the runs:

    fmd2rlbwt --minima-in teralcp.min --minima-out fmd.min IN.fmd OUT

The conversion streams both files. A merged run's prefix and suffix minima
are found among its pieces' stored values and the tops of its pieces after
the first, so no other LCP values are needed; the output's header has flag
bit 2 clear.

## Other output

fmd2rlbwt prints, on stdout: n, the number of sequences (sentinels), the rld
entry count (what `ropebwt3 stat` reports as runs), the run count in
TeraLCP's numbering, the maximal run count, the longest run, and per symbol
its count, rld entries and maximal runs.

## Resources

One pass, one thread. The FMD is memory-mapped and each 64 MiB block's pages
are released once decoded (`--no-release` keeps them), and the outputs are
written through 4 MiB buffers, so resident memory stays small whatever the
FMD's size.

## Building

`make` in this directory. `rld0.c` and `rld0.h` are vendored unchanged from
lh3/ropebwt3 (commit dd8d285, MIT license, `LICENSE.ropebwt3`).
`test_fmd2rlbwt.py ROPEBWT3 ./fmd2rlbwt` round-trips random BWTs through
`ropebwt3 plain2fmd` and checks the heads, the lengths and the minima
conversion.
