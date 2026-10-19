# Query output

`ms query INDEX READS` writes one line per read, in input order, to `-o FILE` or stdout (gzip-compressed when FILE ends in `.gz`). A summary line goes to stderr. Output does not depend on the thread count, interleave depth, block size, mode or ftab.

## Matching statistics (the default, `--report ms`)

```
NAME <tab> MS_0 MS_1 ... MS_{m-1}
```

NAME is the read's name (from FASTA or FASTQ; for one sequence per line, `read<i>`, counting reads from 1), and MS_i is the length of the longest prefix of the read's suffix starting at i that occurs in the text, as space-separated decimals. With `--orient rc` or `qual`, the values are those of the read as given, converted back from the strand matched.

## Positions (`--positions`, tms only)

A third tab-separated field holds, for each MS_i, one text position where that match occurs, space-separated, or `-1` where MS_i is 0. It needs an index with phi.

## SMEMs of one position each (`--report smem-one`, tms only)

A last tab-separated field holds the read's SMEMs (super-maximal exact matches), space-separated, each as `i:L:p`: its start i in the read, its length L, and one text position p where it occurs. `--min-smem-len L` keeps the SMEMs of length L or more. It needs phi.

## Every occurrence (`--report smem-all`, tms only)

In place of the read's line, one line per SMEM, as `ropebwt3 mem -p N` writes them, tab-separated:

| Column | Meaning |
|---|---|
| 1 | the read's name |
| 2 | the SMEM's start in the read (0-based) |
| 3 | its end (exclusive) |
| 4 | c, the count of its occurrences in the text |
| 5 | min(N, c), the number of positions listed (only if N > 0) |
| 6 and on | that many positions, the first ones in BWT row order (only if N > 0) |

N is `--max-smem-positions` (default: every occurrence). It needs phi and phi_inv, and cannot be used with `--positions`.

If the index holds a sequence-ID table (`tms-build --seq-ids`, or `index-add-seq-ids`), each position is `contig:strand:offset`, with a 0-based offset; on the `-` strand it is the offset on the forward strand of the reverse complement's match, as ropebwt3 gives. Otherwise each position is a text position, and the query says so on stderr.

`--min-smem-len L` keeps the SMEMs of length L or more, as ropebwt3's `-l L`. With L = 1 (the default), each read character that occurs nowhere in the text, such as an N, also gets a line with its start, its end and a count of 0, as ropebwt3 gives with `-l 1`. A read with no SMEM and no such character adds no line.

## Summary (stderr)

`query: reads=... bases=... threads=... index_load_s=... ftab_k=... ftab_s=... query_s=... read_s=... parse_s=... format_s=... orient=... flipped=... orient_s=... total_s=... query_ns_per_base=... wall_ns_per_base=...`

query_s, parse_s, format_s and orient_s sum the threads' times; read_s is the time the input lock was held for reading; query_ns_per_base is query_s per base and wall_ns_per_base the elapsed time of the read, query and write phase per base. flipped counts the reads `--orient` matched on the other strand.
