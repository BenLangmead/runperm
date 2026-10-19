# TMSORC01: the reposition oracle file format

A reposition oracle records, for every read of a query, the outcome of each reposition its query makes: the LF interval the match moves to and the match length there. `ms query INDEX READS --oracle-write FILE` records one, in any mode but the oracle modes and without an ftab (which skips some repositions). `--mode oracle --oracle FILE` and `--mode oracle-scan --oracle FILE` replay it. These options are listed by `ms --help-advanced`. Code: `reposition_oracle.hpp`.

## Replay

- `oracle` moves to the recorded interval with the recorded length on every mismatch, without scanning for the candidate runs, reading LCP bounds or computing an LCE. It reads only LF (and phi, for positions or an SMEM report), so it is a lower bound on the time of the other modes.
- `oracle-scan` takes less from the file: the recorded interval says only which side (the run above the current interval, or the one below) holds the candidate, and that side is scanned for it as the other modes do. The recorded length still stands in for the LCE.
- The output of both is that of the mode that recorded the file.
- Replay needs the same reads in the same order, and an index with the same LF intervals: the same build and, for tms, the same layout. It does not depend on the interleave depth, thread count or block size. The ms and tms engines each replay only their own recordings, since when both sides' LCEs reach the match length they choose different sides.

## Layout

All integers are u64, little-endian.

| Field | Count | Meaning |
|---|---|---|
| magic | 8 bytes | `TMSORC01` |
| move_runs | 1 | the index's LF interval count, checked on replay |
| R | 1 | the read count |
| D | 1 | the reposition count |
| per-read counts | R | the repositions of each read, in input order; they sum to D |
| decisions | D pairs | (interval, length) for each reposition, read by read |

In a decision, interval is the LF interval moved to: the nearest run of the read's character above the current interval, whose tail row is taken, if it is above, else the nearest run below, whose head row is taken. length is the match length before the LF step, the LCE capped at the current length.

A file of R reads and D repositions is 32 + 8R + 16D bytes, as the query's `oracle:` summary line reports.
