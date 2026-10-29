# ORBITIDX: the Orbit-ms and Orbit-tms index format

Both `ms build-rlbwt` (Orbit-ms) and `ms tms-build` (Orbit-tms) write this format. The first 4096 bytes of a file describe everything it holds: each move table, its columns in row order, and each column's bit width. The format describes what is stored. Which parts a query loads is the software's choice, so it can change without changing the format. Given that choice, the header alone says exactly what a query will hold in RAM; `ms index-plan INDEX` prints it. Code: `index_format.hpp` and `index_format.cpp`.

All integers are little-endian. The endian check refuses a file written on a big-endian host.

## Layout

```
[0, 4096)        header: fixed fields, section table, info text; zero padded
[4096, ...)      sections, each starting on a 4096-byte boundary, in section-table order
```

The file ends where the last section ends; `ms index-verify` checks this along with every section's CRC-32.

## Fixed fields (bytes 0 to 127)

| offset | type | field |
|---|---|---|
| 0 | char[8] | magic `ORBITIDX` |
| 8 | u32 | format version, 1 |
| 12 | u32 | endian check, 0x01020304 |
| 16 | u32 | header bytes, 4096 |
| 20 | u16 | tool: 1 = ms, 2 = tms |
| 22 | u16 | layout: 0 = none (ms), 1 = full, 2 = psi, 3 = phi |
| 24 | u64 | n, the text length, terminators included |
| 32 | u64 | BWT runs: maximal runs of one character in LF's rows |
| 40 | u16 | sections |
| 42 | u16 | reserved, 0 |
| 44 | u16 | info bytes |
| 48 | u64[4] | tool scalars. ms: max LCP top, max LCP min-sub, spill alignment, spill split bits. tms: 0 |
| 80 | | reserved, 0 |

## Section table (from byte 128, 80 bytes per entry)

| offset | type | field |
|---|---|---|
| 0 | u16 | id: 1 LF, 2 F_STARTS, 3 PHI_LEN, 4 PHI_START, 5 PHI_INV, 6 MS_SPILL, 7 SEQ_IDS |
| 2 | u16 | kind: 1 move table, 2 u64 array, 3 spill arrays, 4 sequence table |
| 4 | u8 | columns (move table; at most 12) |
| 6 | u16 | row bits: the sum of the widths |
| 8 | u64 | offset of the payload, a multiple of 4096 |
| 16 | u64 | payload bytes |
| 24 | u64 | rows (move table), values (u64 array), arrays (spill arrays) or sequences (sequence table) |
| 32 | u64 | domain (move table) |
| 40 | u32 | CRC-32 (IEEE) of the payload |
| 48 | (u8, u8)[12] | (column code, width) per column, in row order |

Column codes: 1 len, 2 start, 3 ptr, 4 off, 5 chr, 16 psi_int, 17 psi_off, 18 phi_int, 19 phi_off, 20 plcp, 21 plcpb, 22 top_lcp, 23 int_min, 32 lcp_top, 33 lcp_min_sub, 34 lcp_spill. A width of 0 means the column takes no bits and reads as 0, as with the LF columns a tms layout omits. Column 0 says whether a move table stores lengths (`len`) or absolute starts (`start`). In a tms LF, `phi_int` holds, at each row that starts a BWT run, the φ interval that starts at that row's text position (0 at other rows), and `phi_off` has width 0, since the φ points queries read all have offset 0. A tms LF may have two more columns after the eight above, `top_lcp` and `int_min`, the LCP summaries of each row that `tms-build --lcp-summary` or `tms-add-summary` adds and `--mode phisum` and `dualsum` read. Without summaries the section leaves them out, entry and payload alike, so it is the LF section of a build that predates them, and a reader loads an eight-column LF section with both at width 0. Indexes built or rewrapped before the summaries therefore load unchanged, and `tms-add-summary` gives one its summaries, rewriting LF alone. Indexes rewrapped from TMSX v8 keep the nonzero `phi_off` width their builds gave it, which loads and answers the same.

Payloads:

- **Move table**: the structure's Orbit `serialize()` bytes. These are its own preamble (magics, versions, split parameters, n, r, rows, widths), then the packed rows, `ceil(rows * row_bits / 8) + 8` bytes. The loader checks the loaded rows, domain and widths against the entry, and the entry's column codes against the ones the code reads.
- **u64 array**: the values. F_STARTS has 7, the first row of each F block after the first.
- **Spill arrays**: per array, a u64 byte count and the bytes.
- **Sequence table**: per sequence, in order, its length as a u64, its name's byte count as a u32, and the name.

## The sequence-ID table

A tms index built from a ropebwt3 FMD may hold the FMD's sequence-ID table, ropebwt3's `.len.gz`, as the SEQ_IDS section, last in the file. Sequence k is the FMD's strings 2k (its forward strand) and 2k + 1 (its reverse complement), so the text is `C_0 $ rc(C_0) $ C_1 $ rc(C_1) $ ...`, and the lengths must give the header's n as 2 * sum(len_k + 1). `query --report smem-all` uses the table to write positions as contig, strand and forward-strand offset, as ropebwt3 does. It is added at build time by `ms tms-build --seq-ids FILE`, which also checks that the BWT has two terminators per sequence. `ms index-add-seq-ids INDEX FILE OUT` adds it to an index already built. It copies the other sections byte for byte, at the same offsets, and checks only n, since counting the terminators would mean loading LF. A table that fits both checks but splits the text elsewhere makes a query fail when a match crosses a sequence's end.

Neither changes the sections an ftab links to, so an ftab built from the index without the table serves it with the table. Earlier `ms` binaries load such an index and ignore the section.

## What this software loads (informative)

This section describes the current `ms` binary, not the format.

- **Orbit-ms:** every query loads every section.
- **Orbit-tms:** a query always loads LF, whole, with every column the layout keeps (`TmsIndex::sections_for`). It then adds:
  - **F_STARTS** if the query walks psi (not in the oracle mode, which computes no LCE);
  - **PHI_START** if it needs text positions (`--positions`, an SMEM report, or `--phi-starts`), or else **PHI_LEN** if it walks phi;
  - **PHI_INV** for `--report smem-all`, and **SEQ_IDS** if the index holds it. In that case it also builds a start table for PHI_START and for PHI_INV. Each has `((domain - 1) >> shift) + 2` entries of `max(1, bit_width(rows - 1))` bits, where `shift` is the smallest with `2^shift * rows >= 4 * domain`, computed from the section's entry.

| query | sections loaded |
|---|---|
| `--mode psi` | LF, F_STARTS |
| `--mode phi` or `phisum` | LF, PHI_LEN |
| `--mode phi` with `--phi-starts`, `--positions` or `--report smem-one` | LF, PHI_START |
| `--mode phiskip`, `dual` or `dualsum` | LF, F_STARTS, PHI_LEN |
| psi with `--positions` or `--report smem-one`; phiskip, dual or dualsum with those or `--phi-starts` | LF, F_STARTS, PHI_START |
| `--mode oracle` or `oracle-scan` (replaying recorded repositions) | LF |
| `--mode oracle` or `oracle-scan` with `--positions` or `--report smem-one` | LF, PHI_START |
| `--report smem-all` | LF, PHI_START, PHI_INV, SEQ_IDS if held, plus F_STARTS outside phi and oracle modes |

`ms index-plan INDEX` prints, for the index at hand, the sections and start tables each tms `query` command line loads and the RAM they take.

Because LF is loaded whole, a phi-mode query on a full-layout index also loads the PSI columns and the union rows it does not read. Build `--layout phi` for phi-only work.

## Companion files

An ftab (FTAB_FORMAT.md), a table of query states by k-mer that `query` takes with `--ftab`, is stored as a file of its own next to its index. Its header links it to the index through copies of the fixed fields and of the section entries its walk read.

## Info text

The info text records provenance as `key=value` lines: `built_by` (the command line) and `built_at` (UTC). A rewrapped index records `rewrapped_from`, `rewrapped_by` and `rewrapped_at` instead.

## Rewrapping earlier indexes

`ms index-rewrap LEGACY NEW` rewrites an MSIX v4 or TMSX v8 index in this format without loading it. It parses each structure's Orbit preamble to find where the structure ends, then copies its bytes into a section. It also counts the BWT runs from LF's character column while copying. Memory use is constant. The output equals a fresh build of the same input with the same options, apart from the info text (see `test_rewrap.sh`).

Each section's bytes are read twice, once to size and checksum it so that the header can go first, and once to copy it.
