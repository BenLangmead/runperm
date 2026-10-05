# ORBITIDX: the Orbit-ms and Orbit-tms index format

Both `ms build`/`build-rlbwt` (Orbit-ms) and `ms tms-build` (Orbit-tms) write this format. The first 4096 bytes of a file describe everything a query loads: which move tables, their columns in row order, and each column's bit width. `ms index-plan INDEX` prints that from the header alone. Code: `index_format.hpp` and `index_format.cpp`.

All integers are little-endian. The endian check refuses a file written on a big-endian host.

## Layout

```
[0, 4096)        header: fixed fields, section table, profile table, info text; zero padded
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
| 42 | u16 | profiles |
| 44 | u16 | info bytes |
| 48 | u64[4] | tool scalars. ms: max LCP top, max LCP min-sub, spill alignment, spill split bits. tms: 0 |
| 80 | | reserved, 0 |

## Section table (from byte 128, 80 bytes per entry)

| offset | type | field |
|---|---|---|
| 0 | u16 | id: 1 LF, 2 F_STARTS, 3 PHI_LEN, 4 PHI_START, 5 PHI_INV, 6 MS_SPILL |
| 2 | u16 | kind: 1 move table, 2 u64 array, 3 spill arrays |
| 4 | u8 | columns (move table; at most 12) |
| 6 | u16 | row bits: the sum of the widths |
| 8 | u64 | offset of the payload, a multiple of 4096 |
| 16 | u64 | payload bytes |
| 24 | u64 | rows (move table), values (u64 array) or arrays (spill arrays) |
| 32 | u64 | domain (move table) |
| 40 | u32 | CRC-32 (IEEE) of the payload |
| 48 | (u8, u8)[12] | (column code, width) per column, in row order |

Column codes: 1 len, 2 start, 3 ptr, 4 off, 5 chr, 16 psi_int, 17 psi_off, 18 phi_int, 19 phi_off, 20 plcp, 21 plcpb, 32 lcp_top, 33 lcp_min_sub, 34 lcp_spill. A width of 0 means the column takes no bits and reads as 0, as with the LF columns a tms layout omits. Column 0 says whether a move table stores lengths (`len`) or absolute starts (`start`).

Payloads:

- **Move table**: the structure's Orbit `serialize()` bytes. These are its own preamble (magics, versions, split parameters, n, r, rows, widths), then the packed rows, `ceil(rows * row_bits / 8) + 8` bytes. The loader checks the loaded rows, domain and widths against the entry, and the entry's column codes against the ones the code reads.
- **u64 array**: the values. F_STARTS has 7, the first row of each F block after the first.
- **Spill arrays**: per array, a u64 byte count and the bytes.

## Profile table (after the section table, 32 bytes per entry)

| offset | type | field |
|---|---|---|
| 0 | char[24] | name, NUL padded |
| 24 | u32 | mask: bit `id` set for each section loaded |
| 28 | u8 | log2 positions per entry of phi's start table built at load, or 255 for none |
| 29 | u8 | the same for phi_inv's start table |

A start table has `((domain - 1) >> shift) + 2` entries of `max(1, bit_width(rows - 1))` bits.

An ms index has one profile, `ms`, which loads LF and MS_SPILL. A tms index has a profile for each set of parts it holds, named by its parts joined with `+`. It always loads LF; `psi` adds F_STARTS; `phi-len` or `phi-starts` adds one form of phi; `phi-inv` adds PHI_INV and builds start tables. `lf` alone names LF with nothing else. `tms-batch` picks its profile as follows:

| query | profile |
|---|---|
| `--mode psi` | `psi` |
| `--mode phi` | `phi-len` |
| `--mode phi --phi-starts`, or `--positions` or `--report smem-one` in phi mode | `phi-starts` |
| `--mode phiskip` or `dual` | `psi+phi-len` |
| psi, phiskip or dual with `--positions` or `--report smem-one` | `psi+phi-starts` |
| `--report smem-all` | `phi-starts+phi-inv`, or `psi+phi-starts+phi-inv` outside phi mode |

`ms index-plan INDEX` prints this table for the index at hand. The loader looks up the profile by name and refuses an index whose profile names other sections than the code expects, so the header's plan is the load's plan.

LF is loaded whole, with every column the layout keeps. On a full-layout index a phi-mode query therefore loads the PSI columns and the union rows it does not read; build `--layout phi` for phi-only work.

## Info text

The info text records provenance as `key=value` lines: `built_by` (the command line) and `built_at` (UTC). A rewrapped index records `rewrapped_from`, `rewrapped_by` and `rewrapped_at` instead.

## Rewrapping earlier indexes

`ms index-rewrap LEGACY NEW` rewrites an MSIX v4 or TMSX v8 index in this format without loading it. It parses each structure's Orbit preamble to find where the structure ends, then copies its bytes into a section. It also counts the BWT runs from LF's character column while copying. Memory use is constant. The output equals a fresh build of the same input with the same options, apart from the info text (see `test_rewrap.sh`).

Each section's bytes are read twice, once to size and checksum it so that the header can go first, and once to copy it.
