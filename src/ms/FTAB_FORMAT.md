# ORBIFTAB: the ftab file format

An ftab is a dense table, indexed by every k-mer over ACGT, of the state an Orbit-ms or Orbit-tms query reaches after reading that k-mer right to left from a cold start. A read whose last k bases occur in the text takes its state from the table and starts k steps in; its statistics in that window are implied (d + 1 at d bases from its right end). Results are byte-identical with and without the table.

`ms ftab-build INDEX FTAB` writes one, as a file of its own next to an ORBITIDX index (INDEX_FORMAT.md); `ms query` loads it with `--ftab FTAB`, or builds one at load with `--build-ftab[=KLEN]` (default KLEN 14) without storing it. A table is used only with the index it was built from: its header carries a link to that index, which `--ftab` and `ms ftab-verify INDEX FTAB` check against the index's header alone. `ms ftab-info FTAB` prints the header. Code: `ftab.hpp` (the table), `ftab_build.hpp` (its walk and the choice of k), `ftab_file.hpp` (this format).

All integers are little-endian. The endian check refuses a file written on a big-endian host.

## Layout

```
[0, 4096)              header: fixed fields, link entries, info text; zero padded
[4096, ...)            table: the packed entries
[next 4096 multiple)   window positions, if the table has them
```

## Fixed fields (bytes 0 to 127)

| offset | type | field |
|---|---|---|
| 0 | char[8] | magic `ORBIFTAB` |
| 8 | u32 | ftab format version, 1 |
| 12 | u32 | endian check, 0x01020304 |
| 16 | u32 | header bytes, 4096 |
| 20 | u8 | k, 1 to 31 |
| 21 | u8 | how k was chosen: 1 = given (`--kmer-len`), 2 = the auto rule (`--auto`) |
| 22 | u8 | flags: bit 0 toehold, bit 1 window positions |
| 23 | u8 | entry bits, the sum of the entry's widths |
| 24 | u8 | LF interval bits |
| 25 | u8 | LF offset bits |
| 26 | u8 | toehold interval bits (0 without a toehold) |
| 27 | u8 | toehold offset bits (0 without a toehold) |
| 28 | u8 | window position bits (0 without window positions) |
| 29 | u8 | reserved, 0 |
| 30 | u16 | link entries |
| 32 | u16 | info bytes |
| 34 | | reserved, 0 |
| 40 | u64 | occurring k-mers |
| 48 | u64 | build time, seconds since the Unix epoch (UTC) |
| 56 | u64 | the index file's bytes at build time (informative) |
| 64 | u16 | index tool (copy): 1 = ms, 2 = tms |
| 66 | u16 | index layout (copy): 0 = none (ms), 1 = full, 2 = psi, 3 = phi |
| 68 | u32 | index format version (copy) |
| 72 | u64 | index n (copy) |
| 80 | u64 | index BWT runs (copy) |
| 88 | u64 | table offset, 4096 |
| 96 | u64 | table bytes |
| 104 | u64 | window offset (0 without window positions) |
| 112 | u64 | window bytes |
| 120 | u32 | CRC-32 (IEEE) of the table |
| 124 | u32 | CRC-32 of the window positions |

Occupancy is the occurring k-mers over 4^k.

## Link entries (from byte 128, 80 bytes each)

Each is the index's section-table entry, in the ORBITIDX encoding, for one section the table's walk read, with the payload offset set to 0 (where a section lies is not part of its content). A table links exactly these sections, among those the index holds:

| index | sections |
|---|---|
| ms | LF, MS_SPILL |
| tms, psi layout | LF, F_STARTS |
| tms, full layout | LF, F_STARTS, PHI_LEN, PHI_START |
| tms, phi layout | LF, PHI_LEN, PHI_START |

A tms table with a toehold links both forms of phi, since its phi points read the same in either and a query loads one or the other.

The link matches an index when the index header's tool, layout, n, BWT runs and format version equal the copies, and each linked section is present with the same kind, rows, domain, payload bytes, row bits, column codes and widths, and CRC-32; and the index holds no other section of the list above. The rest of the index header, notably its info text, does not enter the link, so an index rewrapped from the same content (`ms index-rewrap`) still matches. `ftab-verify` and `query --ftab` name the first field that differs.

## Info text

`key=value` lines: `built_by` (the command line), `built_at` (UTC), `commit` (the commit the `ms` binary was built from; `-dirty` if its source differed), `index_path` (the index's path at build time, informative), `k_rule` (given, or the auto rule with the largest k that fit), and for the auto rule `occupancy`, the occupancy at each of the last few depths the walk counted. Then `threads` and `build_s`.

## Table

Entry x, for the k-mer whose character d places from its right end has code (A 0, C 1, G 2, T 3) in bits 2d and 2d + 1 of x, takes bits [x · E, (x + 1) · E) of the table, E being the entry bits, numbering bits from the least significant bit of the first u64 word. Its fields, from its lowest bit:

1. present, 1 bit: set when the k-mer occurs. The other fields of an absent k-mer are 0.
2. LF interval, then LF offset: the resolved LF position after the walk.
3. With a toehold: phi interval, then phi offset: the resolved phi point of the text position of that LF row, which the tms query carries as its toehold when it walks phi or reports positions.

The table is ceil(4^k · E / 64) + 1 words; the last word is spare, so that a read of an entry never runs past the end.

Tables of ms indexes, and of tms indexes with the psi layout, have no toehold. Tables of tms indexes that hold phi (full and phi layouts) have a toehold and window positions, and serve every tms `query` mode, `--positions` and SMEM reports. A table of a phi layout is filled by a walk that reads repositions' LCEs from phi's PLCP values, since that layout has no psi.

## Window positions

With `--positions`, a tms query reports for each statistic a text position, the one its own walk reached there. Inside the window those come from the walk over each suffix of the read's last k-mer, so the table also holds, for each j from 1 to k - 1 and each j-mer y (keyed as above), the text position the walk over y ends on. Entry (j, y) is at index (4^j - 4) / 3 + y, of the window position bits each, packed as the table is; the k-mer's own position is its toehold's. There are (4^k - 4) / 3 entries in ceil(entries · W / 64) + 1 words, W being the window position bits.

## Choosing k (informative)

`ftab-build --kmer-len KLEN` builds k = KLEN. `--auto` (the default) takes the largest k whose table and window payloads fit in 5% of the index file's bytes; then, while fewer than 80% of that k's k-mers occur, it takes one smaller, stopping at 1. The walk that builds the table at the largest k that fits counts the occurring j-mers at every depth up to it, so the k it steps down to is known from that one walk; the table at that k then takes a second walk.
