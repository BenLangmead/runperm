# AGENTS.md

This file provides guidance to WARP (warp.dev) when working with code in this repository.

## Build Commands

```bash
# Build all test executables (move_test, rlbwt_test, runperm_test) and MS application
make

# Build with debug symbols
make debug

# Clean build artifacts
make clean

# MS application (matching statistics with LCP): build and run tests
make ms && ./src/ms/ms test ../../data

# MS debug build (-O0 -g) for breakpoints
make -C src/ms debug
# Then run under lldb/gdb: lldb src/ms/ms, then (lldb) run test ../../data
```

## Testing

No test framework is used. Tests are assert-based executables that print results to stdout:

```bash
# Build and run all tests
make && ./move_test && ./rlbwt_test && ./runperm_test
make -C src/ms test

# Run individual test
./runperm_test
./src/ms/ms test ../../data
```

Tests use random permutations and verify correctness via assertions. Test files in `src/` contain benchmark-style timing output alongside correctness checks.

## Architecture

This is a **header-only C++17 library** for run-length encoded permutations with move structures. No external dependencies.

### Public Headers (in `include/`)

- `runperm.hpp` - Main entry point for `RunPerm` (with user-defined run data) and `MovePerm` (without user data)
- `rlbwt.hpp` - RLBWT specializations: `RunPermLF`/`MoveLF`, `RunPermFL`/`MoveFL`, `RunPermPhi`/`MovePhi`, `RunPermInvPhi`/`MoveInvPhi`
- `move.hpp` - Low-level `MoveStructure` types (advanced usage)

### MS Application (in `src/ms/`)

The Matching Statistics (MS) application lives in `src/ms/` and is a concrete use of the runperm library. Build with `make ms` or `make -C src/ms`. Subcommands:
- `ms build <input.tsv> <output.idx>` - Build index from movify TSV
- `ms ms <pattern> <index.idx>` - Compute matching statistics
- `ms test [data_dir]` - Run tests
- `ms inspect <index.idx>` - Report bit usage
- `ms lcp-list <index.idx> [output]` - Output LCP columns

### Internal Structure (`include/internal/`)

- `common.hpp` - Type definitions (`ulint`, `uchar`), macros, utility functions (`get_permutation_intervals`, `bwt_to_rlbwt`)
- `ds/` - Data structures: `packed_vector.hpp` (bitpacked storage), `alphabet.hpp` (character mappings)
- `move/` - Move structure implementation: `move_structure.hpp`, `move_table.hpp`, `move_row.hpp`, `move_columns.hpp`, `move_splitting.hpp`
- `runperm/` - Core RunPerm implementation with run column support
- `rlbwt/` - RLBWT-specific implementations (LF/FL/Phi/InvPhi) with specializations in `specializations/`

### Key Design Patterns

**Template parameters control behavior:**
- `IntegratedMoveStructure` - Run data stored alongside move structure (better cache locality) vs separate
- `StoreAbsolutePositions` - Absolute positions (enables idx lookups, more memory) vs interval/offset pairs

**Position type:** Navigation uses a `Position` struct containing either `{interval, offset}` or `{interval, offset, idx}` depending on template params.

**User-defined columns:** `RunPerm` accepts a `RunColsType` enum (must end with `COUNT`) to define custom data stored per run.

### Genetic functionality versus specific BWT functionality

The Move Structure is a generic structure that allows for computation with runny permutations.  A premier example of a useful runny permutation is a BWT of a repetitive text.  BWTs of repetitive texts, once they have been compressed as in the run-length (RL)BWT, are superb foundations for compressed full-text indexes.

One of the text programs -- src/rlbwt_test.cpp -- gives an example of how an RLBWT index can be built atop a move structure, and how to implement important pattern matching queries on top of that, including LF, FL, phi and inverse phi.

### Matching Statistics (MS) application

The MS application in `src/ms/` adds LCP support for one-pass matching statistics. It uses:

- **LCP RunCols**: `lcp_top`, `lcp_min_sub`, `lcp_spill` (spillover for skinny/jumbo rows)
- **Smoovify compression**: gaps replace interior LCP values ≥ boundary max
- **Three LCP query types**: boundary, row-min, range-min (to top or bottom)
- **TSV bridge**: loads movify TSV format (id, off, len, c, sa, lcp)
- **Serialization**: binary index format for build/load

## Citation

If referencing this work: Brown, N. K., & Langmead, B. (2026). arXiv:2602.11029
