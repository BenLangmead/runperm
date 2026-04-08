# RunPerm Design Guide

This document outlines the main design decisions to make when building a new structure on top of the runperm library, along with practical notes on how to achieve each design choice.

---

## 1. Contiguous vs. Separate Address Spaces

### The Decision

Should all (or most) fields in a row live in a single contiguous memory layout, or should some columns reside in separate data structures?

**Contiguous layout** improves cache locality: loading one row brings all its data into the same cache line. This matters when traversal patterns access multiple fields per row (e.g., move structure navigation plus auxiliary data like LCP). **Separate layout** can reduce the size of the “hot” row and may improve performance when only a subset of columns is accessed frequently.

### How to Achieve It

| Choice | Runperm Mechanism | Result |
|--------|-------------------|--------|
| **Contiguous** | `IntegratedMoveStructure = true` | Base columns (LENGTH/START, POINTER, OFFSET, etc.) and RunCols are packed into a single `PackedVector`. Each row is one bit-packed sequence of all columns. |
| **Separate** | `IntegratedMoveStructure = false` (default) | Base structure lives in one `PackedVector`; RunCols live in a separate `PackedVector<RunCols>` (`run_cols_data`). Two distinct memory regions. |

**Usage**:

```cpp
// Contiguous: run data packed with move structure
RunPermLF<MyRunCols, true> index(...);   // integrated

// Separated: run data in separate PackedVector
RunPermLF<MyRunCols, false> index(...);  // default
```

**Trade-off** (from README): Integrated layout improves locality but increases the size of each row fetch, which can slow move queries if only the base permutation fields are needed. Prefer integration when auxiliary data is used on most accesses.

---

## 2. Lengths-Based vs. Starts-Based (Relative vs. Absolute Positions)

### The Decision

How should interval boundaries be represented?

- **Lengths-based (relative)**: Store the *length* of each interval. The start of interval `i` is implied by the sum of lengths of intervals `0..i-1`. Position is `{interval, offset}`.
- **Starts-based (absolute)**: Store the *starting index* of each interval. Length is derived as `start[i+1] - start[i]` (or `domain - start[i]` for the last interval). Position adds `idx` for the absolute position in the permutation.

### How to Achieve It

| Choice | Runperm Mechanism | Column Stored | Position Type |
|--------|-------------------|---------------|---------------|
| **Lengths-based** | `StoreAbsolutePositions = false` (default) | `LENGTH` | `{interval, offset}` |
| **Starts-based** | `StoreAbsolutePositions = true` | `START` | `{interval, offset, idx}` |

Internally, this is implemented via `SwitchColumns<BaseColumns, StoreAbsolutePositions>`, which maps to `MoveCols` (LENGTH, POINTER, OFFSET) or `MoveColsIdx` (START, POINTER, OFFSET).

**Usage**:

```cpp
// Lengths-based (relative)
RunPermLF<MyCols, false, false> index(...);

// Starts-based (absolute) — enables idx lookups
RunPermLF<MyCols, false, true> index(...);
```

**Space and Precision** (from README):
- **Relative**: \( O(r \log \frac{n}{r} + r \log r) \) bits — lengths need fewer bits when runs are short.
- **Absolute**: \( O(r \log n) \) bits — each start is in \( [0, n) \).

**When to use starts-based**: Necessary for `pred`/`succ` by index and for binary search over positions by absolute index. Use with `ExponentialSearch = true` for faster multi-step `next()` when you have `idx`.

---

## 3. Variable-Length Data: Fixed Prefix in Row + Spilled Tail

### The Constraint

`PackedMatrix` and `PackedVector` require **fixed bit-width per column**. There is no built-in type for variable-length per-row data.

### The Pattern: Descriptor in Row, Payload Elsewhere

To support variable-length data while keeping a fixed-length portion in the row:

1. **In-row (RunCols)**: Store a small fixed-width descriptor, for example:
   - `SPILL_OFFSET`: index into a spill buffer where this row’s variable data starts
   - `SPILL_LEN`: number of elements (or bytes) in the spill for this row
   - Optionally a pointer-sized field if indirection is acceptable

2. **Out-of-row (spill)**: A separate structure (e.g. `std::vector<ulint>`, another `PackedMatrix`, or raw byte buffer) holding the variable-length payloads. Access pattern: `spill[get<SPILL_OFFSET>(pos) + local_index]`.

**Example RunCols**:

```cpp
enum class LCPRunCols {
    SPILL_OFFSET,   // start index in spill buffer
    SPILL_LEN,      // number of LCP values for this row
    COUNT
};
```

Construction: scan rows, assign `spill_offset` and `spill_len`, append variable data to the spill buffer. Use `extend_run_data` when the move structure splits runs so that spill ranges can be split accordingly.

**Alternative: Fixed Cap with Overflow**

- Store a fixed number of slots per row (e.g. `LCP_0, ..., LCP_K`).
- If a run has more than K values, use the last slot as an overflow index into the spill, or accept truncation/splitting.

### Interaction with Run Splitting

When `SplitParams` causes runs to be split, `extend_run_data(orig_interval, orig_interval_length, new_offset_from_orig_start, new_length)` is called for each new sub-interval. For spilled data:

- The callback can subdivide the original run’s variable data across the new intervals.
- Ensure `SPILL_OFFSET` and `SPILL_LEN` (or equivalent) are set so each new row points to its slice of the spill.

---

## 4. Table Type: MoveVector vs. MoveTable

Runperm supports two storage backends for the move structure:

| Type | Layout | Use Case |
|------|--------|----------|
| **MoveVector** (default) | `PackedVector` — bit-packed, contiguous rows | Compact, cache-friendly; supports `IntegratedMoveStructure` |
| **MoveTable** | `std::vector<Row>` of structs | Per-field access; **does not support** integrated run data |

**Note**: `IntegratedMoveStructure` is incompatible with `MoveTable`. Use `MoveVector` (the default for RLBWT types) when integrating RunCols with the move structure.

---

## 5. Defining RunCols

RunCols are user-defined fields stored per run. Requirements:

- Enum with `COUNT` as the last value.
- Each field maps to a fixed-width value (bit-width inferred from `max` over the data via `get_run_cols_widths`).

```cpp
enum class MyRunCols {
    FIELD_A,
    FIELD_B,
    COUNT
};

using RunData = std::array<ulint, static_cast<size_t>(MyRunCols::COUNT)>;
std::vector<RunData> run_data(num_runs);
```

Access: `index.get<MyRunCols::FIELD_A>(pos)` or `index.get<MyRunCols::FIELD_A>(interval)`.

---

## 6. Quick Reference

| Design Choice | Parameter | Options |
|---------------|-----------|---------|
| Contiguous vs. separated | `IntegratedMoveStructure` | `true` = contiguous, `false` = separate |
| Lengths vs. starts | `StoreAbsolutePositions` | `false` = lengths-based, `true` = starts-based |
| Variable-length data | Custom | Fixed descriptor in RunCols + spill buffer |
| Table backend | `TableType` | `MoveVector` (default) or `MoveTable` |
