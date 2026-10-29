# tms `--mode phi`: per-row LCP summaries for phi repositioning

Agreed design from the Oct 28, 2026 grilling session. The mode was first `--mode phisum`, beside a `--mode phi` that walked every row; on Oct 29 the summary modes became `phi` and `dual`, the walking ones `phi-walk` and `dual-walk` (advanced), and tms-build writes the summaries into every index with phi.

## Goal

Get tms phi-mode results (with positions, toehold and SMEMs) at close to Orbit-ms speed, without ms's spillover arrays. In the paper, the LCP range-min of a reposition is split into a few expensive local queries (the partial current row) plus per-run queries (whole rows read by the candidate scan), which makes Phi read like Orbit-ms. The walking form stays as `--mode phi-walk`.

## Decomposition

With LCP[i] = lcp(SA[i-1], SA[i]), current row r at offset cur:

- Up candidate (tail of row u): min LCP[top(u+1) .. cur] = whole rows u+1..r-1, plus row r's part top(r)..cur.
- Down candidate (head of row d): min LCP[cur+1 .. top(d)] = row r's part cur+1..bottom(r), whole rows r+1..d-1, plus LCP[top(d)] alone.

## Index

- Two LF data columns, which tms-build writes into every index with phi (`--no-lcp-summary` leaves them out):
  - `TOP_LCP(R)` = LCP[top(R)]
  - `INT_MIN(R)` = min of R's interior LCPs
- Exact width (`bit_width(max)`); saturation is a possible follow-up only if row width hurts.
- Built from `--minima` / `--lcp-bin` with ms_build's per-row cut logic, also when phi itself comes from `--phi` (TeraLCP -ophi), as on the FMD route. `--phi` alone has no interior minima and needs `--no-lcp-summary`; a phi sweep over all n rows could derive the columns from it (as the `phi-shortread` prototype did at load), but is deferred.
- No argmin column.

## Query

1. The candidate scan keeps `sum_up` (min of TOP_LCP and INT_MIN over rows u+1..r-1) and `sum_down` (same over rows r+1..d-1), and reads TOP_LCP(d).
2. Per side, upper bound = summary min (with TOP_LCP(r) on the up side, TOP_LCP(d) on the down side); lower bound = min(upper, INT_MIN(r)). If they agree, or the lower bound reaches the cap, the side is done with no walk.
3. Otherwise it walks phi over the partial piece only, with `phi_min` seeded from the summary:
   - up side: from the toehold, cur - top(r) steps;
   - down side: from the phi point of the next row that starts a BWT run (phi points are stored only at run heads). When LF splitting makes row r+1 continue r's run, the walk also crosses the rest of the run; a counter records how often this happens and the extra steps. If it matters, store phi points at every split row head.
4. The pred/succ early decision, positions and toehold are unchanged.
5. PHISKIP and the walking modes are unchanged. `phi` and `dual` refuse an index without the columns.

## Evaluation baselines

tms `--mode phi-walk`, tms `--mode psi`, Orbit-ms, and the `phi-shortread` prototype (5.4x at h10, 6.1x at h100 over Phi on chr22 short reads, exact outputs). Report index load time separately from query time. Correctness: outputs identical to `--mode phi-walk`.
