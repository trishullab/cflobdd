# CFLOBDDs — Repository Guide for Claude Code

This repository hosts research code on **Context-Free-Language Ordered Binary
Decision Diagrams (CFLOBDDs)** and their weighted variants (WCFLOBDDs), plus
applications: integer multiplication via the Chinese Remainder Theorem (CRT),
ADD-based reference implementations using CUDD, and string compression
experiments combining SEQUITUR with CFLOBDDs/ADDs.

CFLOBDDs are a plug-compatible alternative to BDDs that achieve
**double-exponential** compression in the best case (vs. exponential for BDDs
and SEQUITUR). See [README.md](README.md) and the papers under
[CFLOBDD/docs/](CFLOBDD/docs/):
- TOPLAS 2024 / arXiv 2211.06818 — base CFLOBDDs
- OOPSLA 2024 — Weighted CFLOBDDs
- [CFLOBDD/docs/Multiplication_via_CRT.pdf](CFLOBDD/docs/Multiplication_via_CRT.pdf)
- [CFLOBDD/docs/WPP-CFLOBDD-research-summary.md](CFLOBDD/docs/WPP-CFLOBDD-research-summary.md)
  — current SEQUITUR/WPP compression research planning

---

## Repository layout

| Path | What lives here |
|---|---|
| [CFLOBDD/](CFLOBDD/) | **Core C++ CFLOBDD library and CLI driver.** Self-contained, no CUDD dependency. The main executable (`cflobdd`, or per-width `cflobdd32`..`cflobdd8192`) runs from here. |
| [examples/](examples/) | ADD-based reference implementations using CUDD: integer multiplication ([addMultiplication.cc](examples/addMultiplication.cc)), SEQUITUR→ADD conversion ([sequiturToAdd.cc](examples/sequiturToAdd.cc)), trace→ADD ([traceToAdd.cc](examples/traceToAdd.cc)), and CFLOBDD↔ADD comparison drivers. |
| [cudd-3.0.0/](cudd-3.0.0/) | Stock CUDD 3.0.0 (with `DD_MAX_CACHE_TO_SLOTS_RATIO` reverted to the upstream default of 4 for fair benchmarking). |
| [cudd-big/](cudd-big/), [cudd-complex/](cudd-complex/), [cudd-complex-big/](cudd-complex-big/), [cudd-fourier/](cudd-fourier/), [cudd-addmin/](cudd-addmin/) | CUDD variants used as baselines / specialized weighted ADDs. |
| [sequitur/](sequitur/) | Third-party SEQUITUR implementation (jw/cpp-sequitur), used to produce grammars from input strings/traces. |
| [examples/compression-tests/](examples/compression-tests/) | Synthetic input generators (`gen_repeat`, `gen_luby`) and the test harness for SEQUITUR↔CFLOBDD/ADD compression experiments. |
| [psdd/](psdd/), [ncompress/](ncompress/), [tensor/](tensor/) | Adjacent / experimental code (PSDD, classical compression baselines, tensor work). |
| [vtune/](vtune/), [launch-vtune.bat](launch-vtune.bat) | Intel VTune profiling support; see the VS Code task **"Launch VTune Profiler"**. |

The repo expects a sibling **Boost** install at `../../boost_1_81_0/` on
Windows; on Linux the Makefiles use the system `libboost-all-dev`.

---

## Supported platforms

Two platforms are exercised:

1. **Windows (MSYS2/MinGW UCRT64)** — the Windows branch of
   [CFLOBDD/Makefile](CFLOBDD/Makefile) sets `BOOST_INC` to the sibling
   Boost tree, links `-lpsapi`, and passes `-Wl,--stack,...`. It also
   **pins the compiler** to `C:/msys64/ucrt64/bin/g++.exe` when that path
   exists, falling back to a bare `g++` otherwise — Cygwin's `g++` is
   second on PATH, and mixing the two toolchains produces link errors
   with no obvious cause. `make` itself is Cygwin's
   (`C:\cygwin64\bin\make.exe`), which drives the MSYS2 compiler without
   trouble.
   [examples/Makefile](examples/Makefile) additionally uses
   `-fuse-ld=lld` on Windows to work around a GNU `ld` COMDAT bug when
   linking CFLOBDD + CUDD.
2. **Linux / macOS** (native) — uses system Boost
   (`libboost-all-dev`); no platform-specific link flags.

---

## Building

### Core CFLOBDD library + driver — [CFLOBDD/Makefile](CFLOBDD/Makefile)

```bash
cd CFLOBDD/
make                    # builds ./cflobdd  (default width)
make build-32           # builds ./cflobdd32  (NUM_BITS=32)
make build-64 ... make build-8192   # other widths
make clean
```

- `NUM_BITS=N` is a compile-time macro. Only a small set of object files
  depend on it (`multiplication_crt.o`, `tests_cfl.o`, `cflobdd_c.o`,
  `return_map_specializations.o`); the per-width `build-N` recipes delete
  exactly those before recompiling.
- `WCFLOBDD_SUPPORTED=1` enables the weighted-CFLOBDD sources (~68 of 96
  source files are otherwise excluded). Build it as
  `make WCFLOBDD_SUPPORTED=1 ...`.
- Compiler flags: `-g -O3 -std=c++2a -w -DCFLOBDD_C_EXPORTS -MMD -MP`.
  The `-g` is intentional (debug symbols + `-O3` is correct for VTune
  profiling). Build artifacts include `.d` dependency files.

### ADD/CUDD multiplication baselines — [examples/Makefile](examples/Makefile)

```bash
cd examples/
make multmod-direct     # builds ./addMultiplication using a local cudd_build/
make multmod-32 ... make multmod-4096
make multmod-all
```

- `multmod-direct` bypasses libtool and builds CUDD directly into
  `examples/cudd_build/libcudd.a`. This is the path that works on
  Windows (MSYS2). The libtool-based `multmod` target is preserved for
  Linux/macOS.
- CUDD's C sources are compiled with `-std=gnu99` because they predate
  C23; modern `gcc-14` defaults to C23 which would turn the implicit
  function declarations into errors.

### SEQUITUR — [sequitur/](sequitur/)

```bash
cd sequitur/
make            # or: cmake -S . -B build && cmake --build build
./sequitur <input-file>
```

---

## Running

The CFLOBDD driver dispatches on its first argument; common subcommands
(see [CFLOBDD/main.cpp](CFLOBDD/main.cpp) and
[CFLOBDD/tests_cfl.cpp](CFLOBDD/tests_cfl.cpp)):

```bash
./cflobdd<N> NumsModK <args>          # ProtoCFLOBDDNumsModK / NumsModK tests
./cflobdd<N> MultModK <args>          # (x*y) mod k via ApplyAndReduce
./cflobdd<N> shiftadd      | shiftadd-all      # shift-and-add multiplication
./cflobdd<N> karatsuba     | karatsuba-all     # subtractive Karatsuba
./cflobdd<N> crt-multiply                      # full CRT multiplication
./cflobdd<N> factor                            # factoring relation
```

Multiplication via CRT lives in
[CFLOBDD/multiplication_crt.cpp](CFLOBDD/multiplication_crt.cpp) /
[CFLOBDD/multiplication_crt.h](CFLOBDD/multiplication_crt.h). Garner's
algorithm reconstructs from per-prime residues; the per-prime CFLOBDDs are
built via `ProtoCFLOBDDNumsModK` and `MultModK`.

---

## Code layout (CFLOBDD/)

- Core: [cflobdd_node.h](CFLOBDD/cflobdd_node.h) /
  [cflobdd_node.cpp](CFLOBDD/cflobdd_node.cpp) defines `CFLOBDDInternalNode`
  (level *L*, AConnection + BConnection[]), `CFLOBDDNodeHandle`
  (canonicalized smart pointer), and the cross-product machinery
  (`PairProduct`, `TripleProduct`).
- Templated CFLOBDD class: [cflobdd_t.h](CFLOBDD/cflobdd_t.h);
  per-terminal-type instantiations are split across `cflobdd_int.cpp`,
  `cflobdd_top_node_int.cpp`, the `cflobdd_*_boost.cpp` files, etc.
- Matrix / vector ops: `matrix1234_*` and `vector_*`.
- Weighted variants: `weighted_*.cpp` and `w*.cpp` (gated by
  `WCFLOBDD_SUPPORTED`).
- Tests: `runTests` in [tests_cfl.cpp](CFLOBDD/tests_cfl.cpp), called from
  [main.cpp](CFLOBDD/main.cpp). **Always call `InitModules()` before
  constructing any CFLOBDD and `ClearModules()` before exiting.**
- Namespace: `CFL_OBDD` (note: assignment code uses a separate `SH_OBDD`
  namespace).
- `CFLOBDD_MAX_LEVEL = 13`; widths above the natural max are handled by
  *topmost embedding* (see Conventions below).

---

## Conventions and gotchas

- **Don't mix object files from different toolchains.** Switching between
  MSYS2, Cygwin, and a Linux host (or between substantially different g++
  versions on any one of them) can leave incompatible `.o` files behind.
  Run `make clean` after switching environments — silent link failures
  otherwise. The Windows compiler pin above guards against the MSYS2/Cygwin
  case by accident, but not against a deliberate `make CC=...`.
- **`1ULL << n` overflows for `n ≥ 64`.** When working with `INPUT_TYPE`
  (variably `uint32_t` … `mp::uint4096_t`), use `INPUT_TYPE(1) << n`.
- **Topmost embedding** (multi-width support, see
  [multiplication_crt.cpp](CFLOBDD/multiplication_crt.cpp)):
  `virtualMaxLevel = log2(NUM_BITS)+1`,
  `bottomLevel = CFLOBDDMaxLevel - virtualMaxLevel`,
  `stride = 2^bottomLevel`. Real variables are placed every `stride`
  positions; below `bottomLevel` the structure is a `NoDistinctionNode`
  padding. Use `MkProjection(i * stride)` (not `MkProjection(i)`).
- **Reference counting.** CFLOBDD `NodeHandle`s canonicalize and
  ref-count automatically. CUDD `DD::DD(Cudd*, DdNode*)` auto-refs and
  `~ABDD` auto-derefs — just `return ADD(mgr, result);`, no manual
  `Ref`/`Deref`.
- **CUDD computed-table is direct-mapped.** Collisions evict. The default
  2^18 slots effectively never auto-resizes. For Karatsuba-style workloads
  pass an explicit large initial cache:
  `Cudd(0, 0, numSlots, cacheSize, 0)`. Watch out for caches larger than
  L3 — they slow things down.
- **Don't add `extern` declarations in [hashset.h](CFLOBDD/hashset.h) for
  constants defined in [hashset.cpp](CFLOBDD/hashset.cpp).** `hashset.h`
  *includes* `hashset.cpp`; the constants have internal linkage and adding
  `extern` promotes them to external linkage, causing ODR violations at
  link time.

---

## Documentation pointers

- [README.md](README.md) — top-level project description and citation info.
- [CFLOBDD/docs/](CFLOBDD/docs/) — papers and the SEQUITUR/WPP research
  summary that drives current work on the `sequitur` branch.
- [examples/run.sh](examples/run.sh),
  [examples/run_final_computation.sh](examples/run_final_computation.sh) —
  example invocations of the ADD/CUDD baselines.

---

## License

MIT — see [LICENSE.md](LICENSE.md).
