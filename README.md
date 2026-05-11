# libfci

`libfci` is a compact determinant full-CI solver library developed by Gerald
Knizia, 2010-2012, and packaged as a small C ABI for active-space embedding
workflows.

The original standalone driver accepted FCIDUMP-style Hamiltonians. The library
interface bypasses file I/O and accepts dense active-space one- and
two-electron matrices using CP2K's FCIDUMP conventions.

The exported C entry point is `libfci_solve`. A compatibility wrapper named
`cp2k_fci_solve` is also exported for older CP2K integration experiments.

## Build

```sh
cmake -S . -B build -DCMAKE_INSTALL_PREFIX=/path/to/install
cmake --build build
cmake --install build
```

OpenMP is enabled automatically when available. It can be disabled with
`-DLIBFCI_USE_OPENMP=OFF`.
