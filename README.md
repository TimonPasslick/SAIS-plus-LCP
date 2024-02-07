# SAIS-plus-LCP

This is for a university assignment. Constructs suffix arrays with the SAIS algorithm and LCP arrays with the naive, Kasai and Phi algorithms.

## How to build

I put everything in main.cpp and use no external libraries, so it should be straightforward, on Linux probably:
```
g++ -O3 -g main.cpp -o main
```

If you want to enable exact memory peak measurements for Linux, uncomment `#define LINUX_MEMORY_PEAK`. I didn't compile or test this code snippet because I'm on Windows.

If you want to enable estimated, but platform independent memory peak measurements, set `estimate_memory_peak` to true.

If you want to check for correctness with the naive LCP algorithm, set `check_for_correctness` to true.
