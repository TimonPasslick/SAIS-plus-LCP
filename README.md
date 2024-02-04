# SAIS-plus-LCP

This is for a university assignment.

## How to build

I put everything in main.cpp and use no external libraries, so it should be straightforward.

If you want to enable exact memory peak measurements for Linux, uncomment `#define LINUX_MEMORY_PEAK`. I didn't compile or test this code snippet because I'm on Windows.

If you want to enable estimated, but platform independent memory peak measurements, set `estimate_memory_peak` to true.
