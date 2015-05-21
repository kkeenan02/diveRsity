# TODO

This file lists developments still to be carried out for the `diveRsity` package.

1. Speed-up allelic richness and $F_{IS}$ calculations in `divBasic`
2. Tidy up HWE legacy code in `divBasic`
3. Re-write `divBasic` on top of `rgp` type functions.
4. Allow output from `gpSampler` to be held in memory rather than to be written to file. This allows speedup since other functions like `diffCalc` will not need to read a file from disk.
5. segfault: memory not mapped for large SNP data. Cause seems to be `glbWCcpp`
6. Make sure that all relevant functions return data frames rather than matrices.
7. Add a warning to diffCalc when users specify pairwise = TRUE or bs_pairwise = TRUE, when there are only two pops in the input file.
8. Implement a plink to genepop function.
9. Fix directory bug in snp2gen function