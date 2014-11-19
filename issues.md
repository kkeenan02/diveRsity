# Known issues in the current dev version

This document lists know problems and bugs in the current verison of `diveRsity` on gihub. As issues are fixed they will be deleted from this document.

- When calculating Weir & Cockerham's F-statistics in `diffCalc` for large data sets (>5Mb), the `C++` function `glbWCcpp` throws an error. The cause of the fault is not know yet, and the function works perfectly well for smaller data.