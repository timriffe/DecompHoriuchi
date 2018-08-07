DecompHoriuchi
==============

R implementation of Horiuchi's proposal for continuous decomposition

Superceded by [DemoDecomp](https://github.com/timriffe/DemoDecomp)
========================

`DecompHoriuchi` will no longer be maintained. I changed the name to `DemoDecomp` because the package now contains two general methods, so it didn't make sense to name for just one of them. I also took the chance to modernize the repo a bit with automatic checks. These methods really are applicable to diverse analyses, so I'd like to demonstrate some out-of-the-box applications to give a sense.

The main changes from `DecompHoriuchi` to `DemoDecomp` are that the arguments `rates1` and `rates2` are now `pars1` and `pars2`, respectively. i.e. we're not always decomposing rates, no why not make it more general? Also the only Horiuchi function here is `DecompContinuousOrig()` (the one taking all parameters in a single vector), and it has been renamed to `horiuchi()`. Less typing that way. For now all examples the same. `DecompHoriuchi` will remain here in this state for the foreseeable future, and any new developments will happen in `DemoDecomp`.

Installation
============

To download the most recent version of LifeTable:

1. make sure you have the most recent version of R
2. look at the notes below if you're on Windows or Mac.

Download the [zip ball](https://github.com/timriffe/DecompHoriuchi/zipball/master) or [tar ball](https://github.com/timriffe/DecompHoriuchi/tarball/master), decompress and run `R CMD INSTALL` on it in the terminal command line, or use the **devtools** package to install the development version:

```r
# install.packages("devtools")

library(devtools)
install_github("timriffe/DecompHoriuchi/DecompHoriuchi")
```

**Note**: Windows users need [Rtools](http://cran.r-project.org/bin/windows/Rtools/) to install from github as shown above. Get the most recent version of [R for Windows](http://cran.r-project.org/bin/windows/base/) and download and install the version of Rtools that corresponds to it.

**Note**: Mac users might be required to install the appropriate version [XTools](https://developer.apple.com/xcode/) from the [Apple Developer site](https://developer.apple.com/) in order to install the development version.  You may need to [register as an Apple developer](https://developer.apple.com/programs/register/).  An older version of XTools may also be required.
