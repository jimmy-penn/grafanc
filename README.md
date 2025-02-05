<h1 align="center">

GrafAnc

</h1>

GrafAnc is a software tool that infers genetic ancestry at both continental and subcontinental levels quickly and reliably. GrafAnc source code includes C++ programs and Perl scripts. See [GrafAnc documentation](https://github.com/jimmy-penn/grafanc/tree/master/GrafAncDocumentation.md) for software documentation.

Developer: Yumi (Jimmy) Jin, PhD. Email: [Jimmy.Jin\@Pennmedicine.upenn.edu](mailto:Jimmy.Jin@Pennmedicine.upenn.edu){.email}

### Required libraries and files

GrafAnc requires that htslib be installed on the system. Please refer to the [htslib documentation](https://www.htslib.org) for details.

GrafAnc also requires that file `AncSnpPopAFs.txt` be placed under the `data` directory. File `AncSnpPopAFs.txt.gz` needs to be unzipped:

``` sh
$ cd cpp
$ gunzip data/AncSnpPopAFs.txt.gz
```

### Make C++ binary `grafanc`

If htslib is installed under directory `~/htslib/1.20/`, GrafAnc can be built by running the following command under the root directory:

``` sh
$ make
```

To regenerate the C++ binary after editing the code, execute:

``` sh
$ make clean
$ make
```

If htslib is installed under a different directory, e.g., `/applications/htslib/1.20/`, GrafAnc can be built by executing:

``` sh
$ make HTSLIB="/applications/htslib/1.20"
```

Or alternatively, users can set `HTSLIB` in `Makefile` to the correct htslib directory and execute `make`.

### Run medium tests

Test scripts and test cases are placed under medium_testing directory. Test cases are saved in `test_manifest.txt`. Perl script `test_grafanc.pl` is used for manually running these test cases. 1. Make sure environment variable `PATH` includes current directory `.`. 1. If necessary, set environment variable `GARFPATH` to include the directory where GrafAnc binary and Perl scripts are located. 1. Under `medium_testing` directory, execute: `test_grafanc.pl test_manifest.txt`. 1. If source code is updated, update `test_manifest.txt` to add new test cases or modify existing cases, and execute `test_grafanc.pl test_manifest.txt 1` to update the baseline. 1. Check the baseline files and make sure they are all correct, then execute `test_grafanc.pl test_manifest.txt` (without the second parameter) again. 1. Make sure all test cases pass.
