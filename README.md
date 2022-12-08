# HyPhy standalone analyses

This repository contains various analysis scripts for [HyPhy](https://github.com/veg/hyphy) -- a molecular evolutionary analysis platform. 

```
git clone https://github.com/veg/hyphy-analyses.git hyphy-analyses
```

Unless otherwise stated in directory - specific README file, all analyses assume the availability of the HyPhy develop build. General steps for obtaining one on Linux and OS X systems and running an `analysis-file`, assuming you have 

* C/C++ toolchains, like GCC or Clang
* [CMake](http://cmake.org) 

is as follows

```
git clone https://github.com/veg/hyphy.git hyphy-develop
cd hyphy-develop
git checkout develop
cmake ./
make -j MP

#if you want a system-wide install...
make install 

```

Now, you can run a specific standalone analysis or look for options

```
/path/to/hyphy /path/to/hyphy-analysis/module/analysis.bf --help

...

/path/to/hyphy /path/to/hyphy-analysis/module/analysis.bf --kw1 arg1 --kw2 arg2 ....

```

for example

```
cd /home/sergei/hyphy-analyses/codon-msa
hyphy pre-msa.bf --input example.fas
```


Specific analysis directories have more specific analysis options and examples.


