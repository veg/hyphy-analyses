# HyPhy standalone analyses

This repository contains various analysis scripts for [HyPhy](https://github.com/veg/hyphy) -- a molecular evolutionary analysis platform. Unless otherwise stated in directory - specific README file, all analyses assume the availability of the HyPhy develop build. General steps for obtaining one on Linux and OS X systems and running an `analysis-file`, assuming you have 

* C/C++ toolchains, like GCC or Clang
* [CMake](http://cmake.org) 

is as follows

```
git clone https://github.com/veg/hyphy.git hyphy-develop
cd hyphy-develop
git checkout develop
cmake ./
make -j MP
./HYPHYMP analysis-file --help
```

Specific analysis directories have more specific analysis options and examples.
