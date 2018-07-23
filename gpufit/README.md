# Gpufit for PPE

Modified version of Gpufit to accommodate plane to ray fitting for probabilistic plane extraction (PPE) from 3D laser scans.

Please refer to the Gpufit documententation for further installation details: [github.com/gpufit/Gpufit](https://github.com/gpufit/Gpufit)

## Requirements

Successfully tested using the following:

* CMake 3.11
* CUDA Toolkit 8.0, 9.2
* Matlab R2017b, R2018a
* Visual Studio Community 2017 (Windows only)

Newer and some older versions might work as well.

## Quick start instructions

On Linux (Ubuntu 16.04 LTS):

```bash
cd gpufit/
mkdir build
cd build/
cmake -DCMAKE_BUILD_TYPE=RELEASE -DCUDA_ARCHITECTURES=6.1 ..
make -j8
cd ../Gpufit/matlab/examples/
export LD_PRELOAD=/usr/lib/x86_64-linux-gnu/libstdc++.so.6
matlab
```

On Windows 10 (details can be found [here](https://github.com/gpufit/Gpufit/blob/master/docs/installation.rst#compiler-configuration-via-cmake)):

* Run CMake in a build directory, e.g. `cmake -G "Visual Studio 12 2013 Win64" C:\Sources\ppe\gpufit`
* Open the created Solution file with Visual Studio and build it

In matlab run `planetorays` for an example. Alternatively, load the PPE directory, modify the last line of `extrpln.m` to use gpu fitting and run `extrpln`.
