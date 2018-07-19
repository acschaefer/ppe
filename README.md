# Gpufit_ppe

Modified version of [Gpufit](https://github.com/gpufit/Gpufit) to accommodate plane to ray fitting for probabilistic plane extraction ([PPE](link)) from 3D laser scans.

## Requirements

Successfully tested using the following:

* CMake 3.11
* CUDA Toolkit 8.0
* Matlab R2017b

Newer and some older versions might work as well.

## Quick start instructions

On Linux (Ubuntu 16.04 LTS):

```bash
git clone git@aisgit.informatik.uni-freiburg.de:buescher/gpufit_ppe.git
cd gpufit_ppe/
mkdir build
cd build/
cmake -DCMAKE_BUILD_TYPE=RELEASE -DCUDA_ARCHITECTURES=6.1 ..
make -j8
cd ../Gpufit/matlab/examples/
export LD_PRELOAD=/usr/lib/x86_64-linux-gnu/libstdc++.so.6
matlab
```

In matlab run `planetorays` for an example. Alternatively, load the PPE directory, modify the last line of `extrpln.m` to use gpu fitting and run `extrpln`.
