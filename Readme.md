## Fast Marching Python

This repo is a python wrapper of the fast marching algorithm for computing geodesics on a triangular mesh.
The code is modified from [gproshan](https://github.com/larc/gproshan).


### Note
Fast marching is an approximate algorithm, for exact geodesics, see [gdist](https://github.com/the-virtual-brain/tvb-gdist).

### Dependencies
1. Eigen3
2. OpenMP
3. Pybind11 (already contained in this repo)


### Usage
1. Compile the code:

``` 
    cd fast_marching/fast_marching 
    mkdir build 
    cd build
    cmake -DPYTHON_EXECUTABLE:FILEPATH=python3 ..
    make 
```

2. Run the `test.py` file for demo. The input is `mesh directory` and the output is a `square matrix` recording the pairwise geodesics among all vertices.


