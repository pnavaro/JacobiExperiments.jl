# JacobiExperiments.jl

Just a Julia version of [this cool fortran project](https://github.com/loiseaujc/Jacobi-Experiments)!

It can be noted that the optimized versions have similar performance, but Julia language suffers more from the lack of optimization.

## Instructions for running the Julia program

```bash
git clone https://github.com/pnavaro/JacobiExperiments.jl.git
cd JacobiExperiments.jl/
julia --project
```
               _
```juliarepl
julia> import Pkg; Pkg.instantiate()
julia> import .Threads; Threads.nthreads()
8
julia> include("example/main.jl")
Textbook solver :
    - Number of iterations : 128395
    - l2-norm of the error : 1.4900285278688397e-8
 - Time-to-solution     : 168.993445636
No-copy solver :
    - Number of iterations : 128396
    - l2-norm of the error : 1.4899158919120579e-8
 - Time-to-solution     : 73.2647536
On-the-fly solver :
    - Number of iterations : 129002
    - l2-norm of the error : 1.4232008958046707e-8
 - Time-to-solution     : 30.699903297
Do-concurrent solver :
    - Number of iterations : 129002
    - l2-norm of the error : 1.4232008958046707e-8
 - Time-to-solution     : 20.807622978
```

## Instructions for running the Fortran program 

```
conda install fpm
fpm run --flag "-O3 -mtune=native -march=native"
<WARNING> both openmp and stdlib requested: some functions may not be thread-safe!
Project is up to date
 Textbook solver :
     - Number of iterations :      128395
     - l2-norm of the error :   1.4900285278657732E-008
     - Time-to-solution     :   89.871690439060330
 No-copy solver :
     - Number of iterations :      128396
     - l2-norm of the error :   1.4899158919089789E-008
     - Time-to-solution     :   49.446893895044923
 On-the-fly solver :
     - Number of iterations :      129002
     - l2-norm of the error :   1.4232008958017304E-008
     - Time-to-solution     :   30.556735919788480
 Do-concurrent solver :
     - Number of iterations :      129002
     - l2-norm of the error :   1.4232008958017304E-008
     - Time-to-solution     :   30.745370151475072
```
The concurrent version didi not work. Probably my `gfortran` version is too old
```
gcc version 11.4.0 (Ubuntu 11.4.0-1ubuntu1~22.04)
```
