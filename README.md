# ShellModel.jl (v0.1.0)

[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://SotaYoshita.github.io/ShellModel.jl/stable)
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://SotaYoshita.github.io/ShellModel.jl/dev)
[![Build Status](https://travis-ci.com/SotaYoshita/ShellModel.jl.svg?branch=master)](https://travis-ci.com/SotaYoshita/ShellModel.jl)


<img src="https://github.com/SotaYoshida/ShellModel.jl/blob/main/SMjl_logo.png" width=70%>

**Julia code for nuclear shell-model calculations**  


## Quick Start
1. Download a Julia binary from [Julialang.org](https://julialang.org/)  
or execute the following when using a Linux(-like) environment:
  ```
  $wget https://julialang-s3.julialang.org/bin/linux/x64/1.6/julia-1.6.0-linux-x86_64.tar.gz
  ```

2. Set a path to Julia.
  Then you can run Julia REPL
```
$julia
               _
   _       _ _(_)_     |  Documentation: https://docs.julialang.org
  (_)     | (_) (_)    |
   _ _   _| |_  __ _   |  Type "?" for help, "]?" for Pkg help.
  | | | | | | |/ _` |  |
  | | |_| | | | (_| |  |  Version 1.6.0 (2021-03-24)
 _/ |\__'_|_|_|\__'_|  |  Official https://julialang.org/ release
|__/                   |

julia>
  ```
3. Download ShellModel.jl  
```
$git clone https://github.com/SotaYoshida/ShellModel.jl
```

4. Install other Julia Packages
```
$cd ShellModel.jl
$julia src/package_install.jl
```

5. Run the sample script
```
$julia -t 12 sample_run.jl
```
One can specify the number of execution threads like this.  
(The -t/--threads command line argument requires at least Julia >= 1.5.)  

For more details, please check docs (**not yet provided**).

> Note: 
> The samplecode is to calc.  
> (a) 10 lowest states of 28Si in sd shell  
> (b) 10 lowest states of 28Si with J=0  
> (c) EC estimates of 10 lowest J=0 states of 28Si  
> (d) (b) with the preprocessing  

## Documentation 

in prep.  
For some functions, brief explanations are available in the Julia REPL "help mode".  
One can start the "help mode" by typing "?" at the prompt.

```
julia> include("./src/shellmodel_main.jl")
main (generic function with 5 methods)

julia> ?HbitT1
search: HbitT1

  make bit representation of T=1 (proton-proton&neutron-neutron) interactions for each {m_z}
```


## Features
* Pros (Some are "Pros" of the Julia language though :D)
  * Easy to run
  * Portability (no need to specify "magical" compiler options specific to each environment)
  * Fast (e.g., 10 lowest eigenpairs of 28Si (in full sd shell) can be calculated in ~3 sec. on Mac Mini(2018) & Mac Book Air(2020,M1))
  * One can easily extend the code 
* Cons
  * poorly parallerized (for # of threads >= 12) 
  * greedy (compared to the "on-the-fly" generation of matrix element)
  * one is recommended to include "initial run" for a small system (e.g. Be8 in p shell) in execution scripts  
    This is for JIT compilation of the functions to improve the execution time for a target system.  
    You can also try [PackageCompiler.jl](https://github.com/JuliaLang/PackageCompiler.jl)
  
## What is supported now?  

* To obtain lowest *n* eigenvalues, eigenvectors (in M-scheme)  
  ※ One can specify the total J  
  ※ The current version only supports model spaces with one major shell  
  (i.e., the so-called "0 hbar omega" space such as *sd* shell)  
  ※ interaction fmt: [KSHELL](https://doi.org/10.1016/j.cpc.2019.06.011) fmt (.snt) only

* Eigenvector Continuation (To constuct approximate shell-model wave funcitons for a given effective interaction)
  as efficient emulator & preprocessing for the shell-model
 
* To calculate mu & Q moments and M1, E2 transitions (April 2021 -)

## To be implemented:  

* truncations 
* compatibility to other interaction fmt
* compatibility to wavefunctions by other shell-model codes
* J-scheme
* other operators (Gamow-Teller, etc.)
* MPI parallerization for larger systems
* and more...

Any suggestions, feedbacks, and "Issues&Pull requests" are welcomed.


## Licence  

[MIT License, Copyright (c) 2021 Sota Yoshida](https://github.com/SotaYoshida/ShellModel.jl/blob/main/LICENSE)


## How to cite  

in prep.


