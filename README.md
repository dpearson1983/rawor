# rawor
**RA**ndoms **W**ith**O**ut **R**andoms

This repository will house our rawor code as so called header only libraries for C++, C, and Python. These libraries will
allow one to predict the three point counts involving a random distribution of points without explicity counting the *DDR*,
*DRR* or *RRR* terms. Since the sets of random points used for correlation functions are typically many times larger than
the sets of data points to reduce additional shot noise, predicting these counts should be considerably faster than 
performing the counts directly.

## In what situations will this code be useful?
Currently, the code has only been verified for counts in periodic boxes, such as those found in *N*-body simulations or
various second order lagrangian perturbation theory or fast particle mesh simulations. This means that the code cannot be
used for the analysis of survey data at present, though work to extend the codes capabilities to surveys is currently 
underway.

## Are there any external dependencies?
In an effort to make the use of the code as simple as possible, external dependencies have been largely avoided. The C++ 
implementation only depends on the C++ standard library, meaning that as long as one has a C++ compiler installed on their
system, the code should just work. The only caveat for the C++ code right now is that if you are working in Windows, you must 
pass `_USE_MATH_DEFINES` to the preproccessor as the code uses `M_PI` for the value of &pi; (the forthcomming documentation
will have more details). The C code will also strive to remain within the C standard library.

While the original plan was to have "pure python" port of the library, this was changed to creating a Python 
extension module. In order to wrap the C++ code for Python, we used the Boost Python libraries, which will need
to be installed on your system.

## Are there any operating system dependencies?
Aside from passing `_USE_MATH_DEFINES` to the preproccessor when compiling in Windows, there should be no other operating
system dependencies to worry about. However, we note that as of this writing the C++ code has only been tested on Fedora 
GNU/Linux versions 28 and 29 (g++ (GCC) 8.2.1), and Windows 10 (Microsoft(R) C/C++ Optimizing Compiler Version 19.00.24234.1 for x86). We don't anticipate any difficulties compiling on other operating systems
since only the C++ standard library is used, and we stick to variable types that should be defined the same on all operating
systems, e.g we stick with `int`, and `double`.

## How do I use the library?
For C++, place rawor.hpp in your project directory (i.e. where your code that will use the library is located) then simply add
```c++
#include "rawor.hpp"
```
to your list of include statements to make the rawor class accessible in your code. Create a class object with the needed information for the computations
```c++
rawor ranPredictor(N_data, N_rans, N_shells, V_box, r_max, r_min)
```
then simply get the RRR, DRR, and DDR counts
```c++
std::vector<int> RRR = ranPredictor.getRRR();
std::vector<int> DRR = ranPredictor.getDRR();
std::vector<int> DDR = ranPredictor.getDDR(DD);
```
where `DD` are your pair counts in the `N_shells` bins between `r_min` and `r_max`.

# TO DO
This is a list of all the changes to this repository that we are currently working on.
- Complete work on C header
- ~~Complete work on Python class~~
- ~~Comment the C++ header~~
- Create documentation
