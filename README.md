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

## Why is there only a C++ implementation?
The first, and currently only, implementation is in C++. This is largely due to the author having vastly more experience
coding in C++ than in Python or C. Given the comparitively minor changes that will be needed, the C version should be 
forthcoming. Since the Python version will need to make use of vectorized numpy operations to maintain an acceptable
performance level, the alterations to the algorithm will be more substantial, requiring more time to complete.

## Why release the code as header only libraries?
This code is intended to be used as a small part of a larger piece of software to compute the three point correlation function.
The *DDR* count prediction relies on the user providing an array of *DD* counts from which the number density of *DD* pairs
can be estimated. As such, having the code as a library seemed the natural choice.

As far as using libraries is concerned, header only libraries tend to be easier, as the user simply needs to put a copy of the
header file in their project folder and then include/import the file into their projects main source file. This removes the
often complicated setup process for libraries that must be separately compiled, e.g. no fiddling with configurations or
makefiles.

## Are there any external dependencies?
In an effort to make the use of the code as simple as possible, external dependicies have been largely avoided. The C++ 
implementation only depends on the C++ standard library, meaning that as long as one has a C++ compiler installed on their
system, the code should just work. The only caveat for the C++ code right now is that if you are working in Windows, you must 
pass `_USE_MATH_DEFINES` to the preproccessor as the code uses `M_PI` for the value of &pi; (the forthcomming documentation
will have more details).

The C code will also strive to remain within the C standard library, and the Python code will stick to widely used Python
libraries (in its present state, only numpy is required).

## Are there any operating system dependencies?
Aside from passing `_USE_MATH_DEFINES` to the preproccessor when compiling in Windows, there should be no other operating
system dependencies to worry about. However, we note that as of this writing the C++ code has only been tested on Fedora 
GNU/Linux versions 28 and 29 (g++ (GCC) 8.2.1), and Windows 10 (Microsoft(R) C/C++ Optimizing Compiler Version 19.00.24234.1 for x86). We don't anticipate any difficulties compiling on other operating systems
since only the C++ standard library is used, and we stick to variable types that should be defined the same on all operating
systems, e.g we stick with `int`, and `double`.

# TO DO
This is a list of all the changes to this repository that we are currently working on.
- Complete work on C header
- Complete work on Python class
- ~~Comment the C++ header~~
- Create documentation
