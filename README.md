# HomebrewLib

This is a C/C++ library covering a wide variety of functions.

It has these main features:

  1. several useful special functions including gaussian function, rectangular function, polynomials with variable degrees, etc.

  2. Structs that help construct combinations (in elementary operators) of these functions

  3. Supplementary gsl (GNU Scientific Library) functions, for example,
    a. convert the real vector/matrix to a complex one
    b. perform 'transformation' on a vector by a matrix
    c. Euclidean inner product of a vector
    d. handle gsl vector/matrix's 

There are also some functions handling quantum mechanics, currently:

  1. Perform Wigner Transformation for a equally sampled wavefunction

This library is not intended to be used by others, but sufficient annotation has been added so that you might be able to modify and use them in your own programs.

Many of my programs might use this library, so you might need to install this library first to have my programs compiled. A typical installation procedure is written in compile.sh and install.sh, which helps you compile a static library libhb.a and add this library into your /usr/local/lib directory along with the include files into your /usr/local/include/hb. You might add, for example,

#include <hb/Quantum.h>

to use the library in your program. Modifications are needed if you install this library into other directories.
