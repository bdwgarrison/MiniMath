
#include <iostream>
#include "menus.h"
#include "functions.h"
#include "algorithms.h"

/***************************************************************************************
*********************************** OUTLINE OF PROGRAM *********************************
****************************************************************************************/

/***

@mainpage Numerical Differentiation, Integration, and Interpolation Program
@author Brenden Garrison
@date May 22, 2013
@version 1.0

@brief Methods for approximating solutions to various mathematical problems.

This program runs entirely inside of "menu" objects, which can be found declared in "menus.h"
and defined in "menus.cpp"

The purpose of such implementation is to allow for total recursion of the different stages of the project.
The entire program can be divided into two sections:
    1) Creation of a system of functions.
    2) Applying various algorithms to those functions.

The menu implementation is totally optional and can be commented out entirely, allowing the programmer
to utilize the algorithms in whatever way makes sense to their purpose.

Immediately below is a description of the specific implementation using the menus. Below that are examples
of how one could utilize the algorithms without. (they're fairly intuitive)

@section MAIN

This stage only lasts for a second, displaying a brief warning before loading into the "interface" menu

@section INTERFACE

User is asked to decide how many functions will be in the system.

@section DEFINEFUNCS

Loop scrolls through each position in the system, asking user where to get the function.
    1) Library
    2) Interpolate
    3) Differentiate
    4) Integrate

@section LIBRARY

Asks user to pick from a function of predefined list.
To add a new function to this list, there are four quick steps:
    1) In "functions.h", at the end of the file there is a namespace
        verbosely titled "library_of_predefined_simple_functions"
        declare a raw function:

            double myfunc_f (ntuple);   //< note the _f to distinguish from the "function object"

    2) Then, beneath it, in the same namespace, instantiate the function object:

            extern function myfunc;

    3) Go to "functions.cpp" and again, all the way at the bottom are the
        namespace functions with the alias "flb". Define your raw function:

            double myfunc_f(ntuple x) { return ... }

    4) Similarly, beneath it, define the function object:

            function myfunc( myfunc_f );

    5) Extra: If you want the function to be accessible by the interface:
        In "menus.cpp" find "void mnu::definefuncs::funclib()", this is the menu where the user
        picks from the library. Add in a description with the rest and a new assign block

        else if (selection == mynum) { funcs.assign_func(i, flb::myfunc); }

@section INTERPOLATE

User picks a file to open and read data from. All methods require at least two columns: the first for the
independent variable values and the second for the dependent. One method requires a third column for the
derivative, another just the value of the derivative at the first point and last. User picks from four
methods of interpolation.

    1) Basic Polynomial: This constructs a Lagrangian Polynomial using Newton's Divided Difference method.
                            The function is continuous for all reals, but as soon as the value leaves the
                            domain of the data it was interpolated over, the error grows exponentially.
                            This works better for differentiation techniques over that same domain because
                            it requires evaluating the function at f(x-h) + f(x+h) so if your function is
                            not defined before a small shift outside the domain (as with a piecewise function
                            from another method) then the method will crash.

    2) Derv Polynomial: This constructs a Hermite Polynomial using a modification of Neville's method.
                        The function is again continuous for all reals but the error grows rappidly outside
                        the interpolated domain.

    3) Basic Piecewise: Creates a Natural Cubic Spline, ie: a piecewise function that has a different 3rd
                        degree polynomial on each interval. The function is continuous and differentiable
                        on its domain, but undefined elsewhere. This construction assumes the second
                        derivative is zero at the endpoints.

    4) Derv Piecewise: Creates a Clamped Cubic Spline: Same as above only with the derivative fixed at
                        the endpoints.

@section DIFFERENTIATE

Approximates the derivative at many points, then interpolates them all together with Basic Polynomial.

@section INTEGRATE

Approximates solution to function as an initial value problem differential equation (IVP ODE), then
interpolates the solution together as Derv Polynomial.

@section SOLVEFUNCS

This brings us to stage 2, performing methods on the functions. User is asked wether they'd like to
perform methods on the system as a whole or a single function.


@section SYSTEM

Interprets system as a vector valued function, as each of its members goes from Rn->R1,
making a system of m functions go from Rn->Rm. Here we can select from four methods.

    1) Eval: Simply evaluates each function at an inputted vector and returns the solution as one.
    2) Lin: Solve system as a matrix Ax = b given any b. (only works if all functions are linear)
    3) Nonlin: Solve system F(x) = 0 using "Continuation" method for x.
    4) IVP: Solve system of IVP's using Runge Kutta method on systems.

@section SINGLE

Select one function to perform various numerical methods on.

    1) Eval: Evaluates function at user input.
    2) Root: Pick from various techniques to solve for the roots of the equation.
    3) Diff: Numerically differentiate the function at many points.
    4) Integ: Evaluates definite integral over user input domain.
    5) IVP: Pick from various techniques to solve function as initial value problem.

*/

/**
@remark All examples performed in MAIN.

@example Single Evaluation

	fnc::function myfunc = flb::func6;

	std::cout<< myfunc(1)<<std::endl<< myfunc(1.3)<<std::endl;

@example Differentiation (continues from previous)

	alg::differentiate diff(&myfunc);

	fnc::polynomial derivative = diff.poly();

	std::cout<< derivative(1.5)<<std::endl;

	fnc::pwfunction mypiecewise = diff.piece();

	std::cout<< mypiecewise(2.5)<<std::endl;


@example Integration (continues from previous)

	alg::solve_ivp_fix  ivp( &derivative ); //< here, we integrate our derivative function back again

	fnc::polynomial integral = ivp.integral_poly();

	std::cout<< integral(1.3)<<std::endl;   //< check value with before and see the error


@example Interpolation

    fnc::matrix my_matrix;
    my_matrix.get_nodes();  //< need to specify a file

    alg::interp_basicpoly  interp( my_matrix );

    fnc::polynomial poly = interp.getsoln();

    std::cout<< poly(1.5)<<std::endl;

@example Interpolation (manual entry of file)

    std::ifstream infile;
    infile.open("newtondata.txt");

    fnc::matrix my_matrix(infile);

    alg::interp_basicpoly  interp( my_matrix );

    fnc::polynomial poly = interp.getsoln();

    std::cout<< poly(1.5)<<std::endl;

@example Solve Matrix

    std::ifstream infile;
    infile.open("my_data.txt");

    fnc::matrix my_matrix(infile);

    fnc::ntuple x(2);

    x[0] = 1;
    x[1] = 7;

    alg::solve_linsys solv( my_matrix, x );   //< expects the matrix to be a 2x2

    std::cout<< solv.getsoln();



*/




/********************************************* MAIN **********************************************/

namespace mnu = menus_for_implementing_numerical_analysis_methods;
namespace fnc = functions_and_mathematical_objects_for_numerical_analysis;
namespace flb = fnc::library_of_predefined_simple_functions;
namespace alg = algorithms_for_numerical_analysis_implementation_on_function_objects;

int main()
{






/*
    std::cout<< "This program performs numerical analysis algorithms on either a"<<std::endl
    <<"single or system of functions."<<std::endl
    <<"\"exit\" -- Ends the program if user is at beginning stage (here)."<<std::endl
    <<"\"back\" -- Returns to the previous stage."<<std::endl
    <<"\"help\" -- Displays a discription of all the commands listed at present stage.\n"<<std::endl
    <<std::endl<<std::endl<<std::endl
    <<"         TYPE HELP IF YOU ARE STUCK"<<std::endl;

    ///This program runs completely inside functions to allow repitition of steps.

    mnu:: interface nextstage;
    nextstage.run();
*/
    return 0;
}

