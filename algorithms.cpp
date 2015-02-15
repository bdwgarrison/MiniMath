


#include <iostream>
#include <cmath>
#include <iomanip>
#include <sstream>

#include "functions.h"
#include "algorithms.h"

namespace fnc = functions_and_mathematical_objects_for_numerical_analysis;

namespace algorithms_for_numerical_analysis_implementation_on_function_objects {

/** the helper functions immediately below use a string stream to avoid a system crash **/
void algorithm::setdomain(double* indom)
{
    std::cout<< "On what closed interval to evaluate? Enter lower bound first and upper bound second.\n";
    std::string input;

    std::cin>> input;
    std::stringstream stream(input);
    stream>> indom[0];

    std::cin>> input;
    std::stringstream stream2(input);

    stream2>> indom[1];
}
void algorithm::setdervdomain(double* indom)
{
    std::cout<< "On what closed interval to differentiate? Enter lower bound first and upper bound second.\n";
    std::string input;

    std::cin>> input;
    std::stringstream stream(input);
    stream>> indom[0];

    std::cin>> input;

    std::stringstream stream2(input);

    stream2>> indom[1];
}
void algorithm::setintegdomain(double* indom)
{
    std::cout<< "On what closed interval to integrate? Enter lower bound first and upper bound second.\n";
    std::string input;

    std::cin>> input;
    std::stringstream stream(input);
    stream>> indom[0];

    std::cin>> input;
    std::stringstream stream2(input);

    stream2>> indom[1];
}
void algorithm::setnodecount(long unsigned* innodecount)
{
    std::cout<< "Calculate at how many points within the domain? (equally spaced)\n";
    std::string input;

    std::cin>> input;
    std::stringstream stream(input);
    stream>> *innodecount;
}
void algorithm::setmeshcount(long unsigned* inmeshcount)
{
    std::cout<< "Use how many mesh points while calculating integral?\n";
    std::string input;

    std::cin>> input;
    std::stringstream stream(input);
    stream>> *inmeshcount;
}
void algorithm::setord(short* inord)
{
    std::cout<< "Of what order would you like your approximations?\n";
    std::string input;

    std::cin>> input;
    std::stringstream stream(input);
    stream>> *inord;
}
void algorithm::settol(double* intol)
{
    std::cout<< "Enter the tolerance of error, epsilon: ";
    std::string input;

    std::cin>> input;
    std::stringstream stream(input);
    stream>> *intol;
}
void algorithm::setmaxit(long unsigned * inmaxit)
{
    std::cout<< "Enter the maximum allowed iterations: ";
    std::string input;

    std::cin>> input;
    std::stringstream stream(input);
    stream>> *inmaxit;
}
void algorithm::setinit1(double* init1)
{
    std::cout<< "First approximation in iterative sequence: ";
    std::string input;

    std::cin>> input;
    std::stringstream stream(input);
    stream>> *init1;
}
void algorithm::setinit2(double* init2)
{
    std::cout<< "Second approximation in iterative sequence: ";
    std::string input;

    std::cin>> input;
    std::stringstream stream(input);
    stream>> *init2;
}
void algorithm::setivdep(double* inivdep)
{
    std::cout<< "Initial value of dependent variable: ";
    std::string input;

    std::cin>> input;
    std::stringstream stream(input);
    stream>> *inivdep;
}
void algorithm::setivderv(double* inivderv)
{
    std::cout<< "Initial value of derivative: ";
    std::string input;

    std::cin>> input;
    std::stringstream stream(input);
    stream>> *inivderv;
}
void algorithm::sethmin(double* inhmin)
{
    std::cout<< "Minimum step size: ";
    std::string input;

    std::cin>> input;
    std::stringstream stream(input);
    stream>> *inhmin;
}
void algorithm::sethmax(double* inhmax)
{
    std::cout<< "Maximum step size: ";
    std::string input;

    std::cin>> input;
    std::stringstream stream(input);
    stream>> *inhmax;
}
void algorithm::setstep(double* step)
{
	std::cout<< "Step size: ";
    std::string input;

    std::cin>> input;
    std::stringstream stream(input);
    stream>> *step;
}
void algorithm::setboundary(double* boundary)
{
    std::cout<<"Enter boundary conditions:"<<std::endl;
    std::string input;

    std::cin>> input;
    std::stringstream stream(input);
    stream>> boundary[0];

    std::cin>> input;
    std::stringstream stream2(input);
    stream2>> boundary[0];
}


/*************************************************************************************
************************************ 1D Solvers **************************************
*************************************************************************************/


/** BISECTION METHOD
@brief bisection method finds the root of an equation on an interval such that f has opposite signs at the endpoints.

@example  interval [2,5] such that f(2) = -1 and f(5) = 3

@param f pointer to function object you'd like to solve.
@param ina the lower bound of the interval to perform bisection algorithm on
@param inb the upper bound
@param inmaxit the maximum number of iterations to perform before returning an unsuccessful attempt. typical: 150
@param intol minimum difference between successive iterations to conclude convergence, typical .000000001

@return the root of the function. returns 0 and error message if the method failed to converge.
*/
double bisection::operator() (fnc::function* f, double ina, double inb, long unsigned inmaxit, double intol)
{
    double a = ina, b = inb;
	maxit = inmaxit;
	tol = intol;

    long unsigned n=0;

    if ( (*f)(a) * (*f)(b) == 0)
    {
        if ( (*f)(a)==0 )
            return a;
        if ( (*f)(b)==0 )
            return b;
    }
    else if ( (*f)(a) * (*f)(b) > 0 )
    {
        std::cout<< "Invalid interval specified for bisection1D. Function must have opposite sign at endpoints."<<std::endl;
        return 0;
    }
    else
    {
        while(b-a > tol && n<= maxit)
        {
            n++;
            root = (a+b)/2;
            if ( (*f)(a) * (*f)(root) < 0 ) b = root;
            else if ( (*f)(root) * (*f)(b) < 0) a = root;
            else
                return root;
        }
        if (n > maxit)
        {
            std::cout<< "Method failed to converge in "<< maxit<< "steps."<<std::endl;
            return 0;
        }
        else
            return root;
    }
}

/** get info **/
void bisection::getinfo()
{
    setdomain(domain);
    setmaxit(&maxit);
    settol(&tol);
}
/** user interface version **/
void bisection::operator()( fnc::function* f)
{
	getinfo();
    double a = domain[0], b = domain[1];


    long unsigned n=0;

    if ( (*f)(a) * (*f)(b) == 0)
    {
        if ( (*f)(a)==0 )
        {
            std::cout<< "Root at endpoint, "<<a<<"."<<std::endl;
            return;
        }
        if ( (*f)(b)==0 )
        {
            std::cout<< "Root at endpoint, "<<b<<"."<<std::endl;
            return;
        }
    }
    else if ( (*f)(a) * (*f)(b) > 0 )
    {
        std::cout<< "Invalid interval specified for bisection1D. Function must have opposite sign at endpoints."<<std::endl;
        return;
    }
    else
    {
        while(b-a > tol && n<= maxit)
        {
            n++;
            root = (a+b)/2;
            if ( (*f)(a) * (*f)(root) < 0 ) b = root;
            else if ( (*f)(root) * (*f)(b) < 0) a = root;
            else
            {
                std::cout<< "Function converged to "<< root<<" in "<< n<< " steps."<<std::endl;
                return;
            }
        }
        if (n > maxit)
            std::cout<< "Method failed to converge in "<< maxit<< "steps."<<std::endl;
        else
            std::cout<< "Method converged to "<< root<<" in "<< n<< " steps."<<std::endl;
    }
}






/** NEWTON'S METHOD
@brief Newton's method calculates the root of the function using a tangent line approximation per iteration.
@remark Requires an initial approximation to the root.

@param f pointer to function object you'd like to solve
@param df pointer to that function's derivative
@param inmaxit maximum number of iterations to perform before quitting. typical: 150
@param intol minimum difference between iterations to conclude that the function has converged. typical: .000000001
@param init initial approximation to the root

@return the root of the function. returns 0 and error message if failed to converge

*/

double newton::operator()(fnc::function* f, fnc::function* df, long unsigned inmaxit, double intol, double init)
{
    maxit = inmaxit;
    tol = intol;
    root = init;
    long unsigned n=0;


    while (fabs( (*f)(root)/ (*df)(root) ) > tol && n<=maxit)
    {
        n++;
        root-= (*f)(root)/ (*df)(root);
    }
    if (n< maxit)
        return root;
    else
    {
        std::cout<< "Function failed to converge in "<< maxit<< " steps."<<std::endl;
        return 0;
    }
}

void newton::getinfo()
{
    setmaxit(&maxit);
    settol(&tol);
    setinit1(&root);
}

/** user interface version **/
void newton::operator()(fnc::function* f, fnc::function* df)
{
    getinfo();
    long unsigned n=0;


    while (fabs( (*f)(root)/ (*df)(root) ) > tol && n<=maxit)
    {
        n++;
        root-= (*f)(root)/ (*df)(root);
    }
    if (n< maxit)
        std::cout<< "Function converged to "<< root<<" in "<< n<< " steps."<<std::endl;
    else
        std::cout<< "Function failed to converge in "<< maxit<< " steps."<<std::endl;
}


/** SECANT METHOD
@brief secant method solves for root of function with secant line approximations at each step.
@remark requires two initial values

@param f pointer to function object for evaluating
@param inmaxit maximum iterations allowed before method quits, typical: 150
@param intol minimum difference between successive iterations, typical: .000000001
@param init1 first initial approximation
@param init2 second initial approximation

@return root either the root of the equation or 0 and an error message if method failed to converge.
*/

double secant::operator() (fnc::function* f, long unsigned inmaxit, double intol, double init1, double init2)
{
    maxit = inmaxit;
    tol = intol;
    root = init1;
    init = init2;

    long unsigned n=0;
    while (fabs( (*f)(root)*(root-init)/( (*f)(root) - (*f)(init)) ) > tol && n<=maxit)
    {
        n++;
        root -=  (*f)(root)*(root-init)/( (*f)(root)- (*f)(init));
    }
    if (n<maxit)
    {
        return root;
    }
    else
    {
        std::cout<< "Function failed to converge in "<< maxit<< " steps.\n";
        return 0;
    }
}

void secant::getinfo()
{
    setmaxit(&maxit);
    settol(&tol);
    setinit1(&root);
    setinit2(&init);
}
void secant::operator() (fnc::function* f)
{
    getinfo();

    long unsigned n=0;
    while (fabs( (*f)(root)*(root-init)/( (*f)(root) - (*f)(init)) ) > tol && n<=maxit)
    {
        n++;
        root -=  (*f)(root)*(root-init)/( (*f)(root)- (*f)(init));
    }
    if (n<maxit)
        std::cout<< "Function converged to "<< root<<" in "<< n<< " steps.\n";
    else
        std::cout<< "Function failed to converge in "<< maxit<< " steps.\n";
}


/** RAGULA FALSI
@brief regula falsi uses a combination of bisection and secant method to find the root of a function
@remark requires interval on which f has different signs at the endpoints
@example interval [3,7] such that f(3) = -2 and f(7) = 4

@param f pointer to function object to find root of
@param ina lower bound of interval to perform method on
@param inb upper bound
@param inmaxit maximum allowed iterations before quitting, typical: 150
@param intol minimum difference between successive iterations, typical: .000000001

@return root if success or 0 and error message if failed to converge
*/

double regulafalsi::operator()(fnc::function* f, double ina, double inb, long unsigned inmaxit, double intol)
{
    double p0= ina, p1= inb;
    tol = intol;
    maxit = inmaxit;

    if ( (*f)(p1) * (*f)(p0) > 0)
    {
        std::cout<< "Invalid interval for Regula Falsi method. \nFunction must have opposite signs at each endpoint."<<std::endl;
        return 0;
    }
    else if ( (*f)(p1) * (*f)(p0) == 0)
    {
        if ( (*f)(p1) == 0 ) return p1;
        else return p0;
    }
    else
    {
        int unsigned n=0;
        root = p1 -  (*f)(p1)*(p1-p0)/( (*f)(p1)- (*f)(p0));

        while (fabs( (*f)(p1)*(p1-p0)/( (*f)(p1)- (*f)(p0)) ) > tol && n<maxit)
        {
            n++;
            //root-=  (*f)(p1)*(p1-p0)/( (*f)(p1)- (*f)(p0));
            root = root - (*f)(p1)*(p1-p0)/( (*f)(p1)- (*f)(p0));
            if ( (*f)(root) *  (*f)(p1) < 0) p0=p1;
            p1=root;
        }
        if (n==maxit)
        {
            std::cout<< "Function failed to converge within maximum, "<< maxit<< ", iterations. But reached "<< root<<std::endl;
            return 0;
        }
    return root;
    }
}


void regulafalsi::getinfo()
{
    setdomain(domain);
    setmaxit(&maxit);
    settol(&tol);
}
void regulafalsi::operator()(fnc::function* f)
{
    getinfo();
    double p0= domain[0], p1= domain[1];
    if ( (*f)(p1) * (*f)(p0) > 0)
    {
        std::cout<< "Invalid interval for Regula Falsi method. \nFunction must have opposite signs at each endpoint."<<std::endl;
        return;
    }
    else if ( (*f)(p1) * (*f)(p0) == 0)
    {
        std::cout<< "One starting point is already an intercept."<<std::endl;
        return;
    }
    else
    {
        int unsigned n=0;
        root = p1 -  (*f)(p1)*(p1-p0)/( (*f)(p1)- (*f)(p0));

        while (fabs( (*f)(p1)*(p1-p0)/( (*f)(p1)- (*f)(p0)) ) > tol && n<maxit)
        {
            n++;
            //root-=  (*f)(p1)*(p1-p0)/( (*f)(p1)- (*f)(p0));
            root = root - (*f)(p1)*(p1-p0)/( (*f)(p1)- (*f)(p0));
            if ( (*f)(root) *  (*f)(p1) < 0) p0=p1;
            p1=root;
        }
        if (n==maxit)
        {
            std::cout<< "Function failed to converge within maximum, "<< maxit<< ", iterations. But reached "<< root<<std::endl;
            return;
        }
    std::cout<< "Function converged to "<< root<<" in "<< n<< " steps."<<std::endl;
    return;
    }
}




/**********************************************************************************
******************************** INTERPOLATING METHODS ****************************
**********************************************************************************/

/** INTERPOLATE BASIC POLYNOMIAL

@brief changes a matrix of doubles into a polynomial
@remark requires a matrix of two columns. the first for the independent variable and second for dependent.
        will ignore any extra columns
@example

    fnc::matrix my_data;

    my_data.getnodes(); //< gets entries from a file

    alg::interp_basicpoly interpolated_vals( my_data );

    fnc::polynomial my_poly = interpolated_vals.getsoln();

    std::cout<< my_poly(1.5); //< outputs polynomial evaluated at 1.5

@param mtx matrix to be interpolated

*/
interp_basicpoly:: interp_basicpoly(const fnc::matrix& mtx)
{
    n = mtx.numrows();
	F = fnc::matrix(n,n);

    for (long unsigned i = 0; i<n; i++) F[i][0] = mtx[i][1];

    for (long unsigned i = 1; i<n; i++)
        for (long unsigned j = 1; j<=i; j++)
            F[i][j] = (F[i][j-1] - F[i-1][j-1]) / (mtx[i][0] - mtx[i-j][0]);

	soln = fnc::polynomial(n);

    for (long unsigned i=0; i<n; i++)
    {
        soln.cof(i) = F[i][i];
        soln.ctr(i) = mtx[i][0];
    }
}


/** INTERPOLATE POLYNOMIAL USING DERIVATIVE VALUES
@brief changes matrix into polynomial
@remark requires a matrix of 3 columns. first for the independent variable, second for the dependent, third for the derivative

@example see above
*/
interp_dervpoly::interp_dervpoly(const fnc::matrix& mtx)
{
    n = mtx.numrows();
    n2 = 2*n;
	Q = fnc::matrix(n2,n2);

	soln = fnc::polynomial(n2);

    for (long unsigned i = 0; i < n; i++)
    {
        soln.ctr(2*i) = mtx[i][0];
        soln.ctr(2*i+1) = mtx[i][0];

        Q[2*i][0] = mtx[i][1];
        Q[2*i+1][0] = mtx[i][1];
        Q[2*i+1][1] = mtx[i][2];

        if (i != 0) Q[2*i][1] = ( Q[2*i][0] -  Q[2*i-1][0] )/( soln.ctr(2*i) - soln.ctr(2*i-1) );
    }
    for (long unsigned i = 2; i < n2; i++)
        for (long unsigned j = 2; j<= i; j++)
            Q[i][j] = ( Q[i][j-1] - Q[i-1][j-1] )/( soln.ctr(i) - soln.ctr(i-j) );

    for (long unsigned i = 0; i < n2; i++)
        soln.cof(i) = Q[i][i];
}


/** INTERPOLATE PIECEWISE FUNCTION
@brief changes matrix into piecewise function
@remark requires two columned matrix. disregards more.
@example see above

*/
interp_basicpw::interp_basicpw(const fnc::matrix& mtx)
{
    n = mtx.numrows();
    soln = fnc::pwfunction(n-1);
    h = fnc::ntuple(n-1);

    for (long unsigned i = 0; i<=n-2; i++) h[i] = mtx[i+1][0] - mtx[i][0];

    fnc::ntuple alpha(n);
	alpha[0]=0;

    for (long unsigned i = 1; i<=n-2; i++) alpha[i] = 3/h[i]*( mtx[i+1][1] - mtx[i][1] ) - 3/h[i-1]*( mtx[i][1] - mtx[i-1][1] );

	fnc::ntuple l(n);
	fnc::ntuple u(n-1);
	fnc::ntuple z(n);
	l[0] = 1;
	l[n-1] = 1;



    for (long unsigned i=1; i<=n-2; i++)
    {
        l[i] = 2*(mtx[i+1][0] - mtx[i-1][0]) - h[i-1]*u[i-1];
        u[i] = h[i]/l[i];
        z[i] = ( alpha[i] - h[i-1]*z[i-1] )/l[i];
    }

	fnc::ntuple b(n-1);
	fnc::ntuple c(n);
	fnc::ntuple d(n-1);
	c[0] = 0;
	c[n-1] = 0;

    for (long unsigned i=n-2; i>=1; i--)
    {
        c[i] = z[i] - u[i]*c[i+1];
        b[i] = ( mtx[i+1][1] - mtx[i][1] )/h[i] - h[i]*( c[i+1] + 2*c[i] )/3;
        d[i] = ( c[i+1] - c[i] )/(3*h[i]);
    }


    b[0] = ( mtx[1][1] - mtx[0][1] )/h[0] - h[0]*( c[1] + 2*c[0] )/3;
    d[0] = c[1]/( 3*h[0] );

    for (long unsigned i = 0 ; i<=n-2 ; i++)
    {
        double dom[2] = { mtx[i][0], mtx[i+1][0] };
        double temp[4] = { mtx[i][1], b[i], c[i], d[i] };
        soln[i] = new fnc::polynomial(temp, mtx[i][0], 4, dom);
    }
    std::cout<< mtx[0][0]<<"\n"<< mtx[n-1][0];
    soln.setlobound( mtx[0][0] );
    soln.sethibound( mtx[n-1][0] );

}

/** INTERPOLATE PIECEWISE FUNCTION USING DERIVATIVE VALUES
@brief changes matrix into polynomial
@remark requires value of derivative at first and last nodes. two columns only.
@example see above.
*/

interp_dervpw::interp_dervpw(const fnc::matrix& mtx, double prmval0, double prmvaln)
{
    n = mtx.numrows();
    soln = fnc::pwfunction(n-1);

    h = fnc::ntuple(n-1);
    for (long unsigned i=0; i<=n-2; i++) h[i] = mtx[i+1][0] - mtx[i][0];

    fnc::ntuple alpha(n);

    alpha[0] = 3*( mtx[1][1] - mtx[0][1] )/h[0] - 3*prmval0;
    alpha[n-1] = 3*prmvaln - 3*( mtx[n-1][1] - mtx[n-2][1] )/h[n-2];

    for (long unsigned i=1; i<=n-2; i++) alpha[i] = 3/h[i]*( mtx[i+1][1] - mtx[i][1] ) - 3/h[i-1]*( mtx[i][1] - mtx[i-1][1] );

    fnc::ntuple l(n);
	fnc::ntuple u(n);
	fnc::ntuple z(n);
    l[0] = 2*h[0];
    u[0] = .5;
    z[0] = alpha[0]/l[0];

    for (long unsigned i = 1; i<=n-2; i++)
    {
        l[i] = 2*( mtx[i+1][0] - mtx[i-1][0] ) - h[i-1]*u[i-1];
        u[i] = h[i]/l[i];
        z[i] = ( alpha[i] - h[i-1]*z[i-1] )/l[i];
    }

    l[n-1] = h[n-2]*(2-u[n-2]);
    z[n-1] = ( alpha[n-1] - h[n-2]*z[n-2] )/l[n-1];


    fnc::ntuple b(n-1);
	fnc::ntuple c(n);
	fnc::ntuple d(n-1);

    c[n-1] = z[n-1];

    for (long unsigned i = n-2; i>=1 ; i--)
    {
        c[i] = z[i] - u[i]*c[i+1];
        b[i] = ( mtx[i+1][1] - mtx[i][1] )/h[i] - h[i]*( c[i+1] + 2*c[i] )/3;
        d[i] = ( c[i+1] - c[i] )/(3*h[i]);
    }

    c[0] = z[0] - u[0]*c[1];
    b[0] = ( mtx[1][1] - mtx[0][1] )/h[0] - h[1]*( c[1] + 2*c[0] )/3;
    d[0] = ( c[1] - c[0] )/(3*h[0]);


    for (long unsigned i = 0 ; i<=n-2 ; i++)
    {
        double dom[2] = {mtx[i][0], mtx[i+1][0] };
        double temp[4] = { mtx[i][1], b[i], c[i], d[i] };
        soln[i] = new fnc::polynomial(temp, mtx[i][0], 4, dom);
    }
}


/*****************************************************************************
****************************** LINEAR SYSTEM OF EQs **************************
******************************************************************************/

/** LINEAR SYSTEM MATRIX
@brief produce a matrix from a system of linear equations
@remark will probably crash unless all equations in the system are linear
@example

    fnc::systemfuncs my_system(2);

    my_system.assignfunc(0, flb::linfunc1); //< assuming linfunc1 has been defined in the namespace. see guide in main.
    my_system.assignfunc(1, flb::linfunc2);

    alg::linsys_mtx  lin_trans( funcs );

    fnc::matrix my_mtx = lin_trans.getsoln();

    std::cout<< my_mtx;  //< prints all values of matrix

*/
linsys_mtx::linsys_mtx(fnc::systemfuncs& funcs)
{
    numfuncs = funcs.numfuncs();
    maxdim = funcs.maxdim();
    soln = fnc::matrix( numfuncs , maxdim);

    for ( long unsigned j=0 ; j<maxdim ; j++)
    {
        fnc::elmtuple x(j, maxdim);

        for ( long unsigned i=0 ; i<numfuncs ; i++)
        {
            soln(i,j) = funcs(i, x);
        }
    }
}

/** SOLVE LINEAR SYSTEM
@brief solves a matrix equation Ax = b for x, given a user selected b, using Gaussian Elim with Partial Pivoting
@remark only works on an nXn matrix
@example

    (continuing from the previous example)

    fnc::ntuple b(2);

    b[0] = 1;
    b[1] = 3;

    alg::solve_linsys solve(my_mtx, b);

    fnc::ntuple x = solve.getsoln();

    std::cout<< x;  //< prints the solution


@param inmtx the matrix to solve for
@param column the ntuple b that we're trying to find Ax for.
*/
solve_linsys::solve_linsys(const fnc::matrix& inmtx, const fnc::ntuple& column)
{
    n = inmtx.numrows();

    if (n != inmtx.numcols())
    {
        std::cout<< "ERROR! Cannot use Gauss\' method on "<< n<< "X"<< inmtx.numcols()<<" matrix!"<<std::endl;
        return;
    }

    mtx = inmtx.augment(column);


	long unsigned *nrow = new long unsigned[n];

    for (long unsigned i=0; i<n ; i++)
        nrow[i] = i;

    for (long unsigned i=0; i<n-1; i++)
    {
        fnc::ntuple col = mtx.getcol(i);
        long unsigned p = col.fabsmaxentry(i, n-1);
        if ( col[p] == 0 )
        {
            std::cout<<std::endl<< "Unique solution does not exist."<<std::endl;
            return;
        }
        /** row exchange **/
        if (nrow[i] != nrow[p])
        {
            long unsigned ncopy = nrow[i];
            nrow[i] = nrow[p];
            nrow[p] = ncopy;
        }
        for (long unsigned j=i+1; j<n; j++)
        {
            double m = mtx[nrow[j]][i] / mtx[nrow[i]][i];
            /** to acces vector on other side of equation, use k==n **/
            for (long unsigned k=0; k<=n ; k++)
                mtx[nrow[j]][k] -= m* mtx[nrow[i]][k];
        }
    }
    if ( mtx[nrow[n-1]][n-1] == 0 )
    {
        std::cout<< "Unique solution does not exist."<<std::endl;
		delete[] nrow;
        return;
    }
    /** start backward sub **/
    /** initialize private variable **/
    soln = fnc::ntuple(n);

    /** this whole method assumes that matrix is nXn+1 **/
    soln[n-1] = mtx[nrow[n-1]][n] / mtx[nrow[n-1]][n-1];
    for (long unsigned i=n-2; i>=1; i--)
    {
        double sum = 0;
        for (long unsigned j=i+1; j<n ; j++)
            sum += mtx[nrow[i]][j]*soln[j];

        soln[i] = ( mtx[nrow[i]][n] - sum )/ mtx[nrow[i]][i];
    }

    double sum=0;
    for (long unsigned j=1; j<n; j++)
        sum += mtx[nrow[0]][j]*soln[j];
    soln[0] = ( mtx[nrow[0]][n] - sum )/mtx[nrow[0]][0];

	delete[] nrow;

}

/******************************************************************
****************** NUMERICAL DIFFERENTIATION **********************
*******************************************************************/



double differentiate::threeptbeg(double a, double h)
{
    return ( 1/(2*h)*( -3* (*f)(a)  +  4* (*f)(a + h) - (*f)(a + 2*h) ) );
}

double differentiate::threeptmid(double x, double h)
{
    return ( 1/(2*h) * ( (*f)(x + h) - (*f)(x - h) ) );
}

double differentiate::threeptend(double b, double h)
{
    return ( 1/(2*h) *( (*f)(b-2*h) - 4* (*f)(b - h) + 3* (*f)(b) ) );
}

void differentiate::setvals()
{
    setdervdomain(dom);
    setnodecount(&nodecount);
    setord(&ord);
}

/** NUMERICAL DIFFERENTIATION
@brief approximates the derivative at "innode" points using Richardson Extrapolation and stores the values in a matrix.
@example

    fnc::function my_func = flb::func1;

    alg::differentiate diff(&my_func, 0, 5, 50, 15);

    fnc::polynomial derivative = diff.poly(); //<can return piecewise with diff.piece()

    std::cout<< derivative(1.5);

@param fin pointer to function object to differentiate
@param t0 lower bound of domain to differentiate across
@param tf upper bound
@param innode number of points to differentiate at, typically: (tf-t0)/100 or so
@inord order of approximation you'd like (15 is typical)
*/
differentiate::differentiate(fnc::function* fin, double t0, double tf, long unsigned innode, short inord)
{
    f = fin;
    nodecount = innode;
    ord = inord;

    double a = t0, b = tf, h = ( b - a )/ nodecount;

    /** two column vector to hold indep vals in first and derv vals in sec
        makes array of size nodecount+1 to include right endpoint **/
    soln = fnc::matrix(nodecount+1, 2);

    /** set endpoints **/
    soln[0][0] = a;
    soln[0][1] = threeptbeg(a, h);

    soln[nodecount][0] = b;
    soln[nodecount][1] = threeptend(b, h);

    for (long unsigned i=1 ; i<=nodecount-1 ; i++)
    {
        soln[i][0] = a+i*h;
        soln[i][1] = threeptmid(a + i*h , h);
    }

    for (long unsigned k=0 ; k<nodecount ; k++)
    {

        double** array = new double* [ord];

        array[0] = new double[1];
        array[0][0] = soln[k][1];

        for (short i=1 ; i<ord ; i++)
        {
            array[i] = new double[i+1];

            array[i][0] = threeptmid(a + k*h , h/pow(2.0, i) );

            for (short j=1 ; j<=i ; j++)
                array[i][j] = array[i][j-1] + ( array[i][j-1] - array[i-1][j-1] )/( pow(4.0, j) - 1 );
        }

        /** stores value for derivative into second column **/
        soln[k][1] = array[ord-1][ord-1];

        for (short i=0 ; i<ord ;i++)
            delete[] array[i];
        delete[] array;
    }
}

/** user input **/
differentiate::differentiate(fnc::function* fin)
{
    /*
    f = fin;
    dom[0] = f->lobound();
    dom[1] = f->upbound();
    nodecount = 100;
    ord = 15;
    differentiate derv(f, dom[0], dom[1], nodecount, ord);
    soln = derv.getsoln();
    */

    f = fin;
    setdervdomain(dom);
    setnodecount(&nodecount);
    setord(&ord);
    differentiate derv(f, dom[0], dom[1], nodecount, ord);
    soln = derv.getsoln();

}



/** JACOBIAN MATRIX
@brief approximates the jacobian of a system of functions at a point and stores data as a matrix
@remark process evaluates each function at f(x-h) and f(x+h) so if x is an endpoint of the function's domain,
        then it can cause an error. Also, if encountering lare errors, play with "intol".

@example

        fnc::systemfuncs funcs(2);
        funcs.assign_func(0, flb::func1);
        funcs.assign_func(1, flb::func2);
        fnc::ntuple x(2);
        x[0] = 1;
        x[1] = 3;
        alg::jacobian jac(funcs, x, .0001);
        fnc::matrix my_jac = jac.getsoln();
        std::cout<< my_jac; //< prints the values;


@param funcs the system to evaluate
@param eval the point (ntuple) to evaluate at.
@param inh the tolerance of approximation of the partials. typical: .0001
*/
jacobian::jacobian(fnc::systemfuncs& funcs, const fnc::ntuple& eval, double inh)
{
    /** assumes F: Rn -> Rn **/
    h= inh;
    numfuncs = funcs.numfuncs();
    soln = fnc::matrix(numfuncs, funcs.maxdim() );
    for (long unsigned i=0; i<numfuncs; i++)
        for (long unsigned j=0; j<numfuncs; j++)
            soln[i][j] = partial_fiveptmid (funcs[i] , eval, j);
}

/** approximation to partial **/
double jacobian::partial_fiveptmid(fnc::function* f, fnc::ntuple eval, long unsigned dir)
{
    fnc::elmtuple tup(dir, eval.getsize());
    return ( 1/(12.0*h)* ( (*f)(eval-2*h*tup) - 8* (*f)(eval-h*tup) + 8* (*f)(eval+h*tup) - (*f)(eval+2*h*tup) ));
}

/** user input **/
jacobian::jacobian(fnc::systemfuncs& funcs)
{
    std::cout<< "Enter \"done\" when finished entering coordinates."<<std::endl
        << "Evaluate Jacobian at point: ";
    fnc::ntuple x;
    std::cin>> x;
    setstep(&h);
    numfuncs = funcs.numfuncs();

    jacobian jac(funcs, x, h);
    soln = jac.getsoln();
}



/****************************************************************************
************************* NON LINEAR SYSTEM OF EQs **************************
*****************************************************************************/
/** SOLVE NON-LINEAR SYSTEM
@brief solve a system of non-linear equations, F(x) = 0,  for its root, x, an ntuple
@example

    fnc::systemfuncs funcs(2);
    funcs.assign_funcs(0, flb::func1);
    funcs.assign_funcs(1, flb::func4);
    fnc::ntuple initial_vals(2);
    initial_vals[0] = 1;
    initial_vals[1] = 3;
    alg::solve_nonlinsys  nonlin(funcs, x, .00001);
    std::cout<< nonlin.getsoln();   //< prints the solution

@param funcs the system to solve
@param init initial approximations for the root
@param inmaxit maximum iterations allowed before method quits, typical: 150
@param inh tolerance for jacobian approximation, typical: .0001
*/


solve_nonlinsys::solve_nonlinsys(fnc::systemfuncs& funcs,const fnc::ntuple& init, long unsigned inmaxit, double inh)
{
    numfuncs = funcs.numfuncs();
    maxit = inmaxit;
	jacstep = inh;
    soln = init;

    double h = 1.0/maxit;


    fnc::ntuple b;
    b = -h*funcs(soln);

    for (long unsigned i=0; i<maxit; i++)
    {
		jacobian jac1(funcs, soln , jacstep);
		fnc::matrix A = jac1.getsoln();
		solve_linsys sys1(A, b);
		fnc::ntuple k1 = sys1.getsoln();

        jacobian jac2(funcs, soln + .5*k1, jacstep);
		A = jac2.getsoln();
		solve_linsys sys2(A, b);
		fnc::ntuple k2 = sys2.getsoln();

        jacobian jac3(funcs, soln + .5*k2 , jacstep);
		A = jac3.getsoln();
		solve_linsys sys3(A, b);
		fnc::ntuple k3 = sys3.getsoln();

        jacobian jac4(funcs, soln + k3 , jacstep );
		A = jac4.getsoln();
		solve_linsys sys4(A, b);
		fnc::ntuple k4 = sys4.getsoln();

		soln = soln + (k1 + 2*k2 + 2*k3 + k4)/6;
    }
}

void solve_nonlinsys::getinit()
{
    init = fnc::ntuple(numfuncs);
    for (long unsigned i=0; i<numfuncs ; i++)
    {
        std::cout<< "Initial approximation for x"<<i<<": ";
        std::cin>> init[i];
    }

}

/** user entry of parameters **/
solve_nonlinsys::solve_nonlinsys(fnc::systemfuncs& funcs)
{
    numfuncs = funcs.numfuncs();
    getinit();
	setstep(&jacstep);
	setmaxit(&maxit);

    double h = 1.0/maxit;

    soln = init;

    fnc::ntuple b;
    b = -h*funcs(soln);

    for (long unsigned i=0; i<maxit; i++)
    {
		jacobian jac1(funcs, soln , jacstep);
		fnc::matrix A = jac1.getsoln();
		solve_linsys sys1(A, b);
		fnc::ntuple k1 = sys1.getsoln();

        jacobian jac2(funcs, soln + .5*k1, jacstep);
		A = jac2.getsoln();
		solve_linsys sys2(A, b);
		fnc::ntuple k2 = sys2.getsoln();

        jacobian jac3(funcs, soln + .5*k2 , jacstep);
		A = jac3.getsoln();
		solve_linsys sys3(A, b);
		fnc::ntuple k3 = sys3.getsoln();

        jacobian jac4(funcs, soln + k3 , jacstep );
		A = jac4.getsoln();
		solve_linsys sys4(A, b);
		fnc::ntuple k4 = sys4.getsoln();

		soln = soln + (k1 + 2*k2 + 2*k3 + k4)/6;
    }
}



/****************************** Numerical Integration *****************************/
/** NUMERICAL INTEGRATION
@brief calculate definite integral over small region
@remark after 25 nodes, the speed reduces dramaticaly. 30 nodes will take a while, above 40 is ridiculous
@param fin pointer to function object to take derivative of
@param ina lower bound of interval to integrate on
@param inb upper bound
@param innodes number times to refine approximation
*/

integrate_def::integrate_def(fnc::function* fin, double ina, double inb, double innodes)
{
    f = fin;
    nodecount = innodes;
    double a = ina, b = inb, h = (b - a);

    fnc::ntuple holder(1);
    holder[0] = h/2*( (*f)(a) + (*f)(b) );
    for (long unsigned i=1; i<nodecount; i++ )
    {
        /** fill out the first column with trapezoidal method **/
        fnc::ntuple getter(i+1);
        getter[0] = trapezoidal( holder[0], a, h, i);

        /** use romberg's formula to extrapolate all other entries **/
        for (long unsigned j=1; j<=i; j++)
            getter[j] = getter[j-1] + ( getter[j-1] - holder[j-1] )/( pow(4.0, (double) j) -1);

        holder = getter;
        /** update the step **/
        h /= 2;
    }
    soln = holder[nodecount-1];
}

double integrate_def::trapezoidal (double x, double a, double h, long unsigned i)
{
    long unsigned maxit = pow(2.0, (double) i-1);
    double sum = 0;
    for (long unsigned k=1.0; k<=maxit; k++)
    {
        sum += ( (*f)(a + h*(k - 0.5) ) );
    }
    return (0.5* (x + h*sum));

}
integrate_def::integrate_def(fnc::function* fin)
{
    f = fin;
    setdomain(dom);
    setmeshcount(&nodecount);
    integrate_def integ(fin, dom[0], dom[1], nodecount);
    soln = integ.getsoln();
}


/********************************** Initial Value Problems **************************/

/** SOLVE IVP FIXED STEP / INTEGRATE INDEFINITELY
@brief approximates the solution of an IVP ODE y' = f(x,y) at a bunch of points
@example

    fnc::function myfunc = flb::func1;  //< or any function
    alg::solve_ivp solve( &myfunc, 0, 10, .05, 1000 );
    fnc::polynomial my_poly = solve.integral_poly();
    std::cout<< my_poly(1.5);

@param f pointer to function object to solve
@param ina lower bound on region to solveo ver
@param inb upper bound
@param iniv initial value y(a)
@param innode number of points to approximate at
*/

solve_ivp_fix::solve_ivp_fix(fnc::function* f, double ina, double inb, double iniv, long unsigned innode)
{
    nodecount = innode;
    double t0 = ina, tf = inb, t = t0, w = iniv, h = (tf-t0)/nodecount, k[4];


    /** create nx2 array to hold in 1st col the indep var and 2nd col the apx to dep var **/
    soln = fnc::matrix(nodecount+1, 2);
    soln[0][0] = t0;
    soln[0][1] = w;

    dervval = fnc::ntuple(nodecount+1);
    dervval[0] = (*f)(t,w);

    for (int i = 1; i <= nodecount; i++)
    {
        /** here, since f is multidimensional and thus returns an ntuple,
            we have to derefference to the first element by using the [0]
            at the end **/
        k[0] = h*(*f)(t,w);
        k[1] = h*(*f)(t + h/2 , w + k[0]/2);
        k[2] = h*(*f)(t + h/2 , w + k[1]/2);
        k[3] = h*(*f)(t + h , w + k[2]);

        w = w + (k[0] + 2*k[1] + 2*k[2] + k[3])/6;
        t = t0 + i*h;
        soln[i][0] = t;
        soln[i][1] = w;
        dervval[i] = (*f)(t,w);

    }
}
solve_ivp_fix::solve_ivp_fix(fnc::function* f)
{
    setdomain(domain);
    setivdep(&ivdep);
    setnodecount(&nodecount);
    solve_ivp_fix ivp(f, domain[0], domain[1], ivdep, nodecount);

    soln = ivp.getsoln();
    dervval = ivp.getderv();
}

/** SOLVE IVP ADAPTIVE STEP
@brief uses Runge-Kutta Fehlberg, which adjust the step size depending on the error, to solve y' = f(x,y)
@param ina lower bound for interval to solve over
@param inb upper bound
@param inivdep initial value of dependent variable, y(a)
@param inhmin lower bound the error is allowed to fluctuate to, typical: .00001
@param inhmax upper bound, typical: 1
@param intol minimum difference between iterations to conclude convergence
*/
solve_ivp_adpt::solve_ivp_adpt(fnc::function* f, double ina, double inb, double inivdep, double inhmin, double inhmax, double intol)
{
    ivdep = inivdep;
    hmin = inhmin;
    hmax = inhmax;
    tol = intol;
    double t0= ina, tf= inb, t = t0, w = ivdep, h = hmax, flag = 1, k[6], r, delta;

    long unsigned maxits = (tf - t0)/hmin;

    double c[20] = {3.0/32, 9.0/32,
                    12.0/13, 1932.0/2197, 7200.0/2197, 7296.0/2197,
                    439.0/216, 3680.0/513, 845.0/4104,
                    8.0/27, 3544.0/2565, 1859.0/4104, 11.0/40,
                    1.0/360, 128.0/4275, 2197.0/75240, 2.0/55,
                    25.0/216, 1408.0/2565, 2197.0/4104};


    double** array = new double* [maxits];
    array[0] = new double[2];
    array[0][0] = t;
    array[0][1] = w;


    long unsigned count = 0;

    while (flag == 1)
    {
        //set constants
        /** here, since f is multidimensional and thus returns an ntuple,
            we have to derefference to the first element by using the [0]
            at the end **/
        k[0] = h*(*f)(t,w);
        k[1] = h*(*f)(t + .25*h , w + .25*k[0]);
        k[2] = h*(*f)(t + .375*h , w + c[0]*k[0] + c[1]*k[1]);
        k[3] = h*(*f)(t + c[2]*h , w + c[3]*k[0] - c[4]*k[1] + c[5]*k[2]);
        k[4] = h*(*f)(t + h , w + c[6]*k[0] - 8.0*k[1] + c[7]*k[2] - c[8]*k[3]);
        k[5] = h*(*f)(t + .5*h , w - c[9]*k[0] + 2.0*k[1] - c[10]*k[2] + c[11]*k[3] - c[12]*k[4]);

        //set difference
        r = 1/h * fabs(c[13]*k[0] - c[14]*k[2] - c[15]*k[3] + .02*k[4] + c[16]*k[5]);

        if (r <= tol)   //accept apx
        {
            /** need to make sure this updates immediately **/
            ++count;

            t = t + h;
            w = w + c[17]*k[0] + c[18]*k[2] + c[19]*k[3] - .2*k[4];

            array[count] = new double[2];
            array[count][0] = t;
            array[count][1] = w;
        }

        //verify difference
        delta = .84*pow(tol/r, .25);
        if (delta <= .1) h *= .1;
        else if (delta >= 4) h *= 4;
        else h *= delta;

        //verify h
        if (h > hmax) h = hmax;
        if (t >= tf) flag = 0;
        else if (t + h > tf) h = tf - t;
        else if (h < hmin)
        {
            flag = 0;
            std::cout<< "Minimum h exceeded. Procedure completed unsuccessfully."<<std::endl;
            return;
        }
    }
    //success
    soln = fnc::matrix(array, count+1, 2);
    std::cout<< soln;

    for (long unsigned i=0 ; i<= count ; i++)
        delete[] array[i];
    delete[] array;
}

solve_ivp_adpt::solve_ivp_adpt(fnc::function* f)
{
    setdomain(dom);
    setivdep(&ivdep);
    sethmin(&hmin);
    sethmax(&hmax);
    settol(&tol);
    solve_ivp_adpt adpt(f, dom[0], dom[1], ivdep, hmin, hmax, tol);
    soln = adpt.getsoln();
}

/** SOLVE IVP MULTI STEP
@brief uses a "predictor-corrector" method to calculate the next iteration and then refine it.
@param ina lower bound to calculate solution over
@param inb upper bound
@param inivdep initial value of dependent variable y(a)
@param innode number of points to approximate solution at
*/
solve_ivp_mult::solve_ivp_mult(fnc::function* f, double ina, double inb, double inivdep, double innode )
{
    ivdep = inivdep;
    nodecount = innode;
    double a = ina, b = inb, h = (b - a)/nodecount, k[4], t[4], w[4];
    soln = fnc::matrix(nodecount+1, 2);
    t[0] = a;
    w[0] = ivdep;
    soln(0,0) = a;
    soln(0,1) = ivdep;
    /** uses runge kutta for first 3 iterations **/
    for (short i = 1; i<=3; i++)
    {
        k[0] = h* (*f)( t[i-1] , w[i-1] );
        k[1] = h* (*f)( t[i-1] + .5*h , w[i-1] + .5*k[0] );
        k[2] = h* (*f)( t[i-1] + .5*h , w[i-1] + .5*k[1] );
        k[3] = h* (*f)( t[i-1] + h , w[i-1] + k[2] );

        t[i] = a + i*h;
        soln(i,0) = t[i];

        w[i] = w[i-1] + (k[0] + 2*k[1] + 2*k[2] + k[3])/6.0;
        soln(i,1) = w[i];
    }
    for(long unsigned i = 4; i<=nodecount ; i++)
    {
        soln(i,0) = a + i*h;
        /** predict soln **/
        soln(i,1) = w[3] + h* ( 55* (*f)( t[3] , w[3])
                                    - 59* (*f)( t[2] , w[2])
                                    + 37* (*f)( t[1] , w[1])
                                    - 9* (*f)( t[0] , w[0])  )/24;
        /** correct soln **/
        soln(i,1) = w[3] + h* ( 9* (*f)( soln(i,0) , soln(i,1) )
                                    + 19* (*f)( t[3] , w[3])
                                    - 5* (*f)( t[2] , w[2])
                                    + (*f)( t[1] , w[1])  )/24;
        for (short j=0; j<=2 ; j++)
        {
            t[j] = t[j+1];
            w[j] = w[j+1];
        }
        t[3] = soln(i,0);
        w[3] = soln(i,1);
    }
}

solve_ivp_mult::solve_ivp_mult(fnc::function* f)
{
    setdomain(domain);
    setnodecount(&nodecount);
    setivdep(&ivdep);

    solve_ivp_mult mult(f, domain[0], domain[1], ivdep, nodecount);
    soln = mult.getsoln();
}



/****************************** SYSTEM OF IVP ODE's (Runge Kutta) ********************************/

void solve_ivp_sys::setinitvals()
{
    initvals = fnc::ntuple(numfuncs);
    std::cout<< "Enter the initial value for ";
    for (long unsigned i=0; i<numfuncs; i++)
    {
        std::cout<< "f"<<i<<"(0) = ";
        std::cin>> initvals[i];
    }
}


/** makes use of ntuple arithmetic **/
solve_ivp_sys::solve_ivp_sys(fnc::systemfuncs& funcs, double ina, double inb, double innode, fnc::ntuple init)
{
    numfuncs = funcs.numfuncs();
    initvals = init;
    double a = ina, b = inb, h = (b - a)/nodecount, t = a;

    /** creates matrix to hold a row for every node and a column for the
        approximations of every function but the first column holds the
        values of t **/


    soln = fnc::matrix(nodecount+1, numfuncs+1);
    /** first entry is the tval; **/

    soln[0][0] = t;

    for (long unsigned j = 1; j<=numfuncs ; j++)
        soln[0][j] = initvals[j-1];

    /** creates a matrix to hold the various coefficients
        throughout the calculation **/
    fnc::matrix k(4, numfuncs);


    for (long unsigned i = 1; i<=nodecount; i++)
    {


        /** initialize the ntuple to pass through.
            after each iteration, this gets updated
            to the previous entry of solution **/
        fnc::ntuple base ( soln[i-1], numfuncs+1 );



        for (long unsigned j = 0; j<numfuncs ; j++)
            k[0][j] = h * funcs( j, base );



        /** creates an ntuple for the array k[0], but
            since it only has numfuncs elements and the
            function takes in numfuncs+1 args, have to
            push front**/

        fnc::ntuple tupk0( k[0] , numfuncs);
        tupk0.push_front(h);


        for (long unsigned j = 0; j<numfuncs ; j++)
            k[1][j] = h * funcs( j, base + .5*tupk0 );


        fnc::ntuple tupk1( k[1] , numfuncs );
        tupk1.push_front(h);


        for (long unsigned j = 0; j<numfuncs ; j++)
            k[2][j] = h * funcs( j, base + .5*tupk1);

        fnc::ntuple tupk2( k[2] , numfuncs );
        tupk2.push_front(h);


        for (long unsigned j = 0; j< numfuncs ; j++)
            k[3][j] = h * funcs( j , base + tupk2 );

        t = a+i*h;
        soln(i,0) = t;

        for (long unsigned j = 1 ; j<= numfuncs ; j++)
            soln[i][j] = soln(i-1,j) + ( k(0,j-1) + 2*k(1,j-1) + 2*k(2,j-1) + k(3,j-1) )/6;
    }
}

solve_ivp_sys::solve_ivp_sys(fnc::systemfuncs& funcs)
{
    setdomain(domain);
    setnodecount(&nodecount);
    numfuncs = funcs.numfuncs();
    setinitvals();
    solve_ivp_sys sys(funcs, domain[0], domain[1], nodecount, initvals);
    soln = sys.getsoln();

}

}




