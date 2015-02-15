


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


/********************************** 1D Solvers **************************************/
void bisection::getinfo()
{
    setdomain(domain);
    setmaxit(&maxit);
    settol(&tol);
}

void bisection::operator()( fnc::function* f)
{
	getinfo();
	(*this)(f, domain[0], domain[1], maxit, tol);
	return;
}


/**
@brief bisection method finds the root of an equation on an interval such that f has opposite signs at the endpoints.
	ie: interval [2,5] such that f(2) = -1 and f(5) = 3

@param ina the lower bound of the interval to perform bisection algorithm on
@param inb the upper bound
@param inmaxit the maximum number of iterations to perform before returning an unsuccessful attempt.
@param intol the tolerance for termination of the method. if the difference between the iterations is less than the tolerance,
	then it's assumed to have converge.

@return the root of the function.
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




void newton::getinfo()
{
    setmaxit(&maxit);
    settol(&tol);
    setinit1(&root);
}
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

/******************************* INTERPOLATING METHODS **************************/

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


/** returns hermite polynomial.
    requires x, y, and y' as indval(independent),
    depval(dependent), and prmval(prime) **/
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
}

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


/************************** LINEAR SYSTEM OF EQs **************************/
/** if evaluating a linear function then this produces its matrix **/

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

/** can only take in nXn **/
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

    /** this whole method assumes that matrix is nXn+1 which could
        be a problem with the implementation. will check later **/
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

/************************ NUMERICAL DIFFERENTIATION **********************/

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

/** user input **/
differentiate::differentiate(fnc::function* fin)
{
    f = fin;
    setdervdomain(dom);
    setnodecount(&nodecount);
    setord(&ord);
    differentiate derv(f, dom[0], dom[1], nodecount, ord);
    soln = derv.getsoln();
}

/** manual input **/
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


/******************************** JACOBIAN *********************************/
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

/** manual input **/
jacobian::jacobian(fnc::systemfuncs& funcs, const fnc::ntuple& eval, double inh)
{
    /** assumes F: Rn -> Rn **/
    h= inh;
    numfuncs = funcs.numfuncs();
    soln = fnc::matrix(numfuncs, numfuncs);
    for (long unsigned i=0; i<numfuncs; i++)
        for (long unsigned j=0; j<numfuncs; j++)
            soln[i][j] = partial_fiveptmid (funcs[i] , eval, j);

}



/****************************** NON LINEAR SYSTEM OF EQs **************************/
void solve_nonlinsys::getinit()
{

    //double* vals = new double[numfuncs];
    init = fnc::ntuple(numfuncs);
    for (long unsigned i=0; i<numfuncs ; i++)
    {
        std::cout<< "Initial approximation for x"<<i<<": ";
        //std::cin>> vals[i];
        std::cin>> init[i];
    }
    //init = ntuple(vals, numfuncs);
    //delete[] vals;
}

/** solves F(x) = 0, where x and 0 are vectors and F: Rn-> Rn **/
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

    double a = dom[0], b = dom[1], h = (b - a);

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


/********************************** Initial Value Problems **************************/

/** designed for fnc::functions R2 -> R **/
solve_ivp_fix::solve_ivp_fix(fnc::function* f)
{
    setdomain(domain);
    setivdep(&ivdep);
    setnodecount(&nodecount);

    double t0 = domain[0], tf = domain[1], t = t0, w = ivdep, h = (tf-t0)/nodecount, k[4];


    /** create nx2 array to hold in 1st col the indep var and 2nd col the apx to dep var **/
    /*
    double** array = new double* [nodecount+1];
    for (long unsigned i=0; i<=nodecount ; i++)
        array[i] = new double[2];

    array[0][0] = t0;
    array[0][1] = w;
    */

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
        /*
        array[i][0] = t;
        array[i][1] = w;
        */

    }
/*
    soln = matrix(array, nodecount+1, 2);

    for (long unsigned i=0; i<=nodecount; i++)
        delete[] array[i];
    delete[] array;

    std::cout<< soln;
*/
}

solve_ivp_adpt::solve_ivp_adpt(fnc::function* f)
{
    setdomain(dom);
    setivdep(&ivdep);
    sethmin(&hmin);
    sethmax(&hmax);
    settol(&tol);
    double t0= dom[0], tf= dom[1], t = t0, w = ivdep, h = hmax, flag = 1, k[6], r, delta;

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

/** multi-step method for single IVP **/
solve_ivp_mult::solve_ivp_mult(fnc::function* f)
{
    setdomain(domain);
    setnodecount(&nodecount);
    setivdep(&ivdep);

    double a = domain[0], b = domain[1], h = (b - a)/nodecount, k[4], t[4], w[4];
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
solve_ivp_sys::solve_ivp_sys(fnc::systemfuncs& funcs)
{
    setdomain(domain);
    setnodecount(&nodecount);

    numfuncs = funcs.numfuncs();

    setinitvals();

    double a = domain[0], b = domain[1], h = (b - a)/nodecount, t = a;

    /** creates matrix to hold a row for every node and a column for the
        approximations of every function but the first column holds the
        values of t **/


    //soln = matrix(nodecount, numfuncs+1);
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

}




