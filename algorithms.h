
#ifndef ALGORITHMS_H
#define ALGORITHMS_H




namespace fnc = functions_and_mathematical_objects_for_numerical_analysis;

namespace algorithms_for_numerical_analysis_implementation_on_function_objects{

class algorithm
{
public:

protected:
    void setdomain(double* indom);
    void setdervdomain(double* indom);
    void setintegdomain(double* indom);
    void setnodecount(long unsigned* innodecount);
    void setmeshcount(long unsigned* inmeshcount);
    void setord(short* inord);
    void settol(double* tol);
    void setmaxit(long unsigned* maxit);
    void setinit1(double* init1);
    void setinit2(double* init2);
    void setivdep(double* inivdep);
    void setivderv(double* inivderv);
    void sethmin(double*);
    void sethmax(double*);
	void setstep(double* step);
    void setboundary(double*);
};



/************************************* SOLVERS *****************************************/

class bisection : public algorithm
{
public:
    bisection() {}
    void operator()( fnc::function* f);
	double operator() (fnc::function* f, double ina, double inb, long unsigned inmaxit, double intol);

    void getinfo();
private:
    double domain[2];
    long unsigned maxit;
    double tol;
    //ntuple root;
    double root;
};


class newton : public algorithm
{
public:
    newton() {}
    void operator()(fnc::function* f, fnc::function* df);
    double operator() (fnc::function*f, fnc::function* df, long unsigned inmaxit, double intol, double init);
    void getinfo();
private:
    long unsigned  maxit;
    double tol;
    double root;
};

class secant : public algorithm
{
public:
    secant() {}
    void operator()(fnc::function*);
    double operator()(fnc::function* f, long unsigned, double, double, double);
    void getinfo();
private:
    long unsigned maxit;
    double tol;
    double root;
    double init;
};

class regulafalsi : public algorithm
{
public:
    regulafalsi() {}
    void operator()(fnc::function* f);
    double operator()(fnc::function*, double, double, long unsigned, double);
    void getinfo();
private:
    double domain[2];
    long unsigned maxit;
    double tol;
    double root;
};


/****************************** INTERPOLATION METHODS ***********************************/


class interp_basicpoly
{
public:
    interp_basicpoly(const fnc::matrix&);
    fnc::polynomial& getsoln() {return soln;}
private:
    long unsigned n;
    fnc::matrix F;
    fnc::polynomial soln;
};

class interp_dervpoly
{
public:
    interp_dervpoly(const fnc::matrix&);
    fnc::polynomial& getsoln() {return soln;}
private:
    long unsigned n;
    long unsigned n2;
    fnc::matrix Q;
    fnc::polynomial soln;
};

//pwfunction natcubicspline(const fnc::matrix&);

class interp_basicpw
{
public:
    interp_basicpw(const fnc::matrix&);
    fnc::pwfunction& getsoln() {return soln;}
private:
    long unsigned n;
    fnc::pwfunction soln;

    fnc::ntuple h;
};

class interp_dervpw
{
public:
    interp_dervpw(const fnc::matrix&, double, double);
    fnc::pwfunction& getsoln() {return soln;}
private:
    long unsigned n;
    fnc::pwfunction soln;
    fnc::ntuple h;
};

/********************************* LINEAR METHODS ****************************************/

/** transforms linear function into fnc::matrix **/
class linsys_mtx
{
public:
    linsys_mtx(fnc::systemfuncs&);
    fnc::matrix& getsoln() {return soln;}
private:
    long unsigned numfuncs;
    long unsigned maxdim;
    fnc::matrix soln;
};


class solve_linsys
{
public:
    solve_linsys(const fnc::matrix&, const fnc::ntuple&);
    fnc::ntuple& getsoln() {return soln;}

private:
    long unsigned n;
    fnc::matrix mtx;
    fnc::ntuple soln;
};


/****************************** NUMERICAL DIFFERENTIATION **************************************/

class differentiate : public algorithm
{
public:

    differentiate(fnc::function* fin);
    differentiate(fnc::function* fin, double t0, double tf, long unsigned innode, short inord);

    fnc::matrix& getsoln() {return soln;}

    fnc::polynomial poly() { interp_basicpoly poly(soln); return poly.getsoln();}
    fnc::pwfunction piece() { interp_basicpw piece(soln); return piece.getsoln();}

private:
    void setvals();
    /** function being evaluated **/
    fnc::function* f;

    fnc::matrix soln;

    double dom[2];
    long unsigned nodecount;
    short ord;

    double threeptbeg(double a, double h);
    double threeptmid(double x, double h);
    double threeptend(double b, double h);
};

class jacobian : public algorithm
{
public:
    double partial_fiveptmid(fnc::function*, fnc::ntuple , long unsigned);

    jacobian(fnc::systemfuncs& funcs);
	jacobian(fnc::systemfuncs& funcs, const fnc::ntuple& eval, double inh);
	fnc::matrix& getsoln() {return soln;}

private:

	fnc::matrix soln;
    double h;
    long unsigned numfuncs;
};

/****************************** NON-LINEAR METHODS ************************************/

class solve_nonlinsys : public algorithm
{
public:
    solve_nonlinsys(fnc::systemfuncs&);
    solve_nonlinsys( fnc::systemfuncs& , const fnc::ntuple&, long unsigned, double);

	fnc::ntuple& getsoln() {return soln;}
private:
    void getinit();

	double jacstep;
    fnc::ntuple init;
    long unsigned numfuncs;
    long unsigned maxit;
	fnc::ntuple soln;
};


/****************************** NUMERICAL INTEGRATION (definite) ***************************/

class integrate_def : public algorithm
{
public:
    integrate_def(fnc::function* fin);
    integrate_def(fnc::function* fin, double ina, double inb, double innodes);
    double getsoln() {return soln;}
private:
    fnc::function* f;
    double soln;

    double dom[2];
    long unsigned nodecount;
    double trapezoidal (double x, double a, double h, long unsigned n);
};


/********************************** IVP ODEs *****************************/
/** fixed step **/
class solve_ivp_fix : public algorithm
{
public:
    solve_ivp_fix(fnc::function*);
    solve_ivp_fix(fnc::function* f, double ina, double inb, double iniv, long unsigned innode);

    fnc::matrix& getsoln() {return soln;}

    fnc::ntuple& getderv() {return dervval;}

    fnc::polynomial integral_poly()
    { fnc::matrix mtx = soln.augment(dervval);  interp_dervpoly interp(mtx); return interp.getsoln(); }



    fnc::pwfunction integral_piece()
    { interp_dervpw interp(soln, dervval[0], dervval[nodecount]); return interp.getsoln(); }

private:
    //function* func;

    fnc::matrix soln;
    fnc::ntuple dervval;

    double domain[2];
    double ivdep;
    long unsigned nodecount;
};

/** adaptive step **/
class solve_ivp_adpt : public algorithm
{
public:
    solve_ivp_adpt(fnc::function*);
    solve_ivp_adpt(fnc::function* f, double ina, double inb, double inivdep, double inhmin, double inhmax, double intol);
    fnc::matrix& getsoln() {return soln;}

private:
    fnc::function* f;
    fnc::matrix soln;

    double dom[2];
    double ivdep;
    double hmin;
    double hmax;
    double tol;
};

/** multi-step **/
class solve_ivp_mult : public algorithm
{
public:
    solve_ivp_mult(fnc::function*);
    solve_ivp_mult(fnc::function* f, double ina, double inb, double inivdep, double innode );
    fnc::matrix& getsoln() {return soln;}
private:
    fnc::matrix soln;
    double domain[2];
    long unsigned nodecount;
    double ivdep;

};

/** solve system **/
class solve_ivp_sys : public algorithm
{
public:
    solve_ivp_sys(fnc::systemfuncs&);
    solve_ivp_sys(fnc::systemfuncs& funcs, double ina, double inb, double innode, fnc::ntuple init);
    fnc::matrix& getsoln() {return soln;}
private:
    fnc::matrix soln;
    double domain[2];
    long unsigned nodecount;
    long unsigned numfuncs;
    fnc::ntuple initvals;
    void setinitvals();
};

}

#endif
