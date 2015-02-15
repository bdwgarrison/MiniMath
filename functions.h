#ifndef FUNCTIONS_H
#define FUNCTIONS_H

#include <iostream>
#include <string>
#include <fstream>


/**************************** NTUPLE DECLARATIONS ********************/
namespace functions_and_mathematical_objects_for_numerical_analysis{

class ntuple
{
public:
    /** ntuple random access iterator **/
    class tupiter;

    ntuple();
    virtual ~ntuple() {delete[] vec;}
    ntuple(const ntuple& other);
    ntuple& operator=(const ntuple& other);

    /** empty constructor **/
    ntuple(long unsigned n);
    /** array constructor **/
    ntuple(double* invec, long unsigned n);

    ntuple& operator=(double x);

    /** need to revisit with allocator **/
    void push_front(double);
    void push_back(double);

    /** accessors **/
    long unsigned getsize() const {return size;}
    double& operator[](long unsigned entry);
    double& get(long unsigned i) const {return vec[i];}


    /** returns vector-magnitude **/
    double mag() const;

    /** returns max from first to last entry **/
    long unsigned fabsmaxentry(long unsigned first, long unsigned last) const;
    double fabsmaxval(long unsigned first, long unsigned last) const;

    /** instream operator **/
    friend std::istream& operator>>(std::istream&, ntuple&);

    /** iterator calls
        have no use for these as of yet. **/
    inline tupiter begin();
    inline tupiter end();

    /** incrementors **/
    ntuple operator-();
    ntuple& operator+=(const ntuple&);
    ntuple& operator-=(const ntuple&);
    ntuple& operator+=(double);
    ntuple& operator-=(double);
    ntuple& operator++();
    ntuple& operator--();
    ntuple operator++(int);
    ntuple operator--(int);

private:
    long unsigned size;
    double* vec;
};


/** non-member functions **/
std::ostream& operator<<(std::ostream&, ntuple);
std::istream& operator>>(std::istream&, ntuple&);

const ntuple& smaller(const ntuple& main, const ntuple& other);
const ntuple& bigger(const ntuple& main, const ntuple& other);


/** arithmetic operators **/
ntuple operator+(const ntuple& main, const ntuple& other);
ntuple operator-(const ntuple& main, const ntuple& other);
double operator*(const ntuple& main, const ntuple& other);

/** boolean operators **/
bool operator==(const ntuple& main, const ntuple& other);
bool operator!=(const ntuple& main, const ntuple& other);
bool operator>=(const ntuple& main, const ntuple& other);
bool operator<=(const ntuple& main, const ntuple& other);
bool operator>(const ntuple& main, const ntuple& other);
bool operator<(const ntuple& main, const ntuple& other);


/** arithmetic operators on doubles **/
ntuple operator+(const ntuple& tup, double x);
ntuple operator+(double x, const ntuple& tup);
ntuple operator-(const ntuple& tup, double x);
ntuple operator-(double x, const ntuple& tup);
ntuple operator*(const ntuple& tup, double x);
ntuple operator*(double x, const ntuple& tup);
ntuple operator/(const ntuple& tup, double x);

/** booleans on doubles **/
bool operator==(ntuple tup, double x);
bool operator==(double x, ntuple tup);

bool operator!=(ntuple tup, double x);
bool operator!=(double x, ntuple tup);

bool operator<(ntuple tup, double x);
bool operator<(double x, ntuple tup);

bool operator<=(ntuple tup, double x);
bool operator<=(double x, ntuple tup);

bool operator>(ntuple tup, double x);
bool operator>(double x, ntuple tup);

bool operator>=(ntuple tup, double x);
bool operator>=(double x, ntuple tup);

/************* NTUPLE ITERATOR DECLARATIONS ***************/

class ntuple::tupiter
{
public:
    /** constructors **/
    tupiter();
    ~tupiter();
    tupiter(const tupiter& other) {pos = other.pos;}
    tupiter& operator=(const tupiter& other) {pos = other.pos; return *this;}
    tupiter(double* ptr) {pos = ptr;}

    /** access **/
    inline double& operator*() {return *pos;}
    inline double* operator->() {return pos;}
    inline double& operator[](const int& n) {return pos[n];}

    /** arithmetic **/
    inline tupiter& operator++() {++pos; return *this;}
    inline tupiter& operator--() {--pos; return *this;}
    inline tupiter operator++(int) {tupiter tmp(*this); ++tmp; return tmp;}
    inline tupiter operator--(int) {tupiter tmp(*this); --tmp; return tmp;}
    inline tupiter& operator+=(const int& n) {pos+=n; return *this;}
    inline tupiter& operator-=(const int& n) {pos-=n; return *this;}
    inline tupiter operator+(const int& n) {return tupiter(pos+n);}
    inline tupiter operator-(const int& n) {return tupiter(pos-n);}
    //friend inline tupiter operator+(const int& n, tupiter& it) {return tupiter( n+it[0] );}
    //friend inline tupiter operator-(const int& n, tupiter& it) {return tupiter( n-it[0] );}
    /** bools **/
    inline bool operator==(const tupiter& other) {return pos == other.pos;}
    inline bool operator!=(const tupiter& other) {return pos != other.pos;}
    inline bool operator<=(const tupiter& other) {return pos <= other.pos;}
    inline bool operator>=(const tupiter& other) {return pos >= other.pos;}
    inline bool operator<(const tupiter& other) {return pos < other.pos;}
    inline bool operator>(const tupiter& other) {return pos > other.pos;}

private:
    double* pos;
};



/************ ELEMENTARY NTUPLE DECLARATIONS ****************/

class elmtuple : public ntuple
{
public:
    elmtuple(long unsigned pos, long unsigned insize) : ntuple(insize) { (*this)[pos] = 1; }
};

/*********** QUICK ARGUMENT DECLARATIONS ******************/
/** the purpose of these is to act as a quick constructor for ntuples.
    these will inevitably be sliced, as function arguments expect
    an ntuple and not a qarg, but this isn't a problem because
    qargs hold no specific data.**/
class qarg : public ntuple
{
public:
    qarg(double x) : ntuple(1)
    { (*this)[0] = x;}

    qarg(double x0, double x1) : ntuple(2)
    { (*this)[0] = x0; (*this)[1] = x1;}

    qarg(double x0, double x1, double x2) : ntuple(3)
    { (*this)[0] = x0; (*this)[1] = x1; (*this)[2] = x2;}

    qarg(double x0, double x1, double x2, double x3) : ntuple(4)
    { (*this)[0] = x0; (*this)[1] = x1; (*this)[2] = x2; (*this)[3] = x3;}

    qarg(double x0, double x1, double x2, double x3, double x4) : ntuple(5)
    { (*this)[0] = x0; (*this)[1] = x1; (*this)[2] = x2; (*this)[3] = x3; (*this)[4] = x4;}

    qarg(double x0, double x1, double x2, double x3, double x4, double x5) : ntuple(6)
    { (*this)[0] = x0; (*this)[1] = x1; (*this)[2] = x2; (*this)[3] = x3; (*this)[4] = x4; (*this)[5] = x5;}


};

/****************** FUNCTION DECLARATIONS ******************/
/** FUNCTION OBJECTS
@remark Functions have arguments of up to 6 dimensions built in
@example F(x1,x2,x3,x4,x5,x6) will work
@remark but any higher need to be added in manually
*/
class function
{
public:
    /** BIG 4 **/
    function() {fnd = NULL; domain = new double[2]; domdim=1;}
    virtual ~function() {delete[] domain;}
    function(const function& other);
    function& operator=(const function& other);


    /** f: Rn -> R1 **/
    function(double (*f)(ntuple), long unsigned indim)
    {fnd = f; domdim = indim;}


    /** nd argument **/

    virtual double operator()(ntuple intuple);

    virtual double operator()(double arg)
    {return fnd( qarg(arg));}

    virtual double operator()(double x0, double x1)
    {return fnd( qarg(x0,x1));}

    virtual double operator()(double x0, double x1, double x2)
    {return fnd( qarg(x0,x1,x2));}

    virtual double operator()(double x0, double x1, double x2, double x3)
    {return fnd( qarg(x0,x1,x2,x3));}

    virtual double operator()(double x0, double x1, double x2, double x3, double x4)
    {return fnd( qarg(x0,x1,x2,x3,x4));}

    virtual double operator()(double x0, double x1, double x2, double x3, double x4, double x5)
    {return fnd( qarg(x0,x1,x2,x3,x4,x5));}


    virtual double operator()(double* arg, long unsigned dim) {return fnd( ntuple(arg,dim));}

    virtual void print() {}


    /** domain functions **/
    double* getdom() const {return domain;}
    double lobound() const {return domain[0];}
    double upbound() const {return domain[1];}
    long unsigned getdomdim() const {return domdim;}
    void setlobound(double in) {domain[0]=in;}
    void sethibound(double in) {domain[1]=in;}

    friend std::ostream& operator<<(std::ostream&, function&);

private:

    double (*fnd)(ntuple);


protected:
    double* domain;
    long unsigned domdim;
};


/****************** Piecewise Function Declarations ******************/
/**PIECEWISE FUNCTIONS
@remark these have arguments of up to 6 dimensions coded in, but these are only designed for one varibale
        so they only evaluate the function at the first coordinate
*/
class pwfunction : public function
{
public:
    /** BIG 4 **/
    pwfunction() {pieces = new function*[1]; numfuncs = 0; domdim = 1;}
    ~pwfunction() {delete[] pieces; delete[] domain;}
    pwfunction(const pwfunction&);
    pwfunction& operator=(const pwfunction&);

    /** blank constructor **/
    pwfunction(long unsigned);

    /** used constructor **/
    pwfunction(function** inpieces, long unsigned innum);


    /** index operator **/
    function*& operator[](long unsigned i);
    function* getfunc(long unsigned i) const;

    /** size function **/
    long unsigned size() const {return numfuncs;}

    /** argument operator **/
    virtual double operator()(double arg);
    virtual double operator()(double x0, double x1) {return (*this)(x0);}
    virtual double operator()(double x0, double x1, double x2) {return (*this)(x0);}
    virtual double operator()(double x0, double x1, double x2, double x3) {return (*this)(x0);}
    virtual double operator()(double x0, double x1, double x2, double x3, double x4) {return (*this)(x0);}
    virtual double operator()(double x0, double x1, double x2, double x3, double x4, double x5)
    {return (*this)(x0);}

    /** since multi-dim variable piecewise functions aren't used, passes to single variable **/
    virtual double operator()(ntuple arg) {return (*this)(arg[0]);}



private:
    /** array of function pointers to utilize polymorphism of functions **/
    function** pieces;
    long unsigned numfuncs;
};



/****************** Polynomial Function Declarations ******************/
class polynomial : public function
{
public:
    /** Big 4 **/
    polynomial();
    ~polynomial() {delete[] polycoeffs; delete[] center;}
    polynomial(const polynomial& other);
    polynomial& operator=(const polynomial& other);

    /** poly of insize length with zero for all coeffs and center **/
    polynomial(long unsigned insize);



    /** construct poly from array. incoeffs ptr to first elmt, n is size **/
    polynomial(const double* incoeffs, long unsigned n);

    polynomial(const double* incoeffs, long unsigned n, double* indomain);

    /** construct poly with a center for each term
        this is used for hermite polynomials and newton/lagrange polynomials
        that result from interpolation algorithms.
        Hermite Polynomial:
        H(x) = a0 + a1(x−x0) + a2(x−x0)^2 + a3(x−x0)^2 (x−x1) + a4(x−x0)^2 (x−x1)^2 + · · · +a2n+1 (x − x0 )^2 (x − x1 )^2 · · · (x − xn−1 )^2 (x − xn ).
        where the xi are the entries of the center.

        if one wanted something of a Taylor polynomial, where each term has the same center,
        T(x) = 1 + (x-a1) + (x-a1)^2 + ... + (x-a1)^n
        then center = {a1, a1, a1, ..., a1 }
        **/
    polynomial(const double* incoeffs, const double* incenter, long unsigned n);

    /** constructor poly with coeffs, center, and domain **/
    polynomial(const double* incoeffs, const double* incenter, long unsigned n, const double* indomain);

    /** constructor poly with coeffs, a single center for all terms, and domain **/
    polynomial(const double* incoeffs, double incenter, long unsigned n, const double* indomain);


    /** length & degree operators **/
    long unsigned size() const {return coefflength;}
    long unsigned degree() const {return size()-1;}

    /** argument operator uses coefficients to evaluate**/
    virtual double operator()(double arg);
    virtual double operator()(double x0, double x1) {return (*this)(x0);}
    virtual double operator()(double x0, double x1, double x2) {return (*this)(x0);}
    virtual double operator()(double x0, double x1, double x2, double x3) {return (*this)(x0);}
    virtual double operator()(double x0, double x1, double x2, double x3, double x4) {return (*this)(x0);}
    virtual double operator()(double x0, double x1, double x2, double x3, double x4, double x5)
    {return (*this)(x0);}


    virtual double operator()(ntuple arg) {return (*this)(arg[0]);}

    /** index operators **/
    double& operator[](long unsigned i) const {return polycoeffs[i];}

    double& cof(long unsigned i) {return polycoeffs[i];}
    double& ctr(long unsigned i) {return center[i];}

private:
    double* polycoeffs;
    double* center;
    long unsigned coefflength;
};


/**************************** MATRIX DECLARATIONS **********************/
/** make matrices based on ntuples? **/
class matrix
{
public:
    /** BIG 4 **/
    matrix(){mtx = new double*[1]; mtx[0] = new double[1]; mtx[0][0]=0; rowcount = 1; colcount = 1;}
    ~matrix(){ for (long unsigned i = 0; i<rowcount; i++) delete[] mtx[i];  delete[] mtx; }
    matrix(const matrix& orig);
    matrix& operator=(const matrix& orig);

    /** manual constructor **/
    matrix(double** inmtx, long unsigned inrow, long unsigned incol);
    /** blank constructor **/
    matrix(long unsigned inrow, long unsigned incol);
    /** construct from file **/
    matrix(std::ifstream&);

    /** single row matrix, ie vector **/
    matrix(double* inmtx, long unsigned inrow);

    /** grabs data from file and puts into matrix **/
    void get_nodes();

    matrix augment(const ntuple&) const;

    /** returns row **/
    /** doing [i][j] returns value **/
    double* operator[](long unsigned i) const;    //<accesses rows

    ntuple getrow(long unsigned i) const;
    ntuple getcol(long unsigned i) const;

    /** doing (i,j) returns reference **/
    double& operator()(long unsigned i, long unsigned j) const {return mtx[i][j];}



    void pastecol(long unsigned i, double* incol) const;   //pastes column onto array parameter

    long unsigned numrows() const {return rowcount;}
    long unsigned numcols() const {return colcount;}

private:
    double** mtx;   //matrix of double pointers
    long unsigned rowcount;
    long unsigned colcount;
};

/** non-member functions **/
std::ostream& operator<<(std::ostream&, matrix&);


/************************** System of Functions Declarations ************************/
/** can also be used as multi-dimensional function, F: Rn -> Rm **/
class systemfuncs
{
public:
    /** big 4 **/
    systemfuncs() {system = new function*[1] ; length = 1;}
    ~systemfuncs() {/*for (long unsigned i=0; i<length; i++) delete system[i];*/ delete[] system;}
    systemfuncs(const systemfuncs&);
    systemfuncs& operator=(const systemfuncs&);

    /** construct from double ptr **/
    systemfuncs(function**,long unsigned);

    /** construct blank system of unsigned long length **/
    systemfuncs(long unsigned);


    /** deletes system[entry] and assigns newfunc to system[entry] **/
    void assign_func(long unsigned entry, function newfunc);
    void assign_func(long unsigned entry, polynomial newfunc);
    void assign_func(long unsigned entry, pwfunction newfunc);

    /** length **/
    long unsigned numfuncs() const {return length;}

    /** find maximum dimension **/
    long unsigned maxdim() const;

    /** returns pointer by reference **/
    function*& operator[](long unsigned i) const;
    //function*& getfunc(long unsigned i) const;


    /** evaluate system[entry] at point arg. **/
    double operator()(long unsigned entry, double arg) const { return (( *( (*this)[entry] ) )(arg) ) ; }
    /** nd version **/
    double operator()(long unsigned entry, ntuple arg) const { return (( *( (*this)[entry] ))(arg) ) ; }


    /** treats system of funcs as a function F: Rn -> Rm **/
    ntuple operator()(ntuple arg);


private:
    function** system;
    long unsigned length;
};



namespace library_of_predefined_simple_functions
{
    /// 1 dimensional functions

    double func1_f(ntuple x);

    double func2_f(ntuple x);

    double func3_f(ntuple x);

    double func4_f(ntuple x);

    double func5_f(ntuple x);

    double func6_f(ntuple x);

    double func7_f(ntuple x);

    double func8_f(ntuple x);


    /// n dimensional functions
    double func9_f(ntuple x);

    /// functions being input into system of eqns must have numfuncs+1 arguments?
    double func10_f(ntuple x);

    double func11_f(ntuple x);

    double func12_f(ntuple x);

    double func13_f(ntuple x);

    double func14_f(ntuple x);

    double func15_f(ntuple x);

    double func16_f(ntuple x);

    double func17_f(ntuple x);

    double func18_f(ntuple x);

    double func19_f(ntuple x);

    /// convert all bare functions to class functions

    extern function func1;
    extern function func2;
    extern function func3;
    extern function func4;
    extern function func5;
    extern function func6;
    extern function func7;
    extern function func8;

    extern function func9;
    extern function func10;
    extern function func11;

    extern function func12;
    extern function func13;

    extern function func14;
    extern function func15;
    extern function func16;

    extern function func17;

    extern function func18;
    extern function func19;

}

}

#endif
