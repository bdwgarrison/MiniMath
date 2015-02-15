
#include "functions.h"
#include "algorithms.h"
#include <cmath>
#include <stdexcept>
#include <fstream>
#include <sstream>


/********************************** ntuple DEFINITIONS **************************/

namespace functions_and_mathematical_objects_for_numerical_analysis{

ntuple::tupiter ntuple::begin()
{ return tupiter(vec); }
ntuple::tupiter ntuple::end()
{ return tupiter(vec+size); }

/** constructs 1-tuple of 0**/
ntuple::ntuple()
{
    size = 1;
    vec = new double[1];
}

ntuple::ntuple(const ntuple& other)
{
    size = other.getsize();
    vec = new double[size];
    for(long unsigned i=0; i<size; i++)
        vec[i] = other.get(i);
}

ntuple& ntuple::operator=(const ntuple& other)
{
    size = other.getsize();
    delete[] vec;
    vec = new double[size];
    for (long unsigned i=0; i<size; i++)
        vec[i]=other.get(i);
    return *this;
}

ntuple::ntuple(long unsigned n)
{
    size = n;
    vec = new double[size];
    for (long unsigned i=0; i<size; i++)
        vec[i]=0;
}

/** array constructor **/
ntuple::ntuple(double* invec, long unsigned n)
{
    vec = new double[n];
    size = n;
    for (long unsigned i=0; i<n; i++)
        vec[i] = invec[i];
}

/** deletes all entries and turns into 1-tuple **/
ntuple& ntuple::operator=(double x)
{
    size = 1;
    delete[] vec;
    vec = new double[size];
    vec[0] = x;
    return *this;
}

/** need to revisit with allocators **/
void ntuple::push_front(double x)
{
    double* temp = new double[size];
    for (long unsigned i=0 ; i<size ; i++)
        temp[i] = (*this)[i];
    delete[] vec;
    vec = new double[size+1];
    vec[0] = x;
    for (long unsigned i=0; i<size; i++)
        vec[i+1] = temp[i];
    delete[] temp;
    ++size;
}
void ntuple::push_back(double x)
{
    double* temp = new double[size];
    for (long unsigned i=0 ; i<size ; i++)
        temp[i] = (*this)[i];
    delete[] vec;
    vec = new double[size+1];
    for (long unsigned i=0; i<size; i++)
        vec[i] = temp[i];
    delete[] temp;
    vec[size] = x;
    ++size;
}


double& ntuple::operator[](long unsigned entry)
{
    if (entry >= 0 && entry <getsize())
        return vec[entry];
    else
    {
        std::cout<< std::endl<< entry<< " is not in this ntuple's valid position range: 0 to "<<getsize()-1<<std::endl;
    }
}

double ntuple::mag() const
{
    long unsigned n = getsize();
    double sum=0;
    for (long unsigned i=0; i<n ; i++)
    {
        sum+= ( get(i) * get(i) );
    }
    return pow(sum, .5);
}


/** can use iterators here **/
long unsigned ntuple::fabsmaxentry(long unsigned first, long unsigned last) const
{
    long unsigned max=first;
    for (long unsigned i=first; i<last; i++)
        if ( fabs(get(i)) < fabs(get(i+1)) )
            max = i+1;
    return max;
}

double ntuple::fabsmaxval(long unsigned first, long unsigned last) const
{
    double max=0;
    for (long unsigned i=first; i<last; i++)
    {
        if ( fabs(get(i)) < fabs(get(i+1)) ) max = get(i+1);
        else continue;
    }
    return max;
}

/********* incrementors ***********/
ntuple ntuple::operator-()
{ ntuple tmp(*this); for (long unsigned i=0; i<size; ++i) tmp[i] = -tmp[i]; return tmp;}

ntuple& ntuple::operator+=(const ntuple& tup)
{ for (long unsigned i=0; i<size; ++i) vec[i] += tup.get(i); return *this;}

ntuple& ntuple::operator-=(const ntuple& tup)
{ for (long unsigned i=0; i<size; ++i) vec[i] -= tup.get(i); return *this;}

ntuple& ntuple::operator+=( double x)
{ for (long unsigned i=0; i<size; ++i) vec[i] += x; return *this;}

ntuple& ntuple::operator-=( double x)
{ for (long unsigned i=0; i<size; ++i) vec[i] -= x; return *this;}

ntuple& ntuple::operator++()
{ for (long unsigned i=0; i<size; ++i) ++vec[i]; return *this;}

ntuple& ntuple::operator--()
{ for (long unsigned i=0; i<size; ++i) --vec[i]; return *this;}

ntuple ntuple::operator++(int)
{ ntuple temp(*this); for (long unsigned i=0; i<size; ++i) ++(temp[i]); return temp;}

ntuple ntuple::operator--(int)
{ ntuple temp(*this); for (long unsigned i=0; i<size; ++i) --(temp[i]); return temp;}



std::ostream& operator<<(std::ostream& out, ntuple tup)
{
    if (tup.getsize() > 0)
    {
        long unsigned n = tup.getsize() - 1;
        out<< "(";
        for (long unsigned i = 0; i<n; i++)
        {
            out<< tup[i]<<", ";
        }
        out<< tup[n]<<")";
        return out;
    }
    else
    {
        out<<std::endl<< "ntuple has invalid size."<<std::endl;
        return out;
    }
}

std::istream& operator>>(std::istream& in, ntuple& tup)
{
    long unsigned size = tup.getsize();
    for (long unsigned i=0; i< size; i++)
        std::cin>> tup[i];
	return in;
}

const ntuple& smaller(const ntuple& main, const ntuple& other)
{
    if (main.getsize() < other.getsize())
        return main;
    else return other;
}
const ntuple& bigger(const ntuple& main, const ntuple& other)
{
    if (main.getsize() >= other.getsize())
        return main;
    else return other;
}


/** arithmetic operators **/

ntuple operator+(const ntuple& main, const ntuple& other)
{
    /** only goes up to minimum if trying to add ntuples of differing size **/
    ntuple big = bigger(main, other);
    ntuple small = smaller(main, other);
    long unsigned minsize = small.getsize();
    long unsigned maxsize = big.getsize();

    for (long unsigned i=0; i< minsize; i++)
        big[i] += small[i];
    return big;
}

ntuple operator-(const ntuple& main, const ntuple& other)
{
    ntuple neg = other;
    neg = -neg;
    ntuple big = bigger(main, neg);
    ntuple small = smaller(main, neg);
    long unsigned minsize = small.getsize();
    long unsigned maxsize = big.getsize();

    for (long unsigned i=0; i< minsize; i++)
        big[i] += small[i];
    return big;
}

/** dot product **/
double operator*(const ntuple& main, const ntuple& other)
{
    ntuple big = bigger(main, other);
    ntuple small = smaller(main, other);
    long unsigned minsize = small.getsize();
    double sum = 0;
    for (long unsigned i=0; i< minsize; i++)
        sum += big[i]*small[i];
    return sum;
}

/** arithmetic operators on doubles **/
ntuple operator+(const ntuple& tup, double x)
{
    ntuple temp = tup;
    long unsigned size = temp.getsize();
    for (long unsigned i=0; i< size; i++)
        temp[i] += x;
    return temp;
}

ntuple operator+(double x, const ntuple& tup)
{
    ntuple temp = tup;
    long unsigned size = temp.getsize();
    for (long unsigned i=0; i< size; i++)
        temp[i] += x;
    return temp;
}

ntuple operator-(const ntuple& tup, double x)
{
    ntuple temp = tup;
    long unsigned size = temp.getsize();
    for (long unsigned i=0; i< size; i++)
        temp[i] -= x;
    return temp;
}

ntuple operator-(double x, const ntuple& tup)
{
    ntuple temp = tup;
    temp = -temp;
    long unsigned size = temp.getsize();
    for (long unsigned i=0; i< size; i++)
        temp[i] += x;
    return temp;
}

ntuple operator*(const ntuple& tup, double x)
{
    ntuple temp = tup;
    long unsigned size = temp.getsize();
    for (long unsigned i=0; i< size; i++)
        temp[i] *= x;
    return temp;
}

ntuple operator*(double x, const ntuple& tup)
{
    ntuple temp = tup;
    long unsigned size = temp.getsize();
    for (long unsigned i=0; i< size; i++)
        temp[i] *= x;
    return temp;
}

ntuple operator/(const ntuple& tup, double x)
{
    ntuple temp = tup;
    long unsigned size = temp.getsize();
    for (long unsigned i=0; i< size; i++)
        temp[i] /= x;
    return temp;
}


/** when trying to compare to a double, ntuple assumes you have a 1-tuple **/
bool operator==(ntuple tup, double x) { return tup[0] == x; }
bool operator==(double x, ntuple tup) { return tup[0] == x; }

bool operator<(ntuple tup, double x) { return tup[0] < x; }
bool operator<(double x, ntuple tup) { return x < tup[0]; }

bool operator>(ntuple tup, double x) { return tup[0] > x; }
bool operator>(double x, ntuple tup) { return x > tup[0]; }




/************************* FUNCTION DEFINITIONS *****************************/

function::function(const function& other)
{
    if (other.domain == NULL) domain = NULL;
    else
    {
        domain = new double[2];
        domain[0] = other.domain[0];
        domain[1] = other.domain[1];
    }

    fnd = other.fnd;
    domdim = other.domdim;
}

function& function::operator=(const function& other)
{
    delete[] domain;
    if (other.domain == NULL) domain = NULL;
    else
    {
        domain = new double[2];
        domain[0] = other.domain[0];
        domain[1] = other.domain[1];
    }
    fnd = other.fnd;
    domdim = other.domdim;

    return *this;
}

double function::operator()(ntuple intuple)
{
    /** if ntuple has too many arguments, this will ignore extras **/
    if (intuple.getsize() >= domdim)
        return fnd(intuple);
    /** if ntuple hasn't enough arguments, interprets remainder as 0 **/
    else
    {
        ntuple newtup(domdim);
        newtup += intuple;
        return fnd(newtup);
    }

}


/***************************** Piecewise Function DEFINITIONS ****************************/

/** copy constructor **/
pwfunction::pwfunction(const pwfunction& other)
{
    numfuncs = other.size();
    domdim = other.getdomdim();
    domain = new double[2];
    domain[0] = other.lobound();
    domain[1] = other.upbound();

    /** creates array of function pointers **/
    pieces = new function* [numfuncs];

    for (long unsigned i=0; i<numfuncs ; i++)
        pieces[i] = other.pieces[i];
}

/** assignment operator **/
pwfunction& pwfunction::operator=(const pwfunction& other)
{
    numfuncs = other.size();
    domdim = other.getdomdim();
    domain = new double[2];
    domain[0] = other.lobound();
    domain[1] = other.upbound();

    delete[] pieces;
    /** creates array of function pointers **/
    pieces = new function* [numfuncs];

    for (long unsigned i=0; i<numfuncs ; i++)
        pieces[i] = other.pieces[i];

    return *this;
}

/** blank constructor **/
pwfunction::pwfunction(long unsigned insize)
{
    numfuncs = insize;
    domdim = 1;
    pieces = new function*[insize];
}

pwfunction::pwfunction(function** inpieces, long unsigned innum)
{
    numfuncs = innum;
    domdim = 1;
    pieces = new function* [innum];
    for (long unsigned i=0; i<innum; i++)
        pieces[i] = inpieces[i];
}

function*& pwfunction::operator[](long unsigned n)
{
	if ( n >= 0 && n < size() )
        return pieces[n];
    else
    {
        std::cout<<std::endl<<"Piecewise Function hasn't a function for that index."<<std::endl;
    }
}

function* pwfunction::getfunc(long unsigned n) const
{
	if ( n >= 0 && n < size() )
        return pieces[n];
    else
    {
        std::cout<<std::endl<<"Piecewise Function hasn't a function for that index."<<std::endl;
    }
}


double pwfunction::operator()(double arg)
{
    /** runs through the possibility of evaluating at an endpoint **/
    for (long unsigned i=0; i<numfuncs ; i++)
        if (arg == pieces[i]->lobound() || arg == pieces[i]->upbound())
            return (*pieces[i])(arg) ;
    /** else, function checks in between the domains **/
    long unsigned domcount = 0, correctindex = 0;
    for (long unsigned i=0; i<numfuncs ; i++)
    {
        /** found a domain for which the function can be applied **/
        if (arg > pieces[i]->lobound() && arg < pieces[i]->upbound() )
        {
            domcount++;
            correctindex = i;
        }
        else continue;
    }
    if (domcount >= 2)
    {
        std::cout<< "Piecewise Function is not well defined. At least two of its domain"<<std::endl
        <<"intervals overlap."<<std::endl;
        return 0;
    }
    if (domcount == 0)
    {
        std::cout<< "Argument was not in any of the domains of the piecewise function."<<std::endl;
        return 0;
    }

    /** if we've passed all other flags then this will be the correct output **/
    return (*pieces[correctindex])(arg) ;
}



/************************* POLYNOMIAL DEFINITIONS **************************/
polynomial::polynomial()
{
    polycoeffs = new double[1];
    polycoeffs[0]=0;
    center = new double[1];
    center[0]=0;
    coefflength = 1;
    domdim = 1;
}

/** copy constructor **/
polynomial::polynomial(const polynomial& other)
{
    coefflength = other.size();
    domdim = other.getdomdim();

    polycoeffs = new double[coefflength];
    center = new double[coefflength];
    for (long unsigned i=0 ; i<coefflength ; i++)
    {
        polycoeffs[i] = other.polycoeffs[i];
        center[i] = other.center[i];
    }
}

/** assignment operator **/
polynomial& polynomial::operator=(const polynomial& other)
{
    coefflength = other.size();
    domdim = other.getdomdim();

    delete[] polycoeffs;
    delete[] center;
    polycoeffs = new double[coefflength];
    center = new double[coefflength];
    for (long unsigned i=0 ; i<coefflength ; i++)
    {
        polycoeffs[i] = other.polycoeffs[i];
        center[i] = other.center[i];
    }
    return *this;
}

polynomial::polynomial(long unsigned insize)
{
    coefflength = insize;
    domdim = 1;

    polycoeffs = new double[insize];
    center = new double[insize];
    for (long unsigned i=0; i<insize; i++)
    {
        polycoeffs[i]=0;
        center[i]=0;
    }
}

/** constructor from array without center **/
polynomial::polynomial(const double* incoeffs, long unsigned n)
{
    coefflength = n;
    domdim = 1;

    polycoeffs = new double[n];
    center = new double[n];
    for (long unsigned i = 0 ; i<n ; i++)
    {
        polycoeffs[i] = incoeffs[i];
        center[i] = 0;
    }
}

/** constructor from array with domain array **/
polynomial::polynomial(const double* incoeffs, long unsigned n, double* indomain)
{
    coefflength = n;
    domdim = 1;

    polycoeffs = new double[n];
    center = new double[n];
    for (long unsigned i = 0 ; i<n ; i++)
    {
        polycoeffs[i] = incoeffs[i];
        center[i] = 0;
    }

    domain = new double[2];
    domain[0] = indomain[0];
    domain[1] = indomain[1];
}


/** constructor from array with center **/
polynomial::polynomial(const double* incoeffs, const double* incenter, long unsigned n)
{
    coefflength = n;
    domdim = 1;

    polycoeffs = new double[n];
    center = new double[n];
    for (long unsigned i = 0 ; i<n ; i++)
    {
        polycoeffs[i] = incoeffs[i];
        center[i] = incenter[i];
    }
}

/** constructor from array with center and domain **/
polynomial::polynomial(const double* incoeffs, const double* incenter, long unsigned n, const double* indomain)
{
    coefflength = n;
    domdim = 1;

    polycoeffs = new double[n];
    center = new double[n];
    for (long unsigned i = 0 ; i<n ; i++)
    {
        polycoeffs[i] = incoeffs[i];
        center[i] = incenter[i];
    }

    domain = new double[2];
    domain[0] = indomain[0];
    domain[1] = indomain[1];
}

polynomial::polynomial(const double* incoeffs, double incenter, long unsigned n, const double* indomain)
{
    coefflength = n;
    domdim = 1;

    polycoeffs = new double[n];
    center = new double[n];
    for (long unsigned i = 0 ; i<n ; i++)
    {
        polycoeffs[i] = incoeffs[i];
        center[i] = incenter;
    }

    domain = new double[2];
    domain[0] = indomain[0];
    domain[1] = indomain[1];
}

double polynomial::operator()(double arg)
{
    long unsigned n = size();

    double output = 0;
    /**to evaluate the polynomial, this used Horner's Method which is just Synthetic Division
    that reduces the error by nesting multiplication. **/

    for (long unsigned i = n-1; i>=1; i--)
        output = output*(arg - center[i]) + polycoeffs[i];

    /** evaluate last term because unsigned long is always >=0 **/
    output = output*(arg - center[0]) + polycoeffs[0];

    return output;
}



/******************************* MATRIX DEFINITIONS *******************************/
//copy constructor
matrix::matrix(const matrix& orig)
{
    rowcount = orig.numrows();
    colcount = orig.numcols();

    mtx = new double* [rowcount];

    // Copy the original data into our data
    for (long unsigned i = 0; i < rowcount; i++)
    {
        mtx[i] = new double[colcount];
        for (long unsigned j = 0; j < colcount; j++)
            mtx[i][j] = orig.mtx[i][j];
    }
}

//assignment operator
matrix& matrix::operator=(const matrix& orig)
{
    for (long unsigned i = 0; i<rowcount; i++)
        delete [] mtx[i];
    delete [] mtx;


    rowcount = orig.numrows();
    colcount = orig.numcols();

    // Allocate new space for our copy of the original data
    mtx = new double* [rowcount];

    // Copy the original data into our data
    for (long unsigned i = 0; i < rowcount; i++)
    {
        mtx[i] = new double [colcount];
        for (long unsigned j = 0; j < colcount; j++)
            mtx[i][j] = orig.mtx[i][j];
    }

    return *this;
}



//index operator returns i-th row of matrix
double* matrix::operator[](long unsigned i) const
{
    if (i >= numrows() || i<0 )
    {
        std::cout<<std::endl<<"Tried to access index out of range of matrix."<<std::endl;
    }
    else return mtx[i];
}

//copies i-th column of matrix into incol
void matrix::pastecol(long unsigned i, double* incol) const
{
    if (i < numcols() && i > 0)
    {
        long unsigned m = numrows();

        for (long unsigned j = 0; j< m ; j++)
        {
            incol[j] = mtx[j][i];
        }
    }
    else
        std::cout<<std::endl<<"Tried to access index out of range of matrix."<<std::endl;
}

matrix::matrix(long unsigned inrow, long unsigned incol)
{
    rowcount = inrow;
    colcount = incol;

    mtx = new double* [inrow];
    for (long unsigned i = 0; i < inrow; i++)
    {
        mtx[i] = new double [incol];
        for (long unsigned j = 0; j < incol; j++)
            mtx[i][j] = 0;
    }
}


matrix::matrix(double** inmtx, long unsigned inrow, long unsigned incol)
{

    if (inrow == 0)
        std::cout<<std::endl<<"Trying to construct a matrix with no rows."<<std::endl;

    if (incol == 0)
        std::cout<<std::endl<<"Trying to construct a matrix with no columns."<<std::endl;



    rowcount = inrow;
    colcount = incol;

    mtx = new double* [inrow];

    for (long unsigned i = 0; i < inrow; i++)
    {
        mtx[i] = new double [incol];
        for (long unsigned j = 0; j < incol; j++)
            mtx[i][j] = inmtx[i][j];
    }

}

matrix::matrix(double* inmtx, long unsigned incol)
{
    if (incol == 0)
        std::cout<<std::endl<<"Trying to construct a matrix with no columns."<<std::endl;


    rowcount = 1;
    colcount = incol;

    mtx = new double* [1];
    mtx[0] = new double [incol];
    for (long unsigned j = 0; j < incol; j++)
        mtx[0][j] = inmtx[j];

}

matrix::matrix(std::ifstream& infile)
{
long unsigned m = 0, n = 0;
    std::string unused;
    double temp;


    //find number of rows
    while(getline(infile, unused)) m++;
    infile.clear();
    infile.seekg(0);


    if (m == 0)
    {
        std::cout<<std::endl<<"File is empty."<<std::endl;
        infile.clear();
        infile.seekg(0);
        return;
    }

    //find number of columns
    getline(infile, unused);
    std::stringstream instream(unused);
    while (instream >> temp) n++;
    infile.clear();
    infile.seekg(0);


    if (n == 0)
    {
        std::cout<<std::endl<<"First row of file is empty."<<std::endl;
        infile.clear();
        infile.seekg(0);
        return;
    }

    long unsigned i = 0, j = 0;

    rowcount = m;
    colcount = n;

    mtx = new double* [m];

    while(getline(infile,unused))
    {
        std::stringstream instream(unused);
        mtx[i] = new double[n];
        j=0;
        while(instream>> temp)
        {
            mtx[i][j] = temp;
            j++;
        }
        i++;
    }
}


void matrix::get_nodes()
{
    std::string filename;
    std::cout<< "Enter the path of the file containing the array.\n";
    std::cin>> filename;

    if (filename == "back") return;

    std::ifstream infile;

    infile.open(filename.c_str() );

    if (infile.fail())
    {
        std::cout<<std::endl<<"Failed to open file."<<std::endl;
        infile.clear();
        infile.seekg(0);
        return get_nodes();
    }

    long unsigned m = 0, n = 0;
    std::string unused;
    double temp;

    //find number of rows
    while(getline(infile, unused)) m++;
    infile.clear();
    infile.seekg(0);


    if (m == 0)
    {
        std::cout<<std::endl<<"File is empty."<<std::endl;
        infile.clear();
        infile.seekg(0);
        return get_nodes();
    }

    //find number of columns
    getline(infile, unused);
    std::stringstream instream(unused);
    while (instream >> temp) n++;
    infile.clear();
    infile.seekg(0);


    if (n == 0)
    {
        std::cout<<std::endl<<"First row of file is empty."<<std::endl;
        infile.clear();
        infile.seekg(0);
        return get_nodes();
    }

    long unsigned i = 0, j = 0;

    for (long unsigned i = 0; i<rowcount ; i++)
        delete[] mtx[i];
    delete[] mtx;

    rowcount = m;
    colcount = n;

    mtx = new double* [m];

    while(getline(infile,unused))
    {
        std::stringstream instream(unused);
        mtx[i] = new double[n];
        j=0;
        while(instream>> temp)
        {
            mtx[i][j] = temp;
            j++;
        }
        i++;
    }
    infile.close();
}

matrix matrix::augment(const ntuple& tup) const
{
    long unsigned m = numrows(), n = numcols();
    matrix aug(m,n+1);
    if (tup.getsize() != m)
    {
        std::cout<<std::endl<<"ERROR! Ntuple must have same length as number of rows."<<std::endl;
        return aug;
    }

    for (long unsigned i=0; i<m ; i++)
    {
        for (long unsigned j=0 ; j<n ; j++)
            aug[i][j] = mtx[i][j];
        aug[i][n] = tup.get(i);
    }

    return aug;
}

ntuple matrix::getrow(long unsigned i) const
{
    return ntuple(mtx[i], colcount);
}
ntuple matrix::getcol(long unsigned j) const
{
    ntuple col(rowcount);
    for (long unsigned i=0; i<rowcount; i++)
        col[i] = mtx[i][j];
    return col;
}


std::ostream& operator<<(std::ostream& out, matrix& A)
{
    long unsigned m = A.numrows(), n = A.numcols();
    for (long unsigned i = 0; i<m ; i++)
    {
        for (long unsigned j = 0; j<n ; j++)
            out<< A[i][j]<< " ";
        out<<std::endl;
    }
    return out;
}







/******************************** System of Functions Definitions **********************/

systemfuncs::systemfuncs(const systemfuncs& other)
{
    length = other.numfuncs();
    system = new function* [length];
    for (long unsigned i=0 ; i<length ; i++)
    {
        system[i] = new function;
        system[i] = other[i];
    }
}

systemfuncs& systemfuncs::operator=(const systemfuncs& other)
{
    delete[] system;

    length = other.numfuncs();
    system = new function* [length];
    for (long unsigned i=0 ; i<length ; i++)
    {
        system[i] = new function;
        system[i] = other[i];
    }

    return *this;
}

/** Empty Constructor **/
systemfuncs::systemfuncs(long unsigned n)
{
    length = n;
    system = new function* [n];
    for (long unsigned i=0 ; i<n ; i++)
        system[i] = NULL;
}

/** Base Constructor **/
systemfuncs::systemfuncs(function** infuncs, long unsigned insize)
{
    length = insize;
    system = new function* [insize];

    for (long unsigned i=0; i<insize ; i++)
    {
        if ( infuncs[i] != NULL )
            system[i] = new function(*infuncs[i]);
        else
            system[i] = NULL;
    }
}

/** Parenthetical Assignment **/

void systemfuncs::assign_func(long unsigned entry, function infunc)
{
    delete system[entry];
    system[entry] = dynamic_cast<function*> ( new function( infunc ) );
    *(system[entry]) = infunc;
}
void systemfuncs::assign_func(long unsigned entry, polynomial inpoly)
{
    delete system[entry];
    system[entry] = dynamic_cast<function*> ( new polynomial( inpoly ) );
    *(system[entry]) = inpoly;
}
void systemfuncs::assign_func(long unsigned entry, pwfunction inpw)
{
    delete system[entry];
    system[entry] = dynamic_cast<function*> ( new pwfunction( inpw ) );
    *(system[entry]) = inpw;
}



/** Find Max Dimension of all Functions **/
long unsigned systemfuncs::maxdim() const
{
    long unsigned max = (*this)[0] ->getdomdim();
    for ( long unsigned i=1; i<numfuncs() ; i++ )
        if (max < (*this)[i] ->getdomdim() ) max = (*this)[i]->getdomdim();
    return max;
}


/** index operator **/

function*& systemfuncs::operator[](long unsigned n) const
{
    if ( n >= 0 && n < numfuncs() )
        return system[n];
    else
        std::cout<<std::endl<<"System hasn't a function for that index."<<std::endl;
}

ntuple systemfuncs::operator() (ntuple x)
{
    long unsigned size = (*this).numfuncs();
    ntuple output(size);
    for (long unsigned i=0; i<size; i++)
        output[i] = (*this)(i, x);

    return output;
}


namespace library_of_predefined_simple_functions{

/// 1 dimensional functions
double func1_f(ntuple x)
{return ( -sin(x[0])-2*x[0] );}

double func2_f(ntuple x)
{return (cos(x[0])-x[0]*x[0]);}

double func3_f(ntuple x)
{return (exp(-2*x[0])*sin(x[0])+x[0]*x[0]*x[0]-10*x[0]*x[0]+3*x[0]-1);}

double func4_f(ntuple x)
{return (-2*exp(-2*x[0])*sin(x[0])+exp(x[0])*cos(x[0])+3*x[0]*x[0]-20*x[0]+3);}

double func5_f(ntuple x)
{return (log(x[0])*x[0] - x[0]*sin(x[0]));}

double func6_f(ntuple x)
{return (1-log(x[0]) - sin(x[0]) - x[0]*cos(x[0]));}

double func7_f(ntuple x)
{return (27*x[0]*x[0]*x[0]*x[0]*x[0]-100*x[0]*x[0]*x[0]*x[0]-50*x[0]*x[0]*x[0]+15*x[0]*x[0]+5*x[0]-1);}

double func8_f(ntuple x)
{return (135*x[0]*x[0]*x[0]*x[0]-400*x[0]*x[0]*x[0]-150*x[0]*x[0]+30*x[0]+5);}


/// n dimensional functions
double func9_f(ntuple x)
{return (x[0]+2*x[1]-3*x[2]);}

/// functions being input into system of eqns must have numfuncs+1 arguments?
double func10_f(ntuple x)
{ return ( -4*x[1] + 3*x[2] ); }

double func11_f(ntuple x)
{ return ( -2.4*x[1] + 1.6*x[2] ); }

double func12_f(ntuple x)
{ return ( x[0] + x[1] ); }

double func13_f(ntuple x)
{ return ( x[0] - x[1] ); }

double func14_f(ntuple x)
{ return ( 3*x[0] - cos( x[1]*x[2] ) -.5 ); }

double func15_f(ntuple x)
{ return ( x[0]*x[0] - 81*( x[1] + .1)*( x[1] + .1) + sin( x[2] ) + 1.06 ); }

double func16_f(ntuple x)
{ return ( exp( -x[0]*x[1] ) + 20*x[2] + ( 10*3.1415 - 3)/3 ); }

double func17_f(ntuple x)
{ return ( x[1] - x[0]*x[0] +1 ); }

double func18_f(ntuple x)
{ return ( x[2] ); }

double func19_f(ntuple x)
{ return ( exp(2*x[0])*sin(x[0]) - 2* x[1] + 2*x[2] ); }

function func1(func1_f, 1);
function func2(func2_f, 1);
function func3(func3_f, 1);
function func4(func4_f, 1);
function func5(func5_f, 1);
function func6(func6_f, 1);
function func7(func7_f, 1);
function func8(func8_f, 1);

function func9(func9_f, 3);
function func10(func10_f, 3);
function func11(func11_f, 3);

function func12(func12_f, 2);
function func13(func13_f, 2);

function func14(func14_f, 3);
function func15(func15_f, 3);
function func16(func16_f, 3);

function func17(func17_f, 2);

function func18(func18_f, 3);
function func19(func19_f, 3);

}

}
