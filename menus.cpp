
#include "menus.h"
#include "functions.h"
#include "algorithms.h"
#include <sstream>


namespace mnu = menus_for_implementing_numerical_analysis_methods;
namespace fnc = functions_and_mathematical_objects_for_numerical_analysis;
namespace flb = fnc::library_of_predefined_simple_functions;
namespace alg = algorithms_for_numerical_analysis_implementation_on_function_objects;

void mnu::menu::get_input(std::string& in)
{
    std::cin>>in;
    if ( in == "help" ) return display_help();
    if ( in == "back" ) return;
}

void mnu::menu::bad_input()
{
    std::cout<<std::endl<< "Invalid input."<<std::endl;
}

long unsigned mnu::menu::input_to_num(std::string& in)
{
    std::stringstream stream(in);
    long unsigned num;
    stream>> num;
    return num;
}



/************************************** INTERFACE *******************************/


void mnu::interface::run()
{
    while (input!= "exit")
    {
        std::cout<<std::endl<<"Enter the number of functions on which you'll perform numerical methods: ";

        get_input(input);
        if ( input == "exit")
            break;


        numfuncs = input_to_num(input);
        if (numfuncs <=0)
            std::cout<<std::endl<<"Must enter a non-zero natural number. Type \"help\" for commands."<<std::endl;

        ///empty system of functions is created
        funcs = fnc::systemfuncs(numfuncs);
        ///system of functions is populated
        mnu::definefuncs nextstage;
        nextstage.run(funcs);
    }
    return;
}





/********************************* DEFINE FUNCS ****************************************/


void mnu::definefuncs::prev()
{
    if (i >0 )
        --i;
    else
        std::cout<< "Cannot go further. At first member of system."<<std::endl;
}

void mnu::definefuncs::next()
{
    if (i < numfuncs-1)
        ++i;
    else
        std::cout<< "Cannot go further. At last member of system."<<std::endl;
}


void mnu::definefuncs::run(fnc::systemfuncs& infuncs)
{
    funcs = infuncs;
    numfuncs = funcs.numfuncs();

    display_commands();

    i = 0;
    while (i < numfuncs)
    {
        std::cout<<std::endl<<"Get function "<<i<<" from: ";

        get_input(input);

        if (input == "back")
            return;
        else if (input == "lib")
        {
            funclib();
            ++i;
        }

        else if (input == "interp")
        {
            interpnodes.get_nodes();
            std::cout<< interpnodes;
            interpolate interp;
            interp.run(funcs, i, interpnodes);
            ++i;
        }

        else if (input == "diff")
            diff();

        else if (input == "integ")
            integ();



        /// go back a step
        else if (input == "back")
        {
            mnu::interface prevstage;
            return prevstage.run();
        }

        /// redefine previous function
        else if (input == "prev")
            prev();

        /// move on to next function
        else if (input == "next")
            next();

        else
            bad_input();
    }
    /// Move to next stage
    solvefuncs nextstage;
    nextstage.run(funcs);
	return;
}

/*** SELECT FUNCTION FROM LIBRARY **/

void mnu::definefuncs::funclib()
{
    std::cout<<std::endl<< "Please pick one of the following:"<<std::endl
    <<"1. -sin(x)-2x"<<std::endl
    <<"2. cos(x)-x^2"<<std::endl
    <<"3. exp(-2x)*sin(x)+x^3-10x^2+3x-1"<<std::endl
    <<"4. -2exp(-2x)*sin(x)+exp(x)*cos(x)+3x^2-20x+3"<<std::endl
    <<"5. log(x)*x - x*sin(x)"<<std::endl
    <<"6. 1-log(x) - sin(x) - x*cos(x)"<<std::endl
    <<"7. 27x^5-100x^4-50x^3+15x^2+5x-1"<<std::endl
    <<"8. 135x^4-400x^3-150x^2+30x+5"<<std::endl;

    short selection=0;
    std::string temp;
    std::cin>> temp;

    if (temp == "back") return;


    std::stringstream stream(temp);
    stream>>selection;

    if (selection==1) { funcs.assign_func(i, flb::func1); }
    else if (selection==2) { funcs.assign_func(i, flb::func2); }

    else if (selection==3) { funcs.assign_func(i, flb::func3); }
    else if (selection==4) { funcs.assign_func(i, flb::func4); }
    else if (selection==5) { funcs.assign_func(i, flb::func5); }
    else if (selection==6)  { funcs.assign_func(i, flb::func6); }
    else if (selection==7)  { funcs.assign_func(i, flb::func7); }
    else if (selection==8)  { funcs.assign_func(i, flb::func8); }
    else if (selection==9)  { funcs.assign_func(i, flb::func9); }
    else if (selection==10) { funcs.assign_func(i, flb::func10); }
    else if (selection==11) { funcs.assign_func(i, flb::func11); }
    else if (selection==12) { funcs.assign_func(i, flb::func12); }
    else if (selection==13) { funcs.assign_func(i, flb::func13); }
    else if (selection==14) { funcs.assign_func(i, flb::func14); }
    else if (selection==15) { funcs.assign_func(i, flb::func15); }
    else if (selection==16) { funcs.assign_func(i, flb::func16); }
    else if (selection==17) { funcs.assign_func(i, flb::func17); }
    else if (selection==18) { funcs.assign_func(i, flb::func18); }
    else if (selection==19) { funcs.assign_func(i, flb::func19); }

    else
    {
        std::cout<<std::endl<< "Enter 1-19"<< std::endl;
        return funclib();
    }
    return;
}

/** differentiate a function **/
void mnu::definefuncs::diff()
{
    if (i>0)
    {
        long unsigned n=0;
        if (i>=2)
        {
            std::cout<< "What previously defined function to differentiate? Enter 0 to "<< i-1<<std::endl;
            std::string temp;
            std::cin>> temp;
            n = input_to_num(temp);
        }

		alg::differentiate derv( funcs[n] );
        funcs.assign_func(i, derv.poly() );

        ++i;
    }
    else
        std::cout<<std::endl<<"Define at least one function to approximate derivative of."<<std::endl;
    return;
}


/** integrate a function **/
void mnu::definefuncs::integ()
{
    if (i>0)
    {
        long unsigned n=0;

        if (i>=2)
        {
            std::cout<< "What previously defined function to integrate? Enter 0 to "<< i-1<<std::endl;
            std::string temp;
            std::cin>> temp;
            n = input_to_num(temp);
        }

        alg::solve_ivp_fix ivp( funcs[n] );
        funcs.assign_func(i, ivp.integral_poly());
        ++i;
    }
    else
        std::cout<<std::endl<<"Define at least one function to treat as the integrand."<<std::endl;

    return;
}

/*** INTERPOLATE FUNCTION **/

void mnu::interpolate::run(fnc::systemfuncs& funcs, long unsigned i, const fnc::matrix& mtx)
{
    while (input != "back")
    {
        display_commands();

        get_input(input);

        ///quit
        if (input == "back")
            return;

        else if (input == "poly")
        {

            alg::interp_basicpoly interp( mtx );
            funcs.assign_func(i, interp.getsoln());
            return;
        }

        else if (input == "dpoly")
        {

            alg::interp_dervpoly interp(mtx);
            funcs.assign_func(i, interp.getsoln());

            return;
        }

        else if (input == "piece")
        {
            alg::interp_basicpw interp(mtx);
            funcs.assign_func(i, interp.getsoln());
            return;
        }

        //requires first and last row to have 3 columns
        else if (input == "dpiece")
        {

            double fprm0, fprmn;
            std::cout<< "What is the value of the derivative at the first node?\n";
            std::cin>> fprm0;
            std::cout<< "At the last?\n";
            std::cin>> fprmn;

            alg::interp_dervpw interp(mtx, fprm0, fprmn);
            funcs.assign_func(i, interp.getsoln());
            return;
        }

        else
            bad_input();
    }
}


/*********************** SOLVE FUNCTIONS ********************/

void mnu::solvefuncs::run(fnc::systemfuncs& infuncs)
{
    funcs = infuncs;
    numfuncs = infuncs.numfuncs();

    while ( input != "back" )
    {
        display_commands();
        get_input(input);

        if ( input == "sys" )
            solvesys();

        else if ( input == "single" )
            solvesingle();


        else if (input == "back")
        {
            definefuncs prevstep;
            return prevstep.run(funcs);
        }
        else
            bad_input();
    }
}

/*** solve system of functions **/
void mnu::solvefuncs::solvesys()
{
    n=0;
    std::string temp;
    while ( temp != "back" )
    {
        std::cout<<std::endl<< "\"eval\", \"lin\", \"nonlin\", \"ivp\", \"back\",\"help\""<<std::endl;
        std::cin>> temp;

        if (temp == "help")
        {
            std::cout<< "eval -- Input a point to evaluate the system of functions at."<<std::endl
                << "lin -- If ALL the equations in your system are linear, then this"<<std::endl
                << "        solves the system Ax = b, where A is your system and b you"<<std::endl
                << "        select, for x."<<std::endl
                << "nonlin -- Returns the roots of any non-linear system. F(x) = 0 "<<std::endl
                << "ivp -- Solves if your equations are initial value problem ODE's."<<std::endl;
        }

        else if (temp == "back") return;

        else if (temp == "eval")
        {
            fnc::ntuple x(funcs.maxdim());
            std::cout<< "Evaluate at point: ";
            std::cin>> x;
            std::cout<< "System returned "<< funcs(x)<<std::endl;
        }

        else if (temp == "lin")
        {
            fnc::ntuple b;
            std::cout<< "Enter \"done\" when finished entering coordinates for vector."<<std::endl
                << "Solve matrix Ax=b for b: ";
            std::cin>> b;

            alg::linsys_mtx mtxrep(funcs);
            alg::solve_linsys solve(mtxrep.getsoln(), b);
            std::cout<< solve.getsoln();
        }

        else if (temp == "nonlin")
        {
            alg::solve_nonlinsys sys(funcs);
            std::cout<< "System returned "<<std::endl<< sys.getsoln() <<std::endl;
        }

        else if (temp == "ivp")
        {
            alg::solve_ivp_sys sys(funcs);
            std::cout<< "System returned "<<std::endl<< sys.getsoln();
        }
    }
}

/*** solve single function **/

void mnu::solvefuncs::solvesingle()
{
    std::string temp;
    n = 0;

    while ( temp != "back" )
    {
        if (numfuncs >= 2)
        {
            std::cout<<std::endl<< "Choose function from 0 to "<< numfuncs-1<< " to perform methods on.\n";
            std::cin>> temp;

            if (temp == "back")
                return;

            std::stringstream str(temp);
            str >> n;
            if (n>=0 && n<numfuncs)
            {
                std::cout<<std::endl<<"Chose function "<<n<<"."<<std::endl;
                break;
            }
            else
            {
                std::cout<<std::endl<<"No function at that location."<<std::endl;
                return solvesingle();
            }
        }
    }
    while ( temp != "back" )
    {
        std::cout<<std::endl<<"\"root\", \"eval\", \"diff\", \"integ\", \"ivp\", \"back\", \"help\" "<<std::endl;

        std::cin>> temp;
        if (input == "help")
        {
            std::cout<< "root -- Select an algorithm to find the zeros of your function."<<std::endl
            <<"         (requires function to be f:R1->R1)"<<std::endl
            <<"eval -- Evaluate your function at a given point."<<std::endl
            <<"         (enter \"done\" once all points have been selected)"<<std::endl
            <<"diff -- Numerically differentiate your function on an interval."<<std::endl
            <<"integ -- Numerically integrate your function on an interval."<<std::endl
            <<"ivp -- Select an algorithm to solve your function as an ivp"<<std::endl
            <<"         (requires function to be f:R2->R1)"<<std::endl
            <<"         (will use every function in your system)"<<std::endl;
        }

        /** Evaluate **/
        else if (temp == "eval")
        {
            fnc::ntuple x( funcs[n]->getdomdim() );
            std::cout<< "Evaluate at point: ";
            std::cin>> x;
            std::cout<< "Function returned "<< funcs(n, x) <<std::endl;
        }

        /** Solve **/
        else if (temp == "root")
            root_methods();

        /** Numerically Differentiate **/
        else if (temp == "diff")
        {
            alg::differentiate derv(funcs[n]);
            std::cout<< derv.getsoln();
        }

        /** Numerically Integrate **/
        else if (temp == "integ")
        {
            alg::integrate_def integ( funcs[n] );
            std::cout<< "Integral of function is "<< integ.getsoln()<<std::endl;
        }

        /** IVP **/
        //else if (temp == "ivp")
            //ivp( funcs[n] );

        /** Back **/
        else if (temp == "back")
            return solvesingle();

        else
            bad_input();
    }
}

/*** various root finding methods **/
void mnu::solvefuncs::root_methods()
{
    std::cout<< "\"bis\", \"newt\", \"sec\", \"reg\", \"back\", \"help\""<<std::endl;

    std::string temp;
    std::cin>> temp;

    if (temp == "back")
        return;

    else if (temp == "help")
    {
        std::cout<< "Select one of the following to solve y=f(x) for a root:\n"<<std::endl
        <<"bis -- Bisection Method"<<std::endl
        <<"         (requires an interval on which f has opposite signs at the endpoints)\n"<<std::endl
        <<"newt -- Newton's Method"<<std::endl
        <<"         (requires a function for the first derivative and an initial approximation)\n"<<std::endl
        <<"sec --  Secant Method\n"<<std::endl
        <<"         (requires two initial approximations to the root)"<<std::endl
        <<"reg -- Regula Falsi"<<std::endl
        <<"         (requires an interval on which f has opposite signs at the endpoints)\n"<<std::endl;

        return root_methods();
    }


    //Bisection Method
    else if (temp == "bis")
    {
        alg::bisection bis;
        bis( funcs[n] );
    }
    //Newton's Method
    else if (temp == "newt")
    {
        long unsigned df;
        std::string selection;

        std::cout<< "Which previously defined function will you use as your derivative?"<<std::endl;
        std::cin>> selection;

        if (selection == "back")
            return root_methods();
        else
        {

            std::stringstream stream(selection);
            stream>> df;

            if (df >=0 && df < funcs.numfuncs())
            {
                alg::newton newt;
                newt(funcs[n], funcs[df]);
            }

        }
    }
    //Secant Method
    else if (temp == "sec")
    {
        alg::secant sec;
        sec( funcs[n] );
    }
    //Regula Falsi
    else if (temp == "reg")
    {
        alg::regulafalsi reg;
        reg( funcs[n] );
    }

    else
    {
        bad_input();
        return root_methods();
    }
}

/** various initial value problem methods **/
void mnu::solvefuncs::ivp_methods()
{
    std::string temp;
    while (temp != "back")
    {
        std::cout<< "\"fixed\", \"adapt\", \"multi\", \"back\", \"help\""<<std::endl;
        std::cin>> temp;

        if (temp == "back")
            return;
        else if (temp == "help")
        {
            std::cout<< "Select one of the following to solve y'=f(t,y); y(0)=y0 for the function, y:\n"<<std::endl
            <<"fixed -- Solve with fixed step size method (runge-kutta order 4)."<<std::endl
            <<"adapt -- Solve with step size that adjusts to error (runge-kutta fehlberg)."<<std::endl
            <<"multi -- Solve with multi-step predictor-corrector (adams)"<<std::endl;
        }

        if (temp == "fixed")
        {
            alg::solve_ivp_fix ivp( funcs[n] );
            std::cout<< ivp.getsoln();
        }
        else if (temp == "adapt")
        {
            alg::solve_ivp_adpt ivp( funcs[n] );
            std::cout<< ivp.getsoln();
        }
        else if (temp == "multi")
        {
            alg::solve_ivp_mult ivp( funcs[n] );
            std::cout<< ivp.getsoln();
        }
        else
            bad_input();
    }

}





