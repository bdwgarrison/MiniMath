#ifndef MENUS_H
#define MENUS_H


#include <iostream>
#include "functions.h"


namespace fnc = functions_and_mathematical_objects_for_numerical_analysis;
namespace menus_for_implementing_numerical_analysis_methods{



class menu
{
public:
    menu(std::string incomm, std::string inhelp) {commands = incomm; help = inhelp;}

    void display_commands() {std::cout<<commands<<std::endl;}
    void display_help() {std::cout<<help<<std::endl;}
    void bad_input();


    void get_input(std::string& in);

    long unsigned input_to_num(std::string& in);

    virtual void run(){}

private:

    std::string commands;
    std::string help;
};





class interface : public menu
{
public:
    //friend class definefuncs;

    interface() : menu("\n",
                       "\n\"#\" -- Creates a system of # functions."
                       " \n\"exit\" -- Ends the program. \n") {}

    void run();

private:
    std::string input;

    long unsigned numfuncs;
    fnc::systemfuncs funcs;
};





class definefuncs : public menu
{
public:
    definefuncs() : menu("\n\"lib\", \"interp\", \"diff\", \"integ\", \"prev\", \"next\", \"back\", \"help\", \"done\"\n",
                        "\nlib -- Select function from predefined set of simple functions.\n"
                        "interp -- Interpolate a function from a file containing data.\n"
                        "diff -- Approximate the derivative by numerically differentiating\n"
                        "         at a number of points across an interval, then interpolating.\n"
                        "integ -- Approximate the integral by solving the function as an ivp\n"
                        "back -- Return to previous selection.\n"
                        "prev -- Go back to the function before.\n"
                        "next -- Skip defining this function.\n") {}

    void run(fnc::systemfuncs& infuncs);

    void funclib();

private:
    fnc::systemfuncs funcs;
    long unsigned numfuncs;
    fnc::matrix interpnodes;

    void diff();
    void integ();
    void prev();
    void next();

    long unsigned i;
    std::string input;
};


class interpolate : public menu
{
public:
    interpolate() : menu("\"poly\", \"dpoly\", \"piece\", \"dpiece\", \"back\", \"help\"",
                         "\nSelect your method of interpolation.\n"
                        "    All methods expect a matrix whose first column are the values \n"
                        "    of the independent variable and second are those of the dependent. \n"

                        "poly -- For an input file with two columns. Function is defined \n"
                        "         past endpoints but error increases durastically. \n"
                        "dpoly -- Must have three colums, the third for the derivative. \n"
                        "         This method gains accuracy but still loses control after endpoints. \n"

                        "piece -- Creates a cubic function for every three points, then stores \n"
                        "        them in a piecewise function. Not defined after endpoints. \n"
                        "dpiece -- Works just like piece but with the input of the derivative at \n"
                        "        endpoints for added accuracy. \n") {}
    void run(fnc::systemfuncs&, long unsigned, const fnc::matrix&);

private:
    std::string input;
};


class solvefuncs : public menu
{
public:
    solvefuncs() : menu("\"sys\", \"single\", \"back\", \"help\" ",
                        "\nsys -- perform numerical methods on the system of functions as a whole. \n"
                        "single -- perform numerical methods on one of your functions. \n") {}

    void run(fnc::systemfuncs&);

    void solvesys();
    void solvesingle();

    void root_methods();

    void ivp_methods();
private:
    fnc::systemfuncs funcs;
    long unsigned numfuncs;

    long unsigned n;
    std::string input;
};







}


#endif
