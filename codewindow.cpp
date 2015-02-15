#include "codewindow.h"
#include "menus.h"
#include "functions.h"
#include "algorithms.h"
#include <QTextEdit>
namespace mnu = menus_for_implementing_numerical_analysis_methods;
namespace fnc = functions_and_mathematical_objects_for_numerical_analysis;
namespace flb = fnc::library_of_predefined_simple_functions;
namespace alg = algorithms_for_numerical_analysis_implementation_on_function_objects;

CodeWindow::CodeWindow(QObject *parent) :
    QTextEdit()
{
}
