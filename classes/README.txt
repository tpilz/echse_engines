Contents of folder 'echse_classes':
===================================

This folder contains a library of classes for use with the ECHSE.


Sub-folder 'declaration':
-------------------------

Contains text files named like 'XXX.txt', where 'XXX' is the name of a class.
Each file contains the complete declaration of a class. The declaration
comprises the names and types of the class' data members.

These files represent the input of the code generator (see argument 'files' of
the 'generate' function in the R-package 'codegen').

For each text file, there is a corresponding sub-folder in 'implementation'
containing the C++ source code of the class' internals.


Sub-folder 'implementation':
----------------------------

Contains a larger number of sub-folders, each corresponding to a class
declaration file in folder 'declaration'. In every sub-folder, there are 3 C++
source code files where 'XXX' represents the class name:

  userCode_XXX_simulate.cpp   : Contains the interior (not the interface) of the
                                class' 'simulate' method. The statements in this
                                file are responsible for updating of state
                                variables and output variables.

  userCode_XXX_derivsScal.cpp : Contains the interior (not the interface) of the
                                class' 'derivsScal' method. The statements in
                                this file are responsible for the calculation of
                                the (scalar) state variables' derivatives with
                                respect to time. The file may be empty if the
                                simulation of the scalar state variables does
                                not require a numerical solution.

  userCode_XXX_aux.cpp        : Contains auxiliary code, most likely function
                                definitions, which is to be used within the
                                implementation of 'simulate' and 'derivsScal'.

