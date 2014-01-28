// $File: chp.i $
// $LastChangedDate:  $
// $Rev:  $
// Copyright (c) 2014, Gao Wang <ewanggao@gmail.com>
// GNU General Public License (http://www.gnu.org/licenses/gpl.html)

%module chp

%{
#include "Exception.hpp"
#include "Core.hpp"
#include "chp.hpp"
%}

%include exception.i

%exception
{
    try
    {
        $function
    }
    catch(SEQLinco::IndexError e)
    {
        SWIG_exception(SWIG_IndexError, e.message());
    }
    catch(SEQLinco::ValueError e)
    {
        SWIG_exception(SWIG_ValueError, e.message());
    }
    catch(SEQLinco::SystemError e)
    {
        SWIG_exception(SWIG_SystemError, e.message());
    }
    catch(SEQLinco::RuntimeError e)
    {
        SWIG_exception(SWIG_RuntimeError, e.message());
    }
    catch(...)
    {
        SWIG_exception(SWIG_UnknownError, "Unknown runtime error happened.");
    }
}


%newobject *::clone;

%include "stl.i"
%include "std_vector.i"
%include "std_string.i"
%include "std_map.i"

%template(VecString) std::vector<std::string>; 
%template(VecVecString) std::vector<std::vector<std::string> >; 
%template(VecVecVecString) std::vector<std::vector<std::vector<std::string> > >; 

%include "Exception.hpp"
%include "Core.hpp"
%include "chp.hpp"
