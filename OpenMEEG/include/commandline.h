/*
Project Name : OpenMEEG

© INRIA and ENPC (contributors: Geoffray ADDE, Maureen CLERC, Alexandre
GRAMFORT, Renaud KERIVEN, Jan KYBIC, Perrine LANDREAU, Théodore PAPADOPOULO,
Emmanuel OLIVI
Maureen.Clerc.AT.inria.fr, keriven.AT.certis.enpc.fr,
kybic.AT.fel.cvut.cz, papadop.AT.inria.fr)

The OpenMEEG software is a C++ package for solving the forward/inverse
problems of electroencephalography and magnetoencephalography.

This software is governed by the CeCILL-B license under French law and
abiding by the rules of distribution of free software.  You can  use,
modify and/ or redistribute the software under the terms of the CeCILL-B
license as circulated by CEA, CNRS and INRIA at the following URL
"http://www.cecill.info".

As a counterpart to the access to the source code and  rights to copy,
modify and redistribute granted by the license, users are provided only
with a limited warranty  and the software's authors,  the holders of the
economic rights,  and the successive licensors  have only  limited
liability.

In this respect, the user's attention is drawn to the risks associated
with loading,  using,  modifying and/or developing or reproducing the
software by the user in light of its specific status of free software,
that may mean  that it is complicated to manipulate,  and  that  also
therefore means  that it is reserved for developers  and  experienced
professionals having in-depth computer knowledge. Users are therefore
encouraged to load and test the software's suitability as regards their
requirements in conditions enabling the security of their systems and/or
data to be ensured and,  more generally, to use and operate it in the
same conditions as regards security.

The fact that you are presently reading this means that you have had
knowledge of the CeCILL-B license and that you accept its terms.
*/

#pragma once

#include <iostream>
#include <iomanip>
#include <string>
#include <sstream>

#ifdef WIN32
#pragma warning( disable : 4530)    //MSVC standard library can't be inlined
#pragma warning( disable : 4996)    //MSVC warning C4996: declared deprecated
#pragma warning( disable : 4290)    //MSVC warning C4290: C++ exception specification
#else
#define use_color_terminal
#endif

#include "OpenMEEGConfigure.h"

#ifdef USE_OMP
#include <omp.h>
#endif

#include <filenames.h>

namespace OpenMEEG {

    class CommandLine {
    public:

        CommandLine(const int argc,char* argv[],const std::string& usage): n(argc),args(argv) {
            help = find_argument("-h")!=end() || find_argument("--help")!=end();
            if (help) {
                std::cerr << red << basename(args[0]) << normal;
                if (usage!="")
                    std::cerr << ": " << usage;
                std::cerr << std::endl << std::endl;
            }
        }

        bool help_mode() const { return help; }

        template <typename T>
        T option(const std::string& name,const T defaultvalue,const std::string usage) const {
            char** arg = find_argument(name);
            const T result = (arg==end()) ?  defaultvalue : parse_value(arg+1,defaultvalue);
            if (help)
                std::cerr << "    " << bold << std::left << std::setw(8) << name << normal
                          << " = " << std::left << std::setw(12) << value(result) << purple << usage << normal << std::endl;
            return result;
        }

        bool option(const std::string& name,const bool defaultvalue,const std::string usage) const {
            char** arg = find_argument(name);
            const bool result = (arg==end()) ?  defaultvalue : !defaultvalue;
            if (help)
                std::cerr << "    " << bold << std::left << std::setw(8) << name << normal
                          << " = " << std::left << std::setw(12) << value(result)  << purple << usage << normal << std::endl;
            return result;
        }

    private:

        template <typename T>
        static T value(const T val) { return val; }

        static std::string value(const bool val) { return (val) ? "true" : "false"; }

        static std::string value(const std::string& val) { return '"'+val+'"'; }

        char** end() const { return args+n; }

        char** find_argument(const std::string& name) const {
            for (auto arg = args; arg!=end(); ++arg)
                if (name==*arg)
                    return arg;
            return end();
        }

        template <typename T>
        T parse_value(char* arg[],const T defaultvalue) const {
            if (arg==end())
                return defaultvalue;
            std::istringstream iss(*arg);
            T value = defaultvalue;
            iss >> value;
            return value;
        }

        #if 1
        static constexpr char normal[9]      = { 0x1b, '[', '0', ';', '0', ';', '0', 'm', '\0' };
        static constexpr char red[11]        = { 0x1b, '[', '4', ';', '3', '1', ';', '5', '9', 'm', '\0' };
        static constexpr char bold[5]        = { 0x1b, '[', '1', 'm', '\0' };
        static constexpr char purple[11]     = { 0x1b, '[', '0', ';', '3', '5', ';', '5', '9', 'm', '\0' };
        static constexpr char path_delimiter = '/';
        #else
        static constexpr char normal[1]      = { '\0' };
        static constexpr char red[1]         = { '\0' };
        static constexpr char bold[1]        = { '\0' };
        static constexpr char purple[1]      = { '\0' };
        static constexpr char path_delimiter = '\\';
        #endif

        unsigned n;
        char**   args;
        bool     help;
    };

    inline void
    print_commandline(const int argc,char **argv) {
        std::cout << std::endl << "| ------ " << argv[0] << std::endl;
        for (int i=1;i<argc;++i)
            std::cout << "| " << argv[i] << std::endl;
        std::cout << "| -----------------------" << std::endl;
    }

    inline void print_version(const char* cmd) {
        #ifdef USE_OMP
            std::string omp_support = " using OpenMP\n Executing using " + std::to_string(omp_get_max_threads()) + " threads.";
        #else
            std::string omp_support = "";
        #endif

        std::ostringstream display_info;
        display_info << cmd;
        display_info << " version " << version;
        display_info << " compiled at " << __DATE__ << " " << __TIME__;
        display_info << omp_support;
        std::cout << display_info.str() << std::endl << std::endl;
    }

#if 0
    inline bool option(const char *const name, const int argc, char **argv,
                       const bool defaut, const char *const usage=NULL) 
    {
        const char *s = command_line::option(name, argc, argv, (const char*)NULL);
        const bool res = s?(command_line::strcasecmp(s,"false") && command_line::strcasecmp(s,"off") && command_line::strcasecmp(s,"0")):defaut;
        command_line::option(name, 0, NULL, res?"true":"false", usage);
        return res;
    }

    inline int option(const char *const name, const int argc, char **argv,
                      const int defaut, const char *const usage=NULL) 
    {
        const char *s = command_line::option(name, argc, argv, (const char*)NULL);
        const int res = s?std::atoi(s):defaut;
        char tmp[256];
        std::sprintf(tmp, "%d", res);
        command_line::option(name, 0, NULL, tmp, usage);
        return res;
    }

    inline char option(const char *const name, const int argc, char **argv,
               const char defaut, const char *const usage=NULL) 
    {
        const char *s = command_line::option(name, argc, argv, (const char*)NULL);
        const char res = s?s[0]:defaut;
        char tmp[8];
        tmp[0] = res;
        tmp[1] ='\0';
        command_line::option(name, 0, NULL, tmp, usage);
        return res;
    }

    inline double option(const char *const name, const int argc, char **argv,
             const double defaut, const char *const usage=NULL) 
    {
        const char *s = command_line::option(name, argc, argv, (const char*)NULL);
        const double res = s?command_line::atof(s):defaut;
        char tmp[256];
        std::sprintf(tmp, "%g", res);
        command_line::option(name, 0, NULL, tmp, usage);
        return res;
    }
#endif
}
