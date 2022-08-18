// Project Name: OpenMEEG (http://openmeeg.github.io)
// © INRIA and ENPC under the French open source license CeCILL-B.
// See full copyright notice in the file LICENSE.txt
// If you make a copy of this file, you must either:
// - provide also LICENSE.txt and modify this header to refer to it.
// - replace this header by the LICENSE.txt content.

#pragma once

#include <iostream>
#include <iomanip>
#include <string>
#include <vector>
#include <sstream>
#include <initializer_list>

#ifndef WIN32
#define use_color_terminal
#endif

#include "OpenMEEGConfigure.h"

#ifdef USE_OMP
#include <omp.h>
#endif

#include <filenames.h>

namespace OpenMEEG {

    class CommandLine {

        typedef std::vector<const char*> Strings;

        // Workaround a bug in old gcc compilers which does not allow the conversion of
        // const std::initializer_list<const char* const> to const Strings.

        typedef std::initializer_list<const char* const> List;

        static Strings build_strings(const List& list) {
            Strings strs;
            strs.reserve(list.size());
            for (const auto& item : list)
                strs.push_back(item);
            return strs;
        }

    public:

        CommandLine(const int argc,char* argv[],const std::string& usage=""): n(argc),args(argv) {
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

        char** option(const std::string& option,const Strings& parms,const std::size_t num_mandatory_parms) const {
            char** arg = find_argument(option);
            if (arg==end())
                return nullptr;

            const std::size_t num_parms = num_args(arg);
            if (num_parms<num_mandatory_parms) {
                std::cerr << "\'" << args[0] << "\' option \'" << option << "\' expects at least "
                          << num_mandatory_parms << " arguments (";
                if (parms.size()!=0) {
                    std::cerr << parms[0];
                    for (unsigned i=1; i<parms.size(); ++i)
                        std::cerr << ", " << parms[i];
                }
                std::cerr << ") and you gave only " << num_parms << " arguments." << std::endl;
                exit(1);
            }
            return arg;
        }

        char** option(const std::string& name,const Strings& parms) const { return option(name,parms,parms.size()); }

        char** option(const Strings& options,const Strings& parms) const {
            std::size_t num_mandatory_parms = parms.size();
            for (const auto& parm : parms)
                if (parm[0]=='[')
                    --num_mandatory_parms;

            char** mandatory_args = nullptr;
            for (const char* opt : options) {
                char** arg = option(opt,parms,num_mandatory_parms);
                if (arg!=nullptr) {
                    if (mandatory_args!=nullptr) {
                        std::cerr << "Warning: option " << *(options.begin()) << " provided multiple times!" << std::endl;
                        exit(1);
                    }
                    mandatory_args = arg;
                }
            }
            return mandatory_args;
        }

        // Workaround a bug in old gcc compilers which does not allow the conversion of
        // const std::initializer_list<const char* const> to const Strings.

        char** option(const std::string& name,const List& parms) const { return option(name,build_strings(parms));                   }
        char** option(const List& options,const List& parms)     const { return option(build_strings(options),build_strings(parms)); }

        // End of workaround.

        unsigned num_args(char** argument) const {
            unsigned res = 0;
            for (char** arg=argument+1; arg!=end(); ++arg,++res)
                if ((*arg)[0]=='-')
                    break;
            return res;
        }

        void print() const {
            std::cout << std::endl << "| ------ " << args[0] << std::endl;
            for (unsigned i=1; i<n; ++i)
                std::cout << "| " << args[i] << std::endl;
            std::cout << "| -----------------------" << std::endl;
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

        #ifndef WIN32
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

    // Some command have mutually exclusive options.
    // This helper function ensures that.
    // The counter num_options is counting options in each exclusive group.

    inline void assert_non_conflicting_options(const char* command,const unsigned num_options) {
        if (num_options!=1) {
            std::cerr << "Error: providing mutually exclusive options to " << command << "!" << std::endl;
            exit(1);
        }
    }
}
