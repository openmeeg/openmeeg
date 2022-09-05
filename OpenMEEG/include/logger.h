// Project Name: OpenMEEG (http://openmeeg.github.io)
// © INRIA and ENPC under the French open source license CeCILL-B.
// See full copyright notice in the file LICENSE.txt
// If you make a copy of this file, you must either:
// - provide also LICENSE.txt and modify this header to refer to it.
// - replace this header by the LICENSE.txt content.

#pragma once

#include <iostream>
#include <NullStream.h>

namespace OpenMEEG {

    typedef enum { PROGRESS, INFORMATION, WARNING, ERROR, DEBUG } InfoLevel;

    class Logger {
    public:

        Logger() { }

        void set_info_level(const InfoLevel level) { verbosity = level; }

        bool is_verbose(const InfoLevel level) { return verbosity<=level; }

        static Logger& logger() {
            static Logger logger;
            return logger;
        }

    private:

        InfoLevel verbosity;
    };

    inline std::ostream&
    log_stream(const InfoLevel level) {
        static NullStream<char> nullstream;
        return (Logger::logger().is_verbose(level)) ? std::cout : nullstream;
    }
}
