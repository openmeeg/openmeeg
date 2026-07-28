// Project Name: OpenMEEG (http://openmeeg.github.io)
// © INRIA and ENPC under the French open source license CeCILL-B.
// See full copyright notice in the file LICENSE.txt
// If you make a copy of this file, you must either:
// - provide also LICENSE.txt and modify this header to refer to it.
// - replace this header by the LICENSE.txt content.

#pragma once

#include <atomic>
#include <iostream>
#include <NullStream.h>

namespace OpenMEEG {

    typedef enum { DEBUG, PROGRESS, INFORMATION, WARNING, ERROR } InfoLevel;

    class Logger {
    public:

        //  DEBUG is what the previously uninitialised, zero-initialised singleton got.

        Logger(): verbosity(DEBUG) { }

        void set_info_level(const InfoLevel level) { verbosity = level; }

        InfoLevel get_info_level() const { return verbosity; }

        bool is_verbose(const InfoLevel level) { return verbosity<=level; }

        static Logger& logger() {
            static Logger logger;
            return logger;
        }

    private:

        std::atomic<InfoLevel> verbosity; // set_info_level() can race with logging.
    };

    inline std::ostream&
    log_stream(const InfoLevel level) {
        //  thread_local: callers write to it, so a shared streambuf would race.

        static thread_local NullStream<char> nullstream;
        if (level==WARNING && Logger::logger().is_verbose(level))
            std::cout << "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!" << std::endl
                      << "!!!!!!!!!!! WARNING !!!!!!!!!!!" << std::endl
                      << "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!" << std::endl;
        return (Logger::logger().is_verbose(level)) ? std::cout : nullstream;
    }
}
