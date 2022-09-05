// Project Name: OpenMEEG (http://openmeeg.github.io)
// © INRIA and ENPC under the French open source license CeCILL-B.
// See full copyright notice in the file LICENSE.txt
// If you make a copy of this file, you must either:
// - provide also LICENSE.txt and modify this header to refer to it.
// - replace this header by the LICENSE.txt content.

#pragma once

#ifndef USE_PROGRESSBAR
#include <cmath>
#include <logger.h>
#endif

namespace OpenMEEG {
#ifndef USE_PROGRESSBAR
    class ProgressBar {
    public:

        ProgressBar(const unsigned n,const unsigned sz=20): max_iter(n),bar_size(sz) { }

        void operator++() {
            const unsigned p = std::min(static_cast<unsigned>(floor(static_cast<double>((bar_size+1)*iter++)/max_iter)),bar_size);
            if (p!=pprev && iter>1) {
                log_stream(PROGRESS) << std::string(bar_size+2,'\b') << '[' << std::string(p,'*') << std::string(bar_size-p,'.') << ']';
                pprev = p;
            }
            if (iter>=max_iter)
                log_stream(PROGRESS) << std::endl;
            log_stream(PROGRESS).flush();
        }
        
    private:
        
        unsigned iter = 0;
        unsigned pprev = -1;
        const unsigned max_iter;
        const unsigned bar_size;
    };
#else
    class ProgressBar {
    public:

        ProgressBar(const unsigned,const unsigned=20) { }
        void operator++() { }
    };
#endif
}
