// Project Name: OpenMEEG (http://openmeeg.github.io)
// © INRIA and ENPC under the French open source license CeCILL-B.
// See full copyright notice in the file LICENSE.txt
// If you make a copy of this file, you must either:
// - provide also LICENSE.txt and modify this header to refer to it.
// - replace this header by the LICENSE.txt content.

#pragma once

#include <type_traits>

namespace OpenMEEG {

    namespace Details {

        template<typename R,typename A>
        A argument_helper(R(*)(A));

        template<typename R,typename A>
        A argument_helper(R(*)(A) noexcept);

        template<typename R,typename F,typename A>
        A argument_helper(R(F::*)(A));

        template<typename R,typename F,typename A>
        A argument_helper(R(F::*)(A) const);

        template<typename R,typename F,typename A>
        A argument_helper(R(F::*)(A) noexcept);

        template<typename R,typename F,typename A>
        A argument_helper(R(F::*)(A) const noexcept);

        template <typename F>
        decltype(argument_helper(&F::operator())) argument_helper(F);
    }

    // A class to find the argument of a function, whatever is its implementation (function.
    // member or lambda).

    template <typename T>
    using arg_type = std::remove_reference<decltype(Details::argument_helper(std::declval<T>()))>::type;
}
