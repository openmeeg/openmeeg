// Project Name: OpenMEEG (http://openmeeg.github.io)
// © INRIA and ENPC under the French open source license CeCILL-B.
// See full copyright notice in the file LICENSE.txt
// If you make a copy of this file, you must either:
// - provide also LICENSE.txt and modify this header to refer to it.
// - replace this header by the LICENSE.txt content.

#pragma once

#include <ostream>
#include <streambuf>

namespace OpenMEEG {

    // Null buffer that store nothing.

    template<typename CharType,typename CharTraits=std::char_traits<CharType>>
    class NullBuffer: public std::basic_streambuf<CharType,CharTraits> {

        typedef std::basic_streambuf<CharType,CharTraits> base;

	public:

        using typename base::char_type;
        using typename base::int_type;
        using typename base::traits_type;

		NullBuffer() : std::streambuf() {}

        virtual std::streamsize
        xsputn(char_type const*,std::streamsize n) { return n; }

        virtual int_type
        overflow(int_type c=traits_type::eof()) { return traits_type::not_eof(c); }
    };

    // Null stream that outputs nothing.

    template<typename CharType,typename CharTraits=std::char_traits<CharType>>
    class NullStream: public std::basic_ostream<CharType,CharTraits> {

        typedef std::basic_ostream<CharType,CharTraits> base;

	public:

		NullStream(): base(&nullbuf) {}

	private:

		NullBuffer<CharType,CharTraits> nullbuf;
    };
}
