#include <initializer_list>
#include <iostream>
#include <sstream>

#ifndef _PIXEL_
#define _PIXEL_

struct Pixel {
    unsigned char r = 0;
    unsigned char g = 0;
    unsigned char b = 0;

    Pixel(std::initializer_list<unsigned char> li) {
        auto it = li.begin();
        r = li.size() > 0 ? *(it) : 0;
        g = li.size() > 1 ? *(++it) : 0;
        b = li.size() > 2 ? *(++it) : 0;
    }

    Pixel& operator=(Pixel const& p) {
        r = p.r;
        g = p.g;
        b = p.b;

        return *this;
    }

    Pixel& operator=(std::initializer_list<unsigned char> const& li) {
        auto it = li.begin();
        r = li.size() > 0 ? *(it) : 0;
        g = li.size() > 1 ? *(++it) : 0;
        b = li.size() > 2 ? *(++it) : 0;

        return *this;
    }

    unsigned char& operator[](unsigned idx)
    {
        if ( idx == 0 ) return r;
        if ( idx == 1 ) return g;
        if ( idx == 2 ) return b;

        throw std::invalid_argument( "Invlaid pixel value requested" );
    }

    std::string to_string() {
        std::stringstream ss;
        ss << static_cast<unsigned>(r) << " " << static_cast<unsigned>(g) << " " << static_cast<unsigned>(b);
        return ss.str();
    }
};

#endif