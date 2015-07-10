/**
 * \file logog.hpp Main include file for logog logging functionality.  Include this file to enable logging for your application.
 */

#ifndef __LOGOG_HPP__
#define __LOGOG_HPP__

#include <cstdio>

#define DBUG(...) std::printf( __VA_ARGS__ )
#define INFO(...) std::printf( __VA_ARGS__ )
#define WARN3(...) std::printf( __VA_ARGS__ )
#define WARN2(...) std::printf( __VA_ARGS__ )
#define WARN1(...) std::printf( __VA_ARGS__ )
#define WARN(...) std::printf( __VA_ARGS__ )
#define ERR(...) std::printf( __VA_ARGS__ )
#define ALERT(...) std::printf( __VA_ARGS__ )
#define CRITICAL(...) std::printf( __VA_ARGS__ )
#define EMERGENCY(...) std::printf( __VA_ARGS__ )

#endif // __LOGOG_HPP_
