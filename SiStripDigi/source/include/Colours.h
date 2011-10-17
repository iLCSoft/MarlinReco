#ifndef COLOURS_H
#define COLOURS_H 1

namespace sistrip {

// If a terminal doesn't support any colors, comment following line
#define TERM_COLOR

#ifdef TERM_COLOR

#define DBLUE    "\x1b[36m"
#define DRED     "\x1b[31m"
#define DYELLOW  "\x1b[33m"
#define DGREEN   "\x1b[32m"
#define DUNDERL  "\x1b[4m"
#define ENDCOLOR "\x1b[m"

#else

#define DBLUE    ""
#define DRED     ""
#define DYELLOW  ""
#define DGREEN   ""
#define DUNDERL  ""
#define ENDCOLOR ""

#endif

} // Namespace

#endif // COLOURS_H
