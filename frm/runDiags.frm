#-
off statistics;
#include declarations.h
#include config.inc
.global

#ifndef `diagNum'
#define diagNum "0"
#endif

#define LOOPS "`leftNLoops'"
*#call evaluateDiags(`leftDiagDefFile',`leftMarker',`diagNum',`numPieces',1)
#call evaluateDiags(`leftDiagDefFile',`leftMarker',`diagNum',`numPieces',0)
.clear

#-
off statistics;
#include declarations.h
#include config.inc

#ifndef `diagNum'
#define diagNum "0"
#endif
#define LOOPS "`rightNLoops'"
.global

*#call evaluateDiags(`rightDiagDefFile',`rightMarker',`diagNum',1,0)
#call evaluateDiags(`rightDiagDefFile',`rightMarker',`diagNum',1,1)

.end
