#
#  Creates "DailyDayCentEVI" executable  (trace gas century)
#
# print variable
#  make print-VARIABLE
#
#   Fortran flags:
#   -Aa:    Generate errors for non-ASCII 77 standard features (with all exceptions)
#   -c: suppress loading and produce a .o file for each source file
#   -g: produce debug information for dbx or dbxtool
#   -u: make the default type "undefined"; do not allow implicit typing
#   -C: checks subscripts within array bounds - DO NOT USE THIS
#       OPTION due to FLOW code
#       -D: Compile debug statements (column 1 = D or d)
#       By default these statements are comments.
#   -K: Automatically SAVE all local vars in all subroutines.  All
#       non-initialized vars are set to zero (testing purposes only!).
#   -N: Disallow optimization
#   +T:
#       +FPflag: Specify run-time environment in floating pt. operations
#   +FPV:   Trap on invalid floating point operations
#   +FPZ:   Trap on divide by zero
#   +FPO:   Trap on overflow
#   +FPU:   Trap on underflow
# -ffpe-trap=invalid,zero,overflow : trap IEEE floating exceptions
#       #include <fenv.h>
#       int main(int argc, char **argv) {
#           float a, b, c;
#           feenableexcept(FE_DIVBYZERO | FE_INVALID | FE_OVERFLOW);
#
#   +ppu:   Append _ at end of defns of and refs to externally visible symbols
#       Necessary for linking FORTRAN and C code.
#   +U77:   Include U77 library (BSD 3F extensions)
#       -Xlist  do procedure listings and global program checking.
#               creates *.lst files with debugging information
#               Show errors, listing, and cross reference table
#       -XlistE Show errors
#       -Xlisth Terminate compilation if errors detected
#       -Xlistwar315   suppress "constant declared but never used warning message"
#       -ext_names=underscores - same as +ppu under HP-UX flags
#
# -fcaller-saves
# Enable values to be allocated in registers that will be clobbered by function
# calls, by emitting extra instructions to save and restore the registers around
# such calls. Such allocation is done only when it seems to result in better
# code than would otherwise be produced.
# This option is always enabled by default on certain machines, usually those
# which have no call-preserved registers to use instead.
# Enabled at levels -O2, -O3, -Os.
# OS X clang doesn't recognize!!

#   Modified to link in daily soilwater model. -mdh 10/94
#
# OS X portable executable.  Load the static versions of the fortran libs
# linker flags (OSX)
#  -no_compact_unwind        prevents 'warning: could not create compact unwind' messages
#    This tells the OS X  linker to not produce compact unwind information from the
#    dwarf unwind info in .o files. At runtime the darwin unwinder uses slow dwarf unwind info.
# -nodefaultlibs -l...      static fortran libs for portable OX X executable
#    nodefaultlibs stops linker from using default fortran library.
#    The libraries added are static fortran libs
#     We still end up using OS X's dynamic C libs, but that is what Apple wants
# LDFLAGS    = -Wl,-no_compact_unwind  -nodefaultlibs  -lgfortran -static-libgfortran -lquadmath -lgcc -lSystem

# -fno-underscoring Do not transform names of entities specified in
#                   the Fortran source file by appending underscores to them.
#
# -nodefaultlibs       Do not use the standard system libraries when linking.
# Only the libraries you specify are passed to the linker, and options specifying
# linkage of the system libraries, such as -static-libgcc or -shared-libgcc, are
# ignored. The standard startup files are used normally, unless -nostartfiles is used.
#
# The compiler may generate calls to memcmp, memset, memcpy and memmove. These entries
# are usually resolved by entries in libc. These entry points should be supplied
# through some other mechanism when this option is specified.
#
# NOTE: libquadmath seems to want to use the dynamic version unless you use nodefaults
#       The advice is to specify /usr/local/gfortran/lib/libquadmath.a explicitly byt
#       the GCC manual states the is the same except
#       -l surrounds 'library' with ‘lib’ and ‘.a’ and searches several directories.

 Centx    = DD17centEVI.exe
 L100     = DD17list100.exe
 Sitarc   = DD17sitarchive.exe
 LIBCX    = libDDCEVImingw32.a
 FLOW     = flowstak.f90

# directory for the executables When commented place them in the working directory
 DEST     = /data/paustian/Century/bin
# library directories for source code
Rdir   = 

# optimization level
 OPT        = -O3  # release
# OPT        =  -O0 -g # debug

# linker commands
LDFLAGS    = -static
# force static loads of the fortran/GCC libs Needed for transportable OX X
# LDFLAGS    = -Wl,-no_compact_unwind -nodefaultlibs -lgfortran -static-libgfortran -lquadmath -lgcc -lSystem

# get absolute path to this makefile
mkfilPath := $(abspath $(firstword $(MAKEFILE_LIST)))
Sdir      := $(dir $(mkfilPath))
Mdir      := $(dir $(mkfilPath))# source is the makefile directory
Mdir      := $(Mdir:%Win/=%)#  remove Win64 directory There should be NO code here
# Pdir      = $(realpath $(Mdir)..) # search parent directory for code

# source is the makefile directory
 S         := $(Sdir:%/=%)#    $(patsubst %/,%,$(dir $(mkfile_path)))

SLIB      = $(S) $(Mdir) $(Rdir) $(Pdir)  # $(Pdir)/..
INCLDIR   = $(patsubst %,%/Include,$(SLIB)) $(SLIB)
vpath %.c    $(SLIB)
vpath %.f    $(SLIB)
vpath %.f90  $(SLIB)
vpath %.h    $(INCLDIR)
vpath %.inc  $(INCLDIR)

PATH      := $(PATH):$(VPATH)

INCLUDES   += $(patsubst %,-I%,$(INCLDIR))

SUBDIRS    = # Win64 # 32bit

ILIBS     = -L/usr/lib
LIBS      = -lm
LOPT      = rcsUu

# For GNU compilers gfortran and gcc
 WARNC = -Wall
# normal "all" warnings without type conversion
# WARNF = -Wunused -Waliasing -Wampersand -Wsurprising -Wno-tabs -Wintrinsics-std -Wintrinsic-shadow
 WARNF = -Wunused -Waliasing -Wampersand -Wsurprising -Wno-tabs

# For Debugging with GNU compilers
FFLAGS       = -fno-underscoring ${WARNF}
#FFLAGS      = -fno-underscoring ${WARNF} -fcheck=bounds  -fbacktrace -ffpe-trap=overflow,zero,invalid,denormal #  debug

# For GNU ansi C
# For optimized GNU ansi C:
# Becky's compiler options
#CFLAGS    =  -c -g -ansi -I${SRC} -fbounds-check -fcaller-saves -Wall -lm
CFLAGS    =  ${WARNC} # -ansi

# GCC version.
#COMPGCC  = # use the standard installed GCC
COMPGCC  = /data/paustian/Century/bin/MinGW/bin/x86_64-w64-mingw32-#               unix 64 bit windows
#COMPGCC  = /data/paustian/Century/bin/MinGW2013/bin/x86_64-w64-mingw32-
#COMPGCC  = /usr/local/mingw-w32-bin_i686-darwin_20130414/bin/i686-w64-mingw32-#    OsX 32 bit windows
#COMPGCC  = /usr/local/mingw-w64-bin_i686-darwin_20130622/bin/x86_64-w64-mingw32-#  OsX 32 bit windows
FC     = $(COMPGCC)gfortran
CC     = $(COMPGCC)gcc
CPP    = $(COMPGCC)cpp
LD     = $(COMPGCC)gfortran
AR     = $(COMPGCC)ar


.SUFFIXES : .o  .f  .f90  .a  .c .inc

%.o: %.mod

.f.o:
	${FC} ${OPT} ${FFLAGS} $(INCLUDES) -c $<;
.f90.o:
	${FC} ${OPT} ${FFLAGS} $(INCLUDES) -c $<;
.c.o:
	${CC} ${OPT} ${CFLAGS} $(INCLUDES) -c $<;
.cpp.o:
	${CC} ${OPT} ${CFLAGS} $(INCLUDES) -c $<;
.f.a:
	${FC} ${OPT} ${FFLAGS} $(INCLUDES) -c $<; $(AR) $(LOPT) $@ $%; $(RM) $%
.f90.a:
	${FC} ${OPT} ${FFLAGS} $(INCLUDES) -c $<; $(AR) $(LOPT) $@ $%; $(RM) $%
.c.a:
	${CC} ${OPT} ${CFLAGS} $(INCLUDES) -c $<; $(AR) $(LOPT) $@ $%; $(RM) $%
.cpp.a:
	${CC} ${OPT} ${CFLAGS} $(INCLUDES) -c $<; $(AR) $(LOPT) $@ $%; $(RM) $%

TARARC=`date +%y%m%d%HH`

.PHONY: all
all: $(DEST)/$(Centx)  $(DEST)/$(L100) $(DEST)/$(Sitarc) $(DEST)/DDCtestwin.zip # subdirs

$(DEST)/DD17win.zip:	$(DEST)/$(Centx) $(DEST)/$(L100) $(DEST)/$(Sitarc)
	pushd $(DEST) ;  zip -q  $@ $(Centx) $(L100) $(Sitarc) ; popd

.DELETE_ON_ERROR:
.PHONY: subdirs $(SUBDIRS)
subdirs: $(SUBDIRS)
$(SUBDIRS):
	$(MAKE) -C $@

include $(Mdir)/MakeLib

# make print-VARIABLE
print-%  : ; @echo $* = $($*)

.INTERMEDIATE: %.o

$(DEST)/$(Centx): csa_main.f90 $(LIBCX) cflows.inc const.inc dovars.inc \
                  jday.inc monprd.inc param.inc pheno.inc plot1.inc plot2.inc \
                  plot3.inc timvar.inc wth.inc zztim.inc
	$(FC) $(LDFLAGS) $(OPT) $(FFLAGS) $(INCLUDES) $< $(LIBCX) $(LIBS) -o $@;

# .INTERMEDIATE: message.o
# message.o: message.f90
# 	${FC} ${OPT} ${FFLAGS} $(INCLUDES) -c $<;

$(DEST)/$(L100): list100.f90 message.o outval.inc table.inc \
                 ${LIBCX}(listvar.o) ${LIBCX}(retarg.o) ${LIBCX}(parse.o)
	${FC} $(LDFLAGS) $(OPT) ${FFLAGS} $(INCLUDES) $< message.o $(LIBCX) -o $@;

$(DEST)/$(Sitarc): sitarchive.f90 message.o ${LIBCX}(filext.o) \
                   ${LIBCX}(retarg.o) ${LIBCX}(parse.o) ${LIBCX}(chkdata.o) \
                   ${LIBCX}(writsit.o) ${LIBCX}(sitein.o) ${LIBCX}(listvar.o)
	${FC} $(LDFLAGS) $(OPT) ${FFLAGS} $(INCLUDES) $< message.o $(LIBCX) -o $@;

.PHONY: clean
clean:
	$(RM) $(LIBCX) *.o *.mod core $(DEST)/$(Centx) $(DEST)/$(L100) $(DEST)/$(Sitarc)
