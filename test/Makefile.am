AM_CXXFLAGS = -I$(top_srcdir)/src/lib -DDATADIR="\"$(abs_top_srcdir)/test\"" @AM_CXXFLAGS@
check_PROGRAMS = unittest

EXTRA_DIST = \
A_basic.csv \
x_after_projection.csv \
x_to_project.csv \
x0_basic.csv

unittest_SOURCES = \
test_methods.cpp

unittest_LDADD = -lunittest -larmadillo $(AM_LDFLAGS)
