CXX	= g++
MOSEK = /path/to/Mosek
CPPFLAGS = -I$(MOSEK)/h -I$(MOSEK)/include 			

WARNINGFLAGS = -W -Wall -pedantic -Wcast-qual -Wwrite-strings -Wnon-virtual-dtor -Wcast-align -Wconversion -Wno-char-subscripts -Wpointer-arith -Wundef

LDFLAGS = -Llib -L$(MOSEK)/bin -Wl,-rpath-link,$(MOSEK)/bin '-Wl,-rpath=$(MOSEK)/bin'

maxcut_partial: maxcut_partial.o
		$(CXX) -g maxcut_partial.o $(LDFLAGS) -lm -lfusion64 -lmosek64 -o maxcut_partial
	
maxcut_partial.o: maxcut_partial.cpp
		$(CXX) -std=c++11 $(CPPFLAGS) -O2 -g -c maxcut_partial.cpp -o maxcut_partial.o
