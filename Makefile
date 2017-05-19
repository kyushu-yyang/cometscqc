
NAME   = SimQuench
TARGET = $(NAME).exe
CXX    = g++
VPATH  = src
MAIN   = main.cpp

LIBSRCS = $(wildcard $(VPATH)/*.cpp)
LIBOBJS = $(LIBSRCS:.cpp=.o)

MAINSRC = $(MAIN)
MAINOBJ = $(MAINSRC:.cpp=.o)

SRCS   = $(wildcard $(VPATH)/*.cpp) #$(MAIN)
OBJS   = $(SRCS:.cpp=.o)

ARFLAGS = rv

CXXFLAGS = -Wall -std=c++11 -O3 -I include
CXXLIBS  = 

ROOTFLAGS = $(shell root-config --cflags)
ROOTLIBS  = $(shell root-config --glibs)

BOOSTFLAGS = -I/opt/local/include
BOOSTLIBS  = -L/opt/local/lib #-lboost_python

PYFLAGS = $(shell python-config --include)
PYLIBS  = -L$(shell python-config --exec-prefix)/lib $(shell python-config --libs)

CXXFLAGS += $(ROOTFLAGS) $(BOOSTFLAGS) #$(PYFLAGS)
CXXLIBS  += $(ROOTLIBS) $(BOOSTLIBS) #$(PYLIBS)


.PHONY: all install clean clear
all: lib$(NAME).a

#$(TARGET): $(MAINOBJS) $(NAME).a
#	$(CXX) $(CXXLIBS) $^ -o $@

#$(TARGET): $(OBJS)
#	$(CXX) $(CXXLIBS) $^ -o $@

lib$(NAME).a: $(LIBOBJS)
	$(AR) $(ARFLAGS) $@ $^

.o:.cpp
	$(CXX) $(CXXFLAGS) -c $<


clean:
	$(RM) -rf $(OBJS) $(TARGET) lib$(NAME).a


install:
	mkdir -p bin
	mv -f $(TARGET)

clear:
	$(RM) -r build latex html *.log
