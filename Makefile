#Insert the path to the folder with muparser and json libraries
PACS_ROOT ?= ../../pacs-examples/Examples

CXX       = mpic++
CXXFLAGS ?= -std=c++20
CPPFLAGS ?= -fopenmp -O3 -Wall -pedantic -Iinclude -I${PACS_ROOT}/include

LDFLAGS ?= -L${PACS_ROOT}/lib
LIBS  ?= -lmuparser

DEPEND = make.dep

EXEC = main 
SRCS = $(wildcard *.cpp) $(wildcard src/*.cpp)
OBJS = $(SRCS:.cpp=.o)

.PHONY = all $(EXEC) $(OBJS) clean distclean $(DEPEND)

all: $(DEPEND) $(EXEC)

$(EXEC): $(OBJS)
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) $(LDFLAGS) $^ $(LIBS) -o $@

$(OBJS): %.o: %.cpp
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -c $< -o $@

clean:
	$(RM) *.o

distclean: clean
	$(RM) $(EXEC)
	$(RM) *.csv *.out *.bak *~

$(DEPEND): $(SRCS)
	@$(RM) $(DEPEND)
	@for file in $(SRCS); \
	do \
	  $(CXX) $(CPPFLAGS) $(CXXFLAGS) -MM $${file} >> $(DEPEND); \
	done

-include $(DEPEND)