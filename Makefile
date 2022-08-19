include config.txt

CXX = mpic++
CXXFLAGS = -std=c++17 -O3 -fopenmp -DLMER_LENGTH=8 -DBIDIRECT
INCLUDES = -I$(ZOLTAN_INSTALL)/include -I$(ZOLTAN_HOME)/src/include -I$(PWD)/includes

GET_WINSIZE=1
KVAL=$(shell echo $(ksize)-$(GET_WINSIZE) | bc)

ifeq ($(shell expr $(KVAL) \<= 0), 1)
      KVAL=30
endif

ifeq ($(shell echo $(ksize) \> 64), 1)
$(error ksize is greater than 64)
endif

ifeq ($(shell expr $(KVAL) \<= 31), 1)
  CXXFLAGS += -DWINDW_SIZE=$(KVAL)
else
  CXXFLAGS += -DWINDW_SIZE=$(KVAL) -DEXTEND_KMER
endif 

LDFLAGS = -L$(ZOLTAN_INSTALL)/src -lzoltan
LDFLAGS += -L$(PARMETIS_INSTALL)/lib -lparmetis -lmetis
LDFLAGS += -lm -ldl -lstdc++fs

ENABLE_DEBUG=0
ifeq ($(ENABLE_DEBUG), 1)
  CXXFLAGS+=-fsanitize=address -g
  LDFLAGS+=-fsanitize=address
endif
 
ENABLE_DUMP_DEBUG_DATA=0
ifeq ($(ENABLE_DUMP_DEBUG_DATA), 1)
  CXXFLAGS += -DDEBUG_KSIZE -DDEBUG_WIRE -DDEBUG_IDSET
endif
    
SRCDIR := src
SRCS = $(sort $(wildcard $(SRCDIR)/*.cpp))
OBJSDIR := .objs
OBJS = $(patsubst $(SRCDIR)/%.cpp, $(OBJSDIR)/%.o, $(SRCS))
TARGET = pbucketing
DEPSDIR := .deps
MAKES := Makefile config.txt

DEPGEN = -MT $@ -MMD -MP -MF $(DEPSDIR)/$*.d

.PHONY: all clean distclean
all: $(TARGET)

$(OBJS): $(OBJSDIR)/%.o: $(SRCDIR)/%.cpp $(DEPSDIR)/%.d $(MAKES) | $(DEPSDIR) $(OBJSDIR)
	$(CXX) $(DEPGEN) $(CXXFLAGS) $(INCLUDES) -c -o $@ $<

$(DEPSDIR):
	@mkdir -p $@
$(OBJSDIR):
	@mkdir -p $@

DEPSINCS = $(patsubst $(SRCDIR)/%.cpp, $(DEPSDIR)/%.d, $(SRCS))
$(DEPSINCS):

include $(sort $(wildcard $(DEPSINCS)))

$(TARGET): $(OBJS)
	$(CXX) -o $@ $^ $(LDFLAGS) -fopenmp

clean:
	-rm -rvf $(OBJSDIR)
	-rm -rvf $(DEPSDIR)
	-rm -f $(TARGET)

distclean: clean
