ifndef SRC_ROOT
SRC_ROOT=./../
endif

ifndef SMITHLAB_CPP
SMITHLAB_CPP=./../smithlab_cpp
endif

PROGS = compare-read-mapping-results get_nucleotide_from_chromosome_by_position

CXX = g++
CFLAGS = -Wall -O3 -fPIC -fmessage-length=50 -std=c++0x

TEST_DIR = $(SRC_ROOT)/test

ifeq "$(shell uname)" "Darwin"
CFLAGS += -arch x86_64
endif

SRC = compare-read-mapping-results.cpp \
	  get_nucleotide_from_chromosome_by_position.cpp
OBJ = $(patsubst %.cpp,%.o,$(SRCDB))


all:	$(PROGS)

%.o: %.cpp %.hpp
	$(CXX) $(CFLAGS) -c -o $@ $<  -I$(SMITHLAB_CPP)

makedb : $(OBJDB) \
	$(addprefix $(SMITHLAB_CPP)/, smithlab_os.o smithlab_utils.o OptionParser.o)

bsmapper : $(OBJMAP) \
	$(addprefix $(SMITHLAB_CPP)/, smithlab_os.o smithlab_utils.o OptionParser.o)

%: %.cpp
		$(CXX) $(CFLAGS) -o $@ $^ -I$(SMITHLAB_CPP)
install: all
	@mkdir -p $(SRC_ROOT)/bin
	@install -m 755 $(PROGS) $(SRC_ROOT)/bin

test_%:	%
	@$(TEST_DIR)/$@ $(TEST_DIR)

test:	$(addprefix test_, $(PROGS))

clean:
	@-rm -f $(PROGS) *.o *.so *.a *~

.PHONY: clean
