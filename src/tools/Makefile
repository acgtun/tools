ifndef SRC_ROOT
SRC_ROOT=./../
endif

ifndef ACGTUN_UTIL
ACGTUN_UTIL=./../acgtun_util
endif

PROGS = compare-read-mapping-results \
	get-nucleotide-from-chromosome-by-positoin \
	get-time-memory-from-qsub-output-file

CXX = g++
CFLAGS = -Wall -O3 -fPIC -fmessage-length=50 -std=c++0x

TEST_DIR = $(SRC_ROOT)/test

ifeq "$(shell uname)" "Darwin"
CFLAGS += -arch x86_64
endif

all:	$(PROGS)

%.o: %.cpp %.hpp
	$(CXX) $(CFLAGS) -c -o $@ $<  -I$(ACGTUN_UTIL)

compare-read-mapping-results :  \
	$(addprefix $(ACGTUN_UTIL)/, option.o)

get-nucleotide-from-chromosome-by-positoin : \
	$(addprefix $(ACGTUN_UTIL)/, option.o)

get-time-memory-from-qsub-output-file : \
	$(addprefix $(ACGTUN_UTIL)/, option.o)

%: %.cpp
		$(CXX) $(CFLAGS) -o $@ $^ -I$(ACGTUN_UTIL)
install: all
	@mkdir -p $(SRC_ROOT)/bin
	@install -m 755 $(PROGS) $(SRC_ROOT)/bin

test_%:	%
	@$(TEST_DIR)/$@ $(TEST_DIR)

test:	$(addprefix test_, $(PROGS))

clean:
	@-rm -f $(PROGS) *.o *.so *.a *~

.PHONY: clean
