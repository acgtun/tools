ifndef SRC_ROOT
SRC_ROOT=./../
endif

ifndef ACGTUN_UTIL
ACGTUN_UTIL=./../acgtun_util
endif
ifndef SMITHLAB_CPP
SMITHLAB_CPP=./../smithlab_cpp
endif

PROGS = accuracy-for-single-end \
	accuracy-for-paired-end \
	accuracy-get-best-pairs \
	accuracy-get-best-single \
	mrsFast-count-positions-on-diff-mismatch2 \
	mrsFast-sum-count-positions-on-diff-mismatch3 \
	mrsFast-compress1 \
	mrsFast-get-1000-results4 \
	count-mismatch-for-mapped-reads-SAM


CXX = g++
CFLAGS = -Wall -O3 -fPIC -fmessage-length=50 -std=c++0x

TEST_DIR = $(SRC_ROOT)/test

ifeq "$(shell uname)" "Darwin"
CFLAGS += -arch x86_64
endif

all:	$(PROGS)

%.o: %.cpp %.hpp
	$(CXX) $(CFLAGS) -c -o $@ $<  -I$(ACGTUN_UTIL) -I$(SMITHLAB_CPP)

count-mismatch-for-mapped-reads-SAM :	\
	$(addprefix $(ACGTUN_UTIL)/, option.o) \
	$(addprefix $(SMITHLAB_CPP)/, smithlab_os.o smithlab_utils.o)

accuracy-for-single-end :  \
	$(addprefix $(ACGTUN_UTIL)/, option.o) \
	$(addprefix $(SMITHLAB_CPP)/, smithlab_os.o smithlab_utils.o)

accuracy-for-paired-end :  \
	$(addprefix $(ACGTUN_UTIL)/, option.o) \
	$(addprefix $(SMITHLAB_CPP)/, smithlab_os.o smithlab_utils.o)

accuracy-get-best-pairs :  \
	$(addprefix $(ACGTUN_UTIL)/, option.o) \
	$(addprefix $(SMITHLAB_CPP)/, smithlab_os.o smithlab_utils.o)

accuracy-get-best-single :  \
	$(addprefix $(ACGTUN_UTIL)/, option.o) \
	$(addprefix $(SMITHLAB_CPP)/, smithlab_os.o smithlab_utils.o)

mrsFast-count-positions-on-diff-mismatch2 : \
	$(addprefix $(ACGTUN_UTIL)/, option.o)

mrsFast-sum-count-positions-on-diff-mismatch3 : \
	$(addprefix $(ACGTUN_UTIL)/, option.o)

mrsFast-get-1000-results4 : \
	$(addprefix $(ACGTUN_UTIL)/, option.o)

mrsFast-compress1 : \
	$(addprefix $(ACGTUN_UTIL)/, option.o)


%: %.cpp
		$(CXX) $(CFLAGS) -o $@ $^ -I$(ACGTUN_UTIL) -I$(SMITHLAB_CPP)
		
		
install: all
	@mkdir -p $(SRC_ROOT)/bin
	@install -m 755 $(PROGS) $(SRC_ROOT)/bin

test_%:	%
	@$(TEST_DIR)/$@ $(TEST_DIR)

test:	$(addprefix test_, $(PROGS))

clean:
	@-rm -f $(PROGS) *.o *.so *.a *~

.PHONY: clean
