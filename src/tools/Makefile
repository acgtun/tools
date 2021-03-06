ifndef SRC_ROOT
SRC_ROOT=./../
endif

ifndef ACGTUN_UTIL
ACGTUN_UTIL=./../acgtun_util
endif
ifndef SMITHLAB_CPP
SMITHLAB_CPP=./../smithlab_cpp
endif

PROGS = count-mismatch-for-mapped-reads \
	count-mismatch-reads-for-RMAPBS-format \
	bireads-test.pbs.gen \
	to-mr-bsmap \
	to-cluster-walt \
	trim-reads \
	get-reverse-chromosomes \
	c2t4chroms \
	c2t-reads \
	get-nucleotide-from-chromosome-by-positoin \
	meth_CpG_symmetric \
	meth_comp_two_methylomes \
	get-time-memory-from-qsub-output-file \
	get-time-memory-from-qsub-output-file-pre \
	get-time-memory-from-qsub-error-file \
	generate-test-data 

CXX = g++
CFLAGS = -Wall -O3 -fPIC -fmessage-length=50 -std=c++0x

TEST_DIR = $(SRC_ROOT)/test

ifeq "$(shell uname)" "Darwin"
CFLAGS += -arch x86_64
endif

all:	$(PROGS)

%.o: %.cpp %.hpp
	$(CXX) $(CFLAGS) -c -o $@ $<  -I$(ACGTUN_UTIL) -I$(SMITHLAB_CPP)


count-mismatch-for-mapped-reads : \
	$(addprefix $(ACGTUN_UTIL)/, option.o)

count-mismatch-reads-for-RMAPBS-format : \
	$(addprefix $(ACGTUN_UTIL)/, option.o)

meth_CpG_symmetric : \
	$(addprefix $(ACGTUN_UTIL)/, option.o)

meth_comp_two_methylomes : \
	$(addprefix $(ACGTUN_UTIL)/, option.o)

bireads-test.pbs.gen : \
	$(addprefix $(ACGTUN_UTIL)/, option.o)

to-mr-bsmap : \
	$(addprefix $(ACGTUN_UTIL)/, option.o)

to-cluster-walt : \
	$(addprefix $(SMITHLAB_CPP)/, smithlab_os.o smithlab_utils.o)

trim-reads : \
	$(addprefix $(ACGTUN_UTIL)/, option.o) \
	$(addprefix $(SMITHLAB_CPP)/, smithlab_os.o smithlab_utils.o)

get-reverse-chromosomes : \
	$(addprefix $(ACGTUN_UTIL)/, option.o) \
	$(addprefix $(SMITHLAB_CPP)/, smithlab_os.o smithlab_utils.o)

c2t4chroms : \
	$(addprefix $(ACGTUN_UTIL)/, option.o) \
	$(addprefix $(SMITHLAB_CPP)/, smithlab_os.o smithlab_utils.o)

c2t-reads : \
	$(addprefix $(ACGTUN_UTIL)/, option.o) \
	$(addprefix $(SMITHLAB_CPP)/, smithlab_os.o smithlab_utils.o)

get-nucleotide-from-chromosome-by-positoin : \
	$(addprefix $(ACGTUN_UTIL)/, option.o) \
	$(addprefix $(SMITHLAB_CPP)/, smithlab_os.o smithlab_utils.o)


get-time-memory-from-qsub-output-file : \
	$(addprefix $(ACGTUN_UTIL)/, option.o) \
	$(addprefix $(SMITHLAB_CPP)/, smithlab_os.o smithlab_utils.o)

get-time-memory-from-qsub-output-file-pre : \
	$(addprefix $(ACGTUN_UTIL)/, option.o) \
	$(addprefix $(SMITHLAB_CPP)/, smithlab_os.o smithlab_utils.o)

get-time-memory-from-qsub-error-file : \
	$(addprefix $(ACGTUN_UTIL)/, option.o) \
	$(addprefix $(SMITHLAB_CPP)/, smithlab_os.o smithlab_utils.o)

generate-test-data : \
	$(addprefix $(ACGTUN_UTIL)/, option.o) \
	$(addprefix $(SMITHLAB_CPP)/, smithlab_os.o smithlab_utils.o)

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
