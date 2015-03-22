ifndef SRC_ROOT
SRC_ROOT=./../
endif

ifndef ACGTUN_UTIL
ACGTUN_UTIL=./../acgtun_util
endif
ifndef SMITHLAB_CPP
SMITHLAB_CPP=./../smithlab_cpp
endif

PROGS = compare-read-mapping-results \
	count-mismatch-reads \
	get-reverse-chromosomes \
	c2t4chroms \
	get-nucleotide-from-chromosome-by-positoin \
	get-time-memory-from-qsub-output-file \
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


compare-read-mapping-results :  \
	$(addprefix $(ACGTUN_UTIL)/, option.o) \
	$(addprefix $(SMITHLAB_CPP)/, smithlab_os.o smithlab_utils.o)

count-mismatch-reads : \
	$(addprefix $(ACGTUN_UTIL)/, option.o)

get-reverse-chromosomes : \
	$(addprefix $(ACGTUN_UTIL)/, option.o) \
	$(addprefix $(SMITHLAB_CPP)/, smithlab_os.o smithlab_utils.o)

c2t4chroms : \
	$(addprefix $(ACGTUN_UTIL)/, option.o) \
	$(addprefix $(SMITHLAB_CPP)/, smithlab_os.o smithlab_utils.o)

get-nucleotide-from-chromosome-by-positoin : \
	$(addprefix $(ACGTUN_UTIL)/, option.o) \
	$(addprefix $(SMITHLAB_CPP)/, smithlab_os.o smithlab_utils.o)


get-time-memory-from-qsub-output-file : \
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
