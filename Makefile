# To build gwasify:
#
#   make
#
# run with
#
#   ./build/gwasify

D_COMPILER=ldc2

LDMD=ldmd2

DUB_LIBS =

DFLAGS = -wi -I./source $(DUB_INCLUDE)
RPATH  =
LIBS   =
SRC    = $(wildcard source/gwasify/*.d  source/test/*.d)
IR     = $(wildcard source/gwasify/*.ll source/test/*.ll)
BC     = $(wildcard source/gwasify/*.bc source/test/*.bc)
OBJ    = $(SRC:.d=.o)
OUT    = build/gwasify

debug: DFLAGS += -O0 -g -d-debug $(RPATH) -link-debuglib $(BACKEND_FLAG) -unittest
release: DFLAGS += -O -release $(RPATH)

.PHONY:test

all: debug

build-setup:
	mkdir -p build/

ifeq ($(FORCE_DUPLICATE),1)
  DFLAGS += -d-version=FORCE_DUPLICATE
endif


default debug release profile getIR getBC gperf: $(OUT)

# ---- Compile step
%.o: %.d
	$(D_COMPILER) -lib $(DFLAGS) -c $< -od=$(dir $@) $(BACKEND_FLAG)

# ---- Link step
$(OUT): build-setup $(OBJ)
	$(D_COMPILER) -of=build/gwasify $(DFLAGS)  $(OBJ) $(LIBS   =) $(DUB_LIBS) $(BACKEND_FLAG)

test:
	chmod 755 build/gwasify
	./run_tests.sh

debug-strip: debug

clean:
	rm -rf build/*
	rm -f $(OBJ) $(OUT) trace.{def,log}
