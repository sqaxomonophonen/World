CXX = g++
C99 = gcc -std=c99
LINK = g++
AR = ar
#DEBUG_FLAG=-g
CXXFLAGS = -Wall -fPIC $(DEBUG_FLAG)
#CXXFLAGS += -O2
CXXFLAGS += -O0 -g

CFLAGS = $(CXXFLAGS)
ARFLAGS = -rv
OUT_DIR = ./build
OBJS = $(OUT_DIR)/objs/world/cheaptrick.o $(OUT_DIR)/objs/world/common.o $(OUT_DIR)/objs/world/d4c.o $(OUT_DIR)/objs/world/dio.o $(OUT_DIR)/objs/world/fft.o $(OUT_DIR)/objs/world/harvest.o $(OUT_DIR)/objs/world/matlabfunctions.o $(OUT_DIR)/objs/world/stonemask.o $(OUT_DIR)/objs/world/synthesis.o $(OUT_DIR)/objs/world/synthesisrealtime.o
LIBS =
MKDIR = mkdir -p $(1)
ifeq ($(shell echo "check_quotes"),"check_quotes")
	# Windows
	MKDIR = mkdir $(subst /,\,$(1)) > nul 2>&1 || (exit 0)
endif


CXXFLAGS += $(shell pkg-config --cflags sdl2)
SDL2LIB += $(shell pkg-config --libs sdl2)

all: default test step0analysis step1edit step2synthesis voxpopfix

###############################################################################################################
### Tests
###############################################################################################################
test: $(OUT_DIR)/test $(OUT_DIR)/ctest

step0analysis: $(OUT_DIR)/step0analysis
step1edit: $(OUT_DIR)/step1edit
step2synthesis: $(OUT_DIR)/step2synthesis
voxpopfix: $(OUT_DIR)/voxpopfix

test_OBJS=$(OUT_DIR)/objs/tools/audioio.o $(OUT_DIR)/objs/test/test.o
$(OUT_DIR)/test: $(OUT_DIR)/libworld.a $(test_OBJS)
	$(LINK) $(CXXFLAGS) -o $(OUT_DIR)/test $(test_OBJS) $(OUT_DIR)/libworld.a -lm

step0analysis_OBJS=$(OUT_DIR)/objs/tools/audioio.o $(OUT_DIR)/objs/test/step0analysis.o
$(OUT_DIR)/step0analysis: $(OUT_DIR)/libworld.a $(step0analysis_OBJS)
	$(LINK) $(CXXFLAGS) -o $(OUT_DIR)/step0analysis $(step0analysis_OBJS) $(OUT_DIR)/libworld.a -lm

step1edit_OBJS=$(OUT_DIR)/objs/test/step1edit.o
$(OUT_DIR)/step1edit: $(step1edit_OBJS)
	$(LINK) $(CXXFLAGS) -o $(OUT_DIR)/step1edit $(step1edit_OBJS) $(SDL2LIB) -lm

step2synthesis_OBJS=$(OUT_DIR)/objs/tools/audioio.o $(OUT_DIR)/objs/test/step2synthesis.o
$(OUT_DIR)/step2synthesis: $(OUT_DIR)/libworld.a $(step2synthesis_OBJS)
	$(LINK) $(CXXFLAGS) -o $(OUT_DIR)/step2synthesis $(step2synthesis_OBJS) $(OUT_DIR)/libworld.a -lm

voxpopfix_OBJS=$(OUT_DIR)/objs/tools/audioio.o $(OUT_DIR)/objs/tools/dr_wav.o $(OUT_DIR)/objs/test/voxpopfix.o
$(OUT_DIR)/voxpopfix: $(OUT_DIR)/libworld.a $(voxpopfix_OBJS)
	$(LINK) $(CXXFLAGS) -o $(OUT_DIR)/voxpopfix $(voxpopfix_OBJS) $(OUT_DIR)/libworld.a $(SDL2LIB) -lm



ctest_OBJS=$(OUT_DIR)/objs/tools/audioio.o $(OUT_DIR)/objs/test/ctest.o
$(OUT_DIR)/ctest: $(OUT_DIR)/libworld.a $(ctest_OBJS)
	$(LINK) $(CXXFLAGS) -o $(OUT_DIR)/ctest $(ctest_OBJS) $(OUT_DIR)/libworld.a -lm




$(OUT_DIR)/objs/test/test.o : tools/audioio.h src/world/d4c.h src/world/dio.h src/world/harvest.h src/world/matlabfunctions.h src/world/cheaptrick.h src/world/stonemask.h src/world/synthesis.h src/world/common.h src/world/fft.h src/world/macrodefinitions.h

$(OUT_DIR)/objs/test/ctest.o : tools/audioio.h src/world/d4c.h src/world/dio.h src/world/harvest.h src/world/matlabfunctions.h src/world/cheaptrick.h src/world/stonemask.h src/world/synthesis.h src/world/common.h src/world/fft.h src/world/macrodefinitions.h

$(OUT_DIR)/objs/test/step0analysis.o : tools/audioio.h src/world/d4c.h src/world/dio.h src/world/harvest.h src/world/matlabfunctions.h src/world/cheaptrick.h src/world/stonemask.h src/world/synthesis.h src/world/common.h src/world/fft.h src/world/macrodefinitions.h test/analysis.h

$(OUT_DIR)/objs/test/step1edit.o : test/analysis.h

$(OUT_DIR)/objs/test/step2synthesis.o : tools/audioio.h src/world/d4c.h src/world/dio.h src/world/harvest.h src/world/matlabfunctions.h src/world/cheaptrick.h src/world/stonemask.h src/world/synthesis.h src/world/common.h src/world/fft.h src/world/macrodefinitions.h test/analysis.h

$(OUT_DIR)/objs/test/voxpopfix.o : tools/audioio.h src/world/d4c.h src/world/dio.h src/world/harvest.h src/world/matlabfunctions.h src/world/cheaptrick.h src/world/stonemask.h src/world/synthesis.h src/world/common.h src/world/fft.h src/world/macrodefinitions.h

###############################################################################################################
### Library
###############################################################################################################
default: $(OUT_DIR)/libworld.a

$(OUT_DIR)/libworld.a: $(OBJS)
	$(AR) $(ARFLAGS) $(OUT_DIR)/libworld.a $(OBJS) $(LIBS)
	@echo Done.

$(OUT_DIR)/objs/world/cheaptrick.o : src/world/cheaptrick.h src/world/common.h src/world/constantnumbers.h src/world/matlabfunctions.h src/world/macrodefinitions.h
$(OUT_DIR)/objs/world/common.o : src/world/common.h src/world/constantnumbers.h src/world/matlabfunctions.h src/world/macrodefinitions.h
$(OUT_DIR)/objs/world/d4c.o : src/world/d4c.h src/world/common.h src/world/constantnumbers.h src/world/matlabfunctions.h src/world/macrodefinitions.h
$(OUT_DIR)/objs/world/dio.o : src/world/dio.h src/world/common.h src/world/constantnumbers.h src/world/matlabfunctions.h src/world/macrodefinitions.h
$(OUT_DIR)/objs/world/fft.o : src/world/fft.h src/world/macrodefinitions.h
$(OUT_DIR)/objs/world/harvest.o : src/world/harvest.h src/world/fft.h src/world/common.h src/world/constantnumbers.h src/world/matlabfunctions.h src/world/macrodefinitions.h
$(OUT_DIR)/objs/world/matlabfunctions.o : src/world/constantnumbers.h src/world/matlabfunctions.h src/world/macrodefinitions.h
$(OUT_DIR)/objs/world/stonemask.o : src/world/stonemask.h src/world/fft.h src/world/common.h src/world/constantnumbers.h src/world/matlabfunctions.h src/world/macrodefinitions.h
$(OUT_DIR)/objs/world/synthesis.o : src/world/synthesis.h src/world/common.h src/world/constantnumbers.h src/world/matlabfunctions.h src/world/macrodefinitions.h
$(OUT_DIR)/objs/world/synthesisrealtime.o : src/world/synthesisrealtime.h src/world/common.h src/world/constantnumbers.h src/world/matlabfunctions.h src/world/macrodefinitions.h


###############################################################################################################
### Global rules
###############################################################################################################
$(OUT_DIR)/objs/test/%.o : test/%.c
	$(call MKDIR,$(OUT_DIR)/objs/test)
	$(C99) $(CFLAGS) -Isrc -Itools -o "$@" -c "$<"

$(OUT_DIR)/objs/test/%.o : test/%.cpp
	$(call MKDIR,$(OUT_DIR)/objs/test)
	$(CXX) $(CXXFLAGS) -Isrc -Itools -o "$@" -c "$<"

$(OUT_DIR)/objs/tools/%.o : tools/%.cpp
	$(call MKDIR,$(OUT_DIR)/objs/tools)
	$(CXX) $(CXXFLAGS) -Isrc -o "$@" -c "$<"

$(OUT_DIR)/objs/tools/%.o : tools/%.c
	$(call MKDIR,$(OUT_DIR)/objs/tools)
	$(C99) $(CFLAGS) -o "$@" -c "$<"

$(OUT_DIR)/objs/world/%.o : src/%.cpp
	$(call MKDIR,$(OUT_DIR)/objs/world)
	$(CXX) $(CXXFLAGS) -Isrc -o "$@" -c "$<"

clean:
	@echo 'Removing all temporary binaries... '
	@$(RM) $(OUT_DIR)/libworld.a $(OBJS)
	@$(RM) $(test_OBJS) $(ctest_OBJS) $(OUT_DIR)/test $(OUT_DIR)/ctest
	@echo Done.

clear: clean

.PHONY: clean clear test default
.DELETE_ON_ERRORS:
