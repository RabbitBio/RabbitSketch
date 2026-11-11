INCLUDE=-I../src -I../build/include
# CMake places librabbitsketch_static.a in ../build; shared lib is under ../build/lib
RABBIT_LIBDIRS=-L../build/lib -L../build
CXX=g++
#
# The static library must be rebuilt without LTO after CMake change (see ../CMakeLists.txt).
# If ../build/lib/librabbitsketch_static.a was built with -flto, conda/other g++ will fail
# at link with "bytecode stream ... lto-wrapper". Fix: re-run cmake --build in ../build.
#
# Optional: make USE_LTO=1 only if you rebuilt librabbitsketch_static.a with the same $(CXX)+-flto.
USE_LTO ?= 0
CXXFLAGS=-g1 -O3 -std=c++2a -fopenmp -march=native -DNDEBUG -funroll-loops -fno-lto $(INCLUDE)
ifeq ($(USE_LTO),1)
CXXFLAGS := $(filter-out -fno-lto,$(CXXFLAGS))
CXXFLAGS += -flto
endif

all: rabbitsketch exe_test_HLL exe_test_SetSketch exe_test_Kssd exe_test_MinHash exe_test_OMH exe_test_ProbMinHash exe_test_OneProbMinHash exe_test_FastKMV exe_test_BinDash exe_test_TopLProbMinHash exe_eval_sketch_acc

exe_test_HLL: test_HLL.cpp ../src/HyperLogLog.cpp ../src/InvertedIndex.h
	$(CXX) $(CXXFLAGS) test_HLL.cpp ../src/HyperLogLog.cpp -o exe_test_HLL $(RABBIT_LIBDIRS) $(EXTRA_LDFLAGS) -lz -lrabbitsketch_static

exe_test_SetSketch: test_SetSketch.cpp ../src/SetSketch.cpp ../src/InvertedIndex.h
	$(CXX) $(CXXFLAGS) test_SetSketch.cpp ../src/SetSketch.cpp -o exe_test_SetSketch $(RABBIT_LIBDIRS) $(EXTRA_LDFLAGS) -lz -lrabbitsketch_static

KSSD_HEADERS=../src/Sketch.h ../src/KssdBuild.h ../src/KssdIndex.h ../src/io/KssdFastx.h

exe_test_Kssd: test_Kssd.cpp sub_command.cpp sub_command.h ../src/BinDashIndex.h ../src/InvertedIndex.h $(KSSD_HEADERS) ../build/lib/librabbitsketch_static.a
	$(CXX) $(CXXFLAGS) test_Kssd.cpp sub_command.cpp -o exe_test_Kssd $(RABBIT_LIBDIRS) $(EXTRA_LDFLAGS) -ldeflate -lz -lrabbitsketch_static

exe_test_MinHash: test_MinHash.cpp ../src/InvertedIndex.h
	$(CXX) $(CXXFLAGS) test_MinHash.cpp  -o exe_test_MinHash $(RABBIT_LIBDIRS) $(EXTRA_LDFLAGS) -lz -lrabbitsketch_static

exe_test_OMH: test_OMH.cpp
	$(CXX) $(CXXFLAGS) test_OMH.cpp -o exe_test_OMH $(RABBIT_LIBDIRS) $(EXTRA_LDFLAGS) -lz -lrabbitsketch_static

RABBIT_SRCS=../src/probmh.cpp ../src/fastkmv.cpp ../src/fracminhash.cpp ../src/BinDash.cpp ../src/SetSketch.cpp ../src/shuffle.cpp ../src/HyperLogLog.cpp
RABBIT_CLI_SRCS=main.cpp sub_command.cpp

# The pre-built static library was compiled against a newer libstdc++
# (uses std::__throw_bad_array_new_length, a C++17+ symbol not present
# in older system libstdc++).  If your system g++/libstdc++ is too old, point
# CONDA_LIB at a newer libstdc++ at build time, e.g.
#     make rabbitsketch CONDA_LIB=$CONDA_PREFIX/lib
# By default we use $CONDA_PREFIX/lib if available, otherwise no extra runpath
# (compiler/linker defaults are used).
CONDA_LIB ?= $(if $(CONDA_PREFIX),$(CONDA_PREFIX)/lib,)
ifneq ($(strip $(CONDA_LIB)),)
EXTRA_LDFLAGS=-L$(CONDA_LIB) -Wl,-rpath,$(CONDA_LIB)
else
EXTRA_LDFLAGS=
endif

rabbitsketch: $(RABBIT_CLI_SRCS) sub_command.h $(RABBIT_SRCS) ../src/BinDash.h ../src/BinDashIndex.h ../src/InvertedIndex.h $(KSSD_HEADERS) ../build/lib/librabbitsketch_static.a
	$(CXX) $(CXXFLAGS) $(RABBIT_CLI_SRCS) $(RABBIT_SRCS) -o rabbitsketch $(RABBIT_LIBDIRS) $(EXTRA_LDFLAGS) -ldeflate -lz -lrabbitsketch_static

exe_test_ProbMinHash: test_ProbMinHash.cpp ../src/probmh.cpp ../src/InvertedIndex.h
	$(CXX) $(CXXFLAGS) test_ProbMinHash.cpp ../src/probmh.cpp -o exe_test_ProbMinHash $(RABBIT_LIBDIRS) $(EXTRA_LDFLAGS) -lz -lrabbitsketch_static

exe_test_OneProbMinHash: test_OneProbMinHash.cpp ../src/probmh_op.cpp
	$(CXX) $(CXXFLAGS) test_OneProbMinHash.cpp ../src/probmh_op.cpp -o exe_test_OneProbMinHash $(RABBIT_LIBDIRS) $(EXTRA_LDFLAGS) -lz -lrabbitsketch_static

exe_test_FastKMV: test_FastKMV.cpp ../src/fastkmv.cpp ../src/InvertedIndex.h
	$(CXX) $(CXXFLAGS) test_FastKMV.cpp ../src/fastkmv.cpp -o exe_test_FastKMV $(RABBIT_LIBDIRS) $(EXTRA_LDFLAGS) -lz -lrabbitsketch_static

exe_test_TopLProbMinHash: test_TopLProbMinHash.cpp ../src/probmh.cpp
	$(CXX) $(CXXFLAGS) test_TopLProbMinHash.cpp ../src/probmh.cpp -o exe_test_TopLProbMinHash $(RABBIT_LIBDIRS) $(EXTRA_LDFLAGS) -lz -lrabbitsketch_static

exe_eval_reads: eval_reads.cpp ../src/probmh.cpp ../src/fastkmv.cpp ../src/BinDash.cpp ../src/SetSketch.cpp $(KSSD_HEADERS) ../build/lib/librabbitsketch_static.a
	$(CXX) $(CXXFLAGS) eval_reads.cpp ../src/probmh.cpp ../src/fastkmv.cpp ../src/BinDash.cpp ../src/SetSketch.cpp -o exe_eval_reads $(RABBIT_LIBDIRS) $(EXTRA_LDFLAGS) -lz -lrabbitsketch_static

exe_eval_sketch_acc: eval_sketch_accuracy.cpp ../src/probmh.cpp ../src/probmh_op.cpp ../src/fastkmv.cpp ../src/BinDash.cpp $(KSSD_HEADERS) ../build/lib/librabbitsketch_static.a
	$(CXX) $(CXXFLAGS) eval_sketch_accuracy.cpp ../src/probmh.cpp ../src/probmh_op.cpp ../src/fastkmv.cpp ../src/BinDash.cpp -o exe_eval_sketch_acc $(RABBIT_LIBDIRS) $(EXTRA_LDFLAGS) -lz -lrabbitsketch_static

exe_test_ani_containment: test_ani_containment.cpp ../src/probmh.cpp ../src/fastkmv.cpp
	$(CXX) $(CXXFLAGS) test_ani_containment.cpp ../src/probmh.cpp ../src/fastkmv.cpp -o exe_test_ani_containment $(RABBIT_LIBDIRS) $(EXTRA_LDFLAGS) -lz -lrabbitsketch_static

exe_test_BinDash: test_BinDash.cpp ../src/BinDash.cpp ../src/InvertedIndex.h
	$(CXX) $(CXXFLAGS) test_BinDash.cpp ../src/BinDash.cpp -o exe_test_BinDash $(RABBIT_LIBDIRS) $(EXTRA_LDFLAGS) -lz -lrabbitsketch_static

clean:
	rm -f *.o exe_* rabbitsketch

# KSSD evaluation entry points share the library cache and streaming builder.
exe_eval_param_sweep: eval_param_sweep.cpp $(KSSD_HEADERS) ../build/lib/librabbitsketch_static.a
	$(CXX) $(CXXFLAGS) $< -o $@ $(RABBIT_LIBDIRS) $(EXTRA_LDFLAGS) -lz -lrabbitsketch_static

exe_eval_alisim_sweep: eval_alisim_sweep.cpp $(KSSD_HEADERS) ../build/lib/librabbitsketch_static.a
	$(CXX) $(CXXFLAGS) $< -o $@ $(RABBIT_LIBDIRS) $(EXTRA_LDFLAGS) -lz -lrabbitsketch_static
