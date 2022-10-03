#OTHER_FLAGS=$(ENV_HPC_OPTS) -std=c++17 -O3 -DCLS=$(getconf LEVEL1_DCACHE_LINESIZE) -funroll-loops -ffast-math -mtune=native -march=native -ftree-vectorize -frename-registers -msse3 -mavx -mavx2 -fomit-frame-pointer -ffinite-math-only -fno-signed-zeros -fno-trapping-math -mfma -DNDEBUG -DBOOST_UBLAS_NDEBUG

CXXFLAGS := -std=c++11 -Wall -Werror -ggdb -O0 #-O3 #-DCLS=$(getconf LEVEL1_DCACHE_LINESIZE) -funroll-loops -ffast-math -mtune=native -march=native -ftree-vectorize -frename-registers -msse3 -mavx -mavx2 -fomit-frame-pointer -ffinite-math-only -fno-signed-zeros -fno-trapping-math -mfma -DNDEBUG -DBOOST_UBLAS_NDEBUG
SOURCES := $(wildcard *.cpp)
OBJS := $(SOURCES:%.cpp=%.o)
TARGET := homework_2
.PHONY : all
all : $(TARGET)
$(TARGET) : $(OBJS)
	 $(CXX) $(CXXFLAGS) $^ -o $@
.PHONY : clean
clean :
	rm -rf $(TARGET) $(OBJS)