CXXFLAGS := -std=c++11 -Wall -Werror -ggdb -O0
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