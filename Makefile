# Makefile for VISUALIZING_WORMHOLE

CXX       := g++
CXXFLAGS  := -std=c++17 -O2 -Wall -fopenmp
INCDIR    := src

SRCS      := \
    src/main.cpp \
    src/Camera.cpp \
    src/Geodesic.cpp \
    src/Matrix4x4.cpp \
    src/Vector4.cpp

OBJS      := $(SRCS:.cpp=.o)
TARGET    := wormhole_sim

.PHONY: all plot clean

all: $(TARGET)

$(TARGET): $(OBJS)
	$(CXX) $(CXXFLAGS) -I$(INCDIR) -o $@ $^

src/%.o: src/%.cpp
	$(CXX) $(CXXFLAGS) -I$(INCDIR) -c $< -o $@

plot: $(TARGET)
	@echo "Running Python visualization..."
	python python/plot_results.py

clean:
ifeq ($(OS),Windows_NT)
	del /Q /F wormhole_sim.exe src\*.o
else
	rm -f $(OBJS) $(TARGET)
endif
