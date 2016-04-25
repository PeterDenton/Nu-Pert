CC=g++
CFlags=-c -Wall -O3 -std=c++0x -MMD
Sources=$(wildcard src/*.cpp)
Objects=$(addprefix obj/,$(notdir $(Sources:.cpp=.o)))
py_plots=$(py/*.py)
Executable=main

all: $(Sources) $(Executable)

$(Executable): $(Objects)
	$(CC) $(Objects) -o $@

obj/%.o: src/%.cpp
	@mkdir -p $(@D)
	$(CC) $(CFlags) $< -o $@

-include $(Objects:.o=.d)

clean:
	rm -rf obj/*.o obj/*.d $(Executable)
