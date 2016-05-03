CC=g++
CFlags=-c -Wall -O3 -std=c++0x -MMD
#CoverageFlags=--coverage

Sources=$(wildcard src/*.cpp)
IncludeDir=./include
AllObjects=$(addprefix obj/,$(notdir $(Sources:.cpp=.o)))
Executables=Examples Figures Verify
Objects=$(filter-out $(addprefix obj/,$(Executables:=.o)),$(AllObjects))

all: $(Sources) $(Executables)

$(Executables): $(AllObjects)
	@mkdir -p data
	$(CC) $(CoverageFlags) $(Objects) $(addprefix obj/,$@.o) -o $@

obj/%.o: src/%.cpp
	@mkdir -p $(@D)
	$(CC) $(CFlags) $(CoverageFlags) -I$(IncludeDir) $< -o $@

-include $(AllObjects:.o=.d)

test: $(Executables)
	$(foreach exe,$(Executables),./$(exe);)

clean:
	rm -rf obj/*.o obj/*.d $(Executables)
