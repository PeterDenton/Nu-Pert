ifeq ($(IsTravis),)
	# leave CXX blank for travis, set it otherwise
	CXX=g++
else
	# turn on coverage for travis
	CoverageFlags=--coverage
endif
CFlags=-c -Wall -O3 -std=c++0x -pedantic -MMD

Sources=$(wildcard src/*.cpp)
IncludeDir=./include
AllObjects=$(addprefix obj/,$(notdir $(Sources:.cpp=.o)))
Executables=Examples Figures Verify
Objects=$(filter-out $(addprefix obj/,$(Executables:=.o)),$(AllObjects))

all: $(Sources) $(Executables)

$(Executables): $(AllObjects)
	@mkdir -p data
	$(CXX) $(CoverageFlags) $(Objects) $(addprefix obj/,$@.o) -o $@

obj/%.o: src/%.cpp
	@mkdir -p $(@D)
	$(CXX) $(CFlags) $(CoverageFlags) -I$(IncludeDir) $< -o $@

-include $(AllObjects:.o=.d)

test: $(Executables)
	$(foreach exe,$(Executables),./$(exe);)

clean:
	rm -rf obj/*.o obj/*.d $(Executables)
