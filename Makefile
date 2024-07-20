# Compiler
CXX = g++

# Compiler flags
CXXFLAGS = -g -IC:/SFML-2.6.1/include

# Linker flags
LDFLAGS = -LC:/SFML-2.6.1/lib -lsfml-graphics -lsfml-window -lsfml-system

# Source files for program1
SRCS_XOR_DATA = ANN/XOR_dataset.cpp Dual_Library/dual.cpp
# Object files for program1
OBJS_XOR_DATA = $(SRCS_XOR_DATA:.cpp=.o)

#source files for main -> testing
# SRCS_MAIN = main.cpp
#object files for main -> testing
# OBJS_MAIN = $(SRCS_MAIN:.cpp=.o)

# Source files for program2
SRCS_REG = ANN/Regression.cpp Dual_Library/dual.cpp
# Object files for program2
OBJS_REG = $(SRCS_REG:.cpp=.o)

# Executable names
EXE_XOR_DATA = XOR_dataset
EXE_REG = Regression
# EXE_MAIN = Main

# Default target
all: XOR_dataset Regression

# Build target for XOR_DATA
$(EXE_XOR_DATA): $(OBJS_XOR_DATA)
	$(CXX) $(CXXFLAGS) $^ -o $@ $(LDFLAGS)

# Build target for REG
$(EXE_REG): $(OBJS_REG)
	$(CXX) $(CXXFLAGS) $^ -o $@ $(LDFLAGS)

# # Build target for program2
# $(EXE_MAIN): $(OBJS_MAIN)
# 	$(CXX) $(CXXFLAGS) $^ -o $@ $(LDFLAGS)

# Pattern rule to build object files
%.o: %.cpp
	$(CXX) $(CXXFLAGS) -c $< -o $@

# Run target for program1 (Linux/macOS)
run_XOR: $(EXE_XOR_DATA)
	./$(EXE_XOR_DATA)

# Run target for program2 (Linux/macOS)
run_REG: $(EXE_REG)
	./$(EXE_REG)

# # Run MAIN
# run_MAIN: $(EXE_MAIN)
# 	./$(EXE_MAIN)

# Clean target
clean:
	rm -f $(EXE_XOR_DATA) $(EXE_REG) $(OBJS_XOR_DATA) $(OBJS_REG) $(EXE_MAIN) $(OBJS_MAIN)