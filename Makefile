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

# Source files for program2
SRCS_REG = ANN/Regression.cpp Dual_Library/dual.cpp
# Object files for program2
OBJS_REG = $(SRCS_REG:.cpp=.o)

# Source files for program2
SRCS_GRAD = Others/gradient_tester.cpp Dual_Library/extendedDual.cpp
# Object files for program2
OBJS_GRAD = $(SRCS_GRAD:.cpp=.o)

# Source files for program2
SRCS_DER = Others/ShowcasingDerivative.cpp Dual_Library/extendedDual.cpp
# Object files for program2
OBJS_DER = $(SRCS_DER:.cpp=.o)

# Executable names
EXE_XOR_DATA = XOR_dataset
EXE_REG = Regression
EXE_GRAD = Gradient
EXE_DER = DerivativeShowcase
# EXE_MAIN = Main

# Default target
all: XOR_dataset Regression Gradient Derivative

# Build target for XOR_DATA
$(EXE_XOR_DATA): $(OBJS_XOR_DATA)
	$(CXX) $(CXXFLAGS) $^ -o $@ $(LDFLAGS)

# Build target for REG
$(EXE_REG): $(OBJS_REG)
	$(CXX) $(CXXFLAGS) $^ -o $@ $(LDFLAGS)

# Build target for GRAD
$(EXE_GRAD): $(OBJS_GRAD)
	$(CXX) $(CXXFLAGS) $^ -o $@

$(EXE_DER): $(OBJS_DER)
	$(CXX) $(CXXFLAGS) $^ -o $@

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

run_GRAD: $(EXE_GRAD)
	./$(EXE_GRAD)

run_DER: $(EXE_DER)
	./$(EXE_DER)

# # Run MAIN
# run_MAIN: $(EXE_MAIN)
# 	./$(EXE_MAIN)

# Clean target
clean:
	rm -f $(EXE_XOR_DATA) $(EXE_REG) $(EXE_GRAD) $(EXE_DER) $(OBJS_XOR_DATA) $(OBJS_REG) $(OBJS_GRAD) $(OBJS_DER)