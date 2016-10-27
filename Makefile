CXXFLAGS := --std=c++11 -Wall -Wextra -pedantic -O2 -g -DNDEBUG

BIN = gquery
DISTDIR = build
SRCDIR  = src
SRCFILES = $(wildcard $(SRCDIR)/*.cpp)
OBJFILES = $(patsubst $(SRCDIR)/%.cpp, $(DISTDIR)/%.o, $(SRCFILES))
RM = rm -rf


all: $(BIN)

$(BIN): $(OBJFILES)
	$(CXX) $(CXXFLAGS) $(CPPFLAGS) $^ $(LDLIBS) -o $(BIN)

$(DISTDIR)/%.o: $(SRCDIR)/%.cpp
	$(CXX) $(CXXFLAGS) $(CPPFLAGS) -c $< -o $@

.PHONY: clean distclean

clean:
	$(RM) $(OBJFILES) $(DISTDIR)

distclean: clean
	$(RM) $(BIN)
