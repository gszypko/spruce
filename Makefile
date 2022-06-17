# Based on http:\\make.mad-scientist.net\papers\advanced-auto-dependency-generation
# input argument 
MAINFILE ?= main
FLAG ?= -O3
# detect operating system and update path fixing variable 
ifeq ($(OS),Windows_NT)
define path
$(subst /,\,$1)
endef
endif
# set up source, object, and dependent directories
SRCDIR := source
OBJDIR := objects
DEPDIR := $(OBJDIR)/.deps
SRCS := $(wildcard $(SRCDIR)/*/*.cpp) $(wildcard $(SRCDIR)/*/*/*.cpp) $(MAINFILE).cpp
HEADERS := source $(dir $(SRCS))
OBJS := $(SRCS:%.cpp=$(OBJDIR)/%.o)
EXEC = $(OBJDIR)/run.exe
# define flags for compilation
CXX = g++
CXXFLAGS = -fopenmp -std=c++20 $(FLAG)
CPPFLAGS := $(patsubst %,-I %,$(HEADERS))
LNKFLAGS = -lm -larmadillo -lstdc++fs
DEPFLAGS = -MT $@ -MMD -MP -MF $(DEPDIR)/$*.d

all : init run

init :
	@echo -------
	@echo Headers: $(CPPFLAGS)
	@echo Sources: $(SRCS)
	@echo File: $(MAINFILE).cpp

# set the default goal by making the target of the first rule $(EXEC)
run : $(OBJS)
	@$(CXX) $(CXXFLAGS) $(OBJS) -o $(EXEC) $(LNKFLAGS)
	@echo Complete
	@echo -------

$(DEPDIR) : ; @mkdir $@

# rule for making object files
$(OBJDIR)/%.o : %.cpp $(DEPDIR)/%.d | $(DEPDIR)
	@echo Compiling Object: $@ -- Dependencies: $<
	@if not exist $(dir $@) @mkdir $(call path,$(dir $@))
	@$(CXX) $(DEPFLAGS) $(CXXFLAGS) $(CPPFLAGS) -c -o $@ $<

DEPFILES := $(SRCS:%.cpp=$(DEPDIR)/%.d)
$(DEPFILES):
	@if not exist $(dir $@) @mkdir $(call path,$(dir $@))

include $(wildcard $(DEPFILES))

# option to clean directories at the end
clean :
	if exist $(OBJDIR) rmdir /s /q $(OBJDIR)
	