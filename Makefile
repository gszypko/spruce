# Based on http:\\make.mad-scientist.net\papers\advanced-auto-dependency-generation

# Platform Independency
#	The $(OS) statement below should be updated to include Mac and Linux
#	The <path> command should be updated for each OS
#	An overload for the OS-dependent command line functions such as rmdir, mkdir
#	should be created to ensure syntax is correct

# Notes:
#	This file generates the dependencies and compiles each file within the CPP project.
# 	Files are only recompiled if they or their dependencies have changed.
#	The command <$(wildcard $(SRCDIR)/*.cpp)> must go as deep into source as there are files
#		(i.e., $(wildcard $(SRCDIR)/*/*.cpp)) if necessary and so on
#	Any additional source files not contained within source/ must be manually appended, as main.cpp is
#	Any additional header folders not associated with .cpp files must also be manually appended

# input arguments with default values
file ?= main
flag ?= -O3
version ?= -std=c++20
# define project directories
SRCDIR := source
OBJDIR := objects
DEPDIR := $(OBJDIR)/.deps
# obtain source files, include folders, and object and dependency names from source folder
SRCS := $(wildcard $(SRCDIR)/*/*.cpp) $(wildcard $(SRCDIR)/*/*/*.cpp) $(file).cpp
HEADERS := $(sort $(dir $(SRCS))) source
OBJS := $(SRCS:%.cpp=$(OBJDIR)/%.o)
DEPFILES = $(SRCS:%.cpp=$(DEPDIR)/%.d)
EXEC := $(OBJDIR)/run$(extension)
# define flags for compilation
CXX = g++
DEPFLAGS = -MT $@ -MMD -MP -MF $(DEPDIR)/$*.d
CXXFLAGS = -fopenmp $(version) $(flag)
CPPFLAGS = $(patsubst %,-I %,$(HEADERS))
LNKFLAGS = -lm -lstdc++fs
# detect operating system and define shell commands accordingly
# each of path, makedir, removedir, removefile, and extension must defined for your operatings sytem 

ifeq ($(OS),Windows_NT)
define path
$(subst /,\,$1)
endef
define makedir
if not exist $1 mkdir $(call path,$1)
endef
define removedir
if exist $1 rmdir /s /q $(call path,$1)
endef
define removefile
if exist $1 del /f $(call path,$1)
endef
extension=.exe
else # Linux
define path
$(subst /,/,$1)
endef
define makedir
mkdir -p $(call path,$1)
endef
define removedir
rm -r $(call path,$1)
endef
define removefile
rm -f $(call path,$1)
endef
extension=
endif
# define targets to be completed in order
all : init run
# print compilation data to command line
init :
	@$(call removefile,$(EXEC))
	@echo -------
	@echo File: $(file).cpp
	@echo Flags: $(CXXFLAGS)
	@echo Linker: $(LNKFLAGS)
	@echo Headers: $(CPPFLAGS)
	@echo Sources: $(SRCS)
# define compilation goal
run : $(OBJS)
	@$(CXX) $(CXXFLAGS) $(OBJS) -o $(EXEC) $(LNKFLAGS)
	@echo Complete
# rule for making object files
$(OBJDIR)/%.o : %.cpp $(DEPDIR)/%.d | $(DEPDIR)
	@echo Compiling Object: $@ -- Dependencies: $<
	@$(call makedir,$(dir $@))
	@$(CXX) $(DEPFLAGS) $(CXXFLAGS) $(CPPFLAGS) -c -o $@ $<
# rule for making dependency directories
$(DEPDIR) : ; @$(call makedir,$(dir $@))
# rule for making dependency files
$(DEPFILES):
	@$(call makedir,$(dir $@))
include $(wildcard $(DEPFILES))
# option to clean directories at the end
clean :
	@$(call removedir, $(OBJDIR))