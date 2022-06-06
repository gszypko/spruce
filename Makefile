# Based on http://make.mad-scientist.net/papers/advanced-auto-dependency-generation/
OBJDIR := obj
SRCS := source/solar/solarutils.cpp \
	source/mhd/plasmadomain.cpp \
	source/mhd/derivs.cpp \
	source/mhd/utils.cpp \
	source/mhd/mhd.cpp \
	source/mhd/fileio.cpp \
	source/mhd/grid.cpp \
	source/mhd/evolution.cpp \
	source/equationsets/equationset.cpp \
	source/equationsets/idealmhd.cpp \
	source/modules/solar/thermalconduction.cpp \
	source/modules/solar/localizedheating.cpp \
	source/modules/solar/fieldheating.cpp \
	source/modules/solar/ambientheating.cpp \
	source/modules/solar/radiativelosses.cpp \
	source/modules/modulehandler.cpp \
	source/modules/module.cpp \
	source/modules/sgfilter.cpp \
	main.cpp

OBJS := $(SRCS:%.cpp=$(OBJDIR)/%.o)
DEPDIR := $(OBJDIR)/.deps
DEPFLAGS = -MT $@ -MMD -MP -MF $(DEPDIR)/$*.d
EXEC = run

CXX = g++
CXXFLAGS = -fopenmp -std=c++17
CPPFLAGS = -I./source -I./source/mhd -I./source/equationsets \
	-I./source/modules -I./source/modules/solar -I./source/solar 

COMPILE.c = $(CXX) $(DEPFLAGS) $(CXXFLAGS) $(CPPFLAGS) $(TARGET_ARCH) -c

run: $(OBJS)
	$(CXX) $(CXXFLAGS) $(OBJS) -o $(EXEC)

clean:
	rm -f $(EXEC)
	rm -rf obj/

$(OBJDIR)/%.o : %.cpp $(DEPDIR)/%.d | $(DEPDIR)
	@mkdir -p $(dir $@)
	$(COMPILE.c) -o $@ $<

$(DEPDIR): ; @mkdir -p $@

DEPFILES := $(SRCS:%.cpp=$(DEPDIR)/%.d)
$(DEPFILES):
	@mkdir -p $(dir $@)

include $(wildcard $(DEPFILES))
