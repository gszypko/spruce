# Based on http://make.mad-scientist.net/papers/advanced-auto-dependency-generation/
OBJDIR := obj
SRCS := source/solar/solarutils.cpp \
	source/user-interface/PlasmaSettings.cpp \
	source/user-interface/MhdInp.cpp \
	source/user-interface/ui_utility.cpp \
	source/mhd/plasmadomain.cpp \
	source/mhd/derivs.cpp \
	source/mhd/utils.cpp \
	source/mhd/mhd.cpp \
	source/mhd/fileio.cpp \
	source/mhd/grid.cpp \
	source/mhd/evolution.cpp \
	source/ucnp/ucnputils.cpp \
	source/ucnp/antihelmholtz.cpp \
	source/ucnp/currentloop.cpp \
	source/ucnp/settings.cpp \
	source/modules/solar/thermalconduction.cpp \
	source/modules/solar/localizedheating.cpp \
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
CPPFLAGS = -I./source -I./source/user-interface -I./source/mhd \
	-I./source/modules -I./source/modules/solar -I./source/solar -I./source/ucnp

COMPILE.c = $(CXX) $(DEPFLAGS) $(CXXFLAGS) $(CPPFLAGS) $(TARGET_ARCH) -c

run: $(OBJS)
	$(CXX) $(CXXFLAGS) $(OBJS) -o $(EXEC)

clean:
	rm $(EXEC)
	rm -r obj/

$(OBJDIR)/%.o : %.cpp $(DEPDIR)/%.d | $(DEPDIR)
	@mkdir -p $(dir $@)
	$(COMPILE.c) -o $@ $<

$(DEPDIR): ; @mkdir -p $@

DEPFILES := $(SRCS:%.cpp=$(DEPDIR)/%.d)
$(DEPFILES):
	@mkdir -p $(dir $@)

include $(wildcard $(DEPFILES))
