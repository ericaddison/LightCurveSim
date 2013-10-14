#*************************************
#	File: Makefile
#	Project: Light Curve Simulation (LC)
#	Purpose: makefile for LC simulation
#	Author: Eric Addison
#	Original date: 23Apr13
#*************************************

# This is the makefile for the LC simulation. Relevant calls are
#	make all -- compiles all
#	make run -- compiles and runs
#	make clean -- cleans old object files


# includes options for GSL library functions

# variables
FILES=LC.c LC_funcs.c Roche.c hullLC.c cosmo.c misc_math.c demo.c
CFLAGS=-Wall -std=c99 -g -c -I/usr/local/include -I/usr/local/include/libqhull 
LFLAGS=-L/usr/local/lib -lgsl -lgslcblas -lqhullstatic -lm	# linker flags
SRCDIR=./code/
OBJDIR=./obj/
SOURCES=$(FILES:%.c=$(SRCDIR)%.o)
OBJECTS=$(FILES:%.c=$(OBJDIR)%.o)

# make all
all: LC

# make run, compiles and runs program
run: LC
	@echo "Running LightCurve"
	@./LC

# make prof, include profiling flags
#prof: LFLAGS += -pg 
#prof: CFLAGS += -pg 
#prof: LC

# LC: calls all compile targets
LC: welcome make_dir $(OBJECTS)
	@echo "Compiling LightCurve"
	@gcc $(OBJECTS) $(LFLAGS) -o LC
	
#Object File targets
#note: $@ is replaced by the name of the target
# and $< is the name of the first dependency

$(OBJDIR)LC.o: $(SRCDIR)LC.c
	@echo "Compiling $?"
	@gcc $(CFLAGS) $? -o $@

$(OBJDIR)LC_funcs.o: $(SRCDIR)LC_funcs.c
	@echo "Compiling $?"
	@gcc $(CFLAGS) $? -o $@

$(OBJDIR)Roche.o: $(SRCDIR)Roche.c
	@echo "Compiling $?"
	@gcc $(CFLAGS) $? -o $@

$(OBJDIR)hullLC.o: $(SRCDIR)hullLC.c
	@echo "Compiling $?"
	@gcc $(CFLAGS) $? -o $@

$(OBJDIR)cosmo.o: $(SRCDIR)cosmo.c
	@echo "Compiling $?"
	@gcc $(CFLAGS) $? -o $@

$(OBJDIR)misc_math.o: $(SRCDIR)misc_math.c
	@echo "Compiling $?"
	@gcc $(CFLAGS) $? -o $@

$(OBJDIR)demo.o: $(SRCDIR)demo.c
	@echo "Compiling $?"
	@gcc $(CFLAGS) $? -o $@

make_dir:
	@echo "Checking for object directory"
	@if test -d ./obj; then echo "...Directory Found"; else mkdir ./obj; fi

welcome:
	clear
	@echo "           LightCurve Makefile              "
	@echo "--------------------------------------------"
	@echo

# clean up the object directory if stuff is conflicting
clean:
	rm -rf $(OBJECTS) LC



##### Makefile target template #######
# target: dependency1 dependency2 ƒ
#	rule command1
#	rule command2
#
