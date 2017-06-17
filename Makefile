#  g fortran compiler
   FC = gfortran
#  Compilateur NAG sur mamaku
#  FC = pgf90
# FFLAGS = -fast -r8
# FFLAGS = -fast -r8 -Mprof=lines
# FFLAGS = -r8 -Mprof=lines
# FFLAGS = -g 
# FFLAGS = -g -r8 -C -Miomutex -Minform,inform
# LAPACK = /usr/local/pgi32/linux86/lib-glibc-212
# FFLAGS = -O -s -r8 -check bounds
  FFLAGS = -O
# FFLAGS = -g -s -B111 -m1 -ea
# FFLAGS = -fast
# FFLAGS = -error_limit 2
# LAPACK = /usr/local/absoft/extras/lapack/precompiled

#.MAKEOPTS: -k

 SRC = \
       numerical_recipes.f90 \
       tools.f90 \
       constants.f90 \
       var_VODE.f90 \
       variables_mhd.f90 \
       read_fe_data.f90 \
       chemical_species.f90 \
       chemical_reactions.f90 \
       molecular_cooling.f90 \
       H2.f90 \
       CO.f \
       SiO.f \
       oH2O.f \
       pH2O.f \
       oNH3.f \
       pNH3.f \
       OH.f \
       Atype.f \
       Etype.f \
       line_excit.f90 \
       evolution.f90 \
       integrateur_vode.f \
       energetics.f90 \
       initialize.f90 \
       outputs.f90 \
       mhd_vode.f90 \
       add.f

 OBJ =\
       numerical_recipes.o \
       tools.o \
       constants.o \
       var_VODE.o \
       variables_mhd.o \
       read_fe_data.o \
       chemical_species.o \
       chemical_reactions.o \
       molecular_cooling.o \
       H2.o \
       CO.o \
       SiO.o \
       oH2O.o \
       pH2O.o \
       oNH3.o \
       pNH3.o \
       OH.o \
       Atype.o \
       Etype.o \
       line_excit.o \
       evolution.o \
       integrateur_vode.o \
       energetics.o \
       initialize.o \
       outputs.o \
       mhd_vode.o \
       add.o

.f90.o:;	$(FC) $(FFLAGS) -c $<

mhd_vode: $(OBJ)
	$(FC) $(FFLAGS) $(OBJ) -o mhd_vode -lblas -llapack
	
numerical_recipes.o: numerical_recipes.f90
	$(FC) $(FFLAGS) -c numerical_recipes.f90

tools.o: tools.f90
	$(FC) $(FFLAGS) -c tools.f90

constants.o: constants.f90
	$(FC) $(FFLAGS) -c constants.f90

variables_mhd.o: variables_mhd.f90
	$(FC) $(FFLAGS) -c variables_mhd.f90

read_fe_data.o: read_fe_data.f90
	$(FC) $(FFLAGS) -c read_fe_data.f90

chemical_species.o: chemical_species.f90
	$(FC) $(FFLAGS) -c chemical_species.f90

chemical_reactions.o: chemical_reactions.f90
	$(FC) $(FFLAGS) -c chemical_reactions.f90

molecular_cooling.o: molecular_cooling.f90
	$(FC) $(FFLAGS) -c molecular_cooling.f90

H2.o: H2.f90
	$(FC) $(FFLAGS) -c H2.f90

CO.o: CO.f  
	$(FC) $(FFLAGS) -c CO.f  

SiO.o: SiO.f  
	$(FC) $(FFLAGS) -c SiO.f  

oH2O.o: oH2O.f  
	$(FC) $(FFLAGS) -c oH2O.f  

pH2O.o: pH2O.f  
	$(FC) $(FFLAGS) -c pH2O.f  

oNH3.o: oNH3.f  
	$(FC) $(FFLAGS) -c oNH3.f  

pNH3.o: pNH3.f  
	$(FC) $(FFLAGS) -c pNH3.f  
	
OH.o: OH.f  
	$(FC) $(FFLAGS) -c OH.f  	

Atype.o: Atype.f  
	$(FC) $(FFLAGS) -c Atype.f  

Etype.o: Etype.f  
	$(FC) $(FFLAGS) -c Etype.f  

line_excit.o: line_excit.f90
	$(FC) $(FFLAGS) -c line_excit.f90

evolution.o: evolution.f90
	$(FC) $(FFLAGS) -c evolution.f90

var_VODE.o: var_VODE.f90
	$(FC) $(FFLAGS) -c var_VODE.f90

integrateur_vode.o: integrateur_vode.f  
	$(FC) $(FFLAGS) -c integrateur_vode.f  

energetics.o: energetics.f90
	$(FC) $(FFLAGS) -c energetics.f90

initialize.o: initialize.f90
	$(FC) $(FFLAGS) -c initialize.f90

outputs.o: outputs.f90
	$(FC) $(FFLAGS) -c outputs.f90

mhd_vode.o: mhd_vode.f90
	$(FC) $(FFLAGS) -c mhd_vode.f90

clean:
	rm -f *.o *.mod mhd_vode

print:
	a2ps $(SRC) -f6 -o code.ps

excit_C: excit_C.f90
	$(FC) excit_C.f90 $(FFLAGS) -o excit_C -llapack -lblas
#  Compilateur NAG sur mamaku
#	$(FC) excit_C.f90 $(FFLAGS) -o excit_C -ldxml

