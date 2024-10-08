MOD=../mod/
LOBJ=../obj/
BIN=/home/mjh217/bin/
MV=mv -f
F95 = gfortran


FFLAGS =   -J../mod  -O5
F77FLAGS =  -O3 -e -f77

SHFLAGS= -I/space1/spherical_harmonics/SHTOOLS-4.4/modules
LIB= lapack -L/space1/spherical_harmonics/SHTOOLS-4.4/lib -lSHTOOLS \
        -lfftw3 -lm -s -m64




SRCS =  nrtype.f90            \
        nrutil.f90            \
        module_util.f90       \
	module_function.f90   \
        rotmx2.f              \
        module_maps.f90       \
        gll_library.f90       \
        lagrange_poly.f90     \
        module_sphexp.f90     \
        module_crust.f90      \
        module_ice5g.f90      \
        ran_state.f90         \
        ran1.f90              \
        gasdev.f90



OBJS =  $(LOBJ)nrtype.o            \
        $(LOBJ)nrutil.o            \
        $(LOBJ)module_util.o       \
        $(LOBJ)module_function.o   \
        $(LOBJ)rotmx2.o            \
        $(LOBJ)module_maps.o       \
        $(LOBJ)gll_library.o       \
        $(LOBJ)lagrange_poly.o     \
	$(LOBJ)module_sphexp.o     \
        $(LOBJ)module_crust.o      \
        $(LOBJ)module_ice5g.o      \
        $(LOBJ)ran_state.o         \
        $(LOBJ)ran1.o              \
        $(LOBJ)gasdev.o            




MODS =  $(MOD)nrtype.mod            \
        $(MOD)nrutil.mod            \
        $(MOD)module_util.mod       \
	$(MOD)module_spline.mod     \
        $(MOD)module_function.mod   \
        $(MOD)module_model.mod      \
        $(MOD)module_maps.mod       \
        $(MOD)module_sphexp.mod     \
        $(MOD)module_crust.mod      \
        $(MOD)module_ice5g.mod      \
	$(MOD)ran_state.mod






all:    sphexp_data            \
        sphexp_inv             \
	sphexp_data_id         \
        sphexp_backus_gilbert  \
        sphexp_chi_comp        \
        sphexp_inv_av          \
	sphexp_data_idc        \
        sphexp_data_idlm       \
        sphexp_global_grid     \
        sphexp_global_grid_2   \
        sphexp_inv_short       \
	sphexp_inv_intermediate	\
        sphexp_data_power_spectrum_recovery



sphexp_data: $(BIN)sphexp_data $(OBJS) $(OBJSF77)
$(BIN)sphexp_data: $(OBJS) $(OBJSF77) $(LOBJ)sphexp_data.o
	$(F95) $(FFLAGS) $(LOBJ)sphexp_data.o $(OBJS)  \
        $(OBJSF77) -l $(LIB)  -o $(BIN)sphexp_data

$(LOBJ)sphexp_data.o: sphexp_data.f90 
	$(F95) $(FFLAGS) sphexp_data.f90 -c
	$(MV) sphexp_data.o $(LOBJ)


sphexp_inv: $(BIN)sphexp_inv $(OBJS) $(OBJSF77)
$(BIN)sphexp_inv: $(OBJS) $(OBJSF77) $(LOBJ)sphexp_inv.o
	$(F95) $(FFLAGS) $(LOBJ)sphexp_inv.o $(OBJS)  \
        $(OBJSF77) -l $(LIB)  -o $(BIN)sphexp_inv

$(LOBJ)sphexp_inv.o: sphexp_inv.f90 
	$(F95) $(FFLAGS) sphexp_inv.f90 -c
	$(MV) sphexp_inv.o $(LOBJ)


sphexp_inv_intermediate: $(BIN)sphexp_inv_intermediate $(OBJS) $(OBJSF77)
$(BIN)sphexp_inv_intermediate: $(OBJS) $(OBJSF77) $(LOBJ)sphexp_inv_intermediate.o
	$(F95) $(FFLAGS) $(LOBJ)sphexp_inv_intermediate.o $(OBJS)  \
        $(OBJSF77) -l $(LIB)  -o $(BIN)sphexp_inv_intermediate

$(LOBJ)sphexp_inv_intermediate.o: sphexp_inv_intermediate.f90 
	$(F95) $(FFLAGS) sphexp_inv_intermediate.f90 -c
	$(MV) sphexp_inv_intermediate.o $(LOBJ)


sphexp_inv_short: $(BIN)sphexp_inv_short $(OBJS) $(OBJSF77)
$(BIN)sphexp_inv_short: $(OBJS) $(OBJSF77) $(LOBJ)sphexp_inv_short.o
	$(F95) $(FFLAGS) $(LOBJ)sphexp_inv_short.o $(OBJS)  \
        $(OBJSF77) -l $(LIB)  -o $(BIN)sphexp_inv_short

$(LOBJ)sphexp_inv_short.o: sphexp_inv_short.f90 
	$(F95) $(FFLAGS) sphexp_inv_short.f90 -c
	$(MV) sphexp_inv_short.o $(LOBJ)


sphexp_data_id: $(BIN)sphexp_data_id $(OBJS) $(OBJSF77)
$(BIN)sphexp_data_id: $(OBJS) $(OBJSF77) $(LOBJ)sphexp_data_id.o
	$(F95) $(FFLAGS) $(LOBJ)sphexp_data_id.o $(OBJS)  \
        $(OBJSF77) -l $(LIB)  -o $(BIN)sphexp_data_id

$(LOBJ)sphexp_data_id.o: sphexp_data_id.f90 
	$(F95) $(FFLAGS) sphexp_data_id.f90 -c
	$(MV) sphexp_data_id.o $(LOBJ)


sphexp_backus_gilbert: $(BIN)sphexp_backus_gilbert $(OBJS) $(OBJSF77)
$(BIN)sphexp_backus_gilbert: $(OBJS) $(OBJSF77) $(LOBJ)sphexp_backus_gilbert.o
	$(F95) $(FFLAGS) $(LOBJ)sphexp_backus_gilbert.o $(OBJS)  \
        $(OBJSF77) -l $(LIB)  -o $(BIN)sphexp_backus_gilbert

$(LOBJ)sphexp_backus_gilbert.o: sphexp_backus_gilbert.f90 
	$(F95) $(FFLAGS) sphexp_backus_gilbert.f90 -c
	$(MV) sphexp_backus_gilbert.o $(LOBJ)


sphexp_chi_comp: $(BIN)sphexp_chi_comp $(OBJS) $(OBJSF77)
$(BIN)sphexp_chi_comp: $(OBJS) $(OBJSF77) $(LOBJ)sphexp_chi_comp.o
	$(F95) $(FFLAGS) $(LOBJ)sphexp_chi_comp.o $(OBJS)  \
        $(OBJSF77) -l $(LIB)  -o $(BIN)sphexp_chi_comp

$(LOBJ)sphexp_chi_comp.o: sphexp_chi_comp.f90 
	$(F95) $(FFLAGS) sphexp_chi_comp.f90 -c
	$(MV) sphexp_chi_comp.o $(LOBJ)


sphexp_inv_av: $(BIN)sphexp_inv_av $(OBJS) $(OBJSF77)
$(BIN)sphexp_inv_av: $(OBJS) $(OBJSF77) $(LOBJ)sphexp_inv_av.o
	$(F95) $(FFLAGS) $(LOBJ)sphexp_inv_av.o $(OBJS)  \
        $(OBJSF77) -l $(LIB)  -o $(BIN)sphexp_inv_av

$(LOBJ)sphexp_inv_av.o: sphexp_inv_av.f90 
	$(F95) $(FFLAGS) sphexp_inv_av.f90 -c
	$(MV) sphexp_inv_av.o $(LOBJ)



sphexp_data_idc: $(BIN)sphexp_data_idc $(OBJS) $(OBJSF77)
$(BIN)sphexp_data_idc: $(OBJS) $(OBJSF77) $(LOBJ)sphexp_data_idc.o
	$(F95) $(FFLAGS) $(LOBJ)sphexp_data_idc.o $(OBJS)  \
        $(OBJSF77) -l $(LIB)  -o $(BIN)sphexp_data_idc

$(LOBJ)sphexp_data_idc.o: sphexp_data_idc.f90 
	$(F95) $(FFLAGS) sphexp_data_idc.f90 -c
	$(MV) sphexp_data_idc.o $(LOBJ)




sphexp_data_idlm: $(BIN)sphexp_data_idlm $(OBJS) $(OBJSF77)
$(BIN)sphexp_data_idlm: $(OBJS) $(OBJSF77) $(LOBJ)sphexp_data_idlm.o
	$(F95) $(FFLAGS) $(LOBJ)sphexp_data_idlm.o $(OBJS)  \
        $(OBJSF77) -l $(LIB)  -o $(BIN)sphexp_data_idlm

$(LOBJ)sphexp_data_idlm.o: sphexp_data_idlm.f90 
	$(F95) $(FFLAGS) sphexp_data_idlm.f90 -c
	$(MV) sphexp_data_idlm.o $(LOBJ)


sphexp_global_grid: $(BIN)sphexp_global_grid $(OBJS) $(OBJSF77)
$(BIN)sphexp_global_grid: $(OBJS) $(OBJSF77) $(LOBJ)sphexp_global_grid.o
	$(F95) $(FFLAGS) $(LOBJ)sphexp_global_grid.o $(OBJS)  \
        $(OBJSF77) -l $(LIB)  -o $(BIN)sphexp_global_grid

$(LOBJ)sphexp_global_grid.o: sphexp_global_grid.f90 
	$(F95) $(FFLAGS) sphexp_global_grid.f90 -c
	$(MV) sphexp_global_grid.o $(LOBJ)


sphexp_global_grid_2: $(BIN)sphexp_global_grid_2 $(OBJS) $(OBJSF77)
$(BIN)sphexp_global_grid_2: $(OBJS) $(OBJSF77) $(LOBJ)sphexp_global_grid_2.o
	$(F95) $(FFLAGS) $(LOBJ)sphexp_global_grid_2.o $(OBJS)  \
        $(OBJSF77) -l $(LIB)  -o $(BIN)sphexp_global_grid_2

$(LOBJ)sphexp_global_grid_2.o: sphexp_global_grid_2.f90 
	$(F95) $(FFLAGS) sphexp_global_grid_2.f90 -c
	$(MV) sphexp_global_grid_2.o $(LOBJ)


sphexp_data_power_spectrum_recovery: $(BIN)sphexp_data_power_spectrum_recovery $(OBJS) $(OBJSF77)
$(BIN)sphexp_data_power_spectrum_recovery: $(OBJS) $(OBJSF77) $(LOBJ)sphexp_data_power_spectrum_recovery.o
	$(F95) $(FFLAGS) $(LOBJ)sphexp_data_power_spectrum_recovery.o $(OBJS)  \
        $(OBJSF77) -l $(LIB)  -o $(BIN)sphexp_data_power_spectrum_recovery

$(LOBJ)sphexp_data_power_spectrum_recovery.o: sphexp_data_power_spectrum_recovery.f90 
	$(F95) $(FFLAGS) sphexp_data_power_spectrum_recovery.f90 -c
	$(MV) sphexp_data_power_spectrum_recovery.o $(LOBJ)




$(OBJSF77): $(SRCSF77)	
	$(F95)   $(F77FLAGS) $(SRCSF77) -c 
	$(MV)  *.o $(LOBJ)


$(OBJS): $(SRCS)	
	$(F95)   $(FFLAGS) $(SHFLAGS) $(SRCS) -l $(LIB)    -c
	$(MV)  *.o $(LOBJ)








