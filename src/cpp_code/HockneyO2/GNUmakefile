#our makefile variables.   Good to place commonly changed variables
# at the top of your makefile. Then follow with your rules.
HOME = ./..
DIM=2
include $(HOME)/GNUmakefile
hockneyTest: $(HOME)/GNUmakefile GNUmakefile $(OBJS) HockneyTest$(osuffix) ../lib/libfft$(DIM)D.a ../lib/libfft1D.a
	$(CXX) $(CFLAGS) $(FFTINCFLAGS) $(OBJS) $(LIBS) $(TOOLS) HockneyTest$(osuffix) -o testHockney$(DIM)D.exe

clean:
	rm  *.o *.exe *.d *.vtk

realclean:
	rm  *.o *.exe *.d ../lib/*.a ../fftTools/*.o ../fftTools/*.d ../lib/*

-include $(OBJS:.o=.d)
-include $(EOBJS:.o=.d)
