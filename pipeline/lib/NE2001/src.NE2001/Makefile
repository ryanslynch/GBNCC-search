
OBJ = NE2001.o dmdsm.NE2001.o  density.NE2001.o neLISM.NE2001.o neclumpN.NE2001.o nevoidN.NE2001.o scattering98.o 
CFLAGS = -fPIC -ff2c 
BINDIR = ../../../bin/
FC = gfortran

NE2001: $(OBJ)
	$(FC) $(CFLAGS) $(OBJ) -o $(BINDIR)NE2001


NE2001.o: NE2001.f 

density.NE2001.o: density.NE2001.f

dmdsm.NE2001.o: dmdsm.NE2001.f

neLISM.NE2001.o: neLISM.NE2001.f

neclumpN.NE2001.o: neclumpN.NE2001.f

nevoidN.NE2001.o: nevoidN.NE2001.f


scattering98.o: scattering98.f

clean: 
	'rm' *.o
