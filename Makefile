all : irp

CCC = g++
CCFLAGS = -g -std=gnu++14
LIBS= -lm -lpthread
PATHLIBS=
TARGETDIR=.
CPPFLAGS += -m64 -O3 -fPIC -fno-strict-aliasing -fexceptions -DNDEBUG -DIL_STD


OBJS2 = \
	$(TARGETDIR)/Client.o \
	$(TARGETDIR)/Vehicle.o \
	$(TARGETDIR)/commandline.o \
	$(TARGETDIR)/Genetic.o \
	$(TARGETDIR)/Individual.o \
	$(TARGETDIR)/LocalSearch.o \
	$(TARGETDIR)/main.o \
	$(TARGETDIR)/Node.o \
	$(TARGETDIR)/Rng.o \
	$(TARGETDIR)/Params.o \
	$(TARGETDIR)/Population.o \
	$(TARGETDIR)/Route.o \
	$(TARGETDIR)/Mutations.o \
	$(TARGETDIR)/LinearPiece.o \
	$(TARGETDIR)/PLFunction.o \
	$(TARGETDIR)/LotSizingSolver.o \
	
	
$(TARGETDIR)/irp: $(OBJS2)
	$(CCC)  $(CCFLAGS) $(CPPFLAGS) $(PATHLIBS) -o $(TARGETDIR)/irp $(OBJS2) $(LIBS)


$(TARGETDIR)/LinearPiece.o: LinearPiece.h LinearPiece.cpp
	$(CCC) $(CCFLAGS) $(CPPFLAGS) -c LinearPiece.cpp -o $(TARGETDIR)/LinearPiece.o

$(TARGETDIR)/PLFunction.o: PLFunction.h PLFunction.cpp
	$(CCC) $(CCFLAGS) $(CPPFLAGS) -c PLFunction.cpp -o $(TARGETDIR)/PLFunction.o

$(TARGETDIR)/LotSizingSolver.o: LotSizingSolver.h LotSizingSolver.cpp
	$(CCC) $(CCFLAGS) $(CPPFLAGS) -c LotSizingSolver.cpp -o $(TARGETDIR)/LotSizingSolver.o

$(TARGETDIR)/Client.o: Client.h Client.cpp
	$(CCC) $(CCFLAGS) $(CPPFLAGS) -c Client.cpp -o $(TARGETDIR)/Client.o

$(TARGETDIR)/Vehicle.o: Vehicle.h Vehicle.cpp
	$(CCC) $(CCFLAGS) $(CPPFLAGS) -c Vehicle.cpp -o $(TARGETDIR)/Vehicle.o
	
$(TARGETDIR)/commandline.o: commandline.h commandline.cpp
	$(CCC) $(CCFLAGS) $(CPPFLAGS) -c commandline.cpp -o $(TARGETDIR)/commandline.o
	
$(TARGETDIR)/Genetic.o: Genetic.h Genetic.cpp
	$(CCC) $(CCFLAGS) $(CPPFLAGS) -c Genetic.cpp -o $(TARGETDIR)/Genetic.o

$(TARGETDIR)/Individual.o: Individual.h Individual.cpp
	$(CCC) $(CCFLAGS) $(CPPFLAGS) -c Individual.cpp -o $(TARGETDIR)/Individual.o

$(TARGETDIR)/LocalSearch.o: LocalSearch.h LocalSearch.cpp
	$(CCC) $(CCFLAGS) $(CPPFLAGS) -c LocalSearch.cpp -o $(TARGETDIR)/LocalSearch.o
	
$(TARGETDIR)/main.o: main.cpp
	$(CCC) $(CCFLAGS) $(CPPFLAGS) -c main.cpp -o $(TARGETDIR)/main.o
	
$(TARGETDIR)/Node.o: Node.h Node.cpp
	$(CCC) $(CCFLAGS) $(CPPFLAGS) -c Node.cpp -o $(TARGETDIR)/Node.o

$(TARGETDIR)/SeqData.o: Rng.h Rng.cpp
	$(CCC) $(CCFLAGS) $(CPPFLAGS) -c Rng.cpp -o $(TARGETDIR)/SeqData.o

$(TARGETDIR)/Params.o: Params.h Params.cpp
	$(CCC) $(CCFLAGS) $(CPPFLAGS) -c Params.cpp -o $(TARGETDIR)/Params.o

$(TARGETDIR)/Population.o: Population.h Population.cpp
	$(CCC) $(CCFLAGS) $(CPPFLAGS) -c Population.cpp -o $(TARGETDIR)/Population.o

$(TARGETDIR)/Route.o: Route.h Route.cpp
	$(CCC) $(CCFLAGS) $(CPPFLAGS) -c Route.cpp -o $(TARGETDIR)/Route.o

$(TARGETDIR)/Mutations.o: Mutations.cpp
	$(CCC) $(CCFLAGS) $(CPPFLAGS) -c Mutations.cpp -o $(TARGETDIR)/Mutations.o

clean:
	del *.o $(objects) 



    

   



