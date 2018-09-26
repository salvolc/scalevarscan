#CopyPaste Sachen damit hoffentlich root einigermaßen zufrieden läuft
#--------------------------------------------------------------------
# Root variables (SYSTEM)
ROOTCFLAGS   := $(shell root-config --cflags)
ROOTLIBS     := -lMinuit $(shell root-config --libs)
ROOTTEST	 := -I$(shell root-config --prefix)/test
EXROOT = -L/home/salv/Dokumente/Masterarbeit/MG1/ExRootAnalysis/ -lExRootAnalysis

DELPHES = -L/home/salv/Dokumente/Masterarbeit/Delphes/ -lDelphes

# Root variables (.local) 
# ROOTCFLAGS :=
# ROOTLIBS := 

CXXFLAGS  =  -g -Wall -fPIC -Wno-deprecated -ggdb -O2
CXXFLAGS    += $(ROOTCFLAGS) $(ROOTTEST) 
CXXFLAGS    += -I. -I./include -I/home/salv/Dokumente/Masterarbeit/MG1/ExRootAnalysis/ -I/home/salv/Dokumente/Masterarbeit/MG1/ExRootAnalysis/ExRootAnalysis -I/home/salv/Dokumente/Masterarbeit/Delphes/
LIBS 		= $(ROOTLIBS) -L./ -lTreePlayer -lMinuit $(EXROOT) $(DELPHES)
#--------------------------------------------------------------------
plots = plots/test.pdf
pres = Presentation.pdf
SRC = $(wildcard *.cxx)
#$(info $$plots is [${plots}])
OOBJ = $(SRC:%.cxx=%.o)
#$(info $$OOBJ is [${OOBJ}])

EX = $(wildcard *.cpp)
PRO = $(EX:%.cpp=%)
#$(info $$PRO is [${PRO}])

all: $(PRO) $(OOBJ) $(plots) 
# $(pres)

%.o: %.cxx %.h
	g++ -c $< -o $@ $(CXXFLAGS)

$(PRO): *.cpp $(OOBJ)
	g++ $^ -o $@ -fopenmp $(CXXFLAGS) $(LIBS)

$(plots): $(PRO) truth.py delph.py
	./$(PRO)
	python truth.py >> pyout
	python delph.py >> dpyout

truthl: truth.cpp
	g++ truth.cpp -o truth -fopenmp $(CXXFLAGS) $(LIBS)
	./truth

delphl: delph.cpp
	g++ delph.cpp -o delph -fopenmp $(CXXFLAGS) $(LIBS)
	./delph

truth: truth.cpp truth.py
	g++ truth.cpp -o truth -fopenmp $(CXXFLAGS) $(LIBS)
	./truth
	python truth.py >> pyout

delph: delph.cpp delph.py
	g++ delph.cpp -o delph -fopenmp $(CXXFLAGS) $(LIBS)
	./delph
	python delph.py >> dpyout

tpy: truth.py
	python truth.py >> pyout

dpy: delph.py
	python delph.py >> dpyout

apy: truth.py delph.py
	python truth.py >> pyout
	python delph.py >> dpyout


light: $(PRO)
	./$(PRO)

vertex: 
	g++ vertex.cpp -o vertex -fopenmp $(CXXFLAGS) $(LIBS)
	./vertex

pres:
	cd Vortrag/ && make
	cp Vortrag/build/Presentation.pdf ./

pp:
	python3.5 aufgabe.py
	cd Vortrag/ && make
	cp Vortrag/build/Presentation.pdf ./

#aufgabe: aufgabe.cpp

#$(data): aufgabe
#	./aufgabe

#$(plots): $(data) aufgabe.py
#	python3.5 aufgabe.py


clean:
	@rm ./data/*
	@rm ./plots/*/*

