.SUFFIXES:

# kompilator
CXX = g++
# interpreter pythona
PC = python3
# plik wykonywalny
EXEC = main.out
# flagi kompilacji
CXXFLAGS = -Wall -Wextra -std=c++20
# biblioteki
LIBS = -lm
# pliki źródłowe
SRCS = main.cpp
# pliki z danymi
DATA = psi.dat misc.dat fps.dat V.dat
FPS = $(shell cat fps.dat)

all: evol.mp4

# kompiluj
$(EXEC): $(SRCS)
	$(CXX) $(SRCS) $(CXXFLAGS) $(LIBS) -o $@

# uruchom
$(DATA): $(EXEC)
	./$<

# wykreśl
frames/framesdone.txt: $(DATA) wykres.py
	rm -f frames/frame_*.png
	rm -f frames/framesdone.txt
	$(PC) wykres.py
	touch $@

# animuj
evol.mp4: frames/framesdone.txt fps.dat
	ffmpeg -framerate $(FPS) -i frames/frame_%04d.png -y -loglevel quiet -crf 18 evol.mp4

# posprzątaj
clean: 
	rm -f $(EXEC)
	rm -f $(DATA)
	rm -f frames/frame_*.png
	rm -f frames/framesdone.txt
	rm -f evol.mp4

.PHONY: all clean