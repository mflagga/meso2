.SUFFIXES:

# kompilator
CXX = g++
# interpreter pythona
PC = python3
# plik wykonywalny
EXEC = main.out
# flagi kompilacji
CXXFLAGS = -Wall -std=c++20
# biblioteki
LIBS = -lm
# pliki źródłowe
SRCS = main.cpp
# pliki z danymi
DATA = psi.dat misc.dat fps.dat V.dat Exx.dat
ANIM = evol.mp4 evol.gif
FPS = $(shell cat fps.dat)

all: evol.mp4 evol.gif expect.png

# kompiluj
$(EXEC): $(SRCS)
	$(CXX) $(SRCS) $(CXXFLAGS) $(LIBS) -o $@

# uruchom
$(DATA): $(EXEC)
	./$<

# wykreśl
frames/framesdone.txt expect.png: psi.dat misc.dat V.dat wykres.py Exx.dat
	rm -f frames/frame_*.png
	rm -f frames/framesdone.txt
	$(PC) wykres.py
	touch $@

# animuj
evol.mp4: frames/framesdone.txt fps.dat
	ffmpeg -framerate $(FPS) -i frames/frame_%04d.png -y -loglevel quiet -crf 18 evol.mp4
	
evol.gif: frames/framesdone.txt fps.dat
	ffmpeg -framerate $(FPS) -loglevel quiet -i frames/frame_%04d.png -vf "split[s0][s1];[s0]palettegen=stats_mode=full[p];[s1][p]paletteuse=dither=sierra2_4a" -y frames/evol.gif

# posprzątaj
clean: 
	rm -f $(EXEC)
	rm -f $(DATA)
	rm -f frames/frame_*.png
	rm -f frames/framesdone.txt
	rm -f evol.mp4 expect.png

.PHONY: all clean
