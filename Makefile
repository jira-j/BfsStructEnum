CC = g++
IFLAGS  = -I/usr/local  -I/usr/include/openbabel-2.0
LDFLAGS = -L/usr/lib64/openbabel -lboost_program_options-mt -lcurl

bfsenum: bfsenum.cpp bfsenum.hpp
	$(CC) $(IFLAGS) $(LDFLAGS) -O3 -std=c++0x bfsenum.cpp -o bfsenum -lopenbabel
