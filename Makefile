CC = g++
IFLAGS  = -I/usr/local -I/usr/local/include/openbabel-2.0
LDFLAGS = -L/usr/local/lib -L/usr/local/lib/openbabel/3.0.a1

bfsenum: bfsenum.cpp bfsenum.hpp
	$(CC) $(IFLAGS) $(LDFLAGS) -O3 -std=c++0x bfsenum.cpp -o bfsenum -lopenbabel -lboost_program_options -lcurl
