CC := gcc
CFLAGS := -lgsl -lgslcblas -lm

all : Gaussian_poolOrNot


Gaussian_poolOrNot: GSLfun.c Gaussian_poolOrNot.c 
	$(CC) $(CFLAGS) -o $@ GSLfun.c Gaussian_poolOrNot.c

