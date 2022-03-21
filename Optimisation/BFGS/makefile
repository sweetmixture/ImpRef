a.out: src.o test_f.o
	gcc -std=c17 -o a.out src.o test_f.o -lgslcblas -lgsl

src.o: src.c ./../RosenBrock/test_f.h
	gcc -c src.c
test_f.o: ./../RosenBrock/test_f.c ./../RosenBrock/test_f.h
	gcc -c ./../RosenBrock/test_f.c
