make:
	ifort -c *.f90
	ifort  *.o -o lrz.x
	rm *.o
gen:
	g++ gen.c -o gen.x
	./gen.x > data.txt
clean:
	rm *.o
	rm a.out
tags:  
	@ctags  *.f90 
	@echo tags updated!
