FC  		= ifort
FFLAGS		= -O3 -w
LDFLAGS 	= -lfftw3

default: 	paul_dynamik

paul_dynamik:	main.o   interpolD.o   2dpropagation.o   setpot.o


		${FC} *.o ${FFLAGS} ${LDFLAGS} -o paul_dynamik_pics

main.o:	 main.f90
	${FC} main.f90 ${FFLAGS} -c -o main.o 

interpolD.o:  interpolD.f90
	${FC} interpolD.f90 ${FFLAGS} -c -o interpolD.o

setpot.o:	setpot.f90
	${FC} setpot.f90 ${FFLAGS} -c -o setpot.o 	

2dpropagation.o:  2dpropagation.f90
	${FC} 2dpropagation.f90 ${FFLAGS} -c -o 2dpropagation.o 



clean:
	rm -f *.o 
	rm -f *.mod 
	rm -f paul_dynamik_pics
