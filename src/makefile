objects=freePSI.o \
	EMAlgorithm.o \
	KmerHash.o \
	ConjugateGradientProjection.o 
sources=freePSI.cpp \
	ConjugateGradientProjection.cpp ConjugateGradientProjection.h \
	EMAlgorithm.cpp EMAlgorithm.h \
	KmerHash.cpp KmerHash.h 

boost=/home/jyzhou/lib/cpp/boost_1_61_0
eigen=/home/jyzhou/lib/cpp/eigen-3.3
GLOBALPARS=-O3 -std=c++11 -fopenmp -w -m64 -DNDEBUG

MKL=${MKLROOT}/include
LNKMKL=-Wl,--no-as-needed -L${MKLROOT}/lib/intel64 -lmkl_intel_lp64 -lmkl_core -lmkl_gnu_thread -lpthread -lm -ldl

MARCOS=-D EIGEN_DONT_PARALLELIZE -D EIGEN_USE_MKL_ALL
#MARCOS=-D EIGEN_USE_MKL_ALL #-D EIGEN_RUNTIME_NO_MALLOC

freePSI: ${objects}
	rm -f ../freePSI
	#g++ ${GLOBALPARS} -lgfortran ${LNKMKL} ${objects} -o freePSI
	g++ ${GLOBALPARS} ${LNKMKL} ${objects} -o freePSI
	cp ./freePSI ../

freePSI.o: ${sources}
	g++ ${GLOBALPARS} ${MARCOS} -I ${boost} -I ${eigen} -I ${MKL} -c freePSI.cpp

#BiasCorrection.o: BiasCorrection.cpp BiasCorrection.h \
#    EMAlgorithm.cpp EMAlgorithm.h KmerHash.cpp KmerHash.h
#	g++ ${GLOBALPARS} ${MARCOS} -I ${boost} -I ${eigen} -I ${MKL} -c BiasCorrection.cpp
#
#RosenGradientDescend.o: RosenGradientDescend.cpp RosenGradientDescend.h \
#    EMAlgorithm.cpp EMAlgorithm.h KmerHash.cpp KmerHash.h
#	g++ ${GLOBALPARS} ${MARCOS} -I ${boost} -I ${eigen} -I ${MKL} -c RosenGradientDescend.cpp
#
#Optimization.o: Optimization.cpp Optimization.h EMAlgorithm.cpp EMAlgorithm.h KmerHash.cpp KmerHash.h SLSQP.f
#	g++ ${GLOBALPARS} ${MARCOS} -I ${boost} -I ${eigen} -I ${MKL} -c Optimization.cpp

#SLSQP.o: SLSQP.f
#	gfortran ${GLOBALPARS} -c SLSQP.f

ConjugateGradientProjection.o: ConjugateGradientProjection.cpp ConjugateGradientProjection.h \
    EMAlgorithm.cpp EMAlgorithm.h KmerHash.cpp KmerHash.h
	g++ ${GLOBALPARS} ${MARCOS} -I ${boost} -I ${eigen} -I ${MKL} -c ConjugateGradientProjection.cpp

EMAlgorithm.o: EMAlgorithm.cpp EMAlgorithm.h KmerHash.cpp KmerHash.h SLSQP.f
	g++ ${GLOBALPARS} ${MARCOS} -I ${boost} -I ${eigen} -I ${MKL} -c EMAlgorithm.cpp

KmerHash.o: KmerHash.cpp KmerHash.h
	g++ ${GLOBALPARS} -I ${boost} -c KmerHash.cpp

clean:
	rm -f ${objects} freePSI ../freePSI
