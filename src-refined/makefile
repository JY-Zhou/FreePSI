objects=freePSI.o \
	EMAlgorithm.o \
	KmerHash.o \
	ConjugateGradientProjection.o
sources=freePSI.cpp \
	ConjugateGradientProjection.cpp ConjugateGradientProjection.h \
	EMAlgorithm.cpp EMAlgorithm.h \
	KmerHash.cpp KmerHash.h

GLOBALPARS=-static -O3 -std=c++11 -fopenmp -w -m64 -DNDEBUG

LNKMKL=-Wl,--no-as-needed -lpthread -lm -ldl

MARCOS=-D EIGEN_DONT_PARALLELIZE

freePSI: ${objects}
	rm -f ../freePSI
	g++ ${GLOBALPARS} ${LNKMKL} ${objects} -o freePSI

freePSI.o: ${sources}
	g++ ${GLOBALPARS} ${MARCOS} -I ${BOOST_ROOT} -I ${EIGEN_ROOT} -c freePSI.cpp

ConjugateGradientProjection.o: ConjugateGradientProjection.cpp ConjugateGradientProjection.h \
    EMAlgorithm.cpp EMAlgorithm.h KmerHash.cpp KmerHash.h
	g++ ${GLOBALPARS} ${MARCOS} -I ${BOOST_ROOT} -I ${EIGEN_ROOT} -c ConjugateGradientProjection.cpp

EMAlgorithm.o: EMAlgorithm.cpp EMAlgorithm.h KmerHash.cpp KmerHash.h
	g++ ${GLOBALPARS} ${MARCOS} -I ${BOOST_ROOT} -I ${EIGEN_ROOT} -c EMAlgorithm.cpp

KmerHash.o: KmerHash.cpp KmerHash.h
	g++ ${GLOBALPARS} -I ${BOOST_ROOT} -c KmerHash.cpp

clean:
	rm -f ${objects} freePSI ../freePSI

.PHONY: install
install:
	install -m 555 freePSI freePSI-postProc.py freePSI-summary.py /usr/local/bin
