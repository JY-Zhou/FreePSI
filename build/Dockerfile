FROM centos:7

RUN set -x
RUN yum install -y --quiet wget git
RUN yum groupinstall -y --quiet 'Development Tools'

RUN yum install -y --quiet libmpc-devel mpfr-devel gmp-devel zlib-devel*
RUN wget -q https://ftp.gnu.org/gnu/gcc/gcc-5.4.0/gcc-5.4.0.tar.gz
RUN tar xfz gcc-5.4.0.tar.gz
WORKDIR gcc-5.4.0
RUN ./configure --with-system-zlib --disable-multilib --enable-languages=c,c++
RUN make -j $(nproc)
RUN make install

WORKDIR /

ENV BOOST_VERSION 1.61.0
ENV BOOST_VERSION_STR boost_1_61_0
RUN wget -q https://sourceforge.net/projects/boost/files/boost/$BOOST_VERSION/${BOOST_VERSION_STR}.tar.gz
RUN tar xzf ${BOOST_VERSION_STR}.tar.gz
WORKDIR $BOOST_VERSION_STR
RUN sh bootstrap.sh
RUN ./b2 install -j $(nproc) || true
ENV BOOST_ROOT $(pwd)

WORKDIR /
ENV EIGEN_VERSION 3.2.9
ENV EIGEN_HASH dc6cfdf9bcec
RUN wget -q http://bitbucket.org/eigen/eigen/get/${EIGEN_VERSION}.tar.bz2
RUN tar xjf ${EIGEN_VERSION}.tar.bz2
ENV EIGEN_ROOT $(pwd)/eigen-eigen-${EIGEN_HASH}

ENV PSI_COMMIT=aaf2313aacd8f387efab614c3814f8fa06a87238
RUN git clone https://github.com/mdshw5/FreePSI.git
WORKDIR /FreePSI/src-refined
RUN git pull
RUN git checkout $PSI_COMMIT
RUN make -j $(nproc) BOOST_ROOT=/${BOOST_VERSION_STR} EIGEN_ROOT=/eigen-eigen-${EIGEN_HASH}
RUN make install

WORKDIR /
ENV JELLYFISH_VERSION 2.2.6
RUN wget -q https://github.com/gmarcais/Jellyfish/releases/download/v${JELLYFISH_VERSION}/jellyfish-${JELLYFISH_VERSION}.tar.gz
RUN tar -xzf jellyfish-${JELLYFISH_VERSION}.tar.gz
WORKDIR jellyfish-${JELLYFISH_VERSION}
RUN ./configure
RUN make -j $(nproc)
RUN make install

RUN yum -y install https://centos7.iuscommunity.org/ius-release.rpm
RUN yum -y install python36u python36u-pip python36u-devel
RUN pip3.6 install numpy scipy

WORKDIR /FreePSI/example
RUN pwd
RUN sh run_FreePSI.sh
