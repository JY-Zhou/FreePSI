FROM centos:7

RUN set -x
RUN yum install -y --quiet wget git
RUN yum groupinstall -y --quiet 'Development Tools'

RUN curl https://ftp.gnu.org/gnu/gcc/gcc-5.4.0/gcc-5.4.0.tar.bz2 -O
RUN tar xfj gcc-5.4.0.tar.bz2
RUN yum install -y --quiet gmp-devel mpfr-devel libmpc-devel
RUN mkdir gcc-5.4.0-build
WORKDIR gcc-5.4.0-build
RUN ../gcc-5.4.0/configure --enable-languages=c,c++ --disable-multilib
RUN make -j$(nproc) && make install
WORKDIR /

RUN wget -q https://sourceforge.net/projects/boost/files/boost/1.61.0/boost_1_61_0.tar.gz
RUN tar xzf boost_1_61_0.tar.gz
WORKDIR boost_1_61_0
RUN sh bootstrap.sh
RUN ./b2 install -j $(nproc) || true
ENV BOOST_ROOT $(pwd)

WORKDIR /
RUN wget -q http://bitbucket.org/eigen/eigen/get/3.2.9.tar.bz2
RUN tar xjf 3.2.9.tar.bz2
ENV EIGEN_ROOT $(pwd)/eigen-eigen-dc6cfdf9bcec

ENV PSI_COMMIT=aaf2313aacd8f387efab614c3814f8fa06a87238
RUN git clone https://github.com/mdshw5/FreePSI.git
WORKDIR /FreePSI/src-refined
RUN git pull
RUN git checkout $PSI_COMMIT
RUN make -j $(nproc) BOOST_ROOT=/boost_1_61_0 EIGEN_ROOT=/eigen-eigen-dc6cfdf9bcec
RUN make install

WORKDIR /
RUN wget -q https://github.com/gmarcais/Jellyfish/releases/download/v2.2.10/jellyfish-2.2.10.tar.gz
RUN tar -xzf jellyfish-2.2.10.tar.gz
WORKDIR jellyfish-2.2.10
RUN ./configure
RUN make -j $(nproc)
RUN make install

RUN yum -y install https://centos7.iuscommunity.org/ius-release.rpm
RUN yum -y install python36u python36u-pip python36u-devel
RUN pip3.6 install numpy scipy

WORKDIR /FreePSI/example
RUN pwd
RUN sh run_FreePSI.sh
