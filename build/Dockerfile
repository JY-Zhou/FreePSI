FROM centos:7

RUN set -x
RUN yum install -y --quiet wget git
RUN yum groupinstall -y --quiet 'Development Tools'

RUN wget -q https://sourceforge.net/projects/boost/files/boost/1.61.0/boost_1_61_0.tar.gz
RUN tar xzf boost_1_61_0.tar.gz
WORKDIR boost_1_61_0
RUN sh bootstrap.sh
RUN ./b2
ENV BOOST_ROOT $(pwd)

WORKDIR /
RUN wget -q http://bitbucket.org/eigen/eigen/get/3.2.9.tar.bz2
RUN tar xjf 3.2.9.tar.bz2
ENV EIGEN_ROOT $(pwd)/eigen-eigen-dc6cfdf9bcec

RUN git clone https://github.com/mdshw5/FreePSI.git
WORKDIR FreePSI/src-refined
RUN make
