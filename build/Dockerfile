FROM centos:7

RUN set -x
RUN yum install -y --quiet wget git
RUN yum groupinstall -y --quiet 'Development Tools'

RUN wget https://sourceforge.net/projects/boost/files/boost/1.61.0/boost_1_61_0.tar.gz
RUN tar xzf boost_1_61_0.tar.gz
RUN cd boost_1_61_0
RUN sh bootstrap.sh
RUN ./b2
RUN cd ..
ENV BOOST_ROOT $(pwd)/boost_1_61_0

RUN wget http://bitbucket.org/eigen/eigen/get/3.2.9.tar.bz2
RUN tar xjf 3.2.9.tar.bz2
ENV EIGEN_ROOT $(pwd)/eigen-eigen-dc6cfdf9bcec

RUN git clone https://github.com/mdshw5/FreePSI.git
RUN cd FreePSI/src-refined
RUN make
