FROM python:3.7

ARG ROSETTA_USER=user
ARG ROSETTA_PASS=1234

ARG PYROSETTA_USER=user
ARG PYROSETTA_PASS=1234

ENV ROSETTA_BIN=/flex/rosetta_src_2019.35.60890_bundle/main/source/bin
ENV ROSETTA_DB=/flex/rosetta_src_2019.35.60890_bundle/main/database

RUN mkdir -p /flex/app

# Download and Compile Rosetta
WORKDIR /flex
RUN curl -f -o rosetta.tar.gz  -u ${ROSETTA_USER}:${ROSETTA_PASS} https://www.rosettacommons.org/downloads/academic/3.11/rosetta_src_3.11_bundle.tgz
RUN tar xvf rosetta.tar.gz
RUN rm rosetta.tar.gz

WORKDIR /flex/rosetta_src_2019.35.60890_bundle/main/source

RUN ./scons.py -j 4 mode=release bin

# Download and install PyRosetta
WORKDIR /flex
RUN curl -f -o pyrosetta.tar.bz2 -u ${PYROSETTA_USER}:${PYROSETTA_PASS} https://graylab.jhu.edu/download/PyRosetta4/archive/release/PyRosetta4.Release.python37.linux/PyRosetta4.Release.python37.linux.release-249.tar.bz2
RUN tar xvf pyrosetta.tar.bz2
RUN rm pyrosetta.tar.bz2
RUN pip install PyRosetta4.Release.python37.linux.release-249/setup

COPY requirements.txt /flex/app
COPY src/** /flex/app

WORKDIR /flex/app
RUN pip install -r requirements.txt
RUN chmod +x run_fpb_new.py

ENTRYPOINT [ "/flex/app/run_fpb_new.py" ]