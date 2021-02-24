FROM jupyter/minimal-notebook:latest

USER root

RUN apt-get update && apt-get install -y libgl1-mesa-glx gdb

RUN apt-get -o Dpkg::Options::="--force-confmiss" install --reinstall netbase

USER jovyan

RUN python3 -m pip install gdbgui matplotlib coverage

WORKDIR /home/jovyan

COPY --chown=jovyan:users ./requirements.txt ./requirements.txt

RUN pip install -r ./requirements.txt
