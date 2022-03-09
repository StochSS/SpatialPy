FROM jupyter/minimal-notebook:latest

USER root

RUN apt-get update && apt-get install -y libgl1-mesa-glx gdb

RUN apt-get install -y  libtiff5-dev libjpeg8-dev libopenjp2-7-dev zlib1g-dev \
    libfreetype6-dev liblcms2-dev libwebp-dev tcl8.6-dev tk8.6-dev python3-tk \
        libharfbuzz-dev libfribidi-dev libxcb1-dev

RUN apt-get -o Dpkg::Options::="--force-confmiss" install --reinstall netbase

USER jovyan

RUN python3 -m pip install gdbgui matplotlib coverage

WORKDIR /home/jovyan

COPY --chown=jovyan:users ./requirements.txt ./requirements.txt

RUN pip install -r ./requirements.txt

CMD start-notebook.sh
