#!/usr/bin/env bash
apt-get -y update
apt-get -y install build-essential python-dev
apt-get -y install python-setuptools
apt-get -y install python-matplotlib python-numpy python-scipy
apt-get -y install make
python3 setup.py install
