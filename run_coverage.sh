#!/bin/bash
#pip3 install  coverage
#pip3 install python-libsbml

coverage run --source=spatialpy test/run_unit_tests.py -m develop
coverage run --source=spatialpy test/run_integration_tests.py -m develop
coverage run --source=spatialpy test/run_system_tests.py -m develop
coverage html
