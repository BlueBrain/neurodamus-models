#!/bin/bash

module load unstable neuron python

python -m venv venv_test
. venv_test/bin/activate
pip install nose numpy scipy

sh test/run_all.sh

