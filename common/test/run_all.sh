#!/bin/sh
set -e
which nrniv || (echo "Please load Neuron" && exit 1)

cd "$(dirname "${BASH_SOURCE[0]}")"  # cd to this dir

pushd glusynapse
  nrnivmodl ../../mod/vecevent.mod ../../mod/GluSynapse.mod  # Compile mod
  nosetests -v test_transmission.py
  nosetests -v test_ltpltd.py
popd
