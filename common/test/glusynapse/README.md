How to run tests
================

First we need to compile GluSynapse and VecStim:
$ nrnivmodl ../../common/mod/GluSynapse2019.mod ../../common/mod/vecevent.mod

Then we can use the Nose test suite to run the tests:
$ nosetests -v test\_ltpltd.py
$ nosetests -v test\_transmission.py
Please notice that the tests need to be executed independently to avoid
conflicts (i.e. NEURON global variables overriding).
