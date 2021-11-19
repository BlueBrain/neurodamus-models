# neurodamus-models

This repo is the aggregate Neurodamus repo featuring all simulation models to be built
to run in the BBP supercomputer.

If you are responsible for deployment you will likely deal with this repo directly
to ensure models being built are up to date, namely the shared "common".
Notice that, while neurodamus depends on technical mod files, these have
to be compiled in, together with the models mods, to create a solid "special" Neuron
executable.

Model developers can either use this new structure or the individual model repos.


## Bring models up to date

Models in this aggregate repo are to be stable as a group, and therefore each model 
common should be replaced with the local common.

To avoid doing repetitive work, a script was prepared: update_models.sh
The script takes care of updating submodules and replace common.
