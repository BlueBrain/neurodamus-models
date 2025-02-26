> [!WARNING]
> The Blue Brain Project concluded in December 2024, so development has ceased under the BlueBrain GitHub organization.
> Future development will take place at: https://github.com/openbraininstitute/neurodamus-models

# neurodamus-models

This repo is the aggregate repository featuring all simulation models used within the Blue
Brain Project to run in-silico simulations of mouse brains using the NEURON simulator.

## Installation

### Prerequisites

The models contained here require a compiler, CMake, MPI, and HDF5 to be installed.  On a
Ubuntu system,
```console
sudo apt-get install cmake libopenmpi-dev libhdf5-dev
```
will install the necessary dependencies.  The [NEURON
simulator](https://github.com/neuronsimulator/nrn) and
[Neurodamus](https://github.com/BlueBrain/neurodamus) simulation framework can be
installed via pip:
```console
python -m pip install NEURON-nightly neurodamus
````

Further, an installation of
[libsonatareport](https://github.com/BlueBrain/libsonatareport) is required.  This can be
built with:
```console
git clone https://github.com/BlueBrain/libsonatareport reports/src --recursive --shallow-submodules
cmake \
  -B reports/build \
  -S reports/src \
  -DCMAKE_INSTALL_PREFIX=reports/install \
  -DSONATA_REPORT_ENABLE_SUBMODULES=ON \
  -DSONATA_REPORT_ENABLE_TEST=OFF
cmake --build reports/build
cmake --install reports/build
```

### Compiling the models

The mechanisms for a model can be built as follows, using the `neocortex` model as an
example:
```console
git clone https://github.com/BlueBrain/neurodamus-models.git neurodamus-models/src
export CC=$(which mpicc)
export CXX=$(which mpicxx)
# This can be set directly if the installation location of neurodamus is known
DATADIR=$(python -c "import neurodamus; from pathlib import Path; print(Path(neurodamus.__file__).parent / 'data')")
cmake -B neurodamus-models/build -S neurodamus-models/src \
  -DCMAKE_INSTALL_PREFIX=$PWD/neurodamus-models/install \
  -DCMAKE_INSTALL_RPATH_USE_LINK_PATH=ON \
  -DCMAKE_PREFIX_PATH=$PWD/reports/install \
  -DNEURODAMUS_CORE_DIR=${DATADIR} \
  -DNEURODAMUS_MECHANISMS=neocortex
cmake --build neurodamus-models/build
cmake --install neurodamus-models/build
```
This will also produce a `build_neurodamus.sh` script that can be used to compile future
iterations of the MOD files.

### Testing the models

Following the above step, the compiled mechanisms can be accessed from within NEURON,
e.g., with:
```console
./neurodamus-models/install/bin/special -python -c "from neuron import h; from neurodamus import Neurodamus; h.quit()"
```

For an example on how to run a model with compiled mechanisms on a circuit, see the
[integration test of
Neurodamus](https://github.com/BlueBrain/neurodamus/blob/main/.github/workflows/simulation_test.yml).

## Multi-model builds (MMB project)

For the cases where a combined set of models is desirable a `model_manager.py` tool is provided.

The details on how the models can be combined without clashing (i.e. same `SUFFIX` in MOD files) is to be provided in a configuration
file, where the user specifies a number of source models to be combined with an optional "prefix". Let's refer to the included `model_config.json` for an example:
```
{
    "ALL": {
        "common": {
            "mods": ["common/mod/*.mod"],
            "hocs": ["common/hoc/*.hoc"]
        },
        "neocortex": {
            "prefix": "NCX",
            "mods": [
                "neocortex/mod/v6/*.mod",
                ["neocortex/mod/v5/DetGABAAB.mod", "new_name_for_DetGABAAB.mod"]
            ],
            "hocs": []
        }
    }
}
```

Where
 - `ALL` is the name of this combined model. Resulting set will be stored in `build/ALL`.
 - `common` and `neocortex` are names of two sources of mod/hoc files.
 - For the `common`, we specified "mods" and "hocs" without `prefix` which means they are copied verbatim / without any modifications.
   In our case this is important for existing circuits to keep working with unmodified synapse names.
 - `neocortex` specifies the `NCX` prefix. This means:
    - File names are prefixed with "NCX", unless otherwise specified
    - All .mod files (density channels, Point Processes and Artificial Cells) get their names prefixed
    - All references to these mechanisms in .hoc files are adjusted accordingly.
 - "mods" and "hocs" lists accept a number of entries
    - If the entry is a string then it's the path of the file to be considered (wildcards supported). E.g. `neocortex/mod/v6/*.mod`
    - If the entry is a list, the second argument is taken as the new name.
      E.g. `new_name_for_DetGABAAB.mod` is the new name of `DetGABAAB.mod` (no prefix)

#### Running model manager

```sh
# It's an executable Python application
$ ./model_manager.py
Usage: ./model_manager.py <config_file> [--only-synapses]

$ ./model_manager.py model_config.json
[INFO] Checking 'ALL'
[INFO] Building 'ALL' (only_synapses=False)
[INFO]  - Processing common
[INFO]  - Processing neocortex
[INFO]    > Path neocortex/mod/v6/*.mod
[INFO]    > Path ['neocortex/mod/v5/DetGABAAB.mod', 'new_name_for_DetGABAAB.mod']
[INFO]    > Path neocortex/hoc/AMPANMDAHelper.hoc
[INFO]    > Path neocortex/hoc/GABAABHelper.hoc
[INFO]    > Path neocortex/hoc/GluSynapseHelper.hoc
...

$ ls build/ALL
hoc mod

```

#### Only rename synapses
Notice an optional `--only-synapses` argument. It is intended to apply renaming rules only to synapse files, currently `ProbAMPANMDA_EMS.mod`, `ProbGABAAB_EMS.mod` and `GluSynapse.mod`.

This may be useful in case you are working with existing circuits where you do not want to rename machanisms yet, but still you want to rename and bundle the several synapse types together.

## Acknowledgements

The development of this software was supported by funding to the Blue Brain Project, a
research center of the École polytechnique fédérale de Lausanne (EPFL), from the Swiss
government’s ETH Board of the Swiss Federal Institutes of Technology.

Copyright (c) 2009-2024 Blue Brain Project/EPFL
