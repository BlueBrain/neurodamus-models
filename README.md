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
