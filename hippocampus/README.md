Simulation Models :: hippocampus
--------------------------------

These hoc and mod files implement the hippocampus cell mechanism.


Versions
--------

2019.1 - First import after splitting of neurodamus and models
2021.11 - Start using aggregate neurodamus-models repo, drop common submodule


Common models
-------------

TLDR: Run `./fetch_common.bash` to init and update common.

Most Blue Brain models depend on a set of common mods, among themProbAMPANMDA and ProbGABAAB. 
Previously the current repo would include them as a submodule. However, besides overly
complicating the deployment process, it could lead to situations of outdated versions of these files.

Since 2021.11 all BBP maintained models therefore drop submodules and implement a fetch_common.bash
script which will ensure that you have locally the latest common files.

Therefore, every time you wish to init or fetch updates for common, simply:

 > `./fetch_common.bash`

