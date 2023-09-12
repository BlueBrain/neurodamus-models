#!/bin/env python
#
# Blue Brain Project
#
# This utility combines the several mod and hoc files into combined builds
# according to the config provided as argument
#

import json
import logging
import shutil
import sys
import os
import re
from dataclasses import dataclass
from glob import glob


MOD_PATCH_EXPRS = ("SUFFIX ", "POINT_PROCESS ")  # keep spaces
HOC_NAMES_PREFIX = (
    "ProbAMPANMDA_EMS",
    "ProbGABAAB_EMS",
    "GluSynapse",
)
ONLY_SYNAPSES_PATCH_MOD_FILES = {
    "ProbAMPANMDA_EMS.mod",
    "ProbGABAAB_EMS.mod",
    "GluSynapse.mod",
}


@dataclass
class ModelComponent:
    mods: list
    hocs: list
    prefix: str = None


def merge_model(build_name: str, components: "dict[str, ModelComponent]", only_synapses: False):
    out_mod = f"build/{build_name}/mod"
    out_hoc = f"build/{build_name}/hoc"
    os.makedirs(out_mod)
    os.makedirs(out_hoc)

    for name, component in components.items():
        logging.info(" - Processing %s", name)
        if component.prefix is None:
            for mod_path in component.mods:
                copy_all(mod_path, out_mod)
            for hoc_path in component.hocs:
                copy_all(hoc_path, out_hoc)
        else:
            prefix = component.prefix
            mod_copy_f = copy_patch(
                MOD_PATCH_EXPRS,
                suffix=f"{prefix}_",
                file_prefix=prefix,
                n_times=1,
                subset=ONLY_SYNAPSES_PATCH_MOD_FILES if only_synapses else None
            )
            for mod_path in component.mods:
                logging.info("   > Path %s", mod_path)
                if isinstance(mod_path, list):
                    assert len(mod_path) == 2, "Mod path list must have format [src_file, new_name]"
                    assert os.path.isfile(mod_path[0])
                    mod_copy_f(mod_path[0], out_mod, mod_path[1])
                else:
                    copy_all(mod_path, out_mod, mod_copy_f)

            hoc_copy_f = copy_patch_hoc(prefix)
            for hoc_path in component.hocs:
                logging.info("   > Path %s", hoc_path)
                copy_all(hoc_path, out_hoc, hoc_copy_f)


def copy_all(src, dst, copyfunc=shutil.copy):
    """Copy/process all files in a src dir into a destination dir."""
    for src_pth in glob(src):
        copyfunc(src_pth, dst)


def copy_patch(find_expr, prefix="", suffix="", file_prefix=None, n_times=-1, subset=None):
    """Creates a copy-patcher function according to params

    Args:
        find_expr: The expression(s) to find
        prefix: The prefix to append when `find_expr` is found
        suffix: The suffix to append
    """
    join_path = os.path.join
    filename_prefix = "" if file_prefix is None else file_prefix + "_"
    basename = os.path.basename
    if isinstance(find_expr, str):
        find_expr = (find_expr,)

    def _copy_f(src_f, dst_dir, custom_name=None):
        filename = basename(src_f)
        dst_name = custom_name + os.path.splitext(src_f)[1] if custom_name else filename

        if not custom_name and subset is not None and filename not in subset:
            # In compat mode don't care for overwrites
            dst_f = join_path(dst_dir, dst_name)
            shutil.copy(src_f, dst_f)
            return

        dst_f = join_path(dst_dir, dst_name if custom_name else filename_prefix + dst_name)

        if os.path.exists(dst_f):
            raise Exception("File already exists: " + dst_f)

        with open(src_f) as src:
            full_str = src.read()

        for expr in find_expr:
            if custom_name:
                full_str = re.sub(expr + ".*", expr + custom_name, full_str)
            else:
                new_expr = expr
                if prefix:
                    new_expr = prefix + new_expr
                if suffix:
                    new_expr += suffix
                logging.debug(f"Replacing '{expr}' with '{new_expr}' in [{src_f} -> {dst_f}]")
                full_str = full_str.replace(expr, new_expr, n_times)

        with open(dst_f, "w") as dst:
            dst.write(full_str)

    return _copy_f


def copy_patch_hoc(prefix):
    basename = os.path.basename
    name_re = re.compile(r"(begin|end)template *(.*)")
    name_options = "|".join(HOC_NAMES_PREFIX)
    pp_call_re = re.compile(fr"new *({name_options}) ?\(")

    def copy_f(hoc_f, dst_dir):
        """Rename and adapt the hoc file for the new mods"""
        with open(hoc_f) as src:
            full_str = src.read()

        full_str = name_re.sub(fr"\1template {prefix}_\2", full_str)
        full_str = pp_call_re.sub(fr"{prefix}_\1(", full_str)

        dst_f = os.path.join(dst_dir, prefix + "_" + basename(hoc_f))
        with open(dst_f, "w") as dst:
            dst.write(full_str)

    return copy_f


if __name__ == "__main__":

    if len(sys.argv) == 1 or sys.argv[1] in ("-h", "--help"):
        print(f"Usage: {sys.argv[0]} <config_file> [--only-synapses]")
        exit(-1)

    logging.basicConfig(level=logging.INFO, format="[%(levelname)s] %(message)s")

    filename = sys.argv[1]
    options = sys.argv[2:]
    only_synapses = "--only-synapses" in options

    if not os.path.isfile(filename):
        logging.error("Invalid config file path")
        exit(-1)

    builds = json.load(open(filename))

    for build_name, components in builds.items():
        logging.info("Checking '%s'", build_name)
        components = {name: ModelComponent(**comp) for name, comp in components.items()}
        logging.info(f"Building '{build_name}' (only_synapses={only_synapses})")
        merge_model(build_name, components, only_synapses)

    logging.info("Finished")
