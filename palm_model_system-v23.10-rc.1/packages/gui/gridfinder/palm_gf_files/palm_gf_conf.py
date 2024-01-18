#!/usr/bin/env python
# windows
from configparser import ConfigParser
# linux
#from configparser import configparser
import os

cfg = ConfigParser()


def cfg_write(nx_min, nx_max, ny_min, ny_max, nz_min, nz_max, np_min, np_max, tpn,
              poisfft, switch, temperton, mlt_grid, spectr, dnxny, dpxpy, dpxpy_dev):
    # receives parameters from palm_gf and creates/updates .ini file
    if os.path.exists('./palm_gf_config.ini') is False:
        cfg.add_section("Numerical Grid")
        cfg.add_section("Processor Topology")
        cfg.add_section("Method")
    else:
        cfg.read("palm_gf_config.ini")

    cfg.set("Numerical Grid", "nx min", str(nx_min))
    cfg.set("Numerical Grid", "nx max", str(nx_max))
    cfg.set("Numerical Grid", "ny min", str(ny_min))
    cfg.set("Numerical Grid", "ny max", str(ny_max))
    cfg.set("Numerical Grid", "nz min", str(nz_min))
    cfg.set("Numerical Grid", "nz max", str(nz_max))

    cfg.set("Processor Topology", "np min", str(np_min))
    cfg.set("Processor Topology", "np max", str(np_max))
    cfg.set("Processor Topology", "tasks per node", str(tpn))
    cfg.set("Processor Topology", "ratio nx/ny", str(dnxny))
    cfg.set("Processor Topology", "ratio npex/npey", str(dpxpy))
    cfg.set("Processor Topology", "ratio npex/npey tolerance", str(dpxpy_dev))

    cfg.set("Method", "poisfft", str(poisfft))
    cfg.set("Method", "switch", str(switch))
    cfg.set("Method", "temperton", str(temperton))
    cfg.set("Method", "multigrid", str(mlt_grid))
    cfg.set("Method", "spectra", str(spectr))

    with open("./palm_gf_config.ini", 'w') as cfgfile:
        cfg.write(cfgfile)

    return True


def cfg_read():
    # reads the .ini File and returns all parameters in a list
    cfg.read("palm_gf_config.ini")
    return_list = [cfg.get("Numerical Grid", "nx min"), cfg.get("Numerical Grid", "nx max"),
                   cfg.get("Numerical Grid", "ny min"), cfg.get("Numerical Grid", "ny max"),
                   cfg.get("Numerical Grid", "nz min"), cfg.get("Numerical Grid", "nz max"),
                   cfg.get("Processor Topology", "np min"), cfg.get("Processor Topology", "np max"),
                   cfg.get("Processor Topology", "tasks per node"), cfg.get("Processor Topology", "ratio nx/ny"),
                   cfg.get("Processor Topology", "ratio npex/npey"),
                   cfg.get("Processor Topology", "ratio npex/npey tolerance"), cfg.get("Method", "poisfft"),
                   cfg.get("Method", "switch"), cfg.get("Method", "temperton"), cfg.get("Method", "multigrid"),
                   cfg.get("Method", "spectra")]

    return return_list
