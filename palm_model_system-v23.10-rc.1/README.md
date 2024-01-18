# PALM Model System

This is a collection of software packages that contain PALM and numerous additional tools that support PALM and extend its functionality. PALM is an advanced and modern meteorological model system for atmospheric and oceanic boundary-layer flows. It has been developed as a turbulence-resolving large-eddy simulation (LES) model that is especially designed for performing on massively parallel computer architectures.

The following [packages](packages) are included (in alphabetical order):

-   [Chemistry](packages/chemistry)
    -   **kpp4palm** - A preprocessor that generates the source file, where the gas phase chemistry rate equations are solved within PALM

-   [Dynamic Driver](packages/dynamic_driver)
    -   **inifor** - A mesoscale Interface for Initializing and Forcing PALM using offline nesting
    -   **wrf_interface** - Scripts for processing of WRF and CAMx files to PALM dynamic driver

-   [GUI](packages/gui)
    -   **gridfinder** - A graphical Support tool to find valid PALM grid configurations
    -   **palmrungui** - A graphical user interface for the palmrun script
    -   **watchdog** - A graphical monitoring tool for PALM running remotely on an HPC cluster

-   [Multi Agent System](packages/multi_agent_system)
    -   **agent_preprocessing** - A tool that creates a navigation mesh, which is necessary for the use of the Multi Agent System in PALM

-   [PALM](packages/palm)
    -   **model** - The PALM model source code

-   [Static Driver](packages/static_driver)
    -   **create_basic_static_driver** - A simple python script to create simple static drivers for PALM
    -   **palm_csd** - A tool for creating PIDS conform static drivers from rastered netCDF input

-   [Virtual Measurements](packages/virtual_measurements)
    -   **merge_virtual_measurements** - A script that merges all virtual measurement output files from PALM into one netCDF file
    -   **palm_cvd** - A tool for creating PIDS conform virtual measurement setup files

-   [Visualization](packages/visualization)
    -   **palmplot** - A tool that creates pre-configured plots from PALM netCDF output data

# Installation

First make sure to satisfy the **Software Requirements**. On Debian-based Linux Distributions this can be achieved by the following command:

``` bash
sudo apt-get install gfortran g++ make cmake coreutils libopenmpi-dev openmpi-bin libnetcdff-dev netcdf-bin libfftw3-dev python3-pip python3-pyqt5 flex bison ncl-ncarg
```

Also some additional python dependencies are needed, which can be installed using `pip`. In case you wand to use a virtual environment for these dependencies, please make sure to create one first. Afterwards you can install the python dependencies by executing the following command:

``` bash
python3 -m pip install -r requirements.txt
```

Now the PALM model system can be installed with the following commands (please replace `<install-prefix>` with the desired installation directory):

``` bash
export install_prefix="<install-prefix>"
bash install -p ${install_prefix}
export PATH=${install_prefix}/bin:${PATH}
```

The following optional command permanently adds this installation to your bash environment:

``` bash
echo "export PATH=${install_prefix}/bin:\${PATH}" >> ~/.bashrc
```

Type `bash install -h` to get all available option of the `install` script. During installation the script calls the respective `install` script of all [packages](packages) in this repository and installs them to the chosen `<install-prefix>` directory. Therefore, it is not necessary to manually install any of the [packages](packages).

You can test your installation with the following commands:

``` bash
palmtest --cases urban_environment_restart --cores 4
```

# Usage

After a successful installation the executables for all [packages](packages) have been linked into the directory `<install-prefix>/bin` and a default PALM configuration file can be found at `<install-prefix>/.palm.config.default`. In case you have installed the python dependencies inside a virtual environment, that environment needs to be active whenever you wand to use PALM. For usage of each of the [packages](packages), please refer to their individual documentation. Next, you need to create your first PALM setup in order to start a simulation. To get a simple preconfigured setup and start your first PALM simulation, please execute the following sequence of commands:

``` bash
mkdir -p "${install_prefix}/JOBS/example_cbl/INPUT"
cp "packages/palm/model/tests/cases/example_cbl/INPUT/example_cbl_p3d" "${install_prefix}/JOBS/example_cbl/INPUT/"
cd ${install_prefix}
palmrun -r example_cbl -c default -a "d3#" -X 4 -v -z
```

# License

This repository is free and open-source software licensed under the [GNU GPLv3](LICENSE)

