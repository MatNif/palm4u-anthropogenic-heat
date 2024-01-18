---
title: Get Started
---
# Get Started

A quick tutorial that gets you from downloading PALM all the way to analysing your first successful simulation.

---

!!! warning
    This site is  Work in Progress.

    ToDo:

    - [x] Requirements
    - [x] Download
    - [x] Installation
    - [x] Usage
    - [ ] Analysis (short)

## Software Requirements

PALM is designed to run on Linux.

In order to successfully install PALM, please meet the following software requirements:

- The bash shell (available at /bin/bash).
- A recent Fortran and C++ Compiler (GNU, Intel, Cray, PGI, NEC).
- The build automation tools cmake and make (for library detection and build coordination).
- The Message Passing Interface (MPI) library with MPI-3 support (compiled with the same compiler as PALM).
- A NetCDF library not earlier than 3.6.3 (compiled with the same Compiler as PALM).
- The FFT library FFTW (PALM also comes with a build-in FFT but with less performance).
- Python 3 (needed by the GUI and other helper routines).
- PyQt5 (needed by the GUI).
- The graphic-package NCL from NCAR (needed by the data visualization tool palmplot).
- The FLEX library BISON parser generator (needed by the chemistry tool kpp4palm)


Very important: It is essential that your NetCDF and MPI library has been built with the same Fortran compiler that is used to compile PALM. Furthermore, in case of a NetCDF4 library with parallel I/O support, the NetCDF library needs to be build with the same MPI library as used for compiling PALM.

## Download

Please download PALM by clicking on "Download PALM" in the top right corner of this documentation website. You can choose between several archive formats. Unpack PALM to the directory you would like to work in (please ensure that file permissions are preserved) and follow the installation steps as described in the next section.

## Installation

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

## Usage

After a successful installation the executables for all [packages](packages) have been linked into the directory `<install-prefix>/bin` and a default PALM configuration file can be found at `<install-prefix>/.palm.config.default`. In case you have installed the python dependencies inside a virtual environment, that environment needs to be active whenever you wand to use PALM. For usage of each of the [packages](packages), please refer to their individual documentation. Next, you need to create your first PALM setup in order to start a simulation. To get a simple preconfigured setup and start your first PALM simulation, please execute the following sequence of commands:

``` bash
mkdir -p "${install_prefix}/JOBS/example_cbl/INPUT"
cp "packages/palm/model/tests/cases/example_cbl/INPUT/example_cbl_p3d" "${install_prefix}/JOBS/example_cbl/INPUT/"
cd ${install_prefix}
palmrun -r example_cbl -c default -a "d3#" -X 4 -v -z
```

## Analysis
