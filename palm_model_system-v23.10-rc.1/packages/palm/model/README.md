# PALM Model

This is the PALM model source code repository. PALM is an advanced and modern meteorological model system for atmospheric and oceanic boundary-layer flows. It has been developed as a turbulence-resolving large-eddy simulation (LES) model that is especially designed for performing on massively parallel computer architectures.

PALM comes prepackages with the RRTMG [library](lib).

# Installation

First make sure to satisfy the **Software Requirements**. On Debian-based Linux Distributions this can be achieved by the following command:

``` bash
sudo apt-get install gfortran make cmake libopenmpi-dev openmpi-bin libnetcdff-dev netcdf-bin libfftw3-dev python3-pip
```

Also some additional python dependencies are needed, which can be installed using `pip`. In case you wand to use a virtual environment for these dependencies, please make sure to create one first. Afterwards you can install the python dependencies by executing the following command:

``` bash
python3 -m pip install -r requirements.txt
```

Now the PALM model can be installed with the following commands (please replace `<install-prefix>` with the desired installation directory):

``` bash
export install_prefix="<install-prefix>"
bash install -p ${install_prefix}
export PATH=${install_prefix}/bin:${PATH}
```

The following optional command permanently adds this installation to your bash environment:

``` bash
echo "export PATH=${install_prefix}/bin:\${PATH}" >> ~/.bashrc
```

Type `bash install -h` to get all available option of the `install` script. During installation the script also calls the respective `install` script of the RRTMG and installs it to the chosen `<install-prefix>` directory. Therefore, it is not necessary to manually install the RRTMG.

You can test your installation with the following command:

``` bash
palmtest --cases urban_environment_restart --cores 4
```

# Usage

After a successful installation the PALM executables have been linked into the directory `<install-prefix>/bin` and a default PALM configuration file can be found at `<install-prefix>/.palm.config.default`. In case you have installed the python dependencies inside a virtual environment, that environment needs to be active whenever you wand to use PALM. Next, you need to create your first PALM setup in order to start a simulation. To get a simple preconfigured setup and start your first PALM simulation, please execute the following sequence of commands:

``` bash
mkdir -p "${install_prefix}/JOBS/example_cbl/INPUT"
cp "packages/palm/model/tests/cases/example_cbl/INPUT/example_cbl_p3d" "${install_prefix}/JOBS/example_cbl/INPUT/"
cd ${install_prefix}
palmrun -r example_cbl -c default -a "d3#" -X 4 -v -z
```


# License

This repository is free and open-source software licensed under the [GNU GPLv3](LICENSE)
