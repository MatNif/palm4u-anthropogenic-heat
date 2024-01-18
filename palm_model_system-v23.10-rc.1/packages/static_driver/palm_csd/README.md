# PALM Create Static Driver (palm_csd)

Processing tool for creating PALM Input Data Standard (PIDS) conform static driver files from rastered netCDF input.


# Installation

`palm_csd` is written in Python 3 and can be installed using the provided `install` script. Type `bash install -h` to get all available options. In particular, you can use the `-p <install-prefix>` switch to specify the installation directory and the `-e` switch to create a virtual enviroment for `palm_csd` in `<install-prefix>/.palm_csd_venv`. All required dependencies of `palm_csd` are installed in this virtual environment. This ensures that `palm_csd` uses only tested versions of its dependencies. Both, `conda` and `venv` virtual environments are supported with the former taking precedence over the latter. If you do not use `conda` but python from your distribution, ensure that both, `venv` and `pip` are installed. On Debian-based Linux Distributions this can be achieved by the following command:

``` bash
sudo apt-get install python3-pip python3-venv
```

Now, `palm_csd` can be installed with the following commands (please replace `<install-prefix>` with the desired installation directory and, optionally, add other switches):

``` bash
bash install -p <install-prefix>
export PATH=<install-prefix>/bin:${PATH}
```

The following optional command permanently adds this installation to your bash environment:

``` bash
echo 'export PATH=<install-prefix>/bin:${PATH}' >> ~/.bashrc
```

If you did not create a virtual environment for `palm_csd` as explained above, you need to manually ensure that all dependencies are installed on your system. This can be done using `pip`:

> **WARNING** *Do not* use this command if the `install` script created a virtual environment, i.e. if you used the `install` script *with* the `-e` switch.

``` bash
python3 -m pip install -r requirements.txt
```


# Usage

After a successful installation, the executable `palm_csd` was created in the directory `<install-prefix>/bin`. In the case that a virtual environment as set-up by the `install` script is used, that environment is automatically activated before the program code is executed and automatically deactived afterwards.

You need to create a configuration file and provide all necessary netCDF input files. To start the creation of a PIDS conform static driver file, simply type the following command into your terminal:

```bash
palm_csd <path/to/csd-config>
```

Don't forget to replace `<path/to/csd-config>` with the path to a valid configuration file. An example configuration file can be found in subdirectory [share/](share/csd_default.yml).

If you receive the error `PROJ: internal_proj_create_from_database: Cannot find proj.db`, your system's configuration likely [conflicts with one of the installed python dependencies](https://rasterio.readthedocs.io/en/stable/faq.html#why-can-t-rasterio-find-proj-db-rasterio-from-pypi-versions-1-2-0). In this case, unset the environment variable `PROJ_LIB` or `PROJ_DATA` before running `palm_csd`.  With bash, this could be temporarily achieved in the current terminal with `unset PROJ_LIB` or `unset PROJ_DATA`, respectively.  Note that other software might rely on these environment variables so you might want to close the current terminal after running `palm_csd`.

# License

This repository is free and open-source software licensed under the [GNU GPLv3](LICENSE)
