# palm_cvd

`palm_cvd` is a pre-processor written in python3 to convert sampling coordinates into a so-called "driver file", which contains all information to enable PALM to sample variable at different locations. With this feature, not only single-point measurements can be emulated with PALM, even complex trajectories or tilted profile measurements can be emulated. Please note, the virtual sampling feature of PALM can only be used if the linked NetCDF library is compiled with a parallel output option.

`palm_cvd` is designed to translate sampling coordinates into a driver file, which is a NetCDF file.
Sampling coordinates can be provided in different ways:
- with a geo.json file (for single point sampling only),
- with a list of various NetCDF files (for single- or multiple point sampling, complex profiles and trajectories),
- with a *.csv* or a *.xlsx* file (for complex trajectories only),
- or given by direct specification in a configuration file (for single-point sampling only).

Example files for the different input file formats are located under `/share` or `/data`.


# Installation

First make sure to satisfy the **Software Requirements**. On Debian-based Linux Distributions this can be achieved by the following command:

``` bash
sudo apt-get install python3-pip
```

Also some additional python dependencies are needed, which can be installed using `pip`. In case you wand to use a virtual environment for these dependencies, please make sure to create one first. Afterwards you can install the python dependencies by executing the following command:

``` bash
python3 -m pip install -r requirements.txt
```

Now `palm_cvd` can be installed with the following commands (please replace `<install-prefix>` with the desired installation directory):

``` bash
export install_prefix="<install-prefix>"
bash install -p ${install_prefix}
export PATH=${install_prefix}/bin:${PATH}
```

The following optional command permanently adds this installation to your bash environment:

``` bash
echo 'export PATH=<install-prefix>/bin:${PATH}' >> ~/.bashrc
```

Type `bash install -h` to get all available option of the `install` script.

# Usage

## Options

Type `palm_cvd -h` to get a list of all available option including a description.

Please note, the option *-x*, *-y* and *-z* are optional and depend on static input file availability. In case the PALM setup contains a static input file, the reference coordinate is read from the static file and the options are not required. However, in case no static input file is available, this option is useful to set the PALM reference coordinate. The given sampling coordinates must be given relative to this reference coordinate, else, no sampling output will be created.

## Usage with geo.json input

``` bash
palm_cvd -g data/points.geojson -x 3455249 -y 5424815 -o vm_driver
```

At the moment, input via a `geo.json` format is only realized for single-point measurements. An example input file is provided under `/data/points.geojson`.

With a geo.json input file, the three wind components *u, v, w*, as well as *theta* (potential temperature) and *q* (mixing ratio) will be sampled at the desired locations.

## Usage with configuration file

``` bash
palm_cvd -g share/.cvd.config.default -x 3455249 -y 5424815
```

The configuration file contains different sections, which will be explained in the following.

- `[global]` - In this section, various global attributes can be defined. All of these attributes are optional.
- `[csv_xlsx_data]` - In this section, the input *.csv* or *.xlsx* file (`input_file_coord`) can be described. Further, the variables that shall be measured can be prescribed within the variable `vars_to_be_measured`.
- `[input_from_netcdf]` - In this section, the data path to the NetCDF files can be provided. Please note, `palm_cvd` is designed to follow the [UC]2 data standard. All NetCDF files provided in this standard can be read by `palm_cvd` and sampling coordinates plus various other metadata are processed and written to the driver file.
- `[custom_positions]` - In this section, single-point sampling coordinates can be prescribed by the variable `coordinates<n>`, while *n=1,...,number_positions* is an index in increasing order (starting at 1) to indicate the name the given sampling coordinates. `coordinates<n>` is a list of three values, representing the EUTM, NUMT and z coordinate, while the z-coordinate is the height above the local surface. `number_positions` describes the total number of custom sampling locations given. In `vars_to_be_measured<n>` individual sampling coordinates can be given. Note, the three wind components *u, v, w*, as well as *theta* and the *mixing ratio* will be always sampled.
**Note, sampling variables need to be given in [UC]2 naming convention, which does not necessarily matches the PALM output variables. A list of variable names can be found under <https://palm.muk.uni-hannover.de/trac/wiki/doc/app/virtual_measurement_parameters>.
