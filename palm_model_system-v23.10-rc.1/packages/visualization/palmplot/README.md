# palmplot

NCL-based visualization tool for PALM netCDF output data.


# Installation

First make sure to satisfy the **Software Requirements**. On Debian-based Linux Distributions this can be achieved by the following command:

``` bash
sudo apt-get install ncl-ncarg
```


Now palmplot can be installed with the following commands (please replace `<install-prefix>` with the desired installation directory):

``` bash
bash install -p <install-prefix>
export PATH=<install-prefix>/bin:${PATH}
```

The following optional command permanently adds this installation to your bash environment:

``` bash
echo 'export PATH=<install-prefix>/bin:${PATH}' >> ~/.bashrc
```

Type `bash install -h` to get all available option of the `install` script.


# Usage

After a successful installation, the executable `palmplot` has been linked into the directory `<install-prefix>/bin`. Next, you need to have a set of netCDF output files from a finished PALM simulation ready. To learn how to start creating plots from PALM netCDF output files with `palmplot`, please see the [documentation](https://palm.muk.uni-hannover.de/trac/wiki/doc/app/ncl).


# License

This repository is free and open-source software licensed under the [GNU GPLv3](LICENSE)
