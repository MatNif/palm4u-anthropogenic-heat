# Specific Emission Mode Input Data in LOD 2

For specific emission modes under LOD 2, emissions can be explicitly provided at a cell level. The LOD 2 emissions data are stored in a separate netCDF file in die model `INPUT` directory, and are read automatically. These input files are prefixed with the model name, followed by the keyword `emis`, then by the specific emission mode. For instance, the LOD 2 input for traffic emissions will be called `[model]_emis_traffic`, while that for the all-purpose generic mode will be called `[model]_emis_generic`, etc..


Please note that all emissions are to be expressed in SI units, i.e. kg, mol, m, and s. Further, unless otherwise indicated, reactive gas phase emission species, such as NO<sub>2</sub> or O<sub>3</sub> are to be provided in moles. On  the other hand, tracer or inert species, such as PM<sub>10</sub> or pollen (in the pollen module), are to be specified in kilograms.

## Temporal Emission Profiles

Temporal profiles of emission are stored using the following format:

`YYYY-MM-DD HH:mm:ss +ZZ`

where `YYYY-MM-DD` represents the date in full numeric format, `HH:mm:ss` represents the local time in 24-hour format, `ZZ` is the time zone.  It is generally recommended that the profiles are specified in coordinated universal time (UTC) to maximize reusability.

The user can specify profile with arbitrary start and end times, as well as intervals, in chornological order. Emissions data defined at the time closest but not later than the current model time will be used as input. Should the model start (as defined in the `origin_date_time` option in the `_p3d` file) prior to the earliest defined emission data, the earliest emission data will be used.

## netCDF File Format

Although the netCDF 4 API can be used to create the LOD 2 emisisons input files, storage of the emission data follow the netCDF 3 convention to maintain compatibility with other netCDF files used in the model. In particular, only the time dimension can be used as the `UNLIMITED` dimension, and strings are to be stored as arrays of characters, as opposed to the `STRING` variable type.  User-defined data structures are not used.

### Dimensions

The definition of the LOD 2 emissions variables are predicated on the following mandatory dimensions:

| Dimension name | Value      | Description                           |
|----------------|------------|---------------------------------------|
| `ntime`        | at least 1 | Designated `UNLIMITED` dimension      |
| `field_length` | 64         | Fixed length of all string variables  |
| `nspecies`     | at least 1 | Number of emissions species           |
| `nvsrc`        | at least 1 | Number of volumetric emission sources |

Note that the `field_length` dimension must be set to 64. Other dimensions must be greater than zero.

### Variables

The mandatory variables for storage of emissions can then be defined using the above dimensions:
 
| Variable name | netCDF data type | Dimension(s) | Description |
|---|---|---|---|
| `timestamp` | `NC_CHAR` | (`ntime`, `field_length`) | Individual time stamps for each set of emission data |
| `species` | `NC_CHAR` | (`nspecies`, `field_length`) | Names of individual chemical species |
| `vsrc_i`, `vsrc_j`, `vsrc_k`  | `NC_INT` | `nvsrc` | Cell coordinates (i,j,k) of individual emission sources |
| `vsrc_[species]` | `NC_FLOAT` | (`ntime`, `nvsrc`) | Volumetric emission of chemical species at each source location |

The volumetric emission sources of each species are stored in separate variables with the prefix `vsrc_`, followed by the species name specified in the `species` variable. They are defined for all source cell locations (i,j,k) defined in the varibles `vsrc_i`, `vsrc_j`, and `vsrc_k`. As such, species not emitting at a particular location on a given time stamp must be given a emission value zero. The species indicated in the variable `species` can be of any name, but in general they should appear in the active kinetic mechanism used in the chemistry model. Species names not appeared in said mechanism and corresponding emissions will be ignored. _NOTE_:  All volumetric emission sources are expressed in terms of mol/(m<sup>3</sup>s) for gas-phase species and kg/(m<sup>3</sup>s) for particulate matter (PM).

### Example

The following is the header of an example netCDF emission input file under LOD 2:

```
dimensions:
        ntime = UNLIMITED ; // (1 currently)
        nvsrc = 3 ;
        nspecies = 6 ;
        field_length = 64 ;
variables:
        char timestamp(ntime, field_length) ;
                timestamp:description = "Time stamps" ;
        char species(nspecies, field_length) ;
                species:description = "Emission species" ;
        int vsrc_i(nvsrc) ;
                vsrc_i:description = "i grid indices for volume source location" ;
        int vsrc_j(nvsrc) ;
                vsrc_j:description = "j grid indices for volume source location" ;
        int vsrc_k(nvsrc) ;
                vsrc_k:description = "k grid indices for volume source location" ;
        float vsrc_O3(ntime, nvsrc) ;
                vsrc_O3:description = "volume source values for O3" ;
        float vsrc_NO(ntime, nvsrc) ;
                vsrc_NO:description = "volume source values for NO" ;
        float vsrc_NO2(ntime, nvsrc) ;
                vsrc_NO2:description = "volume source values for NO2" ;
        float vsrc_PM10(ntime, nvsrc) ;
                vsrc_PM10:description = "volume source values for PM10" ;
        float vsrc_PM25(ntime, nvsrc) ;
                vsrc_PM25:description = "volume source values for PM25" ;
        float vsrc_SO2(ntime, nvsrc) ;
                vsrc_SO2:description = "volume source values for SO2" ;
```

While not shown, it should be apparent to the user that the species defined in the netCDF files correspond to the volumetric source variables present.
