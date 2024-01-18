# Anthropogenic Emissions Model

Emissions of gases and/or passive compounds can be provided by the user in different levels of detail (LOD). Three modes were implemented. Two modes require emission data from file. Gridded emission information can be provided in two LODs via a netCDF file (see PALM Input Data Standard document​):

LOD 1 files require gridded annual emission information which are temporally disaggregated by PALM-4U using sector-specific standard time factors

LOD 2 files must contain gridded pre-processed emission information that is already temporally disaggregated. At the moment hourly emission data is required. Other temporal intervals will be possible in later versions. IMPORTANT: In this mode the initial date of the simulation has to coincide with the first day for which emission values are available.

In a third mode traffic emissions are parameterized on the basis of the street type classification provided by  [OpenStreetMap](https://www.openstreetmap.org) and street type specific emissions factors. This option, named ’Parameterized’ mode, is currently only implemented for the traffic sector:

LOD 0 requires a mean surface emission per day (in micromol/m2/d for gases and in kg/m2/d for particulate matter) for each desired species contained in the applied chemical mechanism together with a weight differing between main and side road emissions defined in the PALM Fortran namelist file (see Fortran namelist file and chemistry parameter list). Furthermore, the street type classification provided by Open Street Map must be included in the PALM static driver (see Static input file).

To account for the temporal distribution of the parameterized traffic emissions a diurnal profile based on urban traffic counts is implemented to disaggregate the total emission of the different species to hourly values. Currently, the applied temporal profile is the same for all species and main and side roads.

Please note that LOD 0 is currently not applicable if the horizontal grid exceeds 10 to 15 m at maximum (depending on the width of the streets in the considered area). In this context, it must also be considered, that streets, which are only one grid point wide in the STATIC file will probably be eliminated by the topography filtering during runtime, i.e. emissions from steets which are only one grid point wide will be supressed.

## Input File Specifications

Surface flux data defined under LOD 2 are to be stored in the `_chemistry` file in netCDF format in the `INPUT` directory, which will be read automatically. Although the netCDF 4 API can be used to create the LOD 2 emisisons input files, storage of the emission data follow the netCDF 3 convention to maintain compatibility with other netCDF files used in the model. In particular, only the time dimension can be used as the `UNLIMITED` dimension, and strings are to be stored as arrays of characters, as opposed to the `STRING` variable type.  User-defined data structures are not used.

### Legacy Mode

The `_p3d` chemistry option `emiss_read_legacy_mode` dictates how data are stored in the `_chemistry` file, as well as how the information are loaded during runtime. 

When `emiss_read_legacy_mode` is set to `.TRUE.`, the `_chemistry` file must contain the vertical dimension `z`, and all emission data must also contain storage for the same magntiude of vertical direction representing the initial vertical cell layers. In addition, emission data can only be defined on an hourly basis from the start of the model run, and all emission data will be loaded in memory during initialization.

When `emiss_read_legacy_mode` is set to `.TRUE.`, the `z` dimension will be ignored, and emission data is only defined for the first vertical level. Emission data can be defined at arbitrary time intervals, and are only loaded into program memory at specificed time defined by the time stamp, replacing any previously loaded emission data.

The default value of `emiss_read_legacy_mode` is set to `.FALSE.`.

## Temporal Emission Profiles

Temporal profiles of emission are stored using the following format:

`YYYY-MM-DD HH:mm:ss +ZZ`

where `YYYY-MM-DD` represents the date in full numeric format, `HH:mm:ss` represents the local time in 24-hour format, `ZZ` is the time zone.  It is generally recommended that the profiles are specified in coordinated universal time (UTC) to maximize reusability.

The user can specify profile with arbitrary start and end times, as well as intervals, in chornological order. Emissions data defined at the time closest but not later than the current model time will be used as input. Should the model start (as defined in the `origin_date_time` option in the `_p3d` file) prior to the earliest defined emission data, the earliest emission data will be used.

### Dimensions

The definition of the LOD 2 emissions variables are predicated on the following mandatory dimensions:

| Dimension name | Value      | Description                           |
|----------------|------------|---------------------------------------|
| `time`        | at least 1 | Designated `UNLIMITED` dimension      |
| `field_length` | 64         | Fixed length of all string variables  |
| `nspecies`     | at least 1 | Number of emissions species           |
| `x`        | at least 1 | Number of cells in the west-east direction |
| `y`        | at least 1 | Number of cells in the south-north direction |
| `z`        | at least 1 | Number of cells in the vertical direction (legacy mode only) |
Note that the `field_length` dimension must be set to 64. Other dimensions must be greater than zero. The vertical (`z`) dimension is only required when the `_p3d` option `emiss_read_legacy_mode` is set to `.TRUE.`

### Variables

The mandatory variables for storage of emissions can then be defined using the above dimensions:
 
| Variable name | NetCDF data type | Dimension(s) | Description |
|---|---|---|---|
| `timestamp` | `NC_CHAR` | (`time`, `field_length`) | Individual time stamps for each set of emission data |
| `emission_name` | `NC_CHAR` | (`nspecies`, `field_length`) | Names of individual chemical species |
| `time`  | `NC_INT` | `time` | Emission time index (legacy mode only) |
| `emission_index`  | `NC_FLOAT` | `nspecies` | Emission species index (legacy mode only) |
| `emission_values` | `NC_FLOAT` | (`time`, `z`, `y`, `x`, `nspecies`) | Emission data for pollutant species (`z` dimension only required for legacy mode)  |

Note that the variables `time` and `emission_index` are only needed when the `_p3d` option `emiss_read_legacy_mode` is set to `.TRUE.`, as with the `z` dimension in the variable `emission_values`. The unit of `emission_values` can be set in its attribute. The default is `g/m2/hour` for all species.

### Example

The following is the header of an example netCDF emission input file (`_chemistry` file) under LOD 2 for non-legacy mode:
```
dimensions:
        time = UNLIMITED ; // (3 currently)
        y = 36 ;
        x = 36 ;
        nspecies = 12 ;
        field_length = 64 ;
variables:
        int nspecies(nspecies) ;
        char emission_name(nspecies, field_len) ;
        float emission_index(nspecies) ;
                emission_index:_FillValue = -9999.9f ;
        char timestamp(time, field_len) ;
        int time(time) ;
        float emission_values(time, y, x, nspecies) ;
                emission_values:_FillValue = -9999.9f ;
                emission_values:units = "g/m2/hour" ;
```
