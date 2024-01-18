# Point Source Emission Mode

Pollutant sources found in, for instance, the European Pollutant Release and Transfer Register (E-PRTR) and the Gridding Emission Tool for ArcGIS (GRETA) can be specified using the point source emission mode.  It can be activated by setting the value of the `_p3d` option `emis_pt_source` to `.TRUE.`.  The point source emission mode is only available in LOD 0.

For each point source, annual aggreate emissions for a user-defined species can be specified, which will be converted into a temporally uniform volumetric source distributed over a user-defined vertical distance. The following options can be set, in the `_p3d` file, all of which carry the prefix `emis_pt_source_`:

| `_p3d` option | Default value | Description |
|---|---|---|
| `emis_pt_source_species_names` | None      | Names of chemical species |
| `emis_pt_source_locations_ijk` | None      | Base cell coordinates (i,j,k) of each point source |
| `emis_pt_source_annual_values` | 0         | Annual aggregate of emission for all species at each point source |
| `emis_pt_source_leap_year`     | `.FALSE.` | Leap year correction for annual aggregate distribution |
| `emis_pt_source_k_spread`      | 3         | Number of vertical layers over which the point source is to be distributed |
| `emis_pt_source_k_weights`     | None      | Vertical distribution factors for all point sources |

The name and order of the chemical species specified in `emis_pt_source_speccies_names` can be arbitrary, and dictates the species order in the annual aggregate values defined in `emis_pt_source_source_annual_values`. Species not found in the model active chemical mechanism will be ignored. All point sources emissions must be defined for all species at all source locations. Up to 99 species and 199 point source locations can be specified at any given model. Units for annual aggregate emissions are in mol/year for gas phase species, and kg/year for all PM species. If the model run takes place in a leap year, the option `emis_pt_source_leap_year` must be set to `.TRUE.`.

The point sources can be distributed over a vertical distance of up to five layers using the option `emis_pt_source_k_spread`, with the bottommost location provided by each cell (i,j,k) location defined in `emis_pt_source_locations_ijk`. The manner of vertical distribution can be specified using the option `emis_pt_source_k_weights`. The factors will be automatically normalized so that the sum of all weights up to the value specified by `emis_pt_source_k_spread` will always be unity. When unspecified, the point source will be uniformly distributed across the specified vertical layers.   The same distribution will be applied to all point sources.