# Traffic Emission Mode

The traffic emission mode is available for the user in cases where emissions have been prepared prepared externally. It can be activated by setting the value of the `_p3d` option `emis_traffic` to `.TRUE.`.  The traffic mode is currently available in LOD 2.

Input data for the traffic emission mode are contained in a file called `[model]_emis_traffic` in the model `INPUT` directory, where `[model]` refers to the name of the model. The file is to be stored using netCDF APIs following the format for [specific emission mode input data in LOD 2](./EMISSIONS_LOD2_spec.md).
