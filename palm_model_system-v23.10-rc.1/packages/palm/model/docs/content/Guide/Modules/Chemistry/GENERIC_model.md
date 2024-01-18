# Generic Emission Mode

The generic emission mode is available for the user to accommodate emission modes that have not been implemented, or in cases where emissions from multiple modes have been assembled into a single data set. It can be activated by setting the value of the `_p3d` option `emis_generic` to `.TRUE.`.  The generic mode is only available in LOD 2.

Input data for the generic emission mode are contained in a file called `[model]_emis_generic` in the model `INPUT` directory, where `[model]` refers to the name of the model. The file is to be stored using netCDF APIs following the format for [specific emission mode input data in LOD 2](./EMISSIONS_LOD2_spec.md).