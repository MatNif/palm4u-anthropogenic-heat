#!/usr/bin/env python3
from jinja2 import Template
import argparse
import numpy as np
import os
import xarray as xr

parser = argparse.ArgumentParser(description='Create a xdmf description file for Palm surface output')
parser.add_argument('infile', type=str, help='input file, e.g. <exp>_surf.001.nc')
parser.add_argument('dx', type=float, help='horizontal grid spacing')
parser.add_argument('dy', type=float, help='horizontal grid spacing')
parser.add_argument('dz', type=float, help='vertical grid spacing')
parser.add_argument('-xmf', type=str, help='output file, default is the inputfile with suffix .xmf')
parser.add_argument('-mesh', type=str, help='output file to store the extra mesh data, default is the inputfile with suffix .xmf.mesh')
parser.add_argument('-abs', '--absolute', action='store_true', default=False, help='if the xmf should reference absolute paths, default uses relative paths')
parser.add_argument('-r', '--no-round-time', action='store_true', default=False, help='disable rounding time values according to mean timesteps (default: on)')
parser.add_argument('-v', '--verbose', action='store_true', default=False, help='verbose output')

args = parser.parse_args()

fname_inp = args.infile

if args.xmf is None:
    fname_xmf = os.path.splitext(fname_inp)[0]+'.xmf'
else:
    fname_xmf = args.xmf

if args.mesh is None:
    fname_mesh = fname_xmf+'.mesh'
else:
    fname_mesh = args.mesh

if args.absolute:
    fname_inp = os.path.abspath(fname_inp)
    fname_xmf = os.path.abspath(fname_xmf)
    fname_mesh = os.path.abspath(fname_mesh)

dx, dy, dz = (args.dx, args.dy, args.dz)

rounding_time_vals = not args.no_round_time

verbose = args.verbose

del args

if verbose:
    print(f"Reading from {fname_inp}")
    print(f"Writing to   {fname_xmf} / {fname_mesh}")

D = xr.open_dataset(fname_inp, decode_times=False)

Nfaces = len(D.s)
Ntime = len(D.time)

# time in seconds, used for the time axis
time_sec = ((D.time.data - np.min(D.time.data)) * 1e-9).astype(np.float32)
if rounding_time_vals:
    if len(time_sec) > 1:
        dt = np.mean(time_sec[1:] - time_sec[:-1])
        decimals = np.int(np.round(-np.log10(dt/100.))) # round times to match the avg timestep plus one digit, never more than to the second
        if decimals > 0:
            new = np.round(time_sec, decimals=decimals)
        else:
            new = np.round(time_sec).astype(np.int32)
        if verbose:
            print(f"Rounding time steps:")
            print(f"  from: {time_sec}")
            print(f"  to  : {new}")
        time_sec = new


# first compute all vertex coords of each face
verts = np.empty((Nfaces,4,3))

# preload all geometric vars for better access
# xarray seems not super efficient caching those
xs = D['xs'].load().data[:]
ys = D['ys'].load().data[:]
zs = D['zs'].load().data[:]
zen = D['zenith'].load().data[:]
azi = D['azimuth'].load().data[:]

dx2, dy2, dz2 = [ .5*_ for _ in (dx, dy, dz)]

for i in range(Nfaces):
    center = np.array([xs[i], ys[i], zs[i]])

    if zen[i] == 0:
        verts[i,0,:] = center + np.array([-dx2, -dy2, 0])
        verts[i,1,:] = center + np.array([+dx2, -dy2, 0])
        verts[i,2,:] = center + np.array([+dx2, +dy2, 0])
        verts[i,3,:] = center + np.array([-dx2, +dy2, 0])
    elif zen[i] == 90:
        if (azi[i] == 0) or (azi[i] == 180) :
            verts[i,0,:] = center + np.array([-dx2, 0, -dz2])
            verts[i,1,:] = center + np.array([+dx2, 0, -dz2])
            verts[i,2,:] = center + np.array([+dx2, 0, +dz2])
            verts[i,3,:] = center + np.array([-dx2, 0, +dz2])
        elif (azi[i] == 90) or (azi[i] == 270) :
            verts[i,0,:] = center + np.array([0, -dy2, -dz2])
            verts[i,1,:] = center + np.array([0, +dy2, -dz2])
            verts[i,2,:] = center + np.array([0, +dy2, +dz2])
            verts[i,3,:] = center + np.array([0, -dy2, +dz2])
        else:
            raise ValueError(f"Expected azimuth to be one of cardinal directions but face {i} has {d.azimuth.data}")
    else:
        raise ValueError(f"Expected zenith to be 0 or 90 degrees but face {i} has {d.zenith}")


# vertex indices for each face
verts_of_faces = np.arange(Nfaces*4).reshape((Nfaces,4))

# write vertex mesh data
G = xr.Dataset({
    'verts_of_faces': xr.DataArray(verts_of_faces.astype(np.int32), dims=('N','vof')),
    'verts'         : xr.DataArray(verts.astype(np.float64), dims=('N','vof', 'c')),
    },
    attrs={
        'description': 'This file was generated with a script to describe palm surface variables. This file holds the vertex coordinate information for the faces and is referenced in a corresponding xdmf file',
        'dx': dx,
        'dy': dy,
        'dz': dz,
        })
G.to_netcdf(fname_mesh)

# xmf template for timeseries variables
xmf = '''<?xml version="1.0" encoding="utf-8"?>
 <!DOCTYPE Xdmf SYSTEM "Xdmf.dtd">
 <Xdmf Version="3.0">
 <Domain>
     <Grid CollectionType="Temporal" Name="BuildingsQuads" GridType="Collection">
     {% for time in time_sec %}
         {% set timeloop = loop %}
         <Grid Name="Quads0">
             <Time Value="{{time}}"/>
             <Topology TopologyType="Quadrilateral" NumberOfElements="{{Nfaces}}">  <DataItem ItemType="Uniform" Format="HDF" NumberType="Int" Precision="4" Dimensions="{{Nfaces}} 4">{{fname_mesh}}:/verts_of_faces</DataItem></Topology>
             <Geometry GeometryType="XYZ">  <DataItem ItemType="Uniform" Format="HDF" NumberType="Float" Precision="8" Dimensions="{{Nfaces}} 4 3">{{fname_mesh}}:/verts</DataItem></Geometry>

{% for var in variables %}
<Attribute Center="Cell" Name="{{var.name}}">
    <DataItem Dimensions="{{Nfaces}}" ItemType="HyperSlab">
        <DataItem Dimensions="2 3" Format="XML"> {{timeloop.index0}} 0    1 1    1 {{Nfaces}}</DataItem>
        <DataItem ItemType="Uniform" Format="HDF" NumberType="Float" Precision="{{var.dtype.itemsize}}" Dimensions="{{Ntime}} {{Nfaces}}">{{fname_inp}}:/{{var.name}}</DataItem>
    </DataItem>
</Attribute>
{%endfor%}

         </Grid>
     {%endfor%}
     </Grid>
 </Domain>
</Xdmf>
'''
t=Template(xmf)

xmfout = t.render({
    'variables': [ v for k,v in D.items() if 'time' in v.dims],
    'time_sec': time_sec,
    'Ntime': Ntime,
    'Nfaces': Nfaces,
    'fname_mesh': fname_mesh,
    'fname_inp': fname_inp,
    })

with open(fname_xmf,'w') as fh:
    fh.writelines(xmfout)
