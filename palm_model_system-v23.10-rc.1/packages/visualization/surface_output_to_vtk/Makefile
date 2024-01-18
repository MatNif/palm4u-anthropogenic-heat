#------------------------------------------------------------------------------#
# This file is part of the PALM model system.
#
# PALM is free software: you can redistribute it and/or modify it under the
# terms of the GNU General Public License as published by the Free Software
# Foundation, either version 3 of the License, or (at your option) any later
# version.
#
# PALM is distributed in the hope that it will be useful, but WITHOUT ANY
# WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR
# A PARTICULAR PURPOSE.  See the GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License along with
# PALM. If not, see <http://www.gnu.org/licenses/>.
#
# Copyright 2017-2021 Leibniz Universitaet Hannover
#------------------------------------------------------------------------------#
#
# Authors:
# --------
# @author Matthias Suehring
#
# Description:
# ------------
# This makefile may be used to exclusively compile surface_output_to_vtk
#------------------------------------------------------------------------------#
PROJECT = surface_output_to_vtk
PROJECT_PATH = .
BIN_PATH  = $(PROJECT_PATH)/bin
SRC_PATH  = $(PROJECT_PATH)/src

MODULES = $(SRC_PATH)/posix_interface_mod.mod
SOURCES = $(MODULES:%.mod=%.f90) $(SRC_PATH)/$(PROJECT).f90
OBJECTS = $(SOURCES:%.f90=%.o)

FC      = gfortran
WRNGS   = -Wall -Wextra -Wsurprising -Waliasing -Wcharacter-truncation \
		  -Wline-truncation -fbacktrace -fcheck=all -pedantic \
		  -ffpe-trap=invalid,zero,underflow,overflow
FFLAGS  = -fdefault-real-8 -O3 $(WRNGS) -o $(PROJECT)
INCLUDE =
LIBRARY =


$(PROJECT): $(OBJECTS)
	@echo ""
	$(FC) $(FFLAGS) $(OBJECTS) -I$(SRC_PATH) $(LIBRARY)

$(SRC_PATH)/%.o: $(SRC_PATH)/%.f90
	@echo ""
	$(FC) $(FFLAGS) -c $< -o $@ -J$(SRC_PATH) $(INCLUDE)

$(SRC_PATH)/%.mod: $(SRC_PATH)/%.f90
	@echo ""
	$(FC) $(FFLAGS) -c $< -o $(@:%.mod=%.o) -J$(SRC_PATH) $(INCLUDE)

all:

clean:
	rm -f $(SRC_PATH)/*.o $(SRC_PATH)/*.mod

posix_interface_mod.o:
surface_output_to_vtk.o: posix_interface_mod.o
