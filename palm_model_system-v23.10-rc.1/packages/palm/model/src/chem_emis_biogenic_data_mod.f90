!> @file chem_emis_biogenic_data_mod.f90
!--------------------------------------------------------------------------------------------------!
! This file is part of the PALM model system.
!
! PALM is free software: you can redistribute it and/or modify it under the terms of the GNU General
! Public License as published by the Free Software Foundation, either version 3 of the License, or
! (at your option) any later version.
!
! PALM is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the
! implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General
! Public License for more details.
!
! You should have received a copy of the GNU General Public License along with PALM. If not, see
! <http://www.gnu.org/licenses/>.
!
! Copyright 1997-2021 Leibniz Universitaet Hannover
! Copyright 2017-2021 Karlsruhe Institute of Technology
!--------------------------------------------------------------------------------------------------!
!
! Description:
! ------------
!> Module for  biogenic emission data of VOCs, trees, and plant functional types
!--------------------------------------------------------------------------------------------------!
 MODULE chem_emis_biogenic_data_mod


    USE kinds

    IMPLICIT NONE
!
!-- variables declared before CHAR to define array size/length that are required for CHAR, INTEGER and REAL arrays.
    INTEGER(iwp), PARAMETER ::  n_dcds_blf     = 64  !< max number of deiduous broadleaf trees
    INTEGER(iwp), PARAMETER ::  n_dsoil_layers = 08  !< default number of soil layers
    INTEGER(iwp), PARAMETER ::  n_ebio_class   = 05  !< maximum number of bvoc classes
    INTEGER(iwp), PARAMETER ::  n_ebio_spc     = 30  !< maximum number of bvoc spcs
    INTEGER(iwp), PARAMETER ::  n_evgn_blf     = 04  !< max number of evergreen broadleaf trees
    INTEGER(iwp), PARAMETER ::  n_evgn_nlf     = 16  !< maximum number of evervreen needleleaf trees
    INTEGER(iwp), PARAMETER ::  n_pft          = 18  !< maximum number of vegetation types (ECMWF-IFS)
    INTEGER(iwp), PARAMETER ::  n_root_layers  = 03  !< total number of root layers
    INTEGER(iwp), PARAMETER ::  n_soil_types   = 06  !< total number of soil types
    INTEGER(iwp), PARAMETER ::  n_tree_types   = 86  !< maximum number of tree types
!
    CHARACTER(LEN=4), DIMENSION(1:n_ebio_class) ::  ltab_ebio_class_nam = (/                       &  !< names of bvoc classes
                                                                              'ISOP',              &
                                                                              'MTRP',              &
                                                                              'SQTP',              &
                                                                              'XVOC',              &
                                                                              'OVOC'               &
                                                                              /)
!
    CHARACTER(LEN=6), DIMENSION(1:n_ebio_spc) ::  ltab_ebio_nam = (/                               &  !< list of bvoc spcs in BEM
                                                                    'ISOPRN',                      &
                                                                    'MYRCEN',                      &
                                                                    'SABINE',                      &
                                                                    'LIMONE',                      &
                                                                    '3CAREN',                      &
                                                                    'OCIMEN',                      &
                                                                    'APINEN',                      &
                                                                    'BPINEN',                      &
                                                                    'BCARYO',                      &
                                                                    'ACTLDE',                      &
                                                                    'ETHANL',                      &
                                                                    'FMALDE',                      &
                                                                    'MTHANL',                      &
                                                                    'ACETON',                      &
                                                                    'FMCACI',                      &
                                                                    'ACTACI',                      &
                                                                    '232MBO',                      &
                                                                    'MTHANE',                      &
                                                                    'ETHENE',                      &
                                                                    'HYDCYN',                      &
                                                                    'TOLUEN',                      &
                                                                    'MTHLBR',                      &
                                                                    'MTHLCL',                      &
                                                                    'MTHLIO',                      &
                                                                    'DMTHLS',                      &
                                                                    'ETHANE',                      &
                                                                    'PROPAN',                      &
                                                                    'PROPEN',                      &
                                                                    'BUTANE',                      &
                                                                    'BNZLDE'                       &
                                                                   /)
!
    CHARACTER(26), DIMENSION(1:n_pft), PARAMETER ::  ltab_pft_nam = (/                             &  !< ECMWF-IFC vege types/PFT
                                                                   'bare soil                 ',   &
                                                                   'crops, mixed farming      ',   &
                                                                   'short grass               ',   &
                                                                   'evergreen needleleaf trees',   &
                                                                   'deciduous needleleaf trees',   &
                                                                   'evergreen broadleaf trees ',   &
                                                                   'deciduous broadleaf trees ',   &
                                                                   'tall grass                ',   &
                                                                   'desert                    ',   &
                                                                   'tundra                    ',   &
                                                                   'irrigated crops           ',   &
                                                                   'semidesert                ',   &
                                                                   'ice caps and glaciers     ',   &
                                                                   'bogs and marshes          ',   &
                                                                   'evergreen shrubs          ',   &
                                                                   'deciduous shrubs          ',   &
                                                                   'mixed forest/woodland     ',   &
                                                                   'interrupted forest        '    &
                                                                  /)
!
    CHARACTER(LEN=18), DIMENSION(1:n_tree_types) ::  ltab_tree_nam = (/                            &  !< single tree spcs BERLIN.
                                                                       'Abies             ',       &
                                                                       'Acer              ',       &
                                                                       'Aesculus          ',       &
                                                                       'Ailanthus         ',       &
                                                                       'Alnus             ',       &
                                                                       'Amelanchier       ',       &
                                                                       'Betula            ',       &
                                                                       'Buxus             ',       &
                                                                       'Calocedrus        ',       &
                                                                       'Caragana          ',       &
                                                                       'Carpinus          ',       &
                                                                       'Carya             ',       &
                                                                       'Castanea          ',       &
                                                                       'Catalpa           ',       &
                                                                       'Cedrus            ',       &
                                                                       'Celtis            ',       &
                                                                       'Cercidiphyllum    ',       &
                                                                       'Cercis            ',       &
                                                                       'Chamaecyparis     ',       &
                                                                       'Cladrastis        ',       &
                                                                       'Cornus            ',       &
                                                                       'Corylus           ',       &
                                                                       'Cotinus           ',       &
                                                                       'Crataegus         ',       &
                                                                       'Cryptomeria       ',       &
                                                                       'Cupressocyparis   ',       &
                                                                       'Cupressus         ',       &
                                                                       'Cydonia           ',       &
                                                                       'Davidia           ',       &
                                                                       'Elaeagnus         ',       &
                                                                       'Euodia            ',       &
                                                                       'Euonymus          ',       &
                                                                       'Fagus             ',       &
                                                                       'Fraxinus          ',       &
                                                                       'Ginkgo            ',       &
                                                                       'Gleditsia         ',       &
                                                                       'Gymnocladus       ',       &
                                                                       'Hippophae         ',       &
                                                                       'Ilex              ',       &
                                                                       'Juglans           ',       &
                                                                       'Juniperus         ',       &
                                                                       'Koelreuteria      ',       &
                                                                       'Laburnum          ',       &
                                                                       'Larix             ',       &
                                                                       'Ligustrum         ',       &
                                                                       'Liquidambar       ',       &
                                                                       'Liriodendron      ',       &
                                                                       'Lonicera          ',       &
                                                                       'Magnolia          ',       &
                                                                       'Malus             ',       &
                                                                       'Metasequoia       ',       &
                                                                       'Morus             ',       &
                                                                       'Ostrya            ',       &
                                                                       'Parrotia          ',       &
                                                                       'Paulownia         ',       &
                                                                       'Phellodendron     ',       &
                                                                       'Picea             ',       &
                                                                       'Pinus             ',       &
                                                                       'Platanus          ',       &
                                                                       'Populus           ',       &
                                                                       'Prunus            ',       &
                                                                       'Pseudotsuga       ',       &
                                                                       'Ptelea            ',       &
                                                                       'Pterocaria        ',       &
                                                                       'Pterocarya        ',       &
                                                                       'Pyrus             ',       &
                                                                       'Quercus           ',       &
                                                                       'Rhamnus           ',       &
                                                                       'Rhus              ',       &
                                                                       'Robinia           ',       &
                                                                       'Salix             ',       &
                                                                       'Sambucus          ',       &
                                                                       'Sasa              ',       &
                                                                       'Sequoiadendron    ',       &
                                                                       'Sophora           ',       &
                                                                       'Sorbus            ',       &
                                                                       'Syringa           ',       &
                                                                       'Tamarix           ',       &
                                                                       'Taxodium          ',       &
                                                                       'Taxus             ',       &
                                                                       'Thuja             ',       &
                                                                       'Tilia             ',       &
                                                                       'Tsuga             ',       &
                                                                       'Ulmus             ',       &
                                                                       'Zelkova           ',       &
                                                                       'Zenobia           '        &
                                                                       /)
!
!-- The named indicies below follow numerical order instead of alphabetical order
!-- becuase these are used in the same order in arrays in this data file
!-- and the biogenic module.
!-- bvoc class indecies/IDs
    INTEGER(iwp), PARAMETER ::  class_isop = 1  !< index for isoprene
    INTEGER(iwp), PARAMETER ::  class_mtrp = 2  !< index for monoterpenes
    INTEGER(iwp), PARAMETER ::  class_sqtp = 3  !< index for sesquiterpenes
    INTEGER(iwp), PARAMETER ::  class_xvoc = 4  !< index for oxygenated vocs
    INTEGER(iwp), PARAMETER ::  class_ovoc = 5  !< index for other vocs
!
!-- 30 bvoc spcs IDs
!-- class isoprene
    INTEGER(iwp), PARAMETER ::  isop  = 1   !< isoprene          c5h8 /
!
!-- Class monoterpenes MRTP
    INTEGER(iwp), PARAMETER ::  mtrp1 = 2   !< myrcene          c10h16
    INTEGER(iwp), PARAMETER ::  mtrp2 = 3   !< sabinene         c10h16
    INTEGER(iwp), PARAMETER ::  mtrp3 = 4   !< limonene         c10h16
    INTEGER(iwp), PARAMETER ::  mtrp4 = 5   !< 3carene          c10h16
    INTEGER(iwp), PARAMETER ::  mtrp5 = 6   !< t-b-ocimene      c10h16
    INTEGER(iwp), PARAMETER ::  mtrp6 = 7   !< a-pinene         c10h16
    INTEGER(iwp), PARAMETER ::  mtrp7 = 8   !< b-pinene         c10h16
!
!-- class sesquiterpenes
    INTEGER(iwp), PARAMETER ::  sqtp1  = 9   !< b-caryophyllene  c15h24
!
!-- classs oxygenated BVOCs
    INTEGER(iwp), PARAMETER ::  xvoc1 = 10  !< acetaldehyde      c2h4o
    INTEGER(iwp), PARAMETER ::  xvoc2 = 11  !< ethanol           c2h5oh
    INTEGER(iwp), PARAMETER ::  xvoc3 = 12  !< HCHO'formaldehyde ch2o
    INTEGER(iwp), PARAMETER ::  xvoc4 = 13  !< methanol          ch3oh
    INTEGER(iwp), PARAMETER ::  xvoc5 = 14  !< acetone           c3h6o
    INTEGER(iwp), PARAMETER ::  xvoc6 = 15  !< farmic_acid       ch2o2
    INTEGER(iwp), PARAMETER ::  xvoc7 = 16  !< acetic_acid       ch3cooh
!
!-- class others BVOC<s
    INTEGER(iwp), PARAMETER ::  ovoc1  = 17  !< 232-MBO           c5h10o
    INTEGER(iwp), PARAMETER ::  ovoc2  = 18  !< methane           ch4
    INTEGER(iwp), PARAMETER ::  ovoc3  = 19  !< ethene            c2h4
    INTEGER(iwp), PARAMETER ::  ovoc4  = 20  !< hydrogen_cynide   hcn
    INTEGER(iwp), PARAMETER ::  ovoc5  = 21  !< toluene           c7h8
    INTEGER(iwp), PARAMETER ::  ovoc6  = 22  !< methyl_bromide    ch3br
    INTEGER(iwp), PARAMETER ::  ovoc7  = 23  !< methyl_chloride   ch3cl
    INTEGER(iwp), PARAMETER ::  ovoc8  = 24  !< methyl_iodide     ch3i
    INTEGER(iwp), PARAMETER ::  ovoc9  = 25  !< dimethyl_sulfide  c2h6s
    INTEGER(iwp), PARAMETER ::  ovoc10 = 26  !< ethane            c2h6
    INTEGER(iwp), PARAMETER ::  ovoc11 = 27  !< propane           c3h8
    INTEGER(iwp), PARAMETER ::  ovoc12 = 28  !< propene           c3h6
    INTEGER(iwp), PARAMETER ::  ovoc13 = 29  !< butane            c4h10
    INTEGER(iwp), PARAMETER ::  ovoc14 = 30  !< benzaldehyde      c7h6o
!
!-- Soil types from ECMWF-IFS IDs
    INTEGER(iwp), PARAMETER ::  soil_coarse  = 1  !< soil type coarse
    INTEGER(iwp), PARAMETER ::  soil_medium  = 2  !< soil type medium
    INTEGER(iwp), PARAMETER ::  soil_mfine   = 3  !< soil type med fine
    INTEGER(iwp), PARAMETER ::  soil_fine    = 4  !< soil type fine
    INTEGER(iwp), PARAMETER ::  soil_vfine   = 5  !< soil type very fine
    INTEGER(iwp), PARAMETER ::  soil_organic = 6  !< soil type organic
!
!-- Plant functional types (vegetation types) ECMWF-IFS classification.
    INTEGER(iwp), PARAMETER ::  bare_soil     = 1   !< vegetation type bare soil
    INTEGER(iwp), PARAMETER ::  crops         = 2   !< vegetation type crops
    INTEGER(iwp), PARAMETER ::  short_grass   = 3   !< vegetation type short grass
    INTEGER(iwp), PARAMETER ::  evgn_nlf_tree = 4   !< vegetation type evergreen needleleaf trees
    INTEGER(iwp), PARAMETER ::  dcds_nlf_tree = 5   !< vegetation type deciduous needleleaf trees
    INTEGER(iwp), PARAMETER ::  evgn_blf_tree = 6   !< vegetation type evergreen broadleaf trees
    INTEGER(iwp), PARAMETER ::  dcds_blf_tree = 7   !< vegetation type deciduous broadleaf trees
    INTEGER(iwp), PARAMETER ::  tall_grass    = 8   !< vegetation type tall grass
    INTEGER(iwp), PARAMETER ::  desert        = 9   !< desert
    INTEGER(iwp), PARAMETER ::  tundra        = 10  !< turdra
    INTEGER(iwp), PARAMETER ::  irrig_crops   = 11  !< vegetation type irrigated crops
    INTEGER(iwp), PARAMETER ::  semidesert    = 12  !< semi-desert
    INTEGER(iwp), PARAMETER ::  ice_glaciers  = 13  !< ice caps and glaciers
    INTEGER(iwp), PARAMETER ::  bogs_marshes  = 14  !< vegetation type bogs and marshes
    INTEGER(iwp), PARAMETER ::  evgn_shrubs   = 15  !< vegetation type evergreen shrubs
    INTEGER(iwp), PARAMETER ::  dcds_shrubs   = 16  !< vegetation type deciduous shrubs
    INTEGER(iwp), PARAMETER ::  mixed_forest  = 17  !< vegetation type mixed forest
    INTEGER(iwp), PARAMETER ::  intptd_forest = 18  !< vegetation type interrupted forest
!
!-- List of single trees types in Berlin region
    INTEGER(iwp), PARAMETER ::  abies           = 1   !< pft-evgn_nlf_tree
    INTEGER(iwp), PARAMETER ::  acer            = 2   !< pft-dcds_blf_tree
    INTEGER(iwp), PARAMETER ::  aesculus        = 3   !< pft-dcds_blf_tree
    INTEGER(iwp), PARAMETER ::  ailanthus       = 4   !< pft-dcds_blf_tree
    INTEGER(iwp), PARAMETER ::  alnus           = 5   !< pft-dcds_blf_tree
    INTEGER(iwp), PARAMETER ::  amelanchier     = 6   !< pft-dcds_blf_tree
    INTEGER(iwp), PARAMETER ::  betula          = 7   !< pft-dcds_blf_tree
    INTEGER(iwp), PARAMETER ::  buxus           = 8   !< pft-evgn_blf_tree
    INTEGER(iwp), PARAMETER ::  calocedrus      = 9   !< pft-evgn_nlf_tree
    INTEGER(iwp), PARAMETER ::  caragana        = 10  !< pft-dcds_blf_tree
    INTEGER(iwp), PARAMETER ::  carpinus        = 11  !< pft-dcds_blf_tree
    INTEGER(iwp), PARAMETER ::  carya           = 12  !< pft-dcds_blf_tree
    INTEGER(iwp), PARAMETER ::  castanea        = 13  !< pft-dcds_blf_tree
    INTEGER(iwp), PARAMETER ::  catalpa         = 14  !< pft-dcds_blf_tree
    INTEGER(iwp), PARAMETER ::  cedrus          = 15  !< pft-evgn_nlf_tree
    INTEGER(iwp), PARAMETER ::  celtis          = 16  !< pft-dcds_blf_tree
    INTEGER(iwp), PARAMETER ::  cercidiphyllum  = 17  !< pft-dcds_blf_tree
    INTEGER(iwp), PARAMETER ::  cercis          = 18  !< pft-dcds_blf_tree
    INTEGER(iwp), PARAMETER ::  chamaecyparis   = 19  !< pft-evgn_nlf_tree
    INTEGER(iwp), PARAMETER ::  cladrastis      = 20  !< pft-dcds_blf_tree
    INTEGER(iwp), PARAMETER ::  cornus          = 21  !< pft-dcds_blf_tree
    INTEGER(iwp), PARAMETER ::  corylus         = 22  !< pft-dcds_blf_tree
    INTEGER(iwp), PARAMETER ::  cotinus         = 23  !< pft-dcds_blf_tree
    INTEGER(iwp), PARAMETER ::  crataegus       = 24  !< pft-dcds_blf_tree
    INTEGER(iwp), PARAMETER ::  cryptomeria     = 25  !< pft-evgn_nlf_tree
    INTEGER(iwp), PARAMETER ::  cupressocyparis = 26  !< pft-evgn_nlf_tree
    INTEGER(iwp), PARAMETER ::  cupressus       = 27  !< pft-evgn_nlf_tree
    INTEGER(iwp), PARAMETER ::  cydonia         = 28  !< pft-dcds_blf_tree
    INTEGER(iwp), PARAMETER ::  davidia         = 29  !< pft-dcds_blf_tree
    INTEGER(iwp), PARAMETER ::  elaeagnus       = 30  !< pft-dcds_blf_tree
    INTEGER(iwp), PARAMETER ::  euodia          = 31  !< pft-dcds_blf_tree
    INTEGER(iwp), PARAMETER ::  euonymus        = 32  !< pft-dcds_blf_tree
    INTEGER(iwp), PARAMETER ::  fagus           = 33  !< pft-dcds_blf_tree
    INTEGER(iwp), PARAMETER ::  fraxinus        = 34  !< pft-dcds_blf_tree
    INTEGER(iwp), PARAMETER ::  ginkgo          = 35  !< pft-dcds_blf_tree
    INTEGER(iwp), PARAMETER ::  gleditsia       = 36  !< pft-dcds_blf_tree
    INTEGER(iwp), PARAMETER ::  gymnocladus     = 37  !< pft-dcds_blf_tree
    INTEGER(iwp), PARAMETER ::  hippophae       = 38  !< pft-dcds_blf_tree
    INTEGER(iwp), PARAMETER ::  ilex            = 39  !< pft-evgn_blf_tree
    INTEGER(iwp), PARAMETER ::  juglans         = 40  !< pft-dcds_blf_tree
    INTEGER(iwp), PARAMETER ::  juniperus       = 41  !< pft-evgn_nlf_tree
    INTEGER(iwp), PARAMETER ::  koelreuteria    = 42  !< pft-dcds_blf_tree
    INTEGER(iwp), PARAMETER ::  laburnum        = 43  !< pft-dcds_blf_tree
    INTEGER(iwp), PARAMETER ::  larix           = 44  !< pft-dcds_blf_tree
    INTEGER(iwp), PARAMETER ::  ligustrum       = 45  !< pft-dcds_blf_tree
    INTEGER(iwp), PARAMETER ::  liquidambar     = 46  !< pft-dcds_blf_tree
    INTEGER(iwp), PARAMETER ::  liriodendron    = 47  !< pft-dcds_blf_tree
    INTEGER(iwp), PARAMETER ::  lonicera        = 48  !< pft-dcds_blf_tree
    INTEGER(iwp), PARAMETER ::  magnolia        = 49  !< pft-dcds_blf_tree
    INTEGER(iwp), PARAMETER ::  malus           = 50  !< pft-dcds_blf_tree
    INTEGER(iwp), PARAMETER ::  metasequoia     = 51  !< pft-dcds_blf_tree
    INTEGER(iwp), PARAMETER ::  morus           = 52  !< pft-dcds_blf_tree
    INTEGER(iwp), PARAMETER ::  ostrya          = 53  !< pft-dcds_blf_tree
    INTEGER(iwp), PARAMETER ::  parrotia        = 54  !< pft-dcds_blf_tree
    INTEGER(iwp), PARAMETER ::  paulownia       = 55  !< pft-dcds_blf_tree
    INTEGER(iwp), PARAMETER ::  phellodendron   = 56  !< pft-dcds_blf_tree
    INTEGER(iwp), PARAMETER ::  picea           = 57  !< pft-evgn_nlf_tree
    INTEGER(iwp), PARAMETER ::  pinus           = 58  !< pft-evgn_nlf_tree
    INTEGER(iwp), PARAMETER ::  platanus        = 59  !< pft-dcds_blf_tree
    INTEGER(iwp), PARAMETER ::  populus         = 60  !< pft-dcds_blf_tree
    INTEGER(iwp), PARAMETER ::  prunus          = 61  !< pft-dcds_blf_tree
    INTEGER(iwp), PARAMETER ::  pseudotsuga     = 62  !< pft-evgn_nlf_tree
    INTEGER(iwp), PARAMETER ::  ptelea          = 63  !< pft-dcds_blf_tree
    INTEGER(iwp), PARAMETER ::  Pterocaria      = 64  !< pft-dcds_blf_tree
    INTEGER(iwp), PARAMETER ::  pterocarya      = 65  !< pft-dcds_blf_tree
    INTEGER(iwp), PARAMETER ::  pyrus           = 66  !< pft-dcds_blf_tree
    INTEGER(iwp), PARAMETER ::  quercus         = 67  !< pft-dcds_blf_tree
    INTEGER(iwp), PARAMETER ::  rhamnus         = 68  !< pft-dcds_blf_tree
    INTEGER(iwp), PARAMETER ::  rhus            = 69  !< pft-dcds_blf_tree
    INTEGER(iwp), PARAMETER ::  robinia         = 70  !< pft-dcds_blf_tree
    INTEGER(iwp), PARAMETER ::  salix           = 71  !< pft-dcds_blf_tree
    INTEGER(iwp), PARAMETER ::  sambucus        = 72  !< pft-dcds_blf_tree
    INTEGER(iwp), PARAMETER ::  sasa            = 73  !< pft-evgn_blf_tree
    INTEGER(iwp), PARAMETER ::  sequoiadendron  = 74  !< pft-evgn_nlf_tree
    INTEGER(iwp), PARAMETER ::  sophora         = 75  !< pft-dcds_blf_tree
    INTEGER(iwp), PARAMETER ::  sorbus          = 76  !< pft-dcds_blf_tree
    INTEGER(iwp), PARAMETER ::  syringa         = 77  !< pft-dcds_blf_tree
    INTEGER(iwp), PARAMETER ::  tamarix         = 78  !< pft-dcds_blf_tree
    INTEGER(iwp), PARAMETER ::  taxodium        = 79  !< pft-evgn_nlf_tree
    INTEGER(iwp), PARAMETER ::  taxus           = 80  !< pft-evgn_nlf_tree
    INTEGER(iwp), PARAMETER ::  thuja           = 81  !< pft-evgn_nlf_tree
    INTEGER(iwp), PARAMETER ::  tilia           = 82  !< pft-dcds_blf_tree
    INTEGER(iwp), PARAMETER ::  tsuga           = 83  !< pft-evgn_nlf_tree
    INTEGER(iwp), PARAMETER ::  ulmus           = 84  !< pft-dcds_blf_tree
    INTEGER(iwp), PARAMETER ::  zelkova         = 85  !< pft-dcds_blf_tree
    INTEGER(iwp), PARAMETER ::  zenobia         = 86  !< pft-evgn_blf_tree
!
!-- vegetation and bvoc data arrays
    INTEGER(iwp)  ::  ifill  !< dummy integer to fill array with series
!
    INTEGER(iwp), DIMENSION(1:13), PARAMETER ::  ltab_def_plant = (/                               &  !<  Default plant type (0).
                                                                    1.0_wp,                        &
                                                                    1.0_wp,                        &
                                                                    4.0_wp,                        &
                                                                    12.0_wp,                       &
                                                                    3.0_wp,                        &
                                                                    0.8_wp,                        &
                                                                    0.6_wp,                        &
                                                                    0.025_wp,                      &
                                                                    0.035_wp,                      &
                                                                    0.00035_wp,                    &
                                                                    0.00035_wp,                    &
                                                                    0.00035_wp,                    &
                                                                    0.00035_wp                     &
                                                                   /)
!
!-- data from external model, initialized wiht DATA stament
    INTEGER(iwp), DIMENSION(n_ebio_spc)   ::  ltab_ebio_map  !< mapping individual bvocs to bvoc classes
!
    INTEGER(iwp), DIMENSION (n_ebio_class) ::  ltab_ebio_class = (/                                &  !< 30 bvoc scps in 5 groups.
                                                                   1,                              &
                                                                   2,                              &
                                                                   3,                              &
                                                                   4,                              &
                                                                   5                               &
                                                                  /)
!
    INTEGER(iwp), DIMENSION(1:n_pft), PARAMETER ::  ltab_pft = (/                                  &   !< PFTs class ID
                                                                 -1,                               &
                                                                 -2,                               &
                                                                 -3,                               &
                                                                 -4,                               &
                                                                 -5,                               &
                                                                 -6,                               &
                                                                 -7,                               &
                                                                 -8,                               &
                                                                 -9,                               &
                                                                -10,                               &
                                                                -11,                               &
                                                                -12,                               &
                                                                -13,                               &
                                                                -14,                               &
                                                                -15,                               &
                                                                -16,                               &
                                                                -17,                               &
                                                                -18                                &
                                                               /)
!
!-- specific leaf mass(slm) of 11 PFTs in cm2 g-1
!-- unit coversion to lad (g m-3);   lad(g m-3) =lad(m2 m-3) * lad(g m-2);  lad(g m-2) = 10000/ slm(cm2 g-1) .
!--  specific leaf mass of pfts for unit conv from LAD m2/m3 to g/m3
    INTEGER(iwp), DIMENSION(1:n_pft), PARAMETER ::  ltab_pft_slm = (/                              &  !< specific leaf mass of pfts
                                                                    -9999.0_wp,                    &
                                                                    -9999.0_wp,                    &
                                                                    -9999.0_wp,                    &
                                                                    40.0_wp,                       &
                                                                    100.0_wp,                      &
                                                                    89.0_wp,                       &
                                                                    137.0_wp,                      &
                                                                    -9999.0_wp,                    &
                                                                    -9999.0_wp,                    &
                                                                    -9999.0_wp,                    &
                                                                    -9999.0_wp,                    &
                                                                    -9999.0_wp,                    &
                                                                    -9999.0_wp,                    &
                                                                    -9999.0_wp,                    &
                                                                    71.0_wp,                       &
                                                                    140.0_wp,                      &
                                                                    -9999.0_wp,                    &
                                                                    -9999.0_wp                     &
                                                                   /)
!
    INTEGER(iwp), DIMENSION (1:n_tree_types) ::  ltab_tree = (/ (ifill, ifill = 1, 86) /)  !< Berlin tree IDs
!
    INTEGER(iwp), DIMENSION(n_tree_types), PARAMETER ::  ltab_tree_map = (/                        &  !<  map trees to pft classes
                                                                           evgn_nlf_tree,          &
                                                                           dcds_blf_tree,          &
                                                                           dcds_blf_tree,          &
                                                                           dcds_blf_tree,          &
                                                                           dcds_blf_tree,          &
                                                                           dcds_blf_tree,          &
                                                                           dcds_blf_tree,          &
                                                                           evgn_blf_tree,          &
                                                                           evgn_nlf_tree,          &
                                                                           dcds_blf_tree,          &
                                                                           dcds_blf_tree,          &
                                                                           dcds_blf_tree,          &
                                                                           dcds_blf_tree,          &
                                                                           dcds_blf_tree,          &
                                                                           evgn_nlf_tree,          &
                                                                           dcds_blf_tree,          &
                                                                           dcds_blf_tree,          &
                                                                           dcds_blf_tree,          &
                                                                           evgn_nlf_tree,          &
                                                                           dcds_blf_tree,          &
                                                                           dcds_blf_tree,          &
                                                                           dcds_blf_tree,          &
                                                                           dcds_blf_tree,          &
                                                                           dcds_blf_tree,          &
                                                                           evgn_nlf_tree,          &
                                                                           evgn_nlf_tree,          &
                                                                           evgn_nlf_tree,          &
                                                                           dcds_blf_tree,          &
                                                                           dcds_blf_tree,          &
                                                                           dcds_blf_tree,          &
                                                                           dcds_blf_tree,          &
                                                                           dcds_blf_tree,          &
                                                                           dcds_blf_tree,          &
                                                                           dcds_blf_tree,          &
                                                                           dcds_blf_tree,          &
                                                                           dcds_blf_tree,          &
                                                                           dcds_blf_tree,          &
                                                                           dcds_blf_tree,          &
                                                                           evgn_blf_tree,          &
                                                                           dcds_blf_tree,          &
                                                                           evgn_nlf_tree,          &
                                                                           dcds_blf_tree,          &
                                                                           dcds_blf_tree,          &
                                                                           dcds_blf_tree,          &
                                                                           dcds_blf_tree,          &
                                                                           dcds_blf_tree,          &
                                                                           dcds_blf_tree,          &
                                                                           dcds_blf_tree,          &
                                                                           dcds_blf_tree,          &
                                                                           dcds_blf_tree,          &
                                                                           dcds_blf_tree,          &
                                                                           dcds_blf_tree,          &
                                                                           dcds_blf_tree,          &
                                                                           dcds_blf_tree,          &
                                                                           dcds_blf_tree,          &
                                                                           dcds_blf_tree,          &
                                                                           evgn_nlf_tree,          &
                                                                           evgn_nlf_tree,          &
                                                                           dcds_blf_tree,          &
                                                                           dcds_blf_tree,          &
                                                                           dcds_blf_tree,          &
                                                                           evgn_nlf_tree,          &
                                                                           dcds_blf_tree,          &
                                                                           dcds_blf_tree,          &
                                                                           dcds_blf_tree,          &
                                                                           dcds_blf_tree,          &
                                                                           dcds_blf_tree,          &
                                                                           dcds_blf_tree,          &
                                                                           dcds_blf_tree,          &
                                                                           dcds_blf_tree,          &
                                                                           dcds_blf_tree,          &
                                                                           dcds_blf_tree,          &
                                                                           evgn_blf_tree,          &
                                                                           evgn_nlf_tree,          &
                                                                           dcds_blf_tree,          &
                                                                           dcds_blf_tree,          &
                                                                           dcds_blf_tree,          &
                                                                           dcds_blf_tree,          &
                                                                           evgn_nlf_tree,          &
                                                                           evgn_nlf_tree,          &
                                                                           evgn_nlf_tree,          &
                                                                           dcds_blf_tree,          &
                                                                           evgn_nlf_tree,          &
                                                                           dcds_blf_tree,          &
                                                                           dcds_blf_tree,          &
                                                                           evgn_blf_tree           &
                                                                         /)
!
!-- beta values for bvoc classes
    REAL(wp), PARAMETER ::  beta_aceton = 0.10_wp  !< imperical constant for gamma_LI for acetone
    REAL(wp), PARAMETER ::  beta_isop   = 0.13_wp  !< imperical constant for gamma_LI for isoprene
    REAL(wp), PARAMETER ::  beta_mthanl = 0.08_wp  !< imperical constant for gamma_LI for methanol
    REAL(wp), PARAMETER ::  beta_mtrp   = 0.10_wp  !< imperical constant for gamma_LI for monoterpenes
    REAL(wp), PARAMETER ::  beta_ovoc   = 0.10_wp  !< imperical constant for gamma_LI for other vocs
    REAL(wp), PARAMETER ::  beta_sqtp   = 0.17_wp  !< imperical constant for gamma_LI for sesquiterpenes
    REAL(wp), PARAMETER ::  beta_xvoc   = 0.13_wp  !< imperical constant for gamma_LI for oxygenated vocs
!
!   Data from external model (MEGAN) initialized with DATA stamement
    REAL(wp), DIMENSION(n_ebio_spc)        ::  ltab_ebio_mwt  !< molecular weights for bvoc spcs
    REAL(wp), DIMENSION(n_ebio_spc, n_pft) ::  ltab_ef_pft    !< emission foactor for PFTs
!
!-- imperical coeff (beta) for (gamma_LI) for bvoc spcs !light dependence fraction beta for bvoc spcs
    REAL(wp), DIMENSION(n_ebio_spc) ::  ltab_beta_spc = (/                                         &  !< beta coeff for gamma_LI
                                                          beta_isop ,                              &
                                                          beta_mtrp ,                              &
                                                          beta_mtrp ,                              &
                                                          beta_mtrp ,                              &
                                                          beta_mtrp ,                              &
                                                          beta_mtrp ,                              &
                                                          beta_mtrp ,                              &
                                                          beta_mtrp ,                              &
                                                          beta_sqtp,                               &
                                                          beta_xvoc ,                              &
                                                          beta_xvoc ,                              &
                                                          beta_xvoc ,                              &
                                                          beta_mthanl,                             &
                                                          beta_aceton,                             &
                                                          beta_xvoc ,                              &
                                                          beta_xvoc ,                              &
                                                          beta_ovoc,                               &
                                                          beta_ovoc,                               &
                                                          beta_ovoc,                               &
                                                          beta_ovoc,                               &
                                                          beta_ovoc,                               &
                                                          beta_ovoc,                               &
                                                          beta_ovoc,                               &
                                                          beta_ovoc,                               &
                                                          beta_ovoc,                               &
                                                          beta_ovoc,                               &
                                                          beta_ovoc,                               &
                                                          beta_ovoc,                               &
                                                          beta_ovoc,                               &
                                                          beta_ovoc                                &
                                                         /)
!
    REAL(wp), DIMENSION(1:n_dsoil_layers) ::  ltab_d_soil = (/                                     &  !< default soil layer config
                                                              0.01_wp,                             &
                                                              0.02_wp,                             &
                                                              0.04_wp,                             &
                                                              0.06_wp,                             &
                                                              0.14_wp,                             &
                                                              0.26_wp,                             &
                                                              0.54_wp,                             &
                                                              1.86_wp                              &
                                                             /)
!
    REAL(wp), DIMENSION(n_ebio_class) ::  ltab_ldf = (/                                            &  !< LI depend frac-bvoc class
                                                       1.0_wp,                                     &
                                                       0.4_wp,                                     &
                                                       0.5_wp,                                     &
                                                       0.8_wp,                                     &
                                                       0.5_wp                                      &
                                                      /)
!
    REAL(wp), DIMENSION(n_ebio_spc) ::  ltab_ldf_spc = (/                                          &  !< LI depend frac ldf - bvoc
                                                         1.00_wp,                                  &
                                                         0.60_wp,                                  &
                                                         0.60_wp,                                  &
                                                         0.20_wp,                                  &
                                                         0.20_wp,                                  &
                                                         0.80_wp,                                  &
                                                         0.60_wp,                                  &
                                                         0.20_wp,                                  &
                                                         0.50_wp,                                  &
                                                         0.80_wp,                                  &
                                                         0.80_wp,                                  &
                                                         0.80_wp,                                  &
                                                         0.80_wp,                                  &
                                                         0.2_wp ,                                  &
                                                         0.80_wp,                                  &
                                                         0.80_wp,                                  &
                                                         1.00_wp,                                  &
                                                         0.20_wp,                                  &
                                                         0.20_wp,                                  &
                                                         0.20_wp,                                  &
                                                         0.20_wp,                                  &
                                                         0.20_wp,                                  &
                                                         0.20_wp,                                  &
                                                         0.20_wp,                                  &
                                                         0.20_wp,                                  &
                                                         0.20_wp,                                  &
                                                         0.20_wp,                                  &
                                                         0.20_wp,                                  &
                                                         0.20_wp,                                  &
                                                         0.20_wp                                   &
                                                         /)
!
!-- Root distribution for default soil layer configuration (sum = 1)
    REAL(wp), DIMENSION(0:n_root_layers,1:n_pft), PARAMETER ::  ltab_root_dist = RESHAPE( (/       &  !<  Root distr - soil layers
                                                              1.00_wp, 0.00_wp, 0.00_wp, 0.00_wp,  &
                                                              0.24_wp, 0.41_wp, 0.31_wp, 0.04_wp,  &
                                                              0.35_wp, 0.38_wp, 0.23_wp, 0.04_wp,  &
                                                              0.26_wp, 0.39_wp, 0.29_wp, 0.06_wp,  &
                                                              0.26_wp, 0.38_wp, 0.29_wp, 0.07_wp,  &
                                                              0.24_wp, 0.38_wp, 0.31_wp, 0.07_wp,  &
                                                              0.25_wp, 0.34_wp, 0.27_wp, 0.14_wp,  &
                                                              0.27_wp, 0.27_wp, 0.27_wp, 0.09_wp,  &
                                                              1.00_wp, 0.00_wp, 0.00_wp, 0.00_wp,  &
                                                              0.47_wp, 0.45_wp, 0.08_wp, 0.00_wp,  &
                                                              0.24_wp, 0.41_wp, 0.31_wp, 0.04_wp,  &
                                                              0.17_wp, 0.31_wp, 0.33_wp, 0.19_wp,  &
                                                              0.00_wp, 0.00_wp, 0.00_wp, 0.00_wp,  &
                                                              0.25_wp, 0.34_wp, 0.27_wp, 0.11_wp,  &
                                                              0.23_wp, 0.36_wp, 0.30_wp, 0.11_wp,  &
                                                              0.23_wp, 0.36_wp, 0.30_wp, 0.11_wp,  &
                                                              0.19_wp, 0.35_wp, 0.36_wp, 0.10_wp,  &
                                                              0.19_wp, 0.35_wp, 0.36_wp, 0.10_wp   &
                                                                                 /), (/ 4, 18 /) )
!
    REAL(wp), DIMENSION(0:n_root_layers)       ::  ltab_root_layer_depth =(/                       &  !<  root layer depth of trees
                                                                            0.07_wp,               &
                                                                            0.28_wp,               &
                                                                            1.00_wp,               &
                                                                            2.80_wp                &
                                                                           /)
!
!-- field capacity etc: https://palm.muk.uni-hannover.de/trac/wiki/doc/app/land_surface_parameters#vegetation_type
!-- Field capacity (max water content) in m3/m3 for 6 soil types
    REAL(wp), DIMENSION(1:n_soil_types), PARAMETER :: ltab_qc = (/                                 &  !< field capacity-soil types
                                                                  0.244_wp,                        &
                                                                  0.347_wp,                        &
                                                                  0.383_wp,                        &
                                                                  0.448_wp,                        &
                                                                  0.541_wp,                        &
                                                                  0.663_wp                         &
                                                                 /)
!
!-- soil types for  Wilting point (qw) in m3/m3 for 6 soil types
    REAL(wp), DIMENSION(1:n_soil_types), PARAMETER :: ltab_qw = (/                                 &  !< wiliting point-soil types
                                                                  0.059_wp,                        &
                                                                  0.151_wp,                        &
                                                                  0.131_wp,                        &
                                                                  0.279_wp,                        &
                                                                  0.335_wp,                        &
                                                                  0.267_wp                         &
                                                                 /)
!
!-- individual tree species (roadside) specific leaf mass (slm) in cm2 g-1 for BERLIN.
!-- unit coversion to lad (g m-3);   lad(g m-3) =lad(m2 m-3) * lad(g m-2);  lad(g m-2) = 10000/ slm(cm2 g-1) .
    REAL(wp), DIMENSION(1: n_tree_types), PARAMETER ::  ltab_tree_spc_slm = (/                     &  !< specific leaf mass - trees
                                                                              84.3_wp,             &
                                                                              214.7_wp,            &
                                                                              216.3_wp,            &
                                                                              123.6_wp,            &
                                                                              156.8_wp,            &
                                                                              109.8_wp,            &
                                                                              272.8_wp,            &
                                                                              79.9_wp ,            &
                                                                              82.0_wp,             &
                                                                              -9999.0_wp,          &
                                                                              231.3_wp,            &
                                                                              192.3_wp,            &
                                                                              179.6_wp,            &
                                                                              264.1_wp,            &
                                                                              38.3_wp,             &
                                                                              172.4_wp,            &
                                                                              229.2_wp,            &
                                                                              206.6_wp,            &
                                                                              68.8_wp,             &
                                                                              208.8_wp,            &
                                                                              222.5_wp,            &
                                                                              343.0_wp,            &
                                                                              159.1_wp,            &
                                                                              219.6_wp,            &
                                                                              57.8_wp,             &
                                                                              -9999.0_wp,          &
                                                                              30.6_wp,             &
                                                                              62.1_wp,             &
                                                                              261.2_wp,            &
                                                                              201.6_wp,            &
                                                                              -9999.0_wp,          &
                                                                              132.2_wp,            &
                                                                              268.0_wp,            &
                                                                              144.9_wp,            &
                                                                              139.2_wp,            &
                                                                              177.8_wp,            &
                                                                              222.3_wp,            &
                                                                              130.2_wp,            &
                                                                              65.3_wp,             &
                                                                              240.8_wp,            &
                                                                              39.9_wp,             &
                                                                              132.7_wp,            &
                                                                              242.2_wp,            &
                                                                              146.6_wp,            &
                                                                              81.4_wp,             &
                                                                              229.1_wp,            &
                                                                              223.8_wp,            &
                                                                              136.4_wp,            &
                                                                              269.6_wp,            &
                                                                              190.2_wp,            &
                                                                              208.0_wp,            &
                                                                              323.9_wp,            &
                                                                              234.7_wp,            &
                                                                              337.5_wp,            &
                                                                              136.3_wp,            &
                                                                              210.7_wp,            &
                                                                              46.10_wp,            &
                                                                              47.20_wp,            &
                                                                              159.6_wp,            &
                                                                              102.6_wp,            &
                                                                              126.1_wp,            &
                                                                              162.0_wp,            &
                                                                              -9999.0_wp,          &
                                                                              -9999.0_wp,          &
                                                                              283.6_wp,            &
                                                                              99.8_wp,             &
                                                                              157.8_wp,            &
                                                                              143.8_wp,            &
                                                                              165.5_wp,            &
                                                                              350.3_wp,            &
                                                                              143.6_wp,            &
                                                                              177.0_wp,            &
                                                                              -9999.0_wp,          &
                                                                              -9999.0_wp,          &
                                                                              131.8_wp,            &
                                                                              141.6_wp,            &
                                                                              162.0_wp,            &
                                                                              79.0_wp,             &
                                                                              175.2_wp,            &
                                                                              93.0_wp ,            &
                                                                              77.5_wp,             &
                                                                              260.5_wp,            &
                                                                              97.1_wp,             &
                                                                              199.0_wp,            &
                                                                              148.0_wp,            &
                                                                              -9999.0_wp           &
                                                                              /)
!
!-- Data from extenal models (MEGAN, ECHAM6-HAMMOZ),mapping bvoc spcs to the relevent bvoc classes
!-- Isoprene
    DATA  ltab_ebio_map(isop)  / class_isop /
!
!-- Monoterpenes MRTP
    DATA  ltab_ebio_map(mtrp1) / class_mtrp /
    DATA  ltab_ebio_map(mtrp2) / class_mtrp /
    DATA  ltab_ebio_map(mtrp3) / class_mtrp /
    DATA  ltab_ebio_map(mtrp4) / class_mtrp /
    DATA  ltab_ebio_map(mtrp5) / class_mtrp /
    DATA  ltab_ebio_map(mtrp6) / class_mtrp /
    DATA  ltab_ebio_map(mtrp7) / class_mtrp /
!
!-- SQT
    DATA  ltab_ebio_map(sqtp1) / class_sqtp /
!
!-- Oxygenated VOCs
    DATA  ltab_ebio_map(xvoc1) / class_xvoc /
    DATA  ltab_ebio_map(xvoc2) / class_xvoc /
    DATA  ltab_ebio_map(xvoc3) / class_xvoc /
    DATA  ltab_ebio_map(xvoc4) / class_xvoc /
    DATA  ltab_ebio_map(xvoc5) / class_xvoc /
    DATA  ltab_ebio_map(xvoc6) / class_xvoc /
    DATA  ltab_ebio_map(xvoc7) / class_xvoc /
!
!-- OVOCs
    DATA  ltab_ebio_map(ovoc1)  / class_ovoc /
    DATA  ltab_ebio_map(ovoc2)  / class_ovoc /
    DATA  ltab_ebio_map(ovoc3)  / class_ovoc /
    DATA  ltab_ebio_map(ovoc4)  / class_ovoc /
    DATA  ltab_ebio_map(ovoc5)  / class_ovoc /
    DATA  ltab_ebio_map(ovoc6)  / class_ovoc /
    DATA  ltab_ebio_map(ovoc7)  / class_ovoc /
    DATA  ltab_ebio_map(ovoc8)  / class_ovoc /
    DATA  ltab_ebio_map(ovoc9)  / class_ovoc /
    DATA  ltab_ebio_map(ovoc10) / class_ovoc /
    DATA  ltab_ebio_map(ovoc11) / class_ovoc /
    DATA  ltab_ebio_map(ovoc12) / class_ovoc /
    DATA  ltab_ebio_map(ovoc13) / class_ovoc /
    DATA  ltab_ebio_map(ovoc14) / class_ovoc /
!
!-- defining emission factors (ef) from external model (MEGAN, ECHAM-HAMMOZ) for each of 30 spcs
!-- for all valid pfts from ecmwf-ifs, units  i ug m-2 h-1
    DATA ltab_ef_pft(isop, 1:n_pft)   /-9999.0_wp,     1.0_wp,   200.0_wp,                         &
                                         600.0_wp,     1.0_wp, 10000.0_wp,                         &
                                       10000.0_wp, -9999.0_wp, -9999.0_wp,                         &
                                       -9999.0_wp, -9999.0_wp, -9999.0_wp,                         &
                                       -9999.0_wp, -9999.0_wp,  2000.0_wp,                         &
                                        4000.0_wp, -9999.0_wp, -9999.0_wp /
    DATA ltab_ef_pft(mtrp1, 1:n_pft)  /-9999.0_wp,     0.3_wp,     0.3_wp,                         &
                                          70.0_wp,    60.0_wp,    30.0_wp,                         &
                                          30.0_wp, -9999.0_wp, -9999.0_wp,                         &
                                       -9999.0_wp, -9999.0_wp, -9999.0_wp,                         &
                                       -9999.0_wp, -9999.0_wp,    30.0_wp,                         &
                                          50.0_wp, -9999.0_wp, -9999.0_wp /
    DATA ltab_ef_pft(mtrp2, 1:n_pft)  /-9999.0_wp,     0.7_wp,     0.7_wp,                         &
                                          70.0_wp,     40.0_wp,   50.0_wp,                         &
                                          50.0_wp, -9999.0_wp, -9999.0_wp,                         &
                                       -9999.0_wp, -9999.0_wp, -9999.0_wp,                         &
                                       -9999.0_wp, -9999.0_wp,    50.0_wp,                         &
                                          70.0_wp, -9999.0_wp, -9999.0_wp /
    DATA ltab_ef_pft(mtrp3, 1:n_pft)  /-9999.0_wp,     0.7_wp,     0.7_wp,                         &
                                         100.0_wp,   130.0_wp,    80.0_wp,                         &
                                          80.0_wp, -9999.0_wp, -9999.0_wp,                         &
                                       -9999.0_wp, -9999.0_wp, -9999.0_wp,                         &
                                       -9999.0_wp, -9999.0_wp,    60.0_wp,                         &
                                         100.0_wp, -9999.0_wp, -9999.0_wp /
    DATA ltab_ef_pft(mtrp4, 1:n_pft)  /-9999.0_wp,     0.3_wp,     0.3_wp,                         &
                                         160.0_wp,    80.0_wp,    30.0_wp,                         &
                                          30.0_wp, -9999.0_wp, -9999.0_wp,                         &
                                       -9999.0_wp, -9999.0_wp, -9999.0_wp,                         &
                                       -9999.0_wp, -9999.0_wp,    30.0_wp,                         &
                                         100.0_wp, -9999.0_wp, -9999.0_wp /
    DATA ltab_ef_pft(mtrp5, 1:n_pft)  /-9999.0_wp,     2.0_wp,     2.0_wp,                         &
                                          70.0_wp,    60.0_wp,   120.0_wp,                         &
                                         120.0_wp, -9999.0_wp, -9999.0_wp,                         &
                                       -9999.0_wp, -9999.0_wp, -9999.0_wp,                         &
                                       -9999.0_wp, -9999.0_wp,    90.0_wp,                         &
                                         150.0_wp, -9999.0_wp, -9999.0_wp /
    DATA ltab_ef_pft(mtrp6, 1:n_pft)  /-9999.0_wp,     2.0_wp,     2.0_wp,                         &
                                         500.0_wp,   510.0_wp,   400.0_wp,                         &
                                         400.0_wp, -9999.0_wp, -9999.0_wp,                         &
                                       -9999.0_wp, -9999.0_wp, -9999.0_wp,                         &
                                       -9999.0_wp, -9999.0_wp,   200.0_wp,                         &
                                         300.0_wp, -9999.0_wp, -9999.0_wp /
    DATA ltab_ef_pft(mtrp7, 1:n_pft)  /-9999.0_wp,     1.5_wp,     1.5_wp,                         &
                                         300.0_wp,   200.0_wp,   130.0_wp,                         &
                                         130.0_wp, -9999.0_wp, -9999.0_wp,                         &
                                       -9999.0_wp, -9999.0_wp, -9999.0_wp,                         &
                                       -9999.0_wp, -9999.0_wp,   100.0_wp,                         &
                                         150.0_wp, -9999.0_wp, -9999.0_wp /
    DATA ltab_ef_pft(sqtp1, 1:n_pft)  /-9999.0_wp,     4.0_wp,     1.0_wp,                         &
                                          80.0_wp,    80.0_wp,    40.0_wp,                         &
                                          40.0_wp, -9999.0_wp, -9999.0_wp,                         &
                                       -9999.0_wp, -9999.0_wp, -9999.0_wp,                         &
                                       -9999.0_wp, -9999.0_wp,    50.0_wp,                         &
                                          50.0_wp, -9999.0_wp, -9999.0_wp /
    DATA ltab_ef_pft(xvoc1, 1:n_pft)  /-9999.0_wp,    20.0_wp,    20.0_wp,                         &
                                         200.0_wp,   200.0_wp,   200.0_wp,                         &
                                         200.0_wp, -9999.0_wp, -9999.0_wp,                         &
                                       -9999.0_wp, -9999.0_wp, -9999.0_wp,                         &
                                       -9999.0_wp, -9999.0_wp,   200.0_wp,                         &
                                         200.0_wp, -9999.0_wp, -9999.0_wp /
    DATA ltab_ef_pft(xvoc2, 1:n_pft)  /-9999.0_wp,    20.0_wp,    20.0_wp,                         &
                                         200.0_wp,   200.0_wp,   200.0_wp,                         &
                                         200.0_wp, -9999.0_wp, -9999.0_wp,                         &
                                       -9999.0_wp, -9999.0_wp, -9999.0_wp,                         &
                                       -9999.0_wp, -9999.0_wp,   200.0_wp,                         &
                                         200.0_wp, -9999.0_wp, -9999.0_wp /
    DATA ltab_ef_pft(xvoc3, 1:n_pft)  /-9999.0_wp,    16.0_wp,    16.0_wp,                         &
                                          40.0_wp,    40.0_wp,    40.0_wp,                         &
                                          40.0_wp, -9999.0_wp, -9999.0_wp,                         &
                                       -9999.0_wp, -9999.0_wp, -9999.0_wp,                         &
                                       -9999.0_wp, -9999.0_wp,    40.0_wp,                         &
                                          40.0_wp, -9999.0_wp, -9999.0_wp /
    DATA ltab_ef_pft(xvoc4, 1:n_pft)  /-9999.0_wp,    12.0_wp,    12.0_wp,                         &
                                          30.0_wp,    30.0_wp,    30.0_wp,                         &
                                          30.0_wp, -9999.0_wp, -9999.0_wp,                         &
                                       -9999.0_wp, -9999.0_wp, -9999.0_wp,                         &
                                       -9999.0_wp, -9999.0_wp,    30.0_wp,                         &
                                          30.0_wp, -9999.0_wp, -9999.0_wp /
    DATA ltab_ef_pft(xvoc5, 1:n_pft)  /-9999.0_wp,   900.0_wp,   500.0_wp,                         &
                                         900.0_wp,   900.0_wp,   900.0_wp,                         &
                                         900.0_wp, -9999.0_wp, -9999.0_wp,                         &
                                       -9999.0_wp, -9999.0_wp,- 9999.0_wp,                         &
                                       -9999.0_wp, -9999.0_wp,   900.0_wp,                         &
                                         900.0_wp, -9999.0_wp, -9999.0_wp /
    DATA ltab_ef_pft(xvoc6, 1:n_pft)  /-9999.0_wp,    12.0_wp,    12.0_wp,                         &
                                          30.0_wp,    30.0_wp,    30.0_wp,                         &
                                          30.0_wp, -9999.0_wp, -9999.0_wp,                         &
                                       -9999.0_wp, -9999.0_wp, -9999.0_wp,                         &
                                       -9999.0_wp, -9999.0_wp,    30.0_wp,                         &
                                          30.0_wp, -9999.0_wp, -9999.0_wp /
    DATA ltab_ef_pft(xvoc7, 1:n_pft)  /-9999.0_wp,    80.0_wp,    80.0_wp,                         &
                                         240.0_wp,   240.0_wp,   240.0_wp,                         &
                                         240.0_wp, -9999.0_wp, -9999.0_wp,                         &
                                       -9999.0_wp, -9999.0_wp, -9999.0_wp,                         &
                                       -9999.0_wp, -9999.0_wp,   240.0_wp,                         &
                                         240.0_wp, -9999.0_wp, -9999.0_wp /
    DATA ltab_ef_pft(ovoc1, 1:n_pft)  /-9999.0_wp,     0.0_wp,     0.0_wp,                         &
                                         700.0_wp,     0.0_wp,     0.0_wp,                         &
                                           2.0_wp, -9999.0_wp, -9999.0_wp,                         &
                                       -9999.0_wp, -9999.0_wp, -9999.0_wp,                         &
                                       -9999.0_wp, -9999.0_wp,     0.0_wp,                         &
                                           0.0_wp, -9999.0_wp, -9999.0_wp /
    DATA ltab_ef_pft(ovoc2, 1:n_pft)  /-9999.0_wp,     0.7_wp,     0.7_wp,                         &
                                           0.7_wp,     0.7_wp,     0.7_wp,                         &
                                           0.7_wp, -9999.0_wp, -9999.0_wp,                         &
                                       -9999.0_wp, -9999.0_wp, -9999.0_wp,                         &
                                       -9999.0_wp, -9999.0_wp,     0.7_wp,                         &
                                           0.7_wp, -9999.0_wp, -9999.0_wp /
    DATA ltab_ef_pft(ovoc3, 1:n_pft)  /-9999.0_wp,   174.0_wp,   174.0_wp,                         &
                                         174.0_wp,   174.0_wp,   174.0_wp,                         &
                                         174.0_wp, -9999.0_wp, -9999.0_wp,                         &
                                       -9999.0_wp, -9999.0_wp, -9999.0_wp,                         &
                                       -9999.0_wp, -9999.0_wp,   174.0_wp,                         &
                                         174.0_wp, -9999.0_wp, -9999.0_wp /
    DATA ltab_ef_pft(ovoc4, 1:n_pft)  /-9999.0_wp,     4.5_wp,     4.5_wp,                         &
                                           4.5_wp,     4.5_wp,     4.5_wp,                         &
                                           4.5_wp, -9999.0_wp, -9999.0_wp,                         &
                                       -9999.0_wp, -9999.0_wp, -9999.0_wp,                         &
                                       -9999.0_wp, -9999.0_wp,     4.5_wp,                         &
                                           4.5_wp, -9999.0_wp, -9999.0_wp /
    DATA ltab_ef_pft(ovoc5, 1:n_pft)  /-9999.0_wp,     9.0_wp,     9.0_wp,                         &
                                           9.0_wp,     9.0_wp,     9.0_wp,                         &
                                           9.0_wp, -9999.0_wp, -9999.0_wp,                         &
                                       -9999.0_wp, -9999.0_wp, -9999.0_wp,                         &
                                       -9999.0_wp, -9999.0_wp,     9.0_wp,                         &
                                           9.0_wp, -9999.0_wp, -9999.0_wp /
    DATA ltab_ef_pft(ovoc6, 1:n_pft)  /-9999.0_wp,     2.8_wp,     2.8_wp,                         &
                                           2.8_wp,     2.8_wp,     2.8_wp,                         &
                                           2.8_wp, -9999.0_wp, -9999.0_wp,                         &
                                       -9999.0_wp, -9999.0_wp, -9999.0_wp,                         &
                                       -9999.0_wp, -9999.0_wp,     2.8_wp,                         &
                                           2.8_wp, -9999.0_wp, -9999.0_wp /
    DATA ltab_ef_pft(ovoc7, 1:n_pft)  /-9999.0_wp,     1.4_wp,     1.4_wp,                         &
                                           1.4_wp,     1.4_wp,     1.4_wp,                         &
                                           1.4_wp, -9999.0_wp, -9999.0_wp,                         &
                                       -9999.0_wp, -9999.0_wp, -9999.0_wp,                         &
                                       -9999.0_wp, -9999.0_wp,     1.4_wp,                         &
                                           1.4_wp, -9999.0_wp, -9999.0_wp /
    DATA ltab_ef_pft(ovoc8, 1:n_pft)  /-9999.0_wp,     0.1_wp,     0.1_wp,                         &
                                           0.1_wp,     0.1_wp,     0.1_wp,                         &
                                           0.1_wp, -9999.0_wp, -9999.0_wp,                         &
                                       -9999.0_wp, -9999.0_wp, -9999.0_wp,                         &
                                       -9999.0_wp, -9999.0_wp,     0.1_wp,                         &
                                           0.1_wp, -9999.0_wp, -9999.0_wp /
    DATA ltab_ef_pft(ovoc9, 1:n_pft)  /-9999.0_wp,     0.4_wp,     0.4_wp,                         &
                                           0.4_wp,     0.4_wp,     0.4_wp,                         &
                                           0.4_wp, -9999.0_wp, -9999.0_wp,                         &
                                       -9999.0_wp, -9999.0_wp, -9999.0_wp,                         &
                                       -9999.0_wp, -9999.0_wp,     0.4_wp,                         &
                                           0.4_wp, -9999.0_wp, -9999.0_wp /
    DATA ltab_ef_pft(ovoc10, 1:n_pft) /-9999.0_wp,     1.4_wp,     1.4_wp,                         &
                                           1.4_wp,     1.4_wp,     1.4_wp,                         &
                                           1.4_wp, -9999.0_wp, -9999.0_wp,                         &
                                       -9999.0_wp, -9999.0_wp, -9999.0_wp,                         &
                                       -9999.0_wp, -9999.0_wp,     1.4_wp,                         &
                                           1.4_wp, -9999.0_wp, -9999.0_wp /
    DATA ltab_ef_pft(ovoc11, 1:n_pft) /-9999.0_wp,     0.1_wp,     0.1_wp,                         &
                                           0.1_wp,     0.1_wp,     0.1_wp,                         &
                                           0.1_wp, -9999.0_wp, -9999.0_wp,                         &
                                       -9999.0_wp, -9999.0_wp, -9999.0_wp,                         &
                                       -9999.0_wp, -9999.0_wp,     0.1_wp,                         &
                                           0.1_wp, -9999.0_wp, -9999.0_wp /
    DATA ltab_ef_pft(ovoc12, 1:n_pft) /-9999.0_wp,    33.6_wp,    33.6_wp,                         &
                                          33.6_wp,    33.6_wp,    33.6_wp,                         &
                                          33.6_wp, -9999.0_wp, -9999.0_wp,                         &
                                       -9999.0_wp, -9999.0_wp, -9999.0_wp,                         &
                                       -9999.0_wp, -9999.0_wp,    33.6_wp,                         &
                                          33.6_wp, -9999.0_wp, -9999.0_wp /
    DATA ltab_ef_pft(ovoc13, 1:n_pft) /-9999.0_wp,    67.2_wp,    67.2_wp,                         &
                                          67.2_wp,    67.2_wp,    67.2_wp,                         &
                                          67.2_wp, -9999.0_wp, -9999.0_wp,                         &
                                       -9999.0_wp, -9999.0_wp, -9999.0_wp,                         &
                                       -9999.0_wp, -9999.0_wp,    67.2_wp,                         &
                                          67.2_wp, -9999.0_wp, -9999.0_wp /
    DATA ltab_ef_pft(ovoc14, 1:n_pft) /-9999.0_wp,     0.1_wp,     0.1_wp,                         &
                                           0.1_wp,     0.1_wp,     0.1_wp,                         &
                                           0.1_wp, -9999.0_wp, -9999.0_wp,                         &
                                       -9999.0_wp, -9999.0_wp, -9999.0_wp,                         &
                                       -9999.0_wp, -9999.0_wp,     0.1_wp,                         &
                                           0.1_wp, -9999.0_wp, -9999.0_wp /
!
!-- Molecular weight of biogenic species in g/mol
!-- isoprene
    DATA  ltab_ebio_mwt(isop)  / 68.12_wp /
!
!-- MTP
    DATA  ltab_ebio_mwt(mtrp1) / 136.23_wp /
    DATA  ltab_ebio_mwt(mtrp2) / 136.23_wp /
    DATA  ltab_ebio_mwt(mtrp3) / 136.23_wp /
    DATA  ltab_ebio_mwt(mtrp4) / 136.23_wp /
    DATA  ltab_ebio_mwt(mtrp5) / 136.23_wp /
    DATA  ltab_ebio_mwt(mtrp6) / 136.23_wp /
    DATA  ltab_ebio_mwt(mtrp7) / 136.23_wp /
!
!-- SQT
    DATA  ltab_ebio_mwt(sqtp1) / 204.35_wp /
!
!-- Oxygenated VOCs
    DATA  ltab_ebio_mwt(xvoc1) / 44.05_wp /
    DATA  ltab_ebio_mwt(xvoc2) / 46.07_wp /
    DATA  ltab_ebio_mwt(xvoc3) / 30.03_wp /
    DATA  ltab_ebio_mwt(xvoc4) / 32.04_wp /
    DATA  ltab_ebio_mwt(xvoc5) / 58.08_wp /
    DATA  ltab_ebio_mwt(xvoc6) / 46.03_wp /
    DATA  ltab_ebio_mwt(xvoc7) / 60.05_wp /
!
!-- Other VOC
    DATA  ltab_ebio_mwt(ovoc1)  / 86.13_wp /
    DATA  ltab_ebio_mwt(ovoc2)  / 16.04_wp /
    DATA  ltab_ebio_mwt(ovoc3)  / 28.05_wp /
    DATA  ltab_ebio_mwt(ovoc4)  / 27.03_wp /
    DATA  ltab_ebio_mwt(ovoc5)  / 92.14_wp /
    DATA  ltab_ebio_mwt(ovoc6)  / 94.94_wp /
    DATA  ltab_ebio_mwt(ovoc7)  / 50.49_wp /
    DATA  ltab_ebio_mwt(ovoc8)  / 141.94_wp /
    DATA  ltab_ebio_mwt(ovoc9)  / 62.14_wp /
    DATA  ltab_ebio_mwt(ovoc10) / 30.07_wp /
    DATA  ltab_ebio_mwt(ovoc11) / 44.10_wp /
    DATA  ltab_ebio_mwt(ovoc12) / 42.08_wp /
    DATA  ltab_ebio_mwt(ovoc13) / 58.12_wp /
    DATA  ltab_ebio_mwt(ovoc14) / 106.12_wp /

 END MODULE chem_emis_biogenic_data_mod
