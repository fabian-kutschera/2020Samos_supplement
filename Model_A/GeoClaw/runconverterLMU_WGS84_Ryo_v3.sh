#!/bin/sh
input=Ryo_v3_noWL

dx=0.005
output=${input}_${dx}_tanioka.nc
bathy=gebco_2022_n40.0_s36.0_w23.0_e29.0_converted.nc

#convert to smaller data set and keep .h5 
python /import/freenas-m-05-seissol/kutschera/HIWI/Scripts/Script4Processing/extractDataFromUnstructuredOutput.py /import/freenas-m-05-seissol/kutschera/MAthesis_fkutschera/simulations/Samos_final_Ryo_v3/${input}-surface.xdmf --Data u1 u2 u3 --idt $(seq 0 100) --backend hdf5

#convert u1u2u3 to UVW - works on .h5
python /import/freenas-m-05-seissol/kutschera/HIWI/Scripts/Script4Processing/displacement-converter/convert_u1u2u3_2_UVW.py ${input}_resampled-surface.xdmf

#convert Seissol geometry to geographic coordinate system (lat, lon)
python /import/freenas-m-05-seissol/kutschera/HIWI/Scripts/Script4Processing/displacement-converter/convert_geographic_SeisSol_geom_v0.9.py ${input}_resampledUVW-surface.xdmf '+proj=tmerc +datum=WGS84 +k=0.9996 +lon_0=26.25 +lat_0=37.75'

#convert names in GEBCO file
#bash convert_names.sh

#convert from unstructured seissol format to structured netcdf
/import/freenas-m-05-seissol/SamoaRelated/displacement-converter/build/displacement-converter -i WGS84_${input}_resampledUVW-surface.xdmf -o $output -b $bathy --dx $dx --dy $dx -x26 -w2 -y37 -z 2 --dt 2

#apply tappering on the domain boundary for a smooth transition to 0
python /import/freenas-m-05-seissol/kutschera/HIWI/Scripts/Script4Processing/displacement-converter/tapperNetcdf.py $output --nHann 100

#convert netcdf file to tt3 (here using 4 ranks)
python /import/freenas-m-05-seissol/kutschera/HIWI/Scripts/Script4Processing/displacement-converter/convert_netcdf_tt3.py tapered_$output --MP 4
