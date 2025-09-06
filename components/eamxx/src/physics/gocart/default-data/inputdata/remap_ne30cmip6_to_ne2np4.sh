MAP_FILE="/global/homes/f/fust1/maps/map_ne30cmip6_to_ne2np4_nco.c20250716.nc"
DECK30_DATA_DIR="/global/cfs/cdirs/e3sm/inputdata/atm/cam/chem/trop_mozart_aero/emis/DECK_ne30"
TWODEGREES_DATA_DIR="/global/cfs/cdirs/e3sm/inputdata/atm/cam/chem/trop_mozart_aero/emis/chem_gases/2degrees"
REMAPPED_DATA_DIR="/global/u1/f/fust1/inputdata/ne2np4"

DECK30_NC_FILES=("cmip6_mam4_so2_elev_1x1_2010_clim_c20190821.nc" \
                 "emissions-cmip6_e3sm_SOAG0_elev_2010_clim_1.9x2.5_c20230213.nc" \
                 "cmip6_mam4_bc_a4_elev_1x1_2010_clim_c20190821.nc" \
                 "cmip6_mam4_num_a1_elev_1x1_2010_clim_c20190821.nc" \
                 "cmip6_mam4_num_a2_elev_1x1_2010_clim_c20190821.nc" \
                 "cmip6_mam4_num_a4_elev_1x1_2010_clim_c20190821.nc" \
                 "cmip6_mam4_pom_a4_elev_1x1_2010_clim_c20190821.nc" \
                 "cmip6_mam4_so4_a1_elev_1x1_2010_clim_c20190821.nc" \
                 "cmip6_mam4_so4_a2_elev_1x1_2010_clim_c20190821.nc" \
                 "cmip6_mam4_so2_surf_1x1_2010_clim_c20190821.nc" \
                 "emissions-cmip6_e3sm_SOAG0_surf_2010_clim_1.9x2.5_c20230213.nc" \
                 "cmip6_mam4_bc_a4_surf_1x1_2010_clim_c20190821.nc" \
                 "cmip6_mam4_num_a1_surf_1x1_2010_clim_c20190821.nc" \
                 "cmip6_mam4_num_a2_surf_1x1_2010_clim_c20190821.nc" \
                 "cmip6_mam4_num_a4_surf_1x1_2010_clim_c20190821.nc" \
                 "cmip6_mam4_pom_a4_surf_1x1_2010_clim_c20190821.nc" \
                 "cmip6_mam4_so4_a1_surf_1x1_2010_clim_c20190821.nc" \
                 "cmip6_mam4_so4_a2_surf_1x1_2010_clim_c20190821.nc")
                 
                
TWODEGREES_NC_FILES=("emissions-cmip6_e3sm_MTERP_surface_2010_clim_1.9x2.5_c20230213.nc" \
                     "emissions-cmip6_e3sm_C2H4_surface_2010_clim_1.9x2.5_c20230213.nc" \
                     "emissions-cmip6_e3sm_C2H6_surface_2010_clim_1.9x2.5_c20230213.nc" \
                     "emissions-cmip6_e3sm_C3H8_surface_2010_clim_1.9x2.5_c20230213.nc" \
                     "emissions-cmip6_e3sm_CH2O_surface_2010_clim_1.9x2.5_c20230213.nc" \
                     "emissions-cmip6_e3sm_CH3CHO_surface_2010_clim_1.9x2.5_c20230213.nc" \
                     "emissions-cmip6_e3sm_CH3COCH3_surface_2010_clim_1.9x2.5_c20230213.nc" \
                     "emissions-cmip6_e3sm_CO_surface_2010_clim_1.9x2.5_c20230213.nc" \
                     "emissions_E90_surface_2010_clim_1.9x2.5_c20230213.nc" \
                     "emissions-cmip6_e3sm_ISOP_surface_2010_clim_1.9x2.5_c20230213.nc" \
                     "emissions-cmip6_e3sm_ISOP_surface_2010_clim_1.9x2.5_c20230213.nc" \
                     "emissions-cmip6_e3sm_NO_surface_2010_clim_1.9x2.5_c20230213.nc")

# Load environment
source /global/common/software/e3sm/anaconda_envs/load_latest_e3sm_unified_pm-cpu.sh

# Deck30 data
for FILE in "${DECK30_NC_FILES[@]}"; do
    # Remap files (if they haven't already been remapped)
    if [ -f "${REMAPPED_DATA_DIR}/${FILE}" ]; then
        echo "${FILE} has already been remapped..."
    else
        #echo "ncremap -m ${MAP_FILE} -i ${DATA_DIR}/${FILE} -o ${REMAPPED_DATA_DIR}/${FILE}"
        ncremap -m "${MAP_FILE}" -i "${DECK30_DATA_DIR}/${FILE}" -o "${REMAPPED_DATA_DIR}/${FILE}"
    fi 
done

# 2degrees data
for FILE in "${TWODEGREES_NC_FILES[@]}"; do
    # Remap files (if they haven't already been remapped)
    if [ -f "${REMAPPED_DATA_DIR}/${FILE}" ]; then
        echo "${FILE} has already been remapped..."
    else
        #echo "ncremap -m ${MAP_FILE} -i ${DATA_DIR}/${FILE} -o ${REMAPPED_DATA_DIR}/${FILE}"
        ncremap -m "${MAP_FILE}" -i "${TWODEGREES_DATA_DIR}/${FILE}" -o "${REMAPPED_DATA_DIR}/${FILE}"
    fi 
done

# DMS is only file from emis directory
FILE="DMSflux.2010.1deg_latlon_conserv.POPmonthlyClimFromACES4BGC_c20190220.nc"
DMS_DIR="/global/cfs/cdirs/e3sm/inputdata/atm/cam/chem/trop_mozart_aero/emis" 
MAP_FILE="/global/homes/f/fust1/maps/map_180x360_to_ne2np4_nco.c20250716.nc"
if [ -f "${REMAPPED_DATA_DIR}/${FILE}" ]; then
    echo "${FILE} has already been remapped..."
else
    #echo "ncremap -m ${MAP_FILE} -i ${DATA_DIR}/${FILE} -o ${REMAPPED_DATA_DIR}/${FILE}"
    ncremap -m "${MAP_FILE}" -i "${DMS_DIR}/${FILE}" -o "${REMAPPED_DATA_DIR}/${FILE}"
fi 



