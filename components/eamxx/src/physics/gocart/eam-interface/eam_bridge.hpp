#ifndef GOCART_EAM_BRIDGE_HPP
#define GOCART_EAM_BRIDGE_HPP

#include <string>
#include <cstring>

namespace scream::gocart::eam_bridge
{
    // Initialization of shared EAM datastructures 
    extern "C" void eam_bridge_init_c2f(int         ncol,
                                        int         nlev,
                                        int         nlon,
                                        int         nlat,
                                        const char* eam_nlfile,
                                        int         nl_len);

    // Number of gas tracer species
    extern "C" int gas_pcnst();

    // Get species index
    extern "C" int get_spc_ndx_c2f(const char* spc_name,
                                   int         len);

    // Get name of species at index
    extern "C" void solsym_c2f(int index, char* spc_name);

    // Convert mass mixing ratio to volume mixing ratio
    extern "C" void mmr2vmr_c2f(const double* mmr, double* vmr, const double* mbar, int ncol);

    // Convert volume mixing ratio to mass mixing ratio
    extern "C" void vmr2mmr_c2f(const double* vmr, double* mmr, const double* mbar, int ncol);


    /*---------------------------------------------------------------------------------------------
     * void init_eam_bridge(...)
     *
     * C-side initialization of EAM bridge with cleaner call signature.
     *-------------------------------------------------------------------------------------------*/
    void eam_bridge_init(int         ncol,
                         int         nlev,
                         int         nlon,
                         int         nlat,        
                         std::string eam_nlfile)
    {eam_bridge_init_c2f(ncol, nlev, nlon, nlat, eam_nlfile.c_str(), eam_nlfile.length());}




    /*---------------------------------------------------------------------------------------------
     * int get_spc_ndx(...)
     *
     * Gets species index within contiguous array by name
     * 
     * Note: This function assumes C-style 0-indexing convention while the external Fortran
     *       routine expects 1-indexing convention
     *-------------------------------------------------------------------------------------------*/
    int get_spc_ndx(std::string spc_name)
    {return get_spc_ndx_c2f(spc_name.c_str(), spc_name.length()) - 1;}




    /*---------------------------------------------------------------------------------------------
     * std::string gas_species_name(...)
     *
     * EAMxx-side bridge to `solsym` list of species names. The Fortran can return a C-style
     * string but this should be converted to a C++ string before being used.
     * 
     * Note: This function assumes C-style 0-indexing convention while the external Fortran
     *       routine expects 1-indexing convention
     *-------------------------------------------------------------------------------------------*/
    std::string gas_species_name(int index)
    {
        char spc_name[18];
        solsym_c2f(index+1, spc_name);
        return std::string(spc_name);
    }


    /*---------------------------------------------------------------------------------------------
     * void mmr2vmr(...)
     *
     * EAMxx-side bridge to EAM mass mixing ratio to volume mixing ratio conversion. 
     * Note: This additional function is not particularly necessary as it shadows the external
     *       Fortran subroutine, but does provide slightly more succinct and clearer naming.
     *-------------------------------------------------------------------------------------------*/
    void mmr2vmr(const double* mmr, double* vmr, const double* mbar, int ncol)
    {mmr2vmr_c2f(mmr,vmr,mbar,ncol);}

    /*---------------------------------------------------------------------------------------------
     * void vmr2mmr(...)
     *
     * EAMxx-side bridge to EAM volume mixing ratio to mass mixing ratio conversion. 
     * Note: This additional function is not particularly necessary as it shadows the external
     *       Fortran subroutine, but does provide slightly more succinct and clearer naming.
     *-------------------------------------------------------------------------------------------*/
    void vmr2mmr(const double* vmr, double* mmr, const double* mbar, int ncol)
    {vmr2mmr_c2f(vmr,mmr,mbar,ncol);}

} // namespace gocart::eam_bridge

#endif //GOCART_EAM_BRIDGE_HPP