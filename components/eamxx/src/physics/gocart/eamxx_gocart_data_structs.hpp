/*-------------------------------------------------------------------------------------------------
 * gocart_data_structs.hpp
 *
 * This header file contains definitions of lightweight data structures used internally by the
 * GOCART AtmosphereProcess
 *-----------------------------------------------------------------------------------------------*/
#include <string>
#include <vector>
#include <memory> // shared_ptr

#include "share/eamxx_types.hpp" // Real
#include "share/grid/remap/abstract_remapper.hpp"
#include "share/io/scorpio_input.hpp"

namespace scream::gocart
{
    /*---------------------------------------------------------------------------------------------
     * struct TimeState
     *
     * Organizational structure for date/time used by SurfaceEmission objects during I/O. 
     * Adapted from 'srfEmissTimeState' struct in 'mam' atmosphere process.
     *-------------------------------------------------------------------------------------------*/
    struct TimeState 
    {
        TimeState() = default;
        
        // The current month
        int current_month = -1;
        // Julian Date for the beginning of the month, as defined in
        //           /src/share/util/eamxx_time_stamp.hpp
        // See this file for definition of Julian Date.
        Real t_beg_month;
        // Current simulation Julian Date
        Real t_now;
        // Number of days in the current month, cast as a Real
        Real days_this_month;
    };  // TimeState
}




