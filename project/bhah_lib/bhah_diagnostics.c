#include "BHaH_defines.h"
#include "BHaH_function_prototypes.h"
/**
 * Output BH@H diagnostics
 */
void bhah_diagnostics(BHaH_struct *bhah_struct) { diagnostics(bhah_struct->commondata, bhah_struct->griddata); }
