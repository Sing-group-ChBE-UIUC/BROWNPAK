module m_kmc
    !! Implements charge hopping via Kinetic Monte Carlo.

use m_precision
use m_constants_math
use m_strings
use m_ran_num
use m_globals
use m_interaction, only: ia_calc_forces
use m_logger, only: logger => master_logger

implicit none

private

public :: kmc_init, kmc_hop, kmc_finish

contains
 
!*******************************************************************************
 
subroutine kmc_init()

    end subroutine

!*******************************************************************************
 
subroutine kmc_finish()

    end subroutine

!*******************************************************************************
 
subroutine kmc_hop()

    end subroutine

!******************************************************************************

end module m_kmc
