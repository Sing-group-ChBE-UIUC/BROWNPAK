module m_control_io
    !! Routines for reading and witing control file.

use m_precision
use m_strings
use m_globals

implicit none

contains

!******************************************************************************

subroutine read_control(fn)
    !! Reads simulation parameters from file

    character(len=*), intent(in) :: fn !! Name of parameters file.
    character(len=:), allocatable :: key
    character(len=:), allocatable :: val
    character(len=mxrdln) :: line
    character(len=1) :: cstr = '#' !Comment string
    integer :: fu
    integer :: ios

    open (newunit=fu, file = fn, action = 'read', status = 'old')
    
    do 
        call readline(fu, line, cstr, ios)
        if (ios /= 0) return

        call str_get_keyval(line, key, val)

        if (key=='sim_style') read(val, *) sim_style
        if (key=='use_verlet_tab') read(val, *) use_verlet_tab
        if (key=='rcutoff') rcutoff = str_to_d(val)
        if (key=='tskin') tskin = str_to_d(val)
        if (key=='use_cell_list') read(val, *) use_cell_list
        if (key=='excluded_atoms') excluded_atoms = str_to_i(val)
        if (key=='lvdw') read(val,*) lvdw
        if (key=='lhdia') read(val,*) lhdia
        if (key=='mob_fctr') read(val,*) mob_fctr
        if (key=='lelectrostatics') read(val,*) lelectrostatics

        if (key=='tim_stp')  tim_stp    = str_to_d(val)
        if (key=='nts_md')   nts_md     = int(str_to_d(val))
        if (key=='nts_log')  nts_log    = int(str_to_d(val), ip_long)
        if (key=='nts_dump') nts_dump   = int(str_to_d(val), ip_long)
        if (key=='nts_samp') nts_samp   = int(str_to_d(val), ip_long)
        if (key=='nts_eql')  nts_eql    = int(str_to_d(val), ip_long)
        if (key=='nts_eql_samp') nts_eql_samp = int(str_to_d(val), ip_long)
        if (key=='nts_sim')  nts_sim  = int(str_to_d(val), ip_long)

        if (key == 'fn_cfg')    fn_cfg    = val
        if (key == 'fn_revive') fn_revive = val
        if (key == 'fn_stats')  fn_stats  = val
        if (key == 'fn_traj')   fn_traj   = val

        if (key == 'lrevive') read(val,*) lrevive
        if (key == 'read_seed') read(val,*) read_seed
        if (key == 'write_seed') read(val,*) write_seed
        if (key == 'write_eql_stats') read(val,*) write_eql_stats
        if (key == 'write_traj') read(val,*) write_traj
        if (key == 'traj_frmcmp') read(val,*) traj_frmcmp
        if (key == 'traj_wmpcd') read(val,*) traj_wmpcd
    end do

    close (fu)

    end subroutine

!******************************************************************************

subroutine write_control(fn)
    !! Write simulation parameters to file

    character(len=*), intent(in) :: fn
        !! File name
    integer :: fu

    open(newunit=fu, file=fn, action='write', status='unknown')

    write(fu, '(a,t20,a)'   ) 'sim_style', sim_style
    write(fu, '(a,t20,l1)'  ) 'use_verlet_tab', use_verlet_tab
    write(fu, '(a,t20,g0.6)') 'rcutoff', rcutoff
    write(fu, '(a,t20,g0.6)') 'tskin', tskin
    write(fu, '(a,t20,l1)'  ) 'use_cell_list', use_cell_list
    write(fu, '(a,t20,i0)'  ) 'excluded_atoms', excluded_atoms
    write(fu, '(a,t20,l1)'  ) 'lvdw', lvdw
    write(fu, '(a,t20,l1)'  ) 'lhdia', lhdia
    write(fu, '(a,t20,a)'   ) 'mob_fctr', mob_fctr
    write(fu, '(a,t20,l1)'  ) 'lelectrostatics', lelectrostatics

    write(fu, *)
    write(fu, '(a,t20,g0.6)') 'tim_stp',    tim_stp
    write(fu, '(a,t20,i0)')   'nts_md',     nts_md
    write(fu, '(a,t20,i0)')   'nts_log',    nts_log
    write(fu, '(a,t20,i0)')   'nts_dump',   nts_dump
    write(fu, '(a,t20,i0)')   'nts_samp',   nts_samp
    write(fu, '(a,t20,i0)')   'nts_eql',    nts_eql
    write(fu, '(a,t20,i0)')   'nts_eql_samp', nts_eql_samp
    write(fu, '(a,t20,i0)')   'nts_sim',  nts_sim

    write(fu, *)
    write(fu, '(a,t20,a)') 'fn_cfg', fn_cfg
    write(fu, '(a,t20,a)') 'fn_revive', fn_revive
    write(fu, '(a,t20,a)') 'fn_stats', fn_stats
    write(fu, '(a,t20,a)') 'fn_traj', fn_traj

    write(fu, *)
    write(fu, '(a,t20,l1)') 'lrevive', lrevive
    write(fu, '(a,t20,l1)') 'read_seed', read_seed
    write(fu, '(a,t20,l1)') 'write_seed', write_seed
    write(fu, '(a,t20,l1)') 'write_eql_stats', write_eql_stats
    write(fu, '(a,t20,l1)') 'write_traj', write_traj
    write(fu, '(a,t20,*(i0,2x))') 'traj_frmcmp', traj_frmcmp
    write(fu, '(a,t20,l1)') 'traj_wmpcd', traj_wmpcd

    close(fu)
    
    end subroutine

!******************************************************************************

end module m_control_io
