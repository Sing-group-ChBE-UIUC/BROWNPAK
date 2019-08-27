module m_globals
!! Global variables, primarily dealing with system configuration
!! and simulation execution.

use m_precision
use m_trajectory
use m_simbox

implicit none

!General
integer, parameter :: mxrdln = 1024
    !! Maximum length of character string for input line buffer.
integer, parameter :: mxparam = 12
    !! Maximum number of parameters for bonds, angles, etc.

!Simulation box
type(smbx_t) :: simbox
    !! Simulation box.
integer :: imcon = 0
    !! Flag specifying boundary conditions on the simulation box.  
    !!
    !! * *imcon = 0*: Unbounded domain. While not explicitly enforced this 
    !!    is useful only for a single molecule. The dynamics is performed in the 
    !!    c.o.m. frame of reference.
    !!
    !! * *imcon = 1*: PBC along *x*, *y*, & *z*. There are no restrictions on
    !!    the number of molecules.

!Particle configuration: Atoms
integer :: num_atom_types = 0
    !! Number of *atom_type*s
character(len=8), dimension(:), allocatable :: atom_names
    !! (*num_atom_types*,) array. Name of atoms of each type.
integer, dimension(:), allocatable :: atom_styles
    !! (*num_atom_types*,) array. Style of atoms of each type.
real(rp), dimension(:), allocatable :: atom_mass
    !! (*num_atom_types*,) array. Mass of atoms of each type.
integer :: mpcd_avnc = 0
    !! Average number of MPCD atoms per collision cell
integer :: num_mpcd_atoms = 0
    !!  Number of MPCD atoms
integer :: num_atoms = 0
    !!  Number of atoms (excluding MPCD atoms)
integer :: num_atoms_tot = 0
    !!  Total number of atoms (includes MPCD atoms)
integer , dimension(:), allocatable:: atoms
    !! (*num_atoms*,) array. 
    !!
    !! For atom *i*, its type *at = atoms(i)*, with style
    !! *atom_styles(at)*, name *atom_names(at)*, mass *atom_mass(at)*, charge
    !! *charge(i)*, position *coordinates(:,i)*, velocity *velocities(:,i)*,
    !! orientation (if the style requires) *orientations(:,i)*. The force acting on atom *i*
    !! is *forces(:,i)*.
real(rp), dimension(:), allocatable :: charge
    !! (*num_atoms*,) array.
real(rp), dimension(:,:), allocatable :: coordinates
    !!  (3, *num_atoms_tot*) array
real(rp), dimension(:,:), allocatable :: orientation
    !!  (4, *num_atoms*) array
real(rp), dimension(:,:), allocatable, target :: velocities
    !! (3, *num_atoms_tot*) array. The first *num_atoms* columns stores
    !! velocities of non-MPCD atoms, the rest, i.e *num_atoms+1* to
    !! *num_atoms_tot*, store velocities of MPCD atoms.
real(rp), dimension(:,:), allocatable, target :: forces
    !!  (3, *num_atoms_tot*) array

!Particle configuration: Bonds
integer :: num_bond_types = 0
    !!  Number of *bond_type*s
integer, dimension(:), allocatable :: bond_styles
    !!  (*num_bond_types*,) array.
real(rp), dimension(:,:), allocatable :: bond_params
    !!  (*mxparam*,*num_bond_types*) array.
integer :: num_bonds = 0
    !!  Total number of bonds.
integer, dimension(:,:), allocatable :: bonds
    !! (3, *num_bonds*) array. Bond *i* is of type *bt = bonds(1,i)*,  directed from
    !! atom *bonds(2,i)* to *bonds(3,i)*. Its style is *bond_styles(bt)* with
    !! parameters *bond_params(:,bt)*.

!Particle configuration: Angles
integer :: num_angle_types = 0
    !!  Number of *angle_type*s
integer, dimension(:), allocatable :: angle_styles
    !!  (*num_angle_types*,) array
real(rp), dimension(:,:), allocatable :: angle_params
    !!  (*mxparam*, *num_angle_types*) array
integer :: num_angles = 0
    !!  Number of angles
integer, dimension(:,:), allocatable :: angles
    !! (4, *num_angles*) array. Angle *i* is of type *ant = angles(1,i)*, incident
    !! to atoms *angles(2,i)*, *angles(3,i)*, and *angles(4,i)*. Its style is
    !! *angle_styles(ant)* with parameters *angle_params(:,ant)*.

!Particle configuration: Dihedrals
integer :: num_dihedral_types = 0
    !!  Number of *dihedral_type*s
integer, dimension(:), allocatable :: dihedral_styles
    !!  (*num_dihedral_types*,) array
real(rp), dimension(:,:), allocatable :: dihedral_params
    !!  (*mxparam*, *num_dihedral_types*) array
integer :: num_dihedrals = 0
    !!  Number of dihedrals
integer, dimension(:,:), allocatable :: dihedrals
    !! (5, *num_dihedrals*) array. Dihedral *i* is of type *dt = dihedrals(1,i)*, incident
    !! to atoms *dihedrals(2,i)*, *dihedrals(3,i)*, *dihedrals(4,i)*, and *dihedrals(5,i)*.
    !! Its style is *dihedral_styles(dt)* with parameters *dihedral_params(:,dt)*.

!Particle configuration: Branches
integer :: num_branches = 0
    !! Total number of branches (including the backbone)
integer , dimension(:,:), allocatable:: branches
    !! (3,*num_branches*) array. Branch *i* is tethered to atom *branches(1,i)*,
    !! contains *branches(2,i)* atoms, with the beginning atom index *branches(3,i)*.

!Particle configuration: Molecules
integer :: num_molecule_types = 0
    !!  Number of *molecule_type*s
character(len=8), dimension(:), allocatable :: molecule_names
    !! (*num_molecule_types*,) array
integer, dimension(:), allocatable :: molecule_pop
    !! (*num_molecule_types*,) array
integer :: num_molecules = 0
    !!  Number of molecules
integer, dimension(:,:), allocatable :: molecules
    !! (9,*num_molecules*) array. For molecule *i*, its type *mt = molecules(1,i)*, 
    !! containing *molecules(2,i)* atoms with beginning index *molecules(3,i)*,
    !! *molecules(4,i)* bonds with beginning index *molecules(5,i)*,
    !! *molecules(6,i)* angles with beginning index *molecules(7,i)*, and
    !! *molecules(8,i)* dihedrals with beginning index *molecules(9,i)*.
real(rp), dimension(3) :: molecule_com = 0.0_rp
    !! Center of mass of the molecule. This is used only when imcon == 0, i.e.
    !! for a single molecule without periodic boundaries.

!Particle configuration: Tethers
integer :: num_tether_types = 0
    !!  Number of *tether_type*s
integer, dimension(:), allocatable :: tether_styles
    !!  (*num_tether_types*,) array
real(rp), dimension(:,:), allocatable :: tether_params
    !!  (*mxparam*, *num_tether_types*) array
integer :: num_tethers = 0
    !!  Number of tethers
integer, dimension(:,:), allocatable :: tethers
    !! (2, *num_tethers*) array. Tether *i* is of type *tt = tethers(1,i)*, tethering
    !! atom *tethers(2,i)* to a point *tether_points(:,i)*.
    !! Its style is *tether_styles(tt)* with parameters *tether_params(:,tt)*.
real(rp), dimension(:,:), allocatable :: tether_points
    !!  (3, *num_tethers*) array

!Particle configuration: VDW (pair) interactions
integer :: num_vdw_types = 0
    !!  Number of *vdw_type*s
integer, dimension(:), allocatable :: vdw_styles
    !!  (*num_vdw_types*,) array
real(rp), dimension(:,:), allocatable :: vdw_params
    !!  (*mxparam*, *num_vdw_types*) array
integer, dimension(:,:), allocatable :: vdw_pairs
    !!  (2, *num_vdw_types*) array. Stores atom type of interacting pairs, such
    !! that at_i >= at_j.

!Particle configuration: External force field
integer :: num_externals = 0
    !!  Number of external fields
integer, dimension(:), allocatable :: external_styles
    !!  (*num_external*,) array
real(rp), dimension(:,:), allocatable :: external_params
    !!  (*mxparam*, *num_external*) array

!Particle configuration: Flow field
integer :: flow_style = 0
real(rp), dimension(:), allocatable :: flow_params
    !!  (*mxparam*,) array

!End of configuration related globals

!Variables controlling runtime behavior
integer :: sim_style = 1
    !! '0': Structure relaxation, '1' : Brownian dynamics, '2' : MPCD
logical :: leql = .true.
    !! Is the system equilibrating? {T, F}
logical :: lrevive = .false.
    !! Is this a restart run? {T, F}.
real(rp) :: tim_stp
    !! BD/MPCD time step size
integer(ip_long) :: nts
    !! Counter for BD/MPCD time steps
integer(ip) :: nts_md = 1
    !! Number of MD steps per MPCD step
integer(ip_long) :: nts_log = 1
    !! Interval for logging (in BD time steps)
integer(ip_long) :: nts_dump = 1
    !! Interval for dumping to revive file (in BD time steps)
integer(ip_long) :: nts_samp = 1
    !! Interval for sampling statistics (in BD time steps)
integer(ip_long) :: nts_eql = 0
    !! Number of BD time steps for equilibration
integer(ip_long) :: nts_eql_samp = 1
    !! Sampling interval during equilibration (in BD time steps)
integer(ip_long) :: nts_sim = 0
    !! Total number of BD time steps in production run
logical :: use_verlet_tab = .false.
    !! Use Verlet neighbor table? {T, F}
real(rp) :: rcutoff = 0.0_rp
    !! Cut off for short-ranged interaction. Also used as the radius of the
    !! skin sphere for short-ranged forces
real(rp) :: tskin = 0.0_rp
    !!  Thickness of the skin sphere (same for all)
logical :: use_cell_list = .false.
    !! Use cell list for short-range interactions? {T, F}
!End of Variables controlling runtime behavior

!Variables for I/O
character(len=:), allocatable :: fn_cfg
    !! Name of the file containing the initial configuration
character(len=:), allocatable :: fn_revive
    !! Name of the revive file
character(len=:), allocatable :: fn_traj
    !! Name of the trajectory file
character(len=:), allocatable :: fn_stats
    !! Name of the statistics file
type(trajectory_t) :: traj
    !! Trajectory object
character(len=8) :: job_tag = ''
    !! A tag useful for array jobs, available only as a command line argument
logical :: read_seed = .false.
    !! {T, F}
    !!
    !! Whether to initialize the random number generator by reading a seed from
    !! a file. If `read_seed` == T, the seed will be read from a file
    !! 'random_seed.txt'
logical :: write_seed = .false.
    !! {T, F}
    !!
    !!  Whether to write the random number generator seed. If
    !!  `write_seed` == T the seed will be written to a file named
    !!  'random_seed.txt'
logical :: write_eql_stats = .false.
    !! During equilibration, should the statistics file be written? {T, F}
logical :: write_traj = .false.
    !! Should the trajectory be written to file? {T, F}
integer, dimension(4) :: traj_frmcmp = 0
    !! Control for what gets written to a trajectory frame.
    !! 1: coordinates; 2: velocities; 3: forces; 4: charge
logical :: traj_wmpcd = .false.
    !! Depending on the values in *traj_frmcmp*, whether the corresponding
    !! quantities for the MPCD atoms are written as well.
!End of variables for I/O

!Miscellaneous variables
real(rp), dimension(3,3) :: stress = 0.0_rp
    !! Stress tensor due to non-MPCD atoms
real(rp), dimension(3,3) :: stress_slvnt = 0.0_rp
    !! Stress tensor due to MPCD atoms (solvent)
real(rp) :: energy_kin = 0.0_rp
    !! Kinetic energy
real(rp) :: energy_bond = 0.0_rp
    !! Bond energy
real(rp) :: energy_angle = 0.0_rp
    !! Angle energy
real(rp) :: energy_dihedral = 0.0_rp
    !! Dihedral energy
real(rp) :: energy_vdw = 0.0_rp
    !! vdw interaction energy
real(rp) :: energy_tether = 0.0_rp
    !! Energetic contribution from tethers
real(rp) :: energy_external = 0.0_rp
    !! Energetic contribution from external fields
real(rp) :: energy_tot = 0.0_rp
    !! Total energy
real(rp) :: bndlen = 0.0_rp
    !! Average bond length
real(rp) :: bndlen_min = 0.0_rp
    !! Minimum bond length
real(rp) :: bndlen_max = 0.0_rp
    !! Maximum bond length
integer :: excluded_atoms = 0
    !! Control for excluded atoms in vdw calculation.
    !! 0: No exclusion, 1: exclude 1-ring bonded neighbors,
    !! 2: exclude 2-ring bonded neighbors, 3: exclude 3-ring bonded neighbors.
logical :: lvdw = .true.
    !! Whether to calculate VDW interactions
character(len=4) :: mob_fctr
    !! Factorization method for mobility matrix {'CHOL', 'KRYL'}
logical :: lhdia = .true.
    !! Whether to include hydrodynamic interactions in BD
logical :: lelectrostatics = .false.
    !! Whether to calculate electrostatic interactions

!End of miscellaneous variables

end module m_globals
