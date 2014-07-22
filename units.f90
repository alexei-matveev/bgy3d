module units
  !
  ! Copyright (c) 2013 Alexei Matveev
  !
  use kinds, only: rk
  implicit none
  private

  real (rk), parameter, public :: pi = 4 * atan (1.0_rk)

  !
  ! Working units *are*  angstroms and kcals.  Still when  you wish to
  ! make dimensions  explicit (and  you should) use  combinations like
  ! 0.5 * angstrom, or 1.2 * angstrom**(-1), or 332 kcal * angstrom:
  !
  real (rk), parameter, public :: ANGSTROM = 1
  real (rk), parameter, public :: KCAL = 1

  ! International Steam Table calorie here, see also bgy3d.h:
  real (rk), parameter, public :: MOL = 6.02214129d23
  real (rk), parameter, public :: KJOULE = KCAL / 4.1868d0 ! only for output
  real (rk), parameter, public :: METER = 1.0d10 * ANGSTROM
  real (rk), parameter, public :: KPASCAL = KJOULE / METER**3

  !
  ! The interaction energy of two unit charges separated by 1 A is
  !                  -1
  !   E = 1 * 1 / 1 A   = 0.529 au = 332 kcal [/ mol]
  !
  ! The next parameter appears to have the meaning of this interaction
  ! energy of such two unit charges:
  !
  !   EPSILON0INV = 1 / ε₀
  !
  ! and  is  used  to  covert electrostatic  interaction  energies  to
  ! working units.   It has  to be consistent  with other  force field
  ! parameters defined in bgy3d-solvents.h, notably with Lennard-Jones
  ! parameters σ and ε (FIXME: so maybe it was not a good idea to move
  ! it here). These are the original comments:
  !
  !   You have: e^2/4/pi/epsilon0/angstrom, you want: kcal/avogadro/mol
  !
  !   => 331.84164
  !
  ! Again,  the code  appears  to use  the  International Steam  Table
  ! calorie to define 1/ε₀. Keep this in sync with bgy3d.h.  Fun fact:
  ! some historical  works, e.g.  J.   Aaqvist, J.  Phys.   Chem.  94,
  ! 8021 (1990), even use the exact 332 for electrostatics:
  !
  real (rk), parameter, public :: EPSILON0INV = 331.84164d0 * KCAL * ANGSTROM

  !
  ! Inverse range parameter for  separation of the Coulomb into short-
  ! and long range components. The  inverse of this number 1/α has the
  ! dimension of length. FIXME: with the hardwired parameter like this
  ! the  short-range  Coulomb is  effectively  killed alltogether  for
  ! water-like particles  --- a water-like  particle has a  typical LJ
  ! radius σ = 3.16 A. There are no visible changes in the RDF with or
  ! without short range Coulomb and the charges of the order ±1.
  !
  real (rk), parameter, public :: ALPHA = 1.2d0 * ANGSTROM**(-1)

end module units
