!
! Package supplying the  thermodynamic properties of the Lennard-Jones
! fluid (Credit to http://www.sklogwiki.org).
!
! "The Lennard-Jones fluid: an accurate analytic and
! theoretically-based equation of state", Fluid Phase Equilibria 100
! pp. 1-34 (1994)] Jirí Kolafa, Ivo Nezbeda
! http://dx.doi.org/10.1016/0378-3812(94)80001-4
!
! ALJ(T,rho)...Helmholtz free energy (including the ideal term)
! PLJ(T,rho)...Pressure
! ULJ(T,rho)...Internal energy
!
! A few values for confidence, cp. Fig. 3 of the original paper:
!
!   ULJ (1.3, 0.2) = -1.5371511277396868
!   ULJ (1.3, 0.4) = -2.806016725054152
!   PLJ (1.3, 0.2) = 0.12147936810578708
!   PLJ (1.3, 0.4) = 0.11253334756806482
!
module eos
  use iso_c_binding, only: rk => c_double
  implicit none
  private
  public :: alj, plj, ulj

  real (rk), parameter :: pi = 4 * atan (1.0_rk)

contains

  real (rk) function ALJ (T, rho) bind (c, name="rism_alj")
    !
    ! Helmholtz free energy (including the ideal term)
    !
    implicit none
    real (rk), intent (in), value :: T, rho
    ! *** end of interface ***

    real (rk) :: eta

    eta = PI / 6. * rho * (dC (T))**3
    ALJ =  (log (rho) + betaAHS (eta) &
         + rho * BC (T) / exp (gammaBH (T) * rho**2)) * T &
         + DALJ (T, rho)
  end function ALJ

  real (rk) function ALJres (T, rho)
    !
    ! Helmholtz free energy (without ideal term)
    !
    implicit none
    real (rk), intent (in) :: T, rho
    ! *** end of interface ***

    real (rk) :: eta

    eta = PI / 6. * rho * (dC (T))**3
    ALJres = (betaAHS (eta) &
         + rho * BC (T) / exp (gammaBH (T) * rho**2)) * T &
         + DALJ (T, rho)
  end function ALJres

  real (rk) function PLJ (T, rho) bind (c, name="rism_plj")
    !
    ! Pressure
    !
    implicit none
    real (rk), intent (in), value :: T, rho
    ! *** end of interface ***

    real (rk) :: eta, sum

    eta=PI/6. *rho*(dC(T))**3
    sum=((2.01546797*2+rho*( &
         (-28.17881636)*3+rho*( &
         28.28313847*4+rho* &
         (-10.42402873)*5))) &
         +((-19.58371655)*2+rho*( &
         +75.62340289*3+rho*( &
         (-120.70586598)*4+rho*( &
         +93.92740328*5+rho* &
         (-27.37737354) * 6)))) / sqrt (T) &
         + ((29.34470520*2+rho*( &
         (-112.35356937)*3+rho*( &
         +170.64908980*4+rho*( &
         (-123.06669187)*5+rho* &
         34.42288969*6))))+ &
         ((-13.37031968)*2+rho*( &
         65.38059570*3+rho*( &
         (-115.09233113)*4+rho*( &
         88.91973082*5+rho* &
         (-25.62099890)*6))))/T)/T)*rho**2
    PLJ = ((zHS(eta) &
         + BC(T)/exp(gammaBH(T)*rho**2) &
         *rho*(1-2*gammaBH(T)*rho**2))*T &
         +sum )*rho
  end function PLJ

  real (rk) function ULJ (T, rho) bind (c, name="rism_ulj")
    !
    ! Internal energy
    !
    implicit none
    real (rk), intent (in), value:: T, rho
    ! *** end of interface ***

    real (rk) :: eta, sum, dBHdT, dB2BHdT, d

    dBHdT=dCdT(T)
    dB2BHdT=BCdT(T)
    d=dC(T)
    eta=PI/6. *rho*d**3
    sum= ((2.01546797+rho*( &
         (-28.17881636)+rho*( &
         +28.28313847+rho* &
         (-10.42402873)))) &
         + (-19.58371655*1.5+rho*( &
         75.62340289*1.5+rho*( &
         (-120.70586598)*1.5+rho*( &
         93.92740328*1.5+rho* &
         (-27.37737354) * 1.5)))) / sqrt (T) &
         + ((29.34470520*2+rho*( &
         -112.35356937*2+rho*( &
         170.64908980*2+rho*( &
         -123.06669187*2+rho* &
         34.42288969*2)))) + &
         (-13.37031968*3+rho*( &
         65.38059570*3+rho*( &
         -115.09233113*3+rho*( &
         88.91973082*3+rho* &
         (-25.62099890)*3))))/T)/T) *rho*rho
    ULJ = 3*(zHS(eta)-1)*dBHdT/d &
         +rho*dB2BHdT/exp(gammaBH(T)*rho**2) +sum
  end function ULJ

  real (rk) function zHS (eta)
    implicit none
    real (rk), intent (in) :: eta
    ! *** end of interface ***

    zHS = (1+eta*(1+eta*(1-eta/1.5*(1+eta)))) / (1-eta)**3
  end function zHS

  real (rk) function betaAHS (eta)
    implicit none
    real (rk), intent (in) :: eta
    ! *** end of interface ***

    betaAHS = log (1 - eta) / 0.6 &
         + eta*( (4.0/6*eta-33.0/6)*eta+34.0/6 ) /(1.-eta)**2
  end function betaAHS

  real (rk) function dLJ (T)
    ! hBH diameter
    implicit none
    real (rk), intent (in) :: T
    ! *** end of interface ***

    real (rk) IST
    isT = 1 / sqrt (T)
    dLJ = ((( 0.011117524191338 *isT-0.076383859168060) &
         *isT)*isT+0.000693129033539)/isT+1.080142247540047 &
         + 0.127841935018828 * log (isT)
  end function dLJ

  real (rk) function dC (T)
    implicit none
    real (rk), intent (in) :: T
    ! *** end of interface ***

    real (rk) :: sT

    sT = sqrt (T)
    dC = -0.063920968 * log (T) + 0.011117524 / T &
         -0.076383859/sT+1.080142248+0.000693129*sT
  end function dC

  real (rk) function dCdT (T)
    implicit none
    real (rk), intent (in) :: T
    ! *** end of interface ***

    real (rk) :: sT

    sT = sqrt (T)
    dCdT =   0.063920968*T+0.011117524+(-0.5*0.076383859 &
         -0.5*0.000693129*T)*sT
  end function dCdT

  real (rk) function BC (T)
    implicit none
    real (rk), intent (in) :: T
    ! *** end of interface ***

    real (rk) :: isT

    isT = 1 / sqrt (T)
    BC = (((((-0.58544978*isT+0.43102052)*isT &
         +.87361369)*isT-4.13749995)*isT+2.90616279)*isT &
         -7.02181962)/T+0.02459877
  end function BC

  real (rk) function BCdT (T)
    implicit none
    real (rk), intent (in) :: T
    ! *** end of interface ***

    real (rk) :: isT

    isT = 1 / sqrt (T)
    BCdT = ((((-0.58544978*3.5*isT+0.43102052*3)*isT &
         +0.87361369*2.5)*isT-4.13749995*2)*isT &
         +2.90616279*1.5)*isT-7.02181962
  end function BCdT

  real (rk) function gammaBH (T)
    implicit none
    real (rk), intent (in) :: T
    ! *** end of interface ***

    gammaBH=1.92907278
  end function gammaBH

  real (rk) function DALJ (T, rho)
    implicit none
    real (rk), intent (in) :: T, rho
    ! *** end of interface ***

    DALJ = ((+2.01546797+rho*(-28.17881636 &
         +rho*(+28.28313847+rho*(-10.42402873)))) &
         +(-19.58371655+rho*(75.62340289+rho*((-120.70586598) &
         + rho * (93.92740328 + rho * (-27.37737354))))) / sqrt (T) &
         + ( (29.34470520+rho*((-112.35356937) &
         +rho*(+170.64908980+rho*((-123.06669187) &
         +rho*34.42288969)))) &
         +(-13.37031968+rho*(65.38059570+ &
         rho*((-115.09233113)+rho*(88.91973082 &
         +rho* (-25.62099890)))))/T)/T) *rho*rho
  end function DALJ


  function eos_alj (T, rho) result (A) bind (c)
    use lisp, only: obj, flonum
    implicit none
    type (obj), intent (in), value :: T, rho
    type (obj) :: A
    ! *** end of interface ***

    A = flonum (alj (flonum (T), flonum (rho)))
  end function eos_alj


  function eos_alj_res (T, rho) result (A) bind (c)
    use lisp, only: obj, flonum
    implicit none
    type (obj), intent (in), value :: T, rho
    type (obj) :: A
    ! *** end of interface ***

    A = flonum (aljres (flonum (T), flonum (rho)))
  end function eos_alj_res


  function eos_plj (T, rho) result (P) bind (c)
    use lisp, only: obj, flonum
    implicit none
    type (obj), intent (in), value :: T, rho
    type (obj) :: P
    ! *** end of interface ***

    P = flonum (plj (flonum (T), flonum (rho)))
  end function eos_plj


  function eos_ulj (T, rho) result (U) bind (c)
    use lisp, only: obj, flonum
    implicit none
    type (obj), intent (in), value :: T, rho
    type (obj) :: U
    ! *** end of interface ***

    U = flonum (ulj (flonum (T), flonum (rho)))
  end function eos_ulj

end module eos

