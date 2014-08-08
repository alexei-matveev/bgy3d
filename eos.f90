!
! Package supplying the thermodynamic properties of the
! LENNARD-JONES fluid
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
  private
  public :: alj, plj, ulj

  real (rk), parameter :: pi = 4 * atan (1.0_rk)

contains

  real (rk) FUNCTION ALJ (T, rho) bind (c, name="rism_alj")
    ! Helmholtz free energy (including the ideal term)
    !
    use iso_c_binding
    implicit real (rk) (a-h,o-z)
    real (rk), intent (in), value :: T, rho
    ! *** end of interface ***

    eta = PI/6.*rho * (dC(T))**3
    ALJ =  (log (rho) + betaAHS (eta) &
         &  +rho*BC(T)/exp(gammaBH(T)*rho**2))*T &
         &  +DALJ(T,rho)
    RETURN
  END FUNCTION ALJ
  ! Helmholtz free energy (without ideal term)
  DOUBLE PRECISION FUNCTION ALJres(T,rho)
    implicit double precision (a-h,o-z)
    eta = PI/6. *rho*(dC(T))**3
    ALJres = (betaAHS(eta) &
         & +rho*BC(T)/exp(gammaBH(T)*rho**2))*T &
         & +DALJ(T,rho)
    RETURN
  END FUNCTION ALJres
  ! Pressure
  real (rk) FUNCTION PLJ (T, rho) bind (c, name="rism_plj")
    implicit real (rk) (a-h,o-z)
    real (rk), intent (in), value :: T, rho
    ! *** end of interface ***

    eta=PI/6. *rho*(dC(T))**3
    sum=((2.01546797*2+rho*( &
         & (-28.17881636)*3+rho*( &
         & 28.28313847*4+rho* &
         & (-10.42402873)*5))) &
         & +((-19.58371655)*2+rho*( &
         & +75.62340289*3+rho*( &
         & (-120.70586598)*4+rho*( &
         & +93.92740328*5+rho* &
         & (-27.37737354)*6))))/dsqrt(T) &
         & + ((29.34470520*2+rho*( &
         & (-112.35356937)*3+rho*( &
         & +170.64908980*4+rho*( &
         & (-123.06669187)*5+rho* &
         & 34.42288969*6))))+ &
         & ((-13.37031968)*2+rho*( &
         & 65.38059570*3+rho*( &
         & (-115.09233113)*4+rho*( &
         & 88.91973082*5+rho* &
         & (-25.62099890)*6))))/T)/T)*rho**2
    PLJ = ((zHS(eta) &
         &  + BC(T)/exp(gammaBH(T)*rho**2) &
         &  *rho*(1-2*gammaBH(T)*rho**2))*T &
         &  +sum )*rho
    RETURN
  END FUNCTION PLJ
  ! Internal energy
  real (rk) FUNCTION ULJ (T, rho) bind (c, name="rism_ulj")
    implicit real (rk) (a-h,o-z)
    real (rk), intent (in), value:: T, rho
    ! *** end of interface ***

    dBHdT=dCdT(T)
    dB2BHdT=BCdT(T)
    d=dC(T)
    eta=PI/6. *rho*d**3
    sum= ((2.01546797+rho*( &
         & (-28.17881636)+rho*( &
         & +28.28313847+rho* &
         & (-10.42402873)))) &
         & + (-19.58371655*1.5+rho*( &
         & 75.62340289*1.5+rho*( &
         & (-120.70586598)*1.5+rho*( &
         & 93.92740328*1.5+rho* &
         & (-27.37737354)*1.5))))/dsqrt(T) &
         & + ((29.34470520*2+rho*( &
         & -112.35356937*2+rho*( &
         &  170.64908980*2+rho*( &
         & -123.06669187*2+rho* &
         & 34.42288969*2)))) + &
         & (-13.37031968*3+rho*( &
         &  65.38059570*3+rho*( &
         &  -115.09233113*3+rho*( &
         & 88.91973082*3+rho* &
         & (-25.62099890)*3))))/T)/T) *rho*rho
    ULJ = 3*(zHS(eta)-1)*dBHdT/d &
         & +rho*dB2BHdT/exp(gammaBH(T)*rho**2) +sum
    RETURN
  END FUNCTION ULJ
  DOUBLE PRECISION FUNCTION zHS(eta)
    implicit double precision (a-h,o-z)
    zHS = (1+eta*(1+eta*(1-eta/1.5*(1+eta)))) / (1-eta)**3
    RETURN
  END FUNCTION zHS
  DOUBLE PRECISION FUNCTION betaAHS( eta )
    implicit double precision (a-h,o-z)
    betaAHS = log (1 - eta) / 0.6 &
         &  + eta*( (4.0/6*eta-33.0/6)*eta+34.0/6 ) /(1.-eta)**2
    RETURN
  END FUNCTION betaAHS
  ! hBH diameter
  DOUBLE PRECISION FUNCTION dLJ(T)
    implicit double precision (a-h,o-z)
    DOUBLE PRECISION IST
    isT=1/dsqrt(T)
    dLJ = ((( 0.011117524191338 *isT-0.076383859168060) &
         & *isT)*isT+0.000693129033539)/isT+1.080142247540047 &
         & + 0.127841935018828 * log (isT)
    RETURN
  END FUNCTION dLJ
  DOUBLE PRECISION FUNCTION dC(T)
    implicit double precision (a-h,o-z)
    sT=dsqrt(T)
    dC = -0.063920968 * log (T) + 0.011117524 / T &
         &     -0.076383859/sT+1.080142248+0.000693129*sT
    RETURN
  END FUNCTION dC
  DOUBLE PRECISION FUNCTION dCdT( T)
    implicit double precision (a-h,o-z)
    sT=dsqrt(T)
    dCdT =   0.063920968*T+0.011117524+(-0.5*0.076383859 &
         &   -0.5*0.000693129*T)*sT
    RETURN
  END FUNCTION dCdT
  DOUBLE PRECISION FUNCTION BC( T)
    implicit double precision (a-h,o-z)
    DOUBLE PRECISION isT
    isT=1/dsqrt(T)
    BC = (((((-0.58544978*isT+0.43102052)*isT &
         &  +.87361369)*isT-4.13749995)*isT+2.90616279)*isT &
         &  -7.02181962)/T+0.02459877
    RETURN
  END FUNCTION BC
  DOUBLE PRECISION FUNCTION BCdT( T)
    implicit double precision (a-h,o-z)
    DOUBLE PRECISION iST
    isT=1/dsqrt(T)
    BCdT = ((((-0.58544978*3.5*isT+0.43102052*3)*isT &
         &  +0.87361369*2.5)*isT-4.13749995*2)*isT &
         &  +2.90616279*1.5)*isT-7.02181962
    RETURN
  END FUNCTION BCdT
  DOUBLE PRECISION FUNCTION gammaBH(X)
    implicit double precision (a-h,o-z)
    gammaBH=1.92907278
    RETURN
  END FUNCTION gammaBH
  DOUBLE PRECISION FUNCTION DALJ(T,rho)
    implicit double precision (a-h,o-z)
    DALJ = ((+2.01546797+rho*(-28.17881636 &
         & +rho*(+28.28313847+rho*(-10.42402873)))) &
         & +(-19.58371655+rho*(75.62340289+rho*((-120.70586598) &
         & +rho*(93.92740328+rho*(-27.37737354)))))/dsqrt(T) &
         & + ( (29.34470520+rho*((-112.35356937) &
         & +rho*(+170.64908980+rho*((-123.06669187) &
         & +rho*34.42288969)))) &
         & +(-13.37031968+rho*(65.38059570+ &
         & rho*((-115.09233113)+rho*(88.91973082 &
         & +rho* (-25.62099890)))))/T)/T) *rho*rho
    RETURN
  END FUNCTION DALJ
end module eos

