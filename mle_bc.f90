!!This program solves the modified Lane-Emden equation, as discussed in Ohlmann+16c and Ohlmann16_phdthesis
!!
!!Unlike mle.f90, this program solves the modified LE eq subject to the boundary conditions at the softening radius r = h
!!The BCs lead to the following: 
!! (1) (rho|r=h)/(theta|xi=xi_h)^n - rho_0 = 0
! (2) nnn*(rho|r=h)*(dtheta/dxi|xi=xi_h)/((drho/dr|r=h)*(theta|xi=xi_h)) - alpha = 0
! where xi_h = h/alpha, and where rho_0 and alpha are constants set by the BCs
!
!The parameters nnn, h, rhobar, rho|r=h, and drho/dr|r=h must be specified
!For a chosen h, the parameters rhobar, rho|r=h, and drho/dr|r=h must be determined using the MESA stellar density profile
!The latter procedure is carried out using the code ~/MESA/Projects/make_agb_001/p_bc.pro
!This program p_bc.pro writes the quantities h, m_c, rho|r=h, and drho/dr|r=h to the file ~/Common_envelope/Lane-emden/bc.out
!This boundary conditions file is then read in by this program
!
!Procedure is to solve the equation up to xi = xi_h while iterating over (guessed) values of alpha
!After each iteration, use eq (1) to solve for rho_0. Then check if eq (2) is satisfied.
!If eq (2) is satisfied to the desired precision LHS < eps then exit loop and use the guessed value of alpha
!Otherwise, change value of alpha slightly and repeat. Iterate until converges to a solution.
!*****************************************************
!Fundamental constants and units
module const
  implicit none
!
  real, parameter :: pi=      3.14156295358
  real, parameter :: s_Myr=   1.d6*365.25*24*3600
  real, parameter :: s_Gyr=   1.d9*365.25*24*3600
  real, parameter :: Myr_Gyr= 1.d3
  real, parameter :: g_Msun=  1.989d33
  real, parameter :: cm_Rsun= 6.957d10
  real, parameter :: cgs_G=   6.6743d-8
end module const
!*****************************************************
!Set up the time-stepping
module num_params
  implicit none
!
  integer, parameter :: Nvar=   2 !number of variables 
  integer, parameter :: Nmax=   100000 !max number of points
  integer, parameter :: ialphamax= 100000  !max iterations to find alpha to match BCs
  integer, parameter :: irho0max= 100000  !max iterations to find rho0 to match BCs
  real, parameter ::    eps_alpha=  0.001  !sets tolerance of convergence condition on alpha--smaller eps is more precise
  real, parameter ::    eps2=       0.0005  !sets fractional change in alpha after each iteration
  real, parameter ::    eps_rho0=   0.001  !sets tolerance of convergence condition on rho0--smaller eps is more precise
  !real, parameter ::    dt=     0.001 !0.0001
  real, parameter ::    dt=     1.d-4 !0.0001
  !real ::               tmin=   1.d-6  !must start at finite value to avoid singularity
  real ::               tmin=   1.d-10  !must start at finite value to avoid singularity
  real ::               t  !time variable (equal to xi)
  real ::               first  !for Runge-Kutta routine
end module num_params
!*****************************************************
!Physical parameters
module phys_params
  implicit none
!
  real, parameter ::  nnn= 3. !polytropic index
  ! the following parameters are calculated by MESA/Projects/make_agb_001/p_m.pro
  real, parameter :: alpha_init= 5.87!8.66!8.47!15. !parameter, units of Rsun
  real, parameter :: rho0_init= 0.000514!0.00438!0.0045 !central density, units g/cm^3 
  real :: alpha= alpha_init
  real :: rho0= rho0_init
  real, parameter :: m_c_over_m_h= 0.9446325
end module phys_params
!*****************************************************
!Boundary conditions
module bc
  use const
  use phys_params
!
  implicit none
!
  real :: full_h  !transition radius in units of Rsun (chosen)
  real :: h  !softening length
  real :: m_h  !mass at transition radius in units of Msun (from MESA)
  real :: rho_h  !density at transition radius in units g/cm^3 (from MESA)
  real :: drhodr_h  !density gradient at transition radius in units of g/cm^3/Rsun (from differentiating and smoothing MESA rho)
  real :: rhobar  !mean density of central point particle smoothed over sphere of radius h in units g/cm^3
  real :: xi_h  !dimensionless softening length h/alpha
  real :: m_c  !mass of point particle in units of Msun
contains
  subroutine set_bc
!
    open(11,file="bc.out",status="old")
    read(11,*) full_h, m_h, rho_h, drhodr_h
    close(11)
!
    h = full_h / 2.

    m_c= m_c_over_m_h*m_h
    rhobar= 3.*m_c*g_Msun/4./pi/full_h**3./cm_Rsun**3  !average density within r < h of central gravitating point mass in cgs units
    xi_h= full_h/alpha !dimenless softening length, = r_cutoff/alpha = h/alpha
!
    print*,'h, m_h, rho_h, drhodr_h, rhobar, m_c',h, m_h, rho_h, drhodr_h, rhobar, m_c
!
  end subroutine set_bc
end module bc
!*****************************************************
!Partial differential equations
module equ
  use phys_params
  use bc
  use num_params
!  
  implicit none
!
  real :: xi, theta, eta, dthetadxi
  real :: u, tu, chi, dchidu
  real, dimension(Nvar) :: f
!
contains
  subroutine pde(f,dfdt)
!
  real, dimension(Nvar) :: f, dfdt 
  intent(inout) :: f
  intent(out) :: dfdt
!
! LIST OF VARIABLE NAMES FOR f ARRAY
!
! UNDER FOSA				
! t    = xi
! f(1) = eta				
! f(2) = theta				
!
!use explicit names
  xi=    t
  theta= f(1)
  eta=   f(2)
  dthetadxi= -eta/xi**2
!
  u= xi*alpha/h
!
  if (u >= 0 .and. u < 1.) then
    chi    = (15.*u**3 - 36.*u**2 + 40.)/30.
    dchidu = (15.*u**2 - 24.*u)/10.
  elseif (u >= 1. .and. u < 2.) then 
    chi    = (-5.*u**5 + 36.*u**4 - 90.*u**3 + 80.*u**2 - 2.*u**(-1))/(30.*u**2)
    dchidu = (-5.*u**6 + 24.*u**5 - 30.*u**4 + 2)/(10.*u**4)
  elseif (u >= 2.) then
    chi    = 1./u**3
    dchidu = -3./u**4
  endif


!
  dfdt(1)= -eta/xi**2
  dfdt(2)= xi**2*(theta**nnn +rhobar/rho0*(chi + 1./3.*u*dchidu))
!
  end subroutine pde
end module equ
!*****************************************************
!Timestepping routine (same as Pencil-code)
module timestep 
!
  use num_params
  use equ
!
contains
  subroutine rk(f)
!
  implicit none
!
  real :: gam1, gam2, gam3, zet1, zet2
  real, dimension(Nvar) :: f, dfdt, pdef, ftmp
  real :: ttmp
!
  intent(inout) :: f
!
!  Runge Kutta 3rd order time advance
!  f = f(exact) - 0.0046*dt^3*d3f/dt3
!
    gam1=8./15 !gam1, gam2 AND gam3 ARE THE COEFFICIENTS OF THE TIMESTEPS AT WHICH dfdt IS CALCULATED
    gam2=5./12
    gam3=3./4
    zet1=-17./60
    zet2=-5./12
!
    call pde(f,dfdt)
    pdef=dfdt !pde FUNCTION CALCULATES VALUE OF TIME DERIVATIVE ACCORDING TO P.D.E.s
    f=f+dt*gam1*pdef !FIRST GET INTERMEDIATE VALUE OF f (AT t=t+dt*gam1) USING dfdt AT t_0
    t=t+dt*gam1 !THEN GO TO THAT TIMESTEP
    ftmp=f+dt*zet1*pdef !NOW CALCULATE A TEMPORARY f (AT t=t_0+dt*gam1+dt*zet1) USING dfdt AT t_0
    ttmp=t+dt*zet1 !DEFINE TEMPORARY t AS THAT TIMESTEP (t=t_0+dt*gam1+dt*zet1)
!
    call pde(f,dfdt)
    pdef=dfdt !NOW CALCULATE dfdt AT THE NEW TIMESTEP t_0+dt*gam1
    f=ftmp+dt*gam2*pdef !USE THIS TO GET ANOTHER INTERMEDIATE VALUE OF f (AT t=t_0+dt*gam1+dt*zet1+dt*gam2)
    t=ttmp+dt*gam2 !THEN GO TO THAT TIMESTEP
    ftmp=f+dt*zet2*pdef !NOW CALCULATE A TEMPORARY f (AT t=t_0+dt*gam1+dt*zet1+dt*gam2+dt*zet2) USING dfdt AT t=t_0+dt*gam1
    ttmp=t+dt*zet2 !DEFINE TEMPORARY t AS THAT TIMESTEP (t=t_0+dt*gam1+dt*zet1+dt*gam2+dt*zet2)
!
    call pde(f,dfdt)
    pdef=dfdt !CALCULATE dfdt AT THE NEW TIMESTEP t_0+dt*gam1+dt*zet1+dt*gam2
    f=ftmp+dt*gam3*pdef !USE THIS TO GET THE FINAL VALUE OF f (AT t=t_0+dt*gam1+dt*zet1+dt*gam2+dt*zet2+dt*gam3)
    t=ttmp+dt*gam3 !THEN GO TO THAT TIMESTEP
!
    first=first+1. !COUNTS THE NUMBER OF TIMES RUNGA-KUTTA ROUTINE IS EXCUTED

  end subroutine rk
end module timestep
!*****************************************************
!Initialization routine
module start
!
  use const
  use phys_params
  use num_params 
  use equ
!
  implicit none
!
contains
  subroutine init_start
!
!  INITIALIZE TIMESTEPPING
!
     t= tmin
     first= 0.
!
!  INITIAL CONDITIONS
!
    f(1)= 1.
    f(2)= 0.
!
    !print*,'START VALUES:'
    !print*,'h               =',h
    !print*,'rho_h           =',rho_h
    !print*,'drhodr_h        =',drhodr_h
    !print*,'rhobar          =',rhobar
    print*,'rho0            =',rho0
    print*,'alpha           =',alpha
    !print*,'xi_h            =',xi_h
    !print*,''
!
  end subroutine init_start
end module start
!*****************************************************
module advance_t
  use timestep
  use equ
!
  implicit none 
!
  integer :: it
!
contains
  subroutine eval_rk
    open(20,file= 'out/mle_profile.out',status="replace")
    do it=1,Nmax
      call rk(f) !(f,t,dt,first)
      if (isnan(f(1))) then
        print*,'NANs detected, change time step'
        stop
      endif
!
!      write(* ,*) xi, theta, -eta/xi**2
      write(20,*) xi, theta, dthetadxi, chi, dchidu
      if (theta<=0) then 
        print*,'reached negative theta at it = ',it
        print*,'xi=',xi
        close(20)
        stop
      elseif (xi>=xi_h) then
        !print*,'reached xi >= xi_h at it = ',it
        close(20)
        return
      endif
    enddo
    close(20)    
  end subroutine eval_rk
end module advance_t
!*****************************************************
!Check whether boundary conditions are satisfied
module bc_check
  use num_params
  use phys_params
  use bc
  use equ
!
  implicit none 
!
  logical :: alpha_okay, rho0_okay
  real :: alpha_comp, alpha_sign, rho0_comp, rho0_sign
  real :: F1, F2
!
contains
  subroutine check_alpha_satisfied
!
    F1= theta/dthetadxi
    F2= nnn/alpha*rho_h/drhodr_h
!
    !print*,'abs((F2-F1)/F2)       = ',abs((F2-F1)/F2)
    !print*,'abs((F2-F1)/F2) < eps_alpha = ',abs((F2-F1)/F2) < eps_alpha
!
    alpha_comp= abs((F2-F1)/F2)
    alpha_okay= alpha_comp < eps_alpha
    alpha_sign= (F2-F1)/abs(F2-F1)
    print*, alpha_comp
!
  end subroutine check_alpha_satisfied
!
  subroutine check_rho0_satisfied
!
    F1= rho0*theta**nnn
    F2= rho_h
!
    !print*,'abs((F2-F1)/F2)       = ',abs((F2-F1)/F2)
    !print*,'abs((F2-F1)/F2) < eps_alpha = ',abs((F2-F1)/F2) < eps_alpha
!
    rho0_comp= abs((F2-F1)/F2)
    rho0_okay= rho0_comp < eps_rho0
    rho0_sign= (F2-F1)/abs(F2-F1)
!
  end subroutine check_rho0_satisfied
end module bc_check
!*****************************************************
module alpha_iterate
  use phys_params
  use bc
  use start
  use advance_t
  use bc_check
!
  implicit none 
!
contains
  subroutine find_alpha
!
    integer :: ialpha
!
    do ialpha=1,ialphamax
      !print*,''
      !print*,'alhpa iteration level=',ialpha
      !print*,''
      if (ialpha==ialphamax) then
        stop
        print*,'reached ialphamax without converging'
      endif
      !initialize run
      call init_start
      theta= f(1)
      eta=   f(2)
      dthetadxi= -eta/xi**2
!
      !perform time-stepping
      call eval_rk
      xi=    t
      theta= f(1)
      eta=   f(2)
      dthetadxi= -eta/xi**2
      !check to see if BCs are satisfied
      call check_alpha_satisfied
      if (alpha_okay) then
        print*,''
        print*, '******FOUND BC MATCH at ialpha=',ialpha,'******'
        exit
      elseif (alpha_sign == 1) then  !assign new value to alpha
        alpha= alpha -alpha*eps2*abs(F2/F1)
      else
        alpha= alpha +alpha*eps2*abs(F2/F1)
      endif
      xi_h= full_h/alpha  !assign new value to xi_h to reflect new value of alpha
    enddo
  end subroutine find_alpha
end module alpha_iterate
!*****************************************************
!Main subprogram
module mle
  use phys_params
  use bc
  use start
  use advance_t
  use bc_check
  use alpha_iterate
!
  implicit none 
!
  real :: cpu_time_start, cpu_time_finish
contains
  subroutine mle_run
    integer irho0
!
    call cpu_time(cpu_time_start)
!
!  READ IN BOUNDARY CONDITION PARAMETERS
!
    call set_bc
!
    do irho0=1,irho0max
      call find_alpha
      call check_rho0_satisfied
      if (rho0_okay) then
        print*, '******FOUND BC MATCH FOR BOTH ALPHA AND RHO0*******'
        rho0= rho_h/theta**nnn  !set rho0 to value that is consistent with new alpha (and try again)
        exit
      else 
        print*
        print*, '******RHO0 NOT CONVERGED--SET RHO0 = RHO_H/THETA**NNN AND TRY AGAIN******'
        rho0= rho_h/theta**nnn  !set rho0 to value that is consistent with new alpha (and try again)
      endif
    enddo
    print*,'rho0            =',rho0
    print*,'rho_h/theta**nnn=',rho_h/theta**nnn
    print*,'alpha           =',alpha

    !print*,'rho(h)         ~=',rho0*theta**nnn
    !print*,'drhodr(h)      ~=',nnn*rho0/alpha*theta**(nnn-1.)*dthetadxi
    !print*,'p(h)            =',4*pi*cgs_G*alpha**2*cm_Rsun**2/(nnn+1.)*rho0**2*theta**(nnn+1.)  
    !print*,'dpdr(h)         =',4*pi*cgs_G*alpha*cm_Rsun*rho0**2*theta**nnn*dthetadxi
    print*,''
!
    open (12,file= 'out/nt.out',status="replace")
    write(12,*) it
    close(12)
!
!  WRITE INFO TO FILE
!
    open(10,file= 'out/param.out',status="replace")
    write(10,*) dt, nnn, xi_h, xi, rhobar, rho0, alpha, m_c
    close(10)
!
    call cpu_time(cpu_time_finish)
    print*, "simulation time in seconds: ", (cpu_time_finish -cpu_time_start)
  end subroutine mle_run
end module mle
!*****************************************************
program run
  use mle
!
  implicit none 
!
  call mle_run
!
end program run
!*****************************************************
