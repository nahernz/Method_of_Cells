!  GMC.f90 
!
!  FUNCTIONS:
!  GMC - Entry point of console application.
!

!****************************************************************************
!
!  PROGRAM: GMC
!
!  PURPOSE:  Entry point for the console application.
!
!****************************************************************************   
PROGRAM GMC

    IMPLICIT NONE

    ! Variables
    ! Fiber Material Properties
    REAL, PARAMETER :: E11f = 231e9
    REAL, PARAMETER :: E22f = 15e9
    REAL, PARAMETER :: v12f = 0.14
    REAL, PARAMETER :: G12f = 24e9
    REAL, PARAMETER :: G23f = 5.01e9
    REAL :: v23f = (E22f/(2*G23f))-1
    ! Matrix Material Properties
    REAL, PARAMETER :: Em = 3e9
    REAL, PARAMETER :: vm = 0.36
    ! Build Material Stiffnesses
    REAL :: delf
    REAL, DIMENSION(6, 6) :: Cf, Cm
    REAL :: c1, c2
    
    delf = E22f*v12f**2+E11f*(v23f-1)
    c1 = Em*vm/((1+vm)*(1-2*vm))
    c2 = Em/(2*(1+vm))
    
    ! Body of GMC
    WRITE (*,'(A)') 'Generalized Method of Cells'
    WRITE (*,'(A)') 'written by Michael Kaplan and Rehan Nawaz'
    
    ! Build the stiffness matricees
    
    WRITE (*,*)
    Cm(1:6,1) = (/ c1+2.0*c2, c1, c1, 0.0, 0.0, 0.0 /)
    Cm(1:6,2) = (/ c1, c1+2.0*c2, c1, 0.0, 0.0, 0.0 /)
    Cm(1:6,3) = (/ c1, c1, c1+2.0*c2, 0.0, 0.0, 0.0 /)
    Cm(1:6,4) = (/ 0.0, 0.0, 0.0, c2, 0.0, 0.0 /)
    Cm(1:6,5) = (/ 0.0, 0.0, 0.0, 0.0, c2, 0.0 /)
    Cm(1:6,6) = (/ 0.0, 0.0, 0.0, 0.0, 0.0, c2 /)
    
    
    Cf(1:6,1) = (/ ((E11f**2)*(v23f-1)/delf), (-E11f*E22f*v12f/delf), (-E11f*E22f*v12f/delf), 0.0, 0.0, 0.0 /)
    Cf(1:6,2) = (/ (-E11f*E22f*v12f/delf), (E22f*(E22f*(v12f**2)-E11f)/((1+v23f)*delf)), (-E22f*(E22f*(v12f**2)+E11f*v23f)/((1+v23f)*delf)), 0.0, 0.0, 0.0 /)
    Cf(1:6,3) = (/ (-E11f*E22f*v12f/delf), (-E22f*(E22f*(v12f**2)+E11f*v23f)/((1+v23f)*delf)), (E22f*(E22f*(v12f**2)-E11f)/((1+v23f)*delf)), 0.0, 0.0, 0.0 /)
    Cf(1:6,4) = (/ 0.0, 0.0, 0.0, 2.0*G23f, 0.0, 0.0 /)
    Cf(1:6,5) = (/ 0.0, 0.0, 0.0, 0, 2.0*G12f, 0.0 /)
    Cf(1:6,6) = (/ 0.0, 0.0, 0.0, 0.0, 0.0, 2.0*G12f /)

    
    WRITE (*,*) Cf
    PAUSE
END PROGRAM GMC

