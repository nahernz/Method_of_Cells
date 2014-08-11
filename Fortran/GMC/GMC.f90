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
    
    EXTERNAL DGESV

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
    REAL                  :: delf
    REAL, DIMENSION(6, 6) :: Cf, Cm
    REAL                  :: c1, c2
    ! Geometry Data
    REAL, DIMENSION(7)   :: H, L, y, x
    REAL, DIMENSION(7,7) :: SM
    INTEGER              :: Ng, Nb, g, b
    ! Solving Variables
    DOUBLE PRECISION, DIMENSION(:, :), ALLOCATABLE :: A, K, M
	INTEGER, ALLOCATABLE :: IPIV(:)
    INTEGER :: count = 1, i, NRHS, LDA, LDB, N, INFO
    
    
!-----------------------------------------------------------------------------
!                           CODE STARTS
!-----------------------------------------------------------------------------
    
    ! Body of GMC
    WRITE (*,'(A)') 'Generalized Method of Cells'
    WRITE (*,'(A)') 'written by Michael Kaplan and Rehan Nawaz'
    
    ! Build the stiffness matricees
    
    delf = E22f*v12f**2+E11f*(v23f-1)
    c1 = Em*vm/((1+vm)*(1-2*vm))
    c2 = Em/(2*(1+vm))
    
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
    
    ! Geometry Variables
    ! Square Pack
    
    L = (/ 6.75334e-4, 6.14488e-4, 6.14488e-4, 2.45795e-3,&
           6.14488e-4, 6.14488e-4, 6.75334e-4 /)
    H = L
    SM(1:7,1) = (/ 2,2,2,2,2,2,2 /)
    SM(1:7,2) = (/ 2,2,2,1,2,2,2 /)
    SM(1:7,3) = (/ 2,2,1,1,1,2,2 /)
    SM(1:7,4) = (/ 2,1,1,1,1,1,2 /)
    SM(1:7,5) = (/ 2,2,1,1,1,2,2 /)
    SM(1:7,6) = (/ 2,2,2,1,2,2,2 /)
    SM(1:7,7) = (/ 2,2,2,2,2,2,2 /)

    Ng = SIZE(L)
    Nb = SIZE(H)
    
    x(1) = L(1)/2
    DO g = 2, Ng
        x(g) = SUM(L(1:g-1))+L(g)/2
    END DO
    
    y(1) = H(1)/2
    DO b = 2, Nb
        y(b) = SUM(H(1:b-1))+H(b)/2
    END DO

!------------------------------------------------------------------------------
!                             BUILD SOLUTION
!------------------------------------------------------------------------------
    N = Ng*Nb*6
    ALLOCATE(A(N,N))
    ALLOCATE(K(N,6))
    ALLOCATE(M(N,6))
    ALLOCATE(IPIV(N))
    
    NRHS = 6
    LDA = N
    LDB = N
    
    A = 0.0
    K = 0.0
    
    !e33
    DO b = 1,Nb
        DO g = 1,Ng
            A(count,(g-1)*6+(b-1)*6*Ng+3) = L(g)
        END DO
        K(count,3) = SUM(L)
        count = count + 1
    END DO
    !e22
    DO g = 1,Ng
        DO b = 1,Nb
            A(count,(g-1)*6+(b-1)*6*Ng+2) = H(b)
        END DO
        K(count,2) = sum(H)
        count = count + 1
    END DO
    !e11
    DO i = 1,Nb*Ng
        A(count,(i-1)*6+1) = 1
        K(count,1) = 1
        count = count + 1
    END DO
    !e13
    DO b = 1,Nb
        DO g = 1,Ng
            A(count,(g-1)*6+(b-1)*6*Ng+5) = L(g)
        END DO
    K(count,5) = sum(L)
    count = count + 1
    END DO
    !e12
    DO g = 1,Ng
        DO b = 1,Nb
            A(count,(g-1)*6+(b-1)*6*Ng+6) = H(b)
        END DO
    K(count,6) = sum(H)
    count = count + 1
    END DO
    !e23
DO g = 1,Ng
    DO b = 1,Nb
        A(count,(g-1)*6+(b-1)*6*Ng+4) = H(b)*L(g)
    END DO
END DO

K(count,4) = sum(H)*sum(L)
count = count + 1

!Am
!r22
DO g = 1,Ng
    DO b = 1,(Nb-1)
        IF (SM(b,g) == 1) THEN
            A(count,(g-1)*6+(b-1)*6*Ng+1) = Cf(2,1);
            A(count,(g-1)*6+(b-1)*6*Ng+2) = Cf(2,2);
            A(count,(g-1)*6+(b-1)*6*Ng+3) = Cf(2,3);           
        ELSE
            A(count,(g-1)*6+(b-1)*6*Ng+1) = Cm(2,1);
            A(count,(g-1)*6+(b-1)*6*Ng+2) = Cm(2,2);
            A(count,(g-1)*6+(b-1)*6*Ng+3) = Cm(2,3);
        END IF
        IF (SM(b+1,g) == 1) THEN
            A(count,(g-1)*6+(b)*6*Ng+1) = -Cf(2,1);
            A(count,(g-1)*6+(b)*6*Ng+2) = -Cf(2,2);
            A(count,(g-1)*6+(b)*6*Ng+3) = -Cf(2,3);
       ELSE
            A(count,(g-1)*6+(b)*6*Ng+1) = -Cm(2,1);
            A(count,(g-1)*6+(b)*6*Ng+2) = -Cm(2,2);
            A(count,(g-1)*6+(b)*6*Ng+3) = -Cm(2,3);
        END IF
        count = count + 1;
    END DO
END DO

!r33
DO b = 1,Nb
    DO g = 1,(Ng-1)
        IF (SM(b,g)==1) THEN
            A(count,(g-1)*6+(b-1)*6*Ng+1)=Cf(3,1);
            A(count,(g-1)*6+(b-1)*6*Ng+2)=Cf(3,2);
            A(count,(g-1)*6+(b-1)*6*Ng+3)=Cf(3,3);    
        ELSE
            A(count,(g-1)*6+(b-1)*6*Ng+1)=Cm(3,1);
            A(count,(g-1)*6+(b-1)*6*Ng+2)=Cm(3,2);
            A(count,(g-1)*6+(b-1)*6*Ng+3)=Cm(3,3);
        END IF
        IF (SM(b,g+1)==1) THEN
            A(count,(g)*6+(b-1)*6*Ng+1)=-Cf(3,1);
            A(count,(g)*6+(b-1)*6*Ng+2)=-Cf(3,2);
            A(count,(g)*6+(b-1)*6*Ng+3)=-Cf(3,3);   
        ELSE
            A(count,(g)*6+(b-1)*6*Ng+1)=-Cm(3,1);
            A(count,(g)*6+(b-1)*6*Ng+2)=-Cm(3,2);
            A(count,(g)*6+(b-1)*6*Ng+3)=-Cm(3,3);
        END IF
        count = count + 1;
    END DO
END DO

!r12
DO g = 1,Ng
    DO b = 1,(Nb-1)
        IF (SM(b,g)==1) THEN
            A(count,(g-1)*6+(b-1)*6*Ng+6)=Cf(6,6);            
        ELSE
            A(count,(g-1)*6+(b-1)*6*Ng+6)=Cm(6,6);
        END IF
        IF (SM(b+1,g)==1) THEN
            A(count,(g-1)*6+(b)*6*Ng+6)=-Cf(6,6);
        ELSE
            A(count,(g-1)*6+(b)*6*Ng+6)=-Cm(6,6);
        END IF
        count = count + 1;
    END DO
END DO

!r13
DO b = 1,Nb
    DO g = 1,(Ng-1)
        IF (SM(b,g)==1) THEN
            A(count,(g-1)*6+(b-1)*6*Ng+5)=Cf(5,5);            
        ELSE
            A(count,(g-1)*6+(b-1)*6*Ng+5)=Cm(5,5);
        END IF
        IF (SM(b,g+1)==1) THEN
            A(count,(g)*6+(b-1)*6*Ng+5)=-Cf(5,5);
        ELSE
            A(count,(g)*6+(b-1)*6*Ng+5)=-Cm(5,5);
        END IF
        count = count + 1;
    END DO
END DO

!r23
DO g = 1,Ng
    DO b = 1,(Nb-1)
        IF (SM(b,g)==1) THEN
            A(count,(g-1)*6+(b-1)*6*Ng+4)=Cf(4,4);            
        ELSE
            A(count,(g-1)*6+(b-1)*6*Ng+4)=Cm(4,4);
        END IF
        IF (SM(b+1,g)==1) THEN
            A(count,(g-1)*6+(b)*6*Ng+4)=-Cf(4,4);
        ELSE
            A(count,(g-1)*6+(b)*6*Ng+4)=-Cm(4,4);
        END IF
        count = count + 1;
    END DO
END DO

DO b = 1,Nb
    DO g = 1,(Ng-1)
        IF (count <= Nb*Ng*6) THEN
            IF (SM(b,g)==1) THEN
                A(count,(g-1)*6+(b-1)*6*Ng+4)=Cf(4,4);            
            ELSE 
                A(count,(g-1)*6+(b-1)*6*Ng+4)=Cm(4,4);
            END IF
            IF (SM(b,g+1)==1) THEN
                A(count,(g)*6+(b-1)*6*Ng+4)=-Cf(4,4);
            ELSE
                A(count,(g)*6+(b-1)*6*Ng+4)=-Cm(4,4);
            END IF
            count = count + 1;
        END IF
    END DO
END DO
!------------------------------------------------------------------------------
!                                   SOLVE
!------------------------------------------------------------------------------

!CALL inverse(A_,M,Ng*Nb*6)
M = K
WRITE (*,*) N
WRITE (*,*) NRHS
WRITE (*,*) LDA
WRITE (*,*) LDB
WRITE (*,*) SIZE(A,1), SIZE(A,2)
WRITE (*,*) SIZE(M,1), SIZE(M,2)


call DGESV( N, NRHS, A, LDA, IPIV, M, LDB, INFO )

Write (*,*) INFO


    PAUSE
    END PROGRAM GMC


