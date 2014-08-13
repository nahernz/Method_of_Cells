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
MODULE Cell
    DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: L, H
    DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: SM
END MODULE Cell
    
    
PROGRAM GMC
    USE Cell
    IMPLICIT NONE
    
    EXTERNAL DGESV

    ! Variables
    CHARACTER(32) :: input
    CHARACTER(LEN=9) :: dummy = "input.txt"
    ! Fiber Material Properties
    DOUBLE PRECISION :: E11f, E22f, v12f, v23f, G12f, G23f
    ! Matrix Material Properties
    DOUBLE PRECISION :: Em, vm
    ! Build Material Stiffnesses
    DOUBLE PRECISION                  :: delf
    DOUBLE PRECISION, DIMENSION(6, 6) :: Cf, Cm
    DOUBLE PRECISION                  :: c1, c2
    ! Geometry Data
    INTEGER :: ny, nx
    DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: x, y
    INTEGER              :: Ng, Nb, g, b
    ! Solving Variables
    DOUBLE PRECISION, DIMENSION(:, :), ALLOCATABLE :: A, K, M
	INTEGER, ALLOCATABLE :: IPIV(:)
    INTEGER :: count = 1, i, j, NRHS, LDA, LDB, N, INFO
    ! Output Variables
    DOUBLE PRECISION C(6,6)
    CHARACTER(11), PARAMETER :: frmt = "(6EN13.2E1)"
    
    
!-----------------------------------------------------------------------------
!                           CODE STARTS
!-----------------------------------------------------------------------------
    WRITE (*,*) "Generalized Method of Cells"
    WRITE (*,*) "written by Michael Kaplan and Rehan Nawaz"
    WRITE (*,*)
    
    Call GetArg(1,input)
    IF (LEN_TRIM(input) == 0) THEN
        input = dummy
    END IF
    
!-----------------------------------------------------------------------------
!                            READ INPUT FILE
!-----------------------------------------------------------------------------
    
    CALL GMC_in(input, E11f, E22f, v12f, v23f, G12f, G23f, &
                Em, vm, nx, ny)
    ALLOCATE(x(nx),y(ny))
    
    
    ! Body of GMC
    WRITE (*,'(A)') 'Generalized Method of Cells'
    WRITE (*,'(A)') 'written by Michael Kaplan and Rehan Nawaz'
    
    ! Build the stiffness matricees
    
    delf = E22f*v12f**2+E11f*(v23f-1.0)
    c1 = Em*vm/((1.0+vm)*(1.0-2.0*vm))
    c2 = Em/(2.0*(1.0+vm))

    Cm(1:6,1) = (/ c1+2.0*c2, c1, c1, 0.0d0, 0.0d0, 0.0d0 /)
    Cm(1:6,2) = (/ c1, c1+2.0*c2, c1, 0.0d0, 0.0d0, 0.0d0 /)
    Cm(1:6,3) = (/ c1, c1, c1+2.0*c2, 0.0d0, 0.0d0, 0.0d0 /)
    Cm(1:6,4) = (/ 0.0d0, 0.0d0, 0.0d0, c2, 0.0d0, 0.0d0 /)
    Cm(1:6,5) = (/ 0.0d0, 0.0d0, 0.0d0, 0.0d0, c2, 0.0d0 /)
    Cm(1:6,6) = (/ 0.0d0, 0.0d0, 0.0d0, 0.0d0, 0.0d0, c2 /)
    
    
    Cf(1:6,1) = (/ ((E11f**2)*(v23f-1)/delf), (-E11f*E22f*v12f/delf), (-E11f*E22f*v12f/delf), 0.0d0, 0.0d0, 0.0d0 /)
    Cf(1:6,2) = (/ (-E11f*E22f*v12f/delf), (E22f*(E22f*(v12f**2)-E11f)/((1+v23f)*delf)), (-E22f*(E22f*(v12f**2)+E11f*v23f)/((1+v23f)*delf)), 0.0d0, 0.0d0, 0.0d0 /)
    Cf(1:6,3) = (/ (-E11f*E22f*v12f/delf), (-E22f*(E22f*(v12f**2)+E11f*v23f)/((1+v23f)*delf)), (E22f*(E22f*(v12f**2)-E11f)/((1+v23f)*delf)), 0.0d0, 0.0d0, 0.0d0 /)
    Cf(1:6,4) = (/ 0.0d0, 0.0d0, 0.0d0, 2.0d0*G23f, 0.0d0, 0.0d0 /)
    Cf(1:6,5) = (/ 0.0d0, 0.0d0, 0.0d0, 0, 2.0d0*G12f, 0.0d0 /)
    Cf(1:6,6) = (/ 0.0d0, 0.0d0, 0.0d0, 0.0d0, 0.0d0, 2.0d0*G12f /)
    
    ! Geometry Variables
    ! must change 'nc' to change cell grid size
    ! Four cell
    !L = (/ 6e-3, 4e-3 /)
    !H = L
    !SM(1:2,1) = (/ 1,2 /)
    !SM(1:2,2) = (/ 2,2 /)
    
    ! Square Pack
    !L = (/ 6.75334e-4, 6.14488e-4, 6.14488e-4, 2.45795e-3,&
    !       6.14488e-4, 6.14488e-4, 6.75334e-4 /)
    !H = L
    !SM(1:7,1) = (/ 2,2,2,2,2,2,2 /)
    !SM(1:7,2) = (/ 2,2,2,1,2,2,2 /)
    !SM(1:7,3) = (/ 2,2,1,1,1,2,2 /)
    !SM(1:7,4) = (/ 2,1,1,1,1,1,2 /)
    !SM(1:7,5) = (/ 2,2,1,1,1,2,2 /)
    !SM(1:7,6) = (/ 2,2,2,1,2,2,2 /)
    !SM(1:7,7) = (/ 2,2,2,2,2,2,2 /)

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

M = K
!WRITE (*,*) N
!WRITE (*,*) NRHS
!WRITE (*,*) LDA
!WRITE (*,*) LDB
!WRITE (*,*) SIZE(A,1), SIZE(A,2)
!WRITE (*,*) SIZE(M,1), SIZE(M,2)
!WRITE (*,frmt) A(1,1:6)


call DGESV( N, NRHS, A, LDA, IPIV, M, LDB, INFO )


IF (INFO == 0) THEN
    ! solution found
ELSE IF (INFO < 0) THEN
    WRITE(*,*) "Value ", INFO, " is illegal!"
    STOP
ELSE IF (INFO > 0) THEN
    WRITE(*,*) "The solution is singular at U(", INFO, ",", INFO, ")"
    STOP
END IF

!------------------------------------------------------------------------------
!                            BUILD STIFFNESS
!------------------------------------------------------------------------------
!WRITE (*,*) 
C = 0
i = 1
    DO b=1,Nb
        DO g=1,Ng
            IF (SM(b,g)==1) THEN
                C = C + 1/(sum(H)*sum(L))*H(b)*L(g)*MATMUL(Cf,M((i-1)*6+1:(i)*6,1:6))
            ELSE 
                C = C + 1/(sum(H)*sum(L))*H(b)*L(g)*MATMUL(Cm,M((i-1)*6+1:(i)*6,1:6))
            END IF
            i = i + 1
            
        END DO
    END DO
    !WRITE (*,*)
    !WRITE (*,"(EN12.2E1)") C(1,1) - 2*C(1,2)**2/(C(2,2) + C(2,3))
    
!------------------------------------------------------------------------------
!                         WRITE STIFFNESS TO STIFFNESS.TXT
!------------------------------------------------------------------------------
    WRITE (*,*) "C ="
    WRITE (*,*)
    
    OPEN (unit = 2, file = "Stiffness.txt")
    WRITE (2,*) "C ="
    WRITE (2,*)
    DO i = 1,6
        WRITE (2,frmt) (C(i,j), j = 1,6)
        WRITE (*,frmt) (C(i,j), j = 1,6)
        !WRITE (*,*)
    END DO
    WRITE (*,*)
    
    CLOSE (2)
    PAUSE
    END PROGRAM GMC

! -----------------------------------------------------------------------
!               Generalized Method of Cells Input Reader
! -----------------------------------------------------------------------
!
!     Read the input file specified and return the material properties   
!     and cell architecture to the GMC code. The code is currently only  
!     capable of specifying a transversely isotropic fiber and isotropic
!     matrix. Materials are linear elastic and do not damage. 
      
      SUBROUTINE GMC_in(input, E11f, E22f, v12f, v23f, G12f, G23f, &
                        Em, vm, nx, ny)
      USE Cell
      IMPLICIT NONE
      
      CHARACTER(LEN=*) :: input
      DOUBLE PRECISION :: E11f, E22f, v12f, v23f, G12f, G23f, Em, vm
      INTEGER :: nx, ny, i
      LOGICAL :: there
      
      INQUIRE(file = input, EXIST = there)
      if (there) THEN      
          OPEN(unit = 1, file = input)
          READ(1,*) E11f, E22f, v12f, v23f, G12f, G23f
          READ(1,*) Em, vm
          READ(1,*) nx, ny
      
          ALLOCATE(L(nx)) 
          ALLOCATE(H(ny))
          ALLOCATE(SM(ny,nx))
      
          READ(1,*) L
          READ(1,*) H
      
          DO i = 1,ny
              READ(1,*) SM(1:nx,i)
          END DO
      ELSE
          WRITE (*,*) "Input file """, input(1:LEN_TRIM(input)), """ does not exist"
          WRITE (*,*)
          PAUSE
          STOP
      END IF
      
      END SUBROUTINE GMC_in
