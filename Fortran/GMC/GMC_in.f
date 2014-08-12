c -----------------------------------------------------------------------
c               Generalized Method of Cells Input Reader
c -----------------------------------------------------------------------
c
c     Read the input file specified and return the material properties   
c     and cell architecture to the GMC code. The code is currently only  
c     capable of specifying a transversely isotropic fiber and isotropic
c     matrix. Materials are linear elastic and do not damage. 
      
      SUBROUTINE GMC_in(input, E11f, E22f, v12f, v23f, G12f, G23f, &
                        Em, vm, &
                        nx, ny, &
                        L, H, SM)
      IMPLICIT NONE
      
      CHARACTER(LEN=*) :: input
      DOUBLE PRECISION :: E11f, E22f, v12f, v23f, G12f, G23f, Em, vm
      INTEGER :: nx, ny, i
      DOUBLE PRECISION, DIMENSION(*), ALLOCATABLE :: L, H
      DOUBLE PRECISION, DIMENSION(*,*), ALLOCATABLE :: SM
      
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
      END SUBROUTINE GMC_in
      