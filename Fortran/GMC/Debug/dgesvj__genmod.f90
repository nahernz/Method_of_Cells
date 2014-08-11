        !COMPILER-GENERATED INTERFACE MODULE: Mon Aug 11 10:29:07 2014
        MODULE DGESVJ__genmod
          INTERFACE 
            SUBROUTINE DGESVJ(JOBA,JOBU,JOBV,M,N,A,LDA,SVA,MV,V,LDV,WORK&
     &,LWORK,INFO)
              INTEGER(KIND=4) :: LWORK
              INTEGER(KIND=4) :: LDV
              INTEGER(KIND=4) :: LDA
              INTEGER(KIND=4) :: N
              CHARACTER(LEN=1) :: JOBA
              CHARACTER(LEN=1) :: JOBU
              CHARACTER(LEN=1) :: JOBV
              INTEGER(KIND=4) :: M
              REAL(KIND=8) :: A(LDA,*)
              REAL(KIND=8) :: SVA(N)
              INTEGER(KIND=4) :: MV
              REAL(KIND=8) :: V(LDV,*)
              REAL(KIND=8) :: WORK(LWORK)
              INTEGER(KIND=4) :: INFO
            END SUBROUTINE DGESVJ
          END INTERFACE 
        END MODULE DGESVJ__genmod
