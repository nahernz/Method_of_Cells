        !COMPILER-GENERATED INTERFACE MODULE: Mon Aug 11 10:32:50 2014
        MODULE DGETRF__genmod
          INTERFACE 
            SUBROUTINE DGETRF(M,N,A,LDA,IPIV,INFO)
              INTEGER(KIND=4) :: LDA
              INTEGER(KIND=4) :: M
              INTEGER(KIND=4) :: N
              REAL(KIND=8) :: A(LDA,*)
              INTEGER(KIND=4) :: IPIV(*)
              INTEGER(KIND=4) :: INFO
            END SUBROUTINE DGETRF
          END INTERFACE 
        END MODULE DGETRF__genmod
