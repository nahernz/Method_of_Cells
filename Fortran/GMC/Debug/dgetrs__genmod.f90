        !COMPILER-GENERATED INTERFACE MODULE: Mon Aug 11 10:32:47 2014
        MODULE DGETRS__genmod
          INTERFACE 
            SUBROUTINE DGETRS(TRANS,N,NRHS,A,LDA,IPIV,B,LDB,INFO)
              INTEGER(KIND=4) :: LDB
              INTEGER(KIND=4) :: LDA
              CHARACTER(LEN=1) :: TRANS
              INTEGER(KIND=4) :: N
              INTEGER(KIND=4) :: NRHS
              REAL(KIND=8) :: A(LDA,*)
              INTEGER(KIND=4) :: IPIV(*)
              REAL(KIND=8) :: B(LDB,*)
              INTEGER(KIND=4) :: INFO
            END SUBROUTINE DGETRS
          END INTERFACE 
        END MODULE DGETRS__genmod
