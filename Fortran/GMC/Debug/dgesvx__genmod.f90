        !COMPILER-GENERATED INTERFACE MODULE: Mon Aug 11 10:28:45 2014
        MODULE DGESVX__genmod
          INTERFACE 
            SUBROUTINE DGESVX(FACT,TRANS,N,NRHS,A,LDA,AF,LDAF,IPIV,EQUED&
     &,R,C,B,LDB,X,LDX,RCOND,FERR,BERR,WORK,IWORK,INFO)
              INTEGER(KIND=4) :: LDX
              INTEGER(KIND=4) :: LDB
              INTEGER(KIND=4) :: LDAF
              INTEGER(KIND=4) :: LDA
              CHARACTER(LEN=1) :: FACT
              CHARACTER(LEN=1) :: TRANS
              INTEGER(KIND=4) :: N
              INTEGER(KIND=4) :: NRHS
              REAL(KIND=8) :: A(LDA,*)
              REAL(KIND=8) :: AF(LDAF,*)
              INTEGER(KIND=4) :: IPIV(*)
              CHARACTER(LEN=1) :: EQUED
              REAL(KIND=8) :: R(*)
              REAL(KIND=8) :: C(*)
              REAL(KIND=8) :: B(LDB,*)
              REAL(KIND=8) :: X(LDX,*)
              REAL(KIND=8) :: RCOND
              REAL(KIND=8) :: FERR(*)
              REAL(KIND=8) :: BERR(*)
              REAL(KIND=8) :: WORK(*)
              INTEGER(KIND=4) :: IWORK(*)
              INTEGER(KIND=4) :: INFO
            END SUBROUTINE DGESVX
          END INTERFACE 
        END MODULE DGESVX__genmod
