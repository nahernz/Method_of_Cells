        !COMPILER-GENERATED INTERFACE MODULE: Mon Aug 11 10:29:02 2014
        MODULE DGESVXX__genmod
          INTERFACE 
            SUBROUTINE DGESVXX(FACT,TRANS,N,NRHS,A,LDA,AF,LDAF,IPIV,    &
     &EQUED,R,C,B,LDB,X,LDX,RCOND,RPVGRW,BERR,N_ERR_BNDS,ERR_BNDS_NORM, &
     &ERR_BNDS_COMP,NPARAMS,PARAMS,WORK,IWORK,INFO)
              INTEGER(KIND=4) :: LDX
              INTEGER(KIND=4) :: LDB
              INTEGER(KIND=4) :: LDAF
              INTEGER(KIND=4) :: LDA
              INTEGER(KIND=4) :: NRHS
              CHARACTER(LEN=1) :: FACT
              CHARACTER(LEN=1) :: TRANS
              INTEGER(KIND=4) :: N
              REAL(KIND=8) :: A(LDA,*)
              REAL(KIND=8) :: AF(LDAF,*)
              INTEGER(KIND=4) :: IPIV(*)
              CHARACTER(LEN=1) :: EQUED
              REAL(KIND=8) :: R(*)
              REAL(KIND=8) :: C(*)
              REAL(KIND=8) :: B(LDB,*)
              REAL(KIND=8) :: X(LDX,*)
              REAL(KIND=8) :: RCOND
              REAL(KIND=8) :: RPVGRW
              REAL(KIND=8) :: BERR(*)
              INTEGER(KIND=4) :: N_ERR_BNDS
              REAL(KIND=8) :: ERR_BNDS_NORM(NRHS,*)
              REAL(KIND=8) :: ERR_BNDS_COMP(NRHS,*)
              INTEGER(KIND=4) :: NPARAMS
              REAL(KIND=8) :: PARAMS(*)
              REAL(KIND=8) :: WORK(*)
              INTEGER(KIND=4) :: IWORK(*)
              INTEGER(KIND=4) :: INFO
            END SUBROUTINE DGESVXX
          END INTERFACE 
        END MODULE DGESVXX__genmod
