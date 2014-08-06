MODULE  cell2
    IMPLICIT NONE
    ! Fiber Material Properties
    REAL, PARAMETER :: E11f = 231e9
    REAL, PARAMETER :: E22f = 15e9
    REAL, PARAMETER :: v12f = 0.14
    REAL, PARAMETER :: G12f = 24e9
    REAL, PARAMETER :: G23f = 5.01e9
    !REAL, PARAMETER :: v23f
    ! Matrix Material Properties
    REAL, PARAMETER :: Em = 3e9
    REAL, PARAMETER :: vm = 0.36
    
    END MODULE  cell2