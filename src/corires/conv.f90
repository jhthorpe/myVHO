!------------------------------------------------------------
! conv
!       - module containing using conversion subroutines
!------------------------------------------------------------
MODULE conv

CONTAINS
!------------------------------------------------------------
! conv_cm2MHz
!       - converts cm-1 to MHz
!------------------------------------------------------------
REAL(KIND=8) FUNCTION conv_cm2MHz(val)
  IMPLICIT NONE
  REAL(KIND=8), INTENT(IN) :: val
  conv_cm2MHz = val*29979.2458D0
END FUNCTION conv_cm2MHz
!------------------------------------------------------------

END MODULE conv
!------------------------------------------------------------
