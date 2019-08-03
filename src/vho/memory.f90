!------------------------------------------------------------
! memory
!       - module for managing and analysing memory
!------------------------------------------------------------
MODULE memory
  USE val
  USE linal

CONTAINS
!------------------------------------------------------------
! memory_HObuild 
!       - analyses the memory situation for building the 
!         hamiltonian
!       - three memory scenarios:
!         1) minmem, hamiltonian contructed one vector at a time
!         2) precalc, H one vec, integrals precalced
!         3) inmem, H incore, integrals precalced
!       
!------------------------------------------------------------
! job           : int, jobtype
! mem           : int, memory in MB
! N             : int, size of hamiltonian
! ndim          : int, number of dimensions
! nbas          : int, number of basis in each dimension
! nabs          : 1D int, number of abscissa
! memstat       : int, memory status
! error         : int, exit code

SUBROUTINE memory_HObuild(job,mem,N,ndim,nbas,nabs,memstat,error)
  IMPLICIT NONE
  INTEGER, DIMENSION(0:), INTENT(IN) :: nbas,nabs
  INTEGER(KIND=8), INTENT(IN) :: mem
  INTEGER, INTENT(INOUT) :: memstat,error
  INTEGER, INTENT(IN) :: job,N,ndim
  REAL(KIND=8) :: minmem,premem,inmem,basemem,qw2mb,qmem
  INTEGER :: mbas,mabs,M

  error = 0
  qw2mb = 8.0D0/1000000.0D0
  qmem = mem/qw2mb
  mbas = MAXVAL(nbas)
  mabs = MAXVAL(nabs)
  M = PRODUCT(nabs)
  WRITE(*,*) "Hamiltonian Memory Analysis"
  CALL val_check(qmem,error)
  IF (error .NE. 0) THEN
    WRITE(*,*) "memory_HObuild  : ERROR"
    WRITE(*,*) "qmem was too large to be analysed"
    error = 1
    RETURN
  END IF

  basemem = 1000
  IF (job .EQ. 2) THEN
    ! terms are [minmem], [premem], [inmem]
    ! These are out of date
    ! hamiltonian : N, N, N^2
    ! abscissa    : nabs
    ! weights     : nabs
    ! couplings   : 5*ndim (estimated)
    ! V           : nabs, nabs*ndim 
    ! hermite     : mbas * nabs 
    ! VTint       : ndim*mbas, mbas^2*ndim
    ! Q1-QPint    : ndim*mbas, 6*mbas^2*ndim
    ! norm        : mbas
    WRITE(*,*) "memory_HObuild  : WARNING"
    WRITE(*,*) "Memory analysis for jobtype 2 is probably wrong"
    WRITE(*,*)
    
    minmem = N + 3*mabs + 5*ndim + MAXVAL(nbas)*mabs + &
             MAXVAL(nbas)*ndim + MAXVAL(nbas) + &
             MAXVAL(nbas)*ndim + basemem
    premem = N + 2*mabs + 5*ndim + mabs*ndim + &
                MAXVAL(nbas)*mabs + MAXVAL(nbas)**2.0D0*ndim &
                + MAXVAL(nbas) + 6*MAXVAL(nbas)*ndim + basemem
    inmem = N**2.0D0 + 2*mabs + 5*ndim + mabs*ndim + &
                MAXVAL(nbas)*mabs + MAXVAL(nbas)**2.0D0*ndim &
                + MAXVAL(nbas) + 6*MAXVAL(nbas)*ndim + basemem
  ELSE IF (job .EQ. 1) THEN
    ! terms are [minmem], [premem], [inmem]
    ! hamiltonian : N, N, N^2
    ! FCs         : 5*ndim, 5*ndim, 5*ndim 
    ! Q1int       : 0, mbas*ndim 
    ! Q2int       : 0, 2*mbas*ndim
    ! Q3int       : 0, 2*mbas*ndim
    ! Q4int       : 0, 3*mbas*ndim       
    ! P2int       : 0, 2*mbas*ndim
    minmem = basemem + N + 5*ndim
    premem = basemem + N + 5*ndim + mbas*ndim  + 2*mbas*ndim &
             + 2*mbas*ndim + 3*mbas*ndim + 2*mbas*ndim
    inmem = basemem + N**2 + 5*ndim + mbas*ndim  + 2*mbas*ndim &
             + 2*mbas*ndim + 3*mbas*ndim + 2*mbas*ndim
  ELSE IF (job .EQ. 3) THEN
    !terms are [minmem], [premem], [inmem]
    ! N = product(nbas), M = product(nabs)
    ! hamiltonian : N, N, N^2
    ! Vq          : M, M, M
    ! Hermite     : mbas*mabs*ndim
    ! Heff        : mabs*2*ndim
    ! Norm        : mbas
    minmem = basemem + N + M + mbas*mabs*ndim + mbas
    premem = basemem + N + M + mbas*mabs*ndim + mbas
    inmem = basemem + N**2.0D0 + M + mbas*mabs*ndim + mbas 
  ELSE 
    
    error = 1
    WRITE(*,*) "memory_HO_build  : ERROR"
    WRITE(*,*) "This jobtype is not supported :", job
    RETURN
  END IF

  WRITE(*,'(A20,F12.2)') "Available memory   ", qmem*qw2mb
  WRITE(*,'(A20,F12.2)') "Minimum memory     ", minmem*qw2mb
  WRITE(*,'(A20,F12.2)') "Precalc memory     ", minmem*qw2mb
  WRITE(*,'(A20,F12.2)') "Incore memory      ", inmem*qw2mb
  WRITE(*,*) "Values are in MB"

  WRITE(*,*)
  IF (qmem .LT. minmem) THEN
    WRITE(*,*) "memory_Hbuild  : ERROR"
    WRITE(*,*) "Not enough memory has been given"
    error = 2
    memstat = -1
    RETURN
  ELSE IF (qmem .LT. inmem .AND. qmem .GE. minmem) THEN
    WRITE(*,*) "HO Hamiltonian will be built with minimal memory"
    memstat = 0
  ELSE IF (qmem .GT. premem .AND. qmem .LT. inmem) THEN
    WRITE(*,*) "HO Hamiltonian will be built with precalc memory"
    memstat = 1
  ELSE IF (qmem .GE. inmem) THEN
    WRITE(*,*) "HO Hamiltonian will be built with incore memory"
    memstat = 2
  ELSE
    WRITE(*,*) "memory_Hbuild  : ERROR"
    WRITE(*,*) "Yours truely has somehow coded this badly..."
    error = 3
    RETURN
  END IF
  WRITE(*,*)

END SUBROUTINE memory_HObuild

!------------------------------------------------------------
! memory_Hdiag
!       - figures out the memory situation for the 
!         diagonalization subroutine
!------------------------------------------------------------
! ndim          : int, number of dimensions
! nbas          : 1D int, basis functions per dimension
! enum          : int, number of eigenvalues to calc
! mem           : int*8, memory in MB
! memstat       : int, memory status, output
! lwork         : int, length of working vector allowed
! error         : int, error 

SUBROUTINE memory_Hdiag(ndim,nbas,enum,mem,memstat,lwork,error)
  IMPLICIT NONE
  INTEGER, DIMENSION(0:), INTENT(IN) :: nbas
  INTEGER(KIND=8), INTENT(IN) :: mem
  INTEGER, INTENT(INOUT) :: memstat,error,lwork
  INTEGER, INTENT(IN) :: ndim,enum
  REAL(KIND=8) :: minmem,inmem,basemem,qw2mb,qmem
  INTEGER :: N

  error = 0
  N = PRODUCT(nbas)
  qw2mb = 8.0D0/1000000.0D0

  qmem = mem/qw2mb
  CALL val_check(qmem,error)
  IF (error .NE. 0) THEN
    WRITE(*,*) "memory_Hdiag  : ERROR"
    WRITE(*,*) "qmem was too large to be analysed"
    error = 1
    RETURN
  END IF
 
  !Memory required
  ! in background...
  ! q           : nabs
  ! W           : nabs
  ! in diag routine...
  ! Hij         : N^2
  ! Cij         : N or N*enum 
  ! eval        : enum
  ! WORK        : 3*N, lwork 
  ! IWORK       : 5*N
  
  
  CALL linal_dsyevx_lwork(N,enum,lwork,error)
  IF (error .NE. 0) RETURN
  
  basemem = 1000
  minmem = basemem + N**2 + N + enum + 3*N + 5.0D0*N/2.0D0
  inmem = basemem + N**2 +  N*enum + enum + lwork + 5.0D0*N/2.0D0
  
  WRITE(*,*) "Diagonalization Memory Analysis"
  WRITE(*,'(A20,F12.2)') "Available memory   ", qmem*qw2mb
  WRITE(*,'(A20,F12.2)') "Minimum memory     ", minmem*qw2mb
  WRITE(*,'(A20,F12.2)') "Incore memory      ", inmem*qw2mb
  WRITE(*,*) "Values are in MB"

  WRITE(*,*)
  IF (qmem .LT. minmem) THEN
    WRITE(*,*) "memory_Hdiag  : ERROR"
    WRITE(*,*) "Not enough memory has been given"
    error = 2
    memstat = -1
    RETURN
  ELSE IF (qmem .LT. inmem .AND. qmem .GE. minmem) THEN
    WRITE(*,*) "Hamiltonian will be diagonalized with minimal memory"
    memstat = 0

  ELSE IF (qmem .GE. inmem) THEN
    WRITE(*,*) "Hamiltonian will be diagonalized with incore memory"
    memstat = 2
  ELSE
    WRITE(*,*) "memory_Hdiag  : ERROR"
    WRITE(*,*) "Yours truely has somehow coded this badly..."
    error = 3
    RETURN
  END IF
  WRITE(*,*)

END SUBROUTINE memory_Hdiag

!------------------------------------------------------------


END MODULE memory
