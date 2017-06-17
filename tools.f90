
MODULE MODULE_TOOLS
  !*****************************************************************************
  !** The module 'MODULE_TOOL' contains useful non-numerical procedures       **
  !** used by several subroutines of the MHD shock code                       **
  !*****************************************************************************
  IMPLICIT NONE
  INCLUDE "precision.f90"

  ! variables from "precision.f90" are private
  PRIVATE :: DP, LONG, minus_infinity, plus_infinity

CONTAINS


  FUNCTION GET_FILE_NUMBER() RESULT (num)
    !---------------------------------------------------------------------------
    ! purpose :
    !     gives a file number available for opening, i.e. a free logical unit
    !     example :
    !         n=GET_FILE_NUMBER()
    !         OPEN(n,file='sample.txt')
    ! subroutine/function needed :
    ! input variables :
    ! output variables :
    ! results :
    !     num -> 'LONG' integer : the first free logical unit
    !---------------------------------------------------------------------------
    IMPLICIT NONE
    INTEGER(KIND=LONG) :: num
    LOGICAL :: opened

    opened=.TRUE.
    num=9 ! start at file number 10
    DO WHILE (opened .AND. num < 99)
       num=num+1
       INQUIRE(UNIT=num, opened=opened)
    END DO
    IF (num == 99 .and. opened) STOP "*** WARNING : no file number available ***"
  END FUNCTION GET_FILE_NUMBER




  SUBROUTINE CREATE_ARCHIVE_FILE(shock_type,op_LTE,op,nH,Vs)
    !---------------------------------------------------------------------------
    ! called by :
    !     MHD
    ! purpose :
    !     creates an executable file which allows all the output files to be archived
    !     to the correct directory.
    ! subroutine/function needed :
    ! input variables :
    ! output variables :
    ! results :
    !     file named : archive_mhd.ex
    !---------------------------------------------------------------------------
    IMPLICIT NONE
    CHARACTER(len=1), INTENT(in) :: shock_type
    LOGICAL, INTENT(in) :: op_LTE
    REAL(KIND=DP), INTENT(in) :: op, nH, Vs
    CHARACTER(len=*),PARAMETER :: file_name = 'archive_mhd.ex'
    INTEGER(KIND=LONG) :: file_number
    CHARACTER(len=100) :: path, string_op, string_nH, string_Vs

    ! opening file
    file_number = GET_FILE_NUMBER()
    ! the file must exist and be an executable
    OPEN(file_number, file= file_name, status='REPLACE', action ='WRITE')


    ! define path
    WRITE(path,'("~/DATA_shock",A1,"/")')shock_type
    IF (op_LTE) THEN
       WRITE(string_op,'("op=LTE/")')
    ELSE
       WRITE(string_op,'("op=",ES7.1,"/")') op
    END IF
    WRITE(string_nH,'("nH=",ES7.1,"cm-3/")') nH
    WRITE(string_Vs,'("Vs=",ES7.1,"km.s-1/")') Vs

    ! write title
    WRITE(file_number,'(80("#"))')
    WRITE(file_number,'("##",T79,"##")')
    WRITE(file_number,'("##  execute this file to archive the outputs &
         &created by the MHD shock code",T79,"##")')
!!$    WRITE(file_number,'("##",10X,"warning : do not delete this file, &
!!$         &the code needs it!",T79,"##")')
    WRITE(file_number,'("##",T79,"##")')
    WRITE(file_number,'(80("#"))')
    WRITE(file_number,*)

    ! create directory
    WRITE(file_number,'("### creating the directory ###")')
    WRITE(file_number,'("mkdir ",A)')TRIM(path)
    WRITE(file_number,'("mkdir ",A)')TRIM(path) // TRIM(string_op)
    WRITE(file_number,'("mkdir ",A)')TRIM(path) // TRIM(string_op) &
         // TRIM(string_nH)
    WRITE(file_number,'("mkdir ",A)')TRIM(path) // TRIM(string_op) &
         //TRIM(string_nH) // TRIM(string_Vs)
    WRITE(file_number,*)

    ! mv files to the correct directory
    WRITE(file_number,'("### gzip and mv output files to this directory ###")')
    WRITE(file_number,'("gzip output/*.out ")')
    WRITE(file_number,'("mv output/*.out.gz ")',ADVANCE='NO')
    WRITE(file_number,'(A)',ADVANCE='NO')TRIM(path)
    WRITE(file_number,'(A)',ADVANCE='NO')TRIM(string_op)
    WRITE(file_number,'(A)',ADVANCE='NO')TRIM(string_nH)
    WRITE(file_number,'(A)')TRIM(string_Vs)
    WRITE(file_number,*)

    ! check the contents of the directory
    WRITE(file_number,'("### check the content of the directory ###")')
    WRITE(file_number,'("echo archive directory: ",A)') &
         TRIM(path) // TRIM(string_op) //TRIM(string_nH) // TRIM(string_Vs)
    WRITE(file_number,'("echo ---------------------------------")')
    WRITE(file_number,'("ls ",A)')TRIM(path) // TRIM(string_op) &
         //TRIM(string_nH) // TRIM(string_Vs)

    CLOSE(file_number)

  END SUBROUTINE CREATE_ARCHIVE_FILE





  FUNCTION INTEGRATE_SCALAR(dist_step, value1, value2) RESULT(res)
    !---------------------------------------------------------------------------
    ! purpose :
    !     integrates the function funct using trapezium rule; the independent
    !     variable is the distance (dist_step =distance-distance_old).
    !     Used to calculate column densities and flow times.
    ! subroutine/function needed :
    ! input variables :
    !     funct          -> 'function' to integrate (scalar)
    !     dist_step  -> = distance2 - distance1 (cm)
    !     value1, value2 -> values to integrate at distance1 and distance2
    ! output variables :
    ! results :
    !     res = 0.5_DP*dist_step*(value1+value2)
    !---------------------------------------------------------------------------
    IMPLICIT NONE
    REAL(KIND=DP), INTENT(in) :: dist_step, value1, value2
    REAL(KIND=DP) :: res

    res = 0.5_DP*dist_step*(value1+value2)

  END FUNCTION INTEGRATE_SCALAR



  SUBROUTINE SORT_INCREASING(dim,vect1,vect2,vect3,vect4)
    !---------------------------------------------------------------------------
    ! called by :
    ! purpose :
    !     sorts four vectors in the increasing order of the first vector
    !     method=shell
    ! subroutine/function needed :
    ! input variables :
    !     dim -> dimension of the four vectors, or number of elements to sort
    !     vect1, vect2, vect3, vect4 : vectors to sort
    ! output variables :
    !     vect1, vect2, vect3, vect4 : sorted vectors
    ! results :
    !---------------------------------------------------------------------------
    IMPLICIT NONE
    INTEGER(KIND=LONG), INTENT(in) :: dim
    REAL(KIND=DP),DIMENSION(:),INTENT(inout) :: vect1,vect2,vect3,vect4
    REAL(KIND=DP) :: v1,v2,v3,v4
    INTEGER(KIND=LONG) :: i,j,inc

    ! dim has to be <= dimension of the vectors
    IF (dim > SIZE(vect1)) STOP "*** WARNING, the vector is to small"
    ! the four vectors must have the same dimension
    IF ( SIZE(vect1)==SIZE(vect2) .AND. &
         SIZE(vect2)==SIZE(vect3) .AND. &
         SIZE(vect3)==SIZE(vect4)) THEN
       ! initial value of the increment
       inc=1
       DO WHILE (inc <= dim)
          inc=3*inc+1
       END DO

       ! sub-vectors distants by inc are sorted first and then inc is decreased
       DO WHILE (inc > 1)
          inc=inc/3
          DO i=inc+1,dim
             v1=vect1(i)
             v2=vect2(i)
             v3=vect3(i)
             v4=vect4(i)
             j=i
             DO WHILE ((vect1(j-inc) > v1) .AND. (j > inc))
                vect1(j)=vect1(j-inc)
                vect2(j)=vect2(j-inc)
                vect3(j)=vect3(j-inc)
                vect4(j)=vect4(j-inc)
                j=j-inc
             END DO
             vect1(j)=v1
             vect2(j)=v2
             vect3(j)=v3
             vect4(j)=v4
          END DO
       END DO
    ELSE
       STOP "***, WARNING, the four vectors must have the same dimension in SORT_INCREASING"
    END IF

  END SUBROUTINE SORT_INCREASING

  subroutine JACO
  end subroutine JACO

END MODULE MODULE_TOOLS
