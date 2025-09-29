
!!$ Author:
!!$   Yuanqing Wu, KAUST, Saudi Arabia
!!$
!!$ History:
!!$   2015-2-9 by Yuanqing Wu
!!$
!!$ Support:
!!$   wuyuanq@gmail.com
!!$
!!$ The module provides a subroutine which computes out an optimal scheme to 
!!$ allocate the cells to different processes so that the communication cost
!!$ is the lowest.

module RST_proceAlloc

contains

    subroutine proceAlloc(Np, local_nx, local_ny, local_pncols, local_pnrows)

        implicit none
        integer, intent(in) :: Np ! number of processes
        integer, intent(in) :: local_nx
        integer, intent(in) :: local_ny
        integer, intent(out) :: local_pncols
        integer, intent(out) :: local_pnrows
        integer(kind=8) :: minlen, curlen
        integer :: i, j

        minlen = local_nx*local_ny*4

        do i = 1, Np
            if(mod(Np,i) == 0) then
                j = Np/i
                if((mod(local_nx,i) == 0).and.(mod(local_ny,j) == 0)) then
                    curlen = (i-1)*local_ny + (j-1)*local_nx
                    if(curlen < minlen) then
                        minlen = curlen
                        local_pnrows = j
                        local_pncols = i
                    end if
                end if
            end if
        end do

        if(minlen == local_nx*local_ny*4) then
            print *, 'Error: the number of processes is not right. Please adjust it.'
            print *, 'You must make sure that each process has the same number of cells.'
            stop
        end if

        if((local_nx/local_pncols<2).or.(local_ny/local_pnrows<2)) then
            print *, 'Error: each process must include at least 2 cells in each direction.'
            stop
        end if

    end subroutine proceAlloc

end module RST_proceAlloc