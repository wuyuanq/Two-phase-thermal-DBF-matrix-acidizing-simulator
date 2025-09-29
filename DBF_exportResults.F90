
!!$ Author:
!!$   Yuanqing Wu, DGUT, P.R.China
!!$
!!$ History:
!!$   2025-5-9 by Yuanqing Wu
!!$
!!$ Support:
!!$   wuyuanq@gmail.com

module DBF_exportResults
  
    use DBF_globalData
    use DBF_export2Matlab
    use DBF_export2tecplot
    implicit none

contains

    subroutine exportResults(global_data)

        real(kind=8), dimension(:), pointer, intent(in) :: global_data

        character(len=10) :: chart
        integer :: p_pcol, p_prow
        integer :: p_xlower, p_ylower
        integer :: indexl, indexr, indexd, indexu
        real(kind=8), dimension(:,:), pointer :: g_poro
        real(kind=8), dimension(:,:), pointer :: g_Sw
        real(kind=8), dimension(:,:), pointer :: g_Kxx
        real(kind=8), dimension(:,:), pointer :: g_vxw, g_vxn
        real(kind=8), dimension(:,:), pointer :: g_vyw, g_vyn
        real(kind=8), dimension(:,:), pointer :: g_p
        real(kind=8), dimension(:,:), pointer :: g_Cf
        real(kind=8), dimension(:,:), pointer :: g_Tem
        character(len=40) :: fporotxt, fSwtxt, fKxxtxt, fvxwtxt, fvxntxt, &!
            fvywtxt, fvyntxt, fptxt, fCftxt, fTemtxt
        integer :: pid, ierr
        integer :: i, j, c

        write(chart,'(i10)') t
        fporotxt = trim(adjustl(soludoc))//'/soln_poro_raw_'//trim(adjustl(chart))//'.txt'
        fSwtxt = trim(adjustl(soludoc))//'/soln_Sw_raw_'//trim(adjustl(chart))//'.txt'
        fKxxtxt = trim(adjustl(soludoc))//'/soln_Kxx_raw_'//trim(adjustl(chart))//'.txt'
        fvxwtxt = trim(adjustl(soludoc))//'/soln_vxw_raw_'//trim(adjustl(chart))//'.txt'
        fvxntxt = trim(adjustl(soludoc))//'/soln_vxn_raw_'//trim(adjustl(chart))//'.txt'
        fvywtxt = trim(adjustl(soludoc))//'/soln_vyw_raw_'//trim(adjustl(chart))//'.txt'
        fvyntxt = trim(adjustl(soludoc))//'/soln_vyn_raw_'//trim(adjustl(chart))//'.txt'
        fptxt = trim(adjustl(soludoc))//'/soln_p_raw_'//trim(adjustl(chart))//'.txt'
        fCftxt = trim(adjustl(soludoc))//'/soln_Cf_raw_'//trim(adjustl(chart))//'.txt'
        fTemtxt = trim(adjustl(soludoc))//'/soln_Tem_raw_'//trim(adjustl(chart))//'.txt'

        allocate(g_poro(nx,ny))
        allocate(g_Sw(nx,ny))
        allocate(g_Kxx(nx,ny))
        allocate(g_vxw(nx+1,ny))
        allocate(g_vxn(nx+1,ny))
        allocate(g_vyw(nx,ny+1))
        allocate(g_vyn(nx,ny+1))
        allocate(g_p(nx,ny))
        allocate(g_Cf(nx,ny))
        allocate(g_Tem(nx,ny))

        c = 0
        do pid = 0, nProcs-1

            p_pcol = mod(pid,pncols)+1
            p_prow = pid/pncols+1

            p_xlower = (p_pcol-1)*localncols+1
            p_ylower = (p_prow-1)*localnrows+1

            do j = 1, localnrows
                do i = 1, localncols
                    c = c + 1
                    g_poro(p_xlower+i-1,p_ylower+j-1) = global_data(c)
                end do
            end do

            do j = 1, localnrows
                do i = 1, localncols
                    c = c + 1
                    g_Sw(p_xlower+i-1,p_ylower+j-1) = global_data(c)
                end do
            end do

            do j = 1, localnrows
                do i = 1, localncols
                    c = c + 1
                    g_Kxx(p_xlower+i-1,p_ylower+j-1) = global_data(c)
                end do
            end do

            indexl = 1
            if(p_pcol /= pncols) then
                indexr = localncols
            else
                indexr = localncols + 1
            end if
            indexd = 1
            indexu = localnrows

            do j = indexd, indexu
                do i = indexl, indexr
                    c = c + 1
                    g_vxw(p_xlower+i-1,p_ylower+j-1) = global_data(c)
                end do
            end do

            do j = indexd, indexu
                do i = indexl, indexr
                    c = c + 1
                    g_vxn(p_xlower+i-1,p_ylower+j-1) = global_data(c)
                end do
            end do

            indexl = 1
            indexr = localncols
            indexd = 1
            if(p_prow /= pnrows) then
                indexu = localnrows
            else
                indexu = localnrows + 1
            end if

            do j = indexd, indexu
                do i = indexl, indexr
                    c = c + 1
                    g_vyw(p_xlower+i-1,p_ylower+j-1) = global_data(c)
                end do
            end do

            do j = indexd, indexu
                do i = indexl, indexr
                    c = c + 1
                    g_vyn(p_xlower+i-1,p_ylower+j-1) = global_data(c)
                end do
            end do

            do j = 1, localnrows
                do i = 1, localncols
                    c = c + 1
                    g_p(p_xlower+i-1,p_ylower+j-1) = global_data(c)
                end do
            end do

            do j = 1, localnrows
                do i = 1, localncols
                    c = c + 1
                    g_Cf(p_xlower+i-1,p_ylower+j-1) = global_data(c)
                end do
            end do

            do j = 1, localnrows
                do i = 1, localncols
                    c = c + 1
                    g_Tem(p_xlower+i-1,p_ylower+j-1) = global_data(c)
                end do
            end do

        end do

        ! export results to the raw txt files
        open(unit=10, file=fporotxt, status='replace', iostat=ierr)
        if(ierr /= 0) then
            print *, 'open file ', fporotxt, ' error. ', ierr
            stop
        end if
        do j = 1, ny
            do i = 1, nx
                write(10, fmt='(es24.16)', iostat=ierr) g_poro(i,j)
                if(ierr /= 0) then
                    print *, 'write file ', fporotxt, ' error. ', ierr
                    stop
                end if
            end do
        end do
        close(10)

        open(unit=10, file=fSwtxt, status='replace', iostat=ierr)
        if(ierr /= 0) then
            print *, 'open file ', fSwtxt, ' error. ', ierr
            stop
        end if
        do j = 1, ny
            do i = 1, nx
                write(10, fmt='(es24.16)', iostat=ierr) g_Sw(i,j)
                if(ierr /= 0) then
                    print *, 'write file ', fSwtxt, ' error. ', ierr
                    stop
                end if
            end do
        end do
        close(10)

        open(unit=10, file=fKxxtxt, status='replace', iostat=ierr)
        if(ierr /= 0) then
            print *, 'open file ', fKxxtxt, ' error. ', ierr
            stop
        end if
        do j = 1, ny
            do i = 1, nx
                write(10, fmt='(es24.16)', iostat=ierr) g_Kxx(i,j)
                if(ierr /= 0) then
                    print *, 'write file ', fKxxtxt, ' error. ', ierr
                    stop
                end if
            end do
        end do
        close(10)

        open(unit=10, file=fvxwtxt, status='replace', iostat=ierr)
        if(ierr /= 0) then
            print *, 'open file ', fvxwtxt, ' error. ', ierr
            stop
        end if
        do j = 1, ny
            do i = 1, nx+1
                write(10, fmt='(es24.16)', iostat=ierr) g_vxw(i,j)
                if(ierr /= 0) then
                    print *, 'write file ', fvxwtxt, ' error. ', ierr
                    stop
                end if
            end do
        end do
        close(10)

        open(unit=10, file=fvxntxt, status='replace', iostat=ierr)
        if(ierr /= 0) then
            print *, 'open file ', fvxntxt, ' error. ', ierr
            stop
        end if
        do j = 1, ny
            do i = 1, nx+1
                write(10, fmt='(es24.16)', iostat=ierr) g_vxn(i,j)
                if(ierr /= 0) then
                    print *, 'write file ', fvxntxt, ' error. ', ierr
                    stop
                end if
            end do
        end do
        close(10)

        open(unit=10, file=fvywtxt, status='replace', iostat=ierr)
        if(ierr /= 0) then
            print *, 'open file ', fvywtxt, ' error. ', ierr
            stop
        end if
        do j = 1, ny+1
            do i = 1, nx
                write(10, fmt='(es24.16)', iostat=ierr) g_vyw(i,j)
                if(ierr /= 0) then
                    print *, 'write file ', fvywtxt, ' error. ', ierr
                    stop
                end if
            end do
        end do
        close(10)

        open(unit=10, file=fvyntxt, status='replace', iostat=ierr)
        if(ierr /= 0) then
            print *, 'open file ', fvyntxt, ' error. ', ierr
            stop
        end if
        do j = 1, ny+1
            do i = 1, nx
                write(10, fmt='(es24.16)', iostat=ierr) g_vyn(i,j)
                if(ierr /= 0) then
                    print *, 'write file ', fvyntxt, ' error. ', ierr
                    stop
                end if
            end do
        end do
        close(10)

        open(unit=10, file=fptxt, status='replace', iostat=ierr)
        if(ierr /= 0) then
            print *, 'open file ', fptxt, ' error. ', ierr
            stop
        end if
        do j = 1, ny
            do i = 1, nx
                write(10, fmt='(es24.16)', iostat=ierr) g_p(i,j)
                if(ierr /= 0) then
                    print *, 'write file ', fptxt, ' error. ', ierr
                    stop
                end if
            end do
        end do
        close(10)

        open(unit=10, file=fCftxt, status='replace', iostat=ierr)
        if(ierr /= 0) then
            print *, 'open file ', fCftxt, ' error. ', ierr
            stop
        end if
        do j = 1, ny
            do i = 1, nx
                write(10, fmt='(es24.16e3)', iostat=ierr) g_Cf(i,j)
                if(ierr /= 0) then
                    print *, 'write file ', fCftxt, ' error. ', ierr
                    stop
                end if
            end do
        end do
        close(10)

        open(unit=10, file=fTemtxt, status='replace', iostat=ierr)
        if(ierr /= 0) then
            print *, 'open file ', fTemtxt, ' error. ', ierr
            stop
        end if
        do j = 1, ny
            do i = 1, nx
                write(10, fmt='(es24.16e3)', iostat=ierr) g_Tem(i,j)
                if(ierr /= 0) then
                    print *, 'write file ', fTemtxt, ' error. ', ierr
                    stop
                end if
            end do
        end do
        close(10)

        ! export results to Matlab file
        call export2Matlab()

        ! export results to tecplot file
        call export2tecplot(g_poro, g_Sw, g_Kxx, g_vxw, g_vxn, g_vyw, g_vyn, g_p, g_Cf, g_Tem)

        deallocate(g_poro)
        deallocate(g_Sw)
        deallocate(g_Kxx)
        deallocate(g_vxw)
        deallocate(g_vxn)
        deallocate(g_vyw)
        deallocate(g_vyn)
        deallocate(g_p)
        deallocate(g_Cf)
        deallocate(g_Tem)

    end subroutine exportResults

end module DBF_exportResults
