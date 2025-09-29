
!!$ Author:
!!$   Yuanqing Wu, DGUT, P.R.China
!!$
!!$ History:
!!$   2025-5-9 by Yuanqing Wu
!!$
!!$ Support:
!!$   wuyuanq@gmail.com

module DBF_export2tecplot

    use DBF_globalData
    implicit none

contains

    subroutine exportPointCenteredIJ(data, varName, zoneName, title, fileName)

        real(kind=8), dimension(:,:), pointer, intent(in) :: data
        character(len=*), intent(in) :: varName, zoneName, title, fileName

        integer :: ix, iy, ierr

        open(unit=10, file=fileName, status="replace", action="write", iostat=ierr)
        if(ierr /= 0) then
            print *, "Failed to open '", fileName, "'"
            stop
        end if
        write(10, fmt=*) 'TITLE = "', trim(adjustl(title)), '"'
        write(10, fmt=*) 'VARIABLES = "X", "Y", "', trim(adjustl(varName)), '"'
        write(10,*) 'ZONE T="', trim(adjustl(zoneName)), '" DATAPACKING=POINT, I=', nx+1, ', J=', ny+1
        do iy = 1, ny+1
            do ix = 1, nx+1
                write(10,*) xs(ix), ys(iy), data(ix, iy)
            end do
        end do
        close(10)

    end subroutine exportPointCenteredIJ

    subroutine exportPointCenteredIJ_vars(dats, varNames, zoneName, title, fileName)

        real(kind=8), dimension(:,:,:), pointer, intent(in) :: dats
        character(len=*), dimension(:), intent(in) :: varNames
        character(len=*), intent(in) :: zoneName, title, fileName

        integer :: ix, iy, iVar, nVar, ierr

        nVar = size(varNames)
        open(unit=10, file=fileName, status="replace", action="write", iostat=ierr)
        if(ierr /= 0) then
            print *, "Failed to open '", fileName, "'" 
            stop
        end if
        write(10,fmt='(a,a,a)') 'TITLE = "', trim(adjustl(title)), '"'
        write(10,fmt='(a)', advance='no') 'VARIABLES = "X", "Y"'
        do iVar = 1, nVar 
            write(10,fmt='(a,a,a)', advance='no') ', "', trim(adjustl(varNames(iVar))), '"'
        end do 
        write(10,*)
        nVar = size(dats, 3)  
        write(10,*) 'ZONE T="', trim(adjustl(zoneName)), '" DATAPACKING=POINT, I=', nx+1, ', J=', ny+1
        do iy = 1, ny+1
            do ix = 1, nx+1
                write(10, fmt='(g15.5,g15.5)', advance='no') xs(ix), ys(iy)
                do iVar = 1, nVar
                    write(10, fmt='(g15.5)', advance='no') dats(ix, iy, iVar)
                end do
                write(10,*)
            end do
        end do
        close(10)

    end subroutine exportPointCenteredIJ_vars

    subroutine xEdgCtr2nodCtr(xEdgCtrDat, nodCtrDat)

        real(kind=8), dimension(:,:), pointer, intent(in) :: xEdgCtrDat
        real(kind=8), dimension(:,:), pointer, intent(in out) :: nodCtrDat

        integer :: iy

        nodCtrDat(:,:) = 0
        do iy = 0, 1
            nodCtrDat(:,1+iy:ny+iy) = nodCtrDat(:,1+iy:ny+iy) + xEdgCtrDat
        end do
        nodCtrDat(:,2:ny) = nodCtrDat(:,2:ny) / 2.D0

    end subroutine xEdgCtr2nodCtr

    subroutine yEdgCtr2nodCtr(yEdgCtrDat, nodCtrDat)

        real(kind=8), dimension(:,:), pointer, intent(in) :: yEdgCtrDat
        real(kind=8), dimension(:,:), pointer, intent(in out) :: nodCtrDat

        integer :: ix

        nodCtrDat(:,:) = 0
        do ix = 0, 1
            nodCtrDat(1+ix:nx+ix,:) = nodCtrDat(1+ix:nx+ix,:) + yEdgCtrDat
        end do
        nodCtrDat(2:nx,:) = nodCtrDat(2:nx,:) / 2.D0

    end subroutine yEdgCtr2nodCtr

    subroutine cellCtr2nodCtr(cellCtrDat, nodCtrDat)

        real(kind=8), dimension(:,:), pointer, intent(in) :: cellCtrDat
        real(kind=8), dimension(:,:), pointer, intent(in out) :: nodCtrDat

        integer :: ix, iy

        nodCtrDat(:,:) = 0
        do ix = 0, 1
            do iy = 0, 1
                nodCtrDat(1+ix:nx+ix,1+iy:ny+iy) = nodCtrDat(1+ix:nx+ix,1+iy:ny+iy) + cellCtrDat
            end do
        end do
        nodCtrDat(2:nx,:) = nodCtrDat(2:nx,:) / 2.D0
        nodCtrDat(:,2:ny) = nodCtrDat(:,2:ny) / 2.D0

    end subroutine cellCtr2nodCtr

    subroutine export2tecplot(g_poro, g_Sw, g_Kxx, g_vxw, g_vxn, g_vyw, g_vyn, g_p, g_Cf, g_Tem)

        real(kind=8), dimension(:,:), pointer, intent(in) :: g_poro
        real(kind=8), dimension(:,:), pointer, intent(in) :: g_Sw
        real(kind=8), dimension(:,:), pointer, intent(in) :: g_Kxx
        real(kind=8), dimension(:,:), pointer, intent(in) :: g_vxw, g_vxn
        real(kind=8), dimension(:,:), pointer, intent(in) :: g_vyw, g_vyn
        real(kind=8), dimension(:,:), pointer, intent(in) :: g_p
        real(kind=8), dimension(:,:), pointer, intent(in) :: g_Cf
        real(kind=8), dimension(:,:), pointer, intent(in) :: g_Tem

        character(len=10) :: chart
        real(kind=8), dimension(:,:), pointer :: poro_a
        real(kind=8), dimension(:,:), pointer :: Sw_a
        real(kind=8), dimension(:,:), pointer :: Kxx_a
        real(kind=8), dimension(:,:), pointer :: vxw_a, vxn_a
        real(kind=8), dimension(:,:), pointer :: vyw_a, vyn_a
        real(kind=8), dimension(:,:), pointer :: p_a
        real(kind=8), dimension(:,:), pointer :: Cf_a
        real(kind=8), dimension(:,:), pointer :: Tem_a
        real(kind=8), dimension(:,:,:), pointer :: vars_a

        allocate(poro_a(nx+1,ny+1))
        allocate(Sw_a(nx+1,ny+1))
        allocate(Kxx_a(nx+1,ny+1))
        allocate(vxw_a(nx+1,ny+1))
        allocate(vxn_a(nx+1,ny+1))
        allocate(vyw_a(nx+1,ny+1))
        allocate(vyn_a(nx+1,ny+1))
        allocate(p_a(nx+1,ny+1))
        allocate(Cf_a(nx+1,ny+1))
        allocate(Tem_a(nx+1,ny+1))
        allocate(vars_a(nx+1,ny+1,10))

        call cellCtr2nodCtr(g_poro, poro_a)
        call cellCtr2nodCtr(g_Sw, Sw_a)
        call cellCtr2nodCtr(g_Kxx, Kxx_a)
        call xEdgCtr2nodCtr(g_vxw, vxw_a)
        call xEdgCtr2nodCtr(g_vxn, vxn_a)
        call yEdgCtr2nodCtr(g_vyw, vyw_a)
        call yEdgCtr2nodCtr(g_vyn, vyn_a)
        call cellCtr2nodCtr(g_p, p_a)
        call cellCtr2nodCtr(g_Cf, Cf_a)
        call cellCtr2nodCtr(g_Tem, Tem_a)

        vars_a(:,:,1) = poro_a
        vars_a(:,:,2) = Sw_a
        vars_a(:,:,3) = Kxx_a
        vars_a(:,:,4) = vxw_a
        vars_a(:,:,5) = vxn_a
        vars_a(:,:,6) = vyw_a
        vars_a(:,:,7) = vyn_a
        vars_a(:,:,8) = p_a
        vars_a(:,:,9) = Cf_a
        vars_a(:,:,10) = Tem_a

        write(chart,'(i10)') t

        call exportPointCenteredIJ_vars(vars_a, (/"poro", "Sw  ", "Kxx ", "vxw ", "vxn ", "vyw ", "vyn ", "p   ", "Cf  ", "Tem "/), &!
            "zone1", "result from two-phase flow", trim(adjustl(soludoc))//"/out_twoPhase_"// &!
            trim(adjustl(chart))//".plt")

        deallocate(poro_a)
        deallocate(Sw_a)
        deallocate(Kxx_a)
        deallocate(vxw_a)
        deallocate(vxn_a)
        deallocate(vyw_a)
        deallocate(vyn_a)
        deallocate(p_a)
        deallocate(Cf_a)
        deallocate(Tem_a)
        deallocate(vars_a)

    end subroutine export2tecplot

end module DBF_export2tecplot

