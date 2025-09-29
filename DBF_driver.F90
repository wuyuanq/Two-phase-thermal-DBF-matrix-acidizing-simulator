
!!$ Author:
!!$   Yuanqing Wu, DGUT, P.R.China
!!$
!!$ History:
!!$   2025-5-9 by Yuanqing Wu
!!$
!!$ Support:
!!$   wuyuanq@gmail.com

module DBF_driver

    use DBF_globalData
    use DBF_resi
    use DBF_constructMat
    use DBF_exportResults

    implicit none
    include 'mpif.h'

contains

    subroutine genRandomNum()

        real(kind=8) :: rand
        real(kind=8), dimension(:), allocatable :: random
        integer :: ierr, i

        open(unit=10, file=trim(adjustl(FRANDOMTXT)), status='replace', iostat=ierr)
        if(ierr /= 0) then
            print *, 'open file error. ', ierr
            stop
        end if

        allocate(random(RANDOMSIZE))
        i = 1
        do while (i <= RANDOMSIZE)
            call random_number(rand)
            rand = 6.D-2*rand + 1.5D-1
            random(i) = rand
            i = i + 1
        end do
        write(10, fmt="(f8.6)") random(:)

        deallocate(random)
        close(10)

    end subroutine genRandomNum

    subroutine initialize()

        integer :: xmomSize, ymomSize, conSize, CfSize, SwSize, TemSize
        ! The arrays that will communicate between different functions must use pointer type instead of allocatalbe
        ! type. Pointer type can make sure that the subscripts of the arrays keep the same in the calling and called
        ! functions. However, if the subscripts of the arrays begin with 1, using allocatable type is also OK.
        integer :: indexl, indexr, indexd, indexu
        integer :: global_ind_b, global_ind_e
        logical :: alive
        integer :: i, j, n, c, ierr

        call MPI_Init(ierr)
        call MPI_Comm_size(MPI_COMM_WORLD, nProcs, ierr)
        call MPI_Comm_rank(MPI_COMM_WORLD, myid, ierr)

#if defined(MUMPS) || defined(HYPRE)
        call MPI_BUFFER_ATTACH(buffer, buffer_size, ierr)
#endif

        call MPI_Barrier(MPI_COMM_WORLD, ierr)
        timestart = MPI_Wtime()
        solvertime = 0.D0

        localncols = nx/pncols
        localnrows = ny/pnrows

        pcol = mod(myid,pncols)+1
        prow = myid/pncols+1

        xlower = (pcol-1)*localncols+1 
        xupper = pcol*localncols
        ylower = (prow-1)*localnrows+1
        yupper = prow*localnrows

        call index_convert_local_global(myid, 11, 1, 1, ilower_vp)
        call index_convert_local_global(myid, 3, localncols, localnrows, iupper_vp)
        local_x_size_vp = iupper_vp - ilower_vp + 1

        if((nProcs>1).and.(myid==0)) then
            allocate(slave_vp_data_size(nProcs-1))
            do i = 1, nProcs-1
                call index_convert_local_global(i, 11, 1, 1, global_ind_b)
                call index_convert_local_global(i, 3, localncols, localnrows, global_ind_e)
                slave_vp_data_size(i) = global_ind_e - global_ind_b + 1
            end do
        end if

        call index_convert_local_global(myid, 4, 1, 1, ilower_Cf)
        call index_convert_local_global(myid, 4, localncols, localnrows, iupper_Cf)
        local_x_size_Cf = iupper_Cf - ilower_Cf + 1

        call index_convert_local_global(myid, 5, 1, 1, ilower_Sw)
        call index_convert_local_global(myid, 5, localncols, localnrows, iupper_Sw)
        local_x_size_Sw = iupper_Sw - ilower_Sw + 1

        call index_convert_local_global(myid, 6, 1, 1, ilower_Tem)
        call index_convert_local_global(myid, 6, localncols, localnrows, iupper_Tem)
        local_x_size_Tem = iupper_Tem - ilower_Tem + 1

        t = 2

        ! initialize dm, kc, ks
        allocate(dm(0:localncols+1, 0:localnrows+1))
        allocate(kc(1:localncols, 1:localnrows))
        allocate(ks(1:localncols, 1:localnrows))

        ! initialize hx, hy
        if(pcol /= 1) then
            indexl = -1
        else
            indexl = 1
        end if
        if(pcol /= pncols) then
            indexr = localncols + 1
        else
            indexr = localncols
        end if
        allocate(hx(indexl:indexr))
        do i = indexl, indexr
            hx(i) = xs(xlower+i) - xs(xlower+i-1)
        end do

        if(prow /= 1) then
            indexd = -1
        else
            indexd = 1
        end if
        if(prow /= pnrows) then
            indexu = localnrows + 1
        else
            indexu = localnrows
        end if
        allocate(hy(indexd:indexu))
        do j = indexd, indexu
            hy(j) = ys(ylower+j) - ys(ylower+j-1)
        end do

        ! initialize poro
        allocate(poro(indexl:indexr,indexd:indexu))
        allocate(poro_old(indexl:indexr,indexd:indexu))
        do j = indexd, indexu
            do i = indexl, indexr
                poro(i,j) = poroInit(xlower+i-1,ylower+j-1)
                poro_old(i,j) = poro(i,j)
            end do
        end do

        ! initialize Sw
        allocate(Sw(indexl:indexr,indexd:indexu))
        allocate(Sw_old(indexl:indexr,indexd:indexu))
        do j = indexd, indexu
            do i = indexl, indexr
                Sw(i,j) = SwInit(xlower+i-1,ylower+j-1)
                Sw_old(i,j) = Sw(i,j)
            end do
        end do

        ! initialize poroEdge, SwEdge
        if(pcol /= 1) then
            indexl = 0
        else
            indexl = 1
        end if
        indexr = localncols + 1
        if(prow /= 1) then
            indexd = 0
        else
            indexd = 1
        end if
        if(prow /= pnrows) then
            indexu = localnrows + 1
        else
            indexu = localnrows
        end if

        allocate(poroEdgeX(indexl:indexr,indexd:indexu))
        allocate(poroEdgeX_old(indexl:indexr,indexd:indexu))
        allocate(poroEdgeXInit(indexl:indexr,indexd:indexu))

        allocate(SwEdgeX(indexl:indexr,indexd:indexu))
        allocate(SwEdgeX_old(indexl:indexr,indexd:indexu))

        if(pcol /= 1) then
            indexl = 0
        else
            indexl = 1
        end if
        if(pcol /= pncols) then
            indexr = localncols + 1
        else
            indexr = localncols
        end if
        if(prow /= 1) then
            indexd = 0
        else
            indexd = 1
        end if
        indexu = localnrows + 1

        allocate(poroEdgeY(indexl:indexr,indexd:indexu))
        allocate(poroEdgeY_old(indexl:indexr,indexd:indexu))
        allocate(poroEdgeYInit(indexl:indexr,indexd:indexu))

        allocate(SwEdgeY(indexl:indexr,indexd:indexu))
        allocate(SwEdgeY_old(indexl:indexr,indexd:indexu))

        ! initialize K
        if(pcol /= 1) then
            indexl = 0
        else
            indexl = 1
        end if
        indexr = localncols
        if(prow /= 1) then
            indexd = 0
        else
            indexd = 1
        end if
        indexu = localnrows
        allocate(Kxx(indexl:indexr,indexd:indexu))
        allocate(Kyy(indexl:indexr,indexd:indexu))
        do j = indexd, indexu
            do i = indexl, indexr
                Kxx(i,j) = KxxInit(xlower+i-1,ylower+j-1)
                Kyy(i,j) = KyyInit(xlower+i-1,ylower+j-1)
            end do
        end do

        ! initialize KEdge
        if(pcol /= pncols) then
            indexr = localncols
        else
            indexr = localncols + 1
        end if
        if(prow /= pnrows) then
            indexu = localnrows
        else
            indexu = localnrows + 1
        end if

        allocate(KxxEdge(1:indexr,1:localnrows))
        allocate(KyyEdge(1:localncols,1:indexu))

        ! initialize av
        indexl = 1
        indexr = localncols
        indexd = 1
        indexu = localnrows
        allocate(av(indexl:indexr,indexd:indexu))
        do j = indexd, indexu
            do i = indexl, indexr
                av(i,j) = avInit(xlower+i-1,ylower+j-1)
            end do
        end do

        ! initialize vxw, vyw, vxn, vyn
        indexl = 1
        indexr = localncols + 1
        if(prow /= 1) then
            indexd = 0
        else
            indexd = 1
        end if
        if(prow /= pnrows) then
            indexu = localnrows + 1
        else
            indexu = localnrows
        end if
        allocate(vxw(indexl:indexr,indexd:indexu))
        vxw(:,:) = 0.D0
        allocate(vxn(indexl:indexr,indexd:indexu))
        vxn(:,:) = 0.D0
        if(pcol == 1) then
            do j = 1, localnrows
                if(isDiriX0_p(ylower+j-1) == 0) then
                    vxw(1,j) = vxwBdryX0(ylower+j-1)
                    vxn(1,j) = vxnBdryX0(ylower+j-1)
                end if
            end do
        end if
        if(pcol == pncols) then
            do j = 1, localnrows
                if(isDiriX1_p(ylower+j-1) == 0) then
                    vxw(localncols+1,j) = vxwBdryX1(ylower+j-1)
                    vxn(localncols+1,j) = vxnBdryX1(ylower+j-1)
                end if
            end do
        end if

        if(pcol /= 1) then
            indexl = 0
        else
            indexl = 1
        end if
        if(pcol /= pncols) then
            indexr = localncols + 1
        else
            indexr = localncols
        end if
        indexd = 1
        indexu = localnrows + 1
        allocate(vyw(indexl:indexr,indexd:indexu))
        vyw(:,:) = 0.D0
        allocate(vyn(indexl:indexr,indexd:indexu))
        vyn(:,:) = 0.D0
        if(prow == 1) then
            do i = 1, localncols
                if(isDiriY0_p(xlower+i-1) == 0) then
                    vyw(i,1) = vywBdryY0(xlower+i-1)
                    vyn(i,1) = vynBdryY0(xlower+i-1)
                end if
            end do
        end if
        if(prow == pnrows) then
            do i = 1, localncols
                if(isDiriY1_p(xlower+i-1) == 0) then
                    vyw(i,localnrows+1) = vywBdryY1(xlower+i-1)
                    vyn(i,localnrows+1) = vynBdryY1(xlower+i-1)
                end if
            end do
        end if

        ! initialize p
        allocate(p(1:localncols,1:localnrows))
        p(1:localncols,1:localnrows) = pInit(xlower:xupper,ylower:yupper)

        ! initialize Cf
        allocate(Cf(1:localncols,1:localnrows))
        Cf(1:localncols,1:localnrows) = CfInit(xlower:xupper,ylower:yupper)

        ! initialize Tem
        allocate(Tem(1:localncols,1:localnrows))
        Tem(1:localncols,1:localnrows) = TemInit(xlower:xupper,ylower:yupper)

        ! compute the size of the equations
        if(pcol /= pncols) then
            xmomSize = localncols*localnrows
        else
            xmomSize = (localncols+1)*localnrows
        end if
        if(prow /= pnrows) then
            ymomSize = localncols*localnrows
        else
            ymomSize = localncols*(localnrows+1)
        end if
        conSize = localncols*localnrows
        CfSize = localncols*localnrows
        SwSize = localncols*localnrows
        TemSize = localncols*localnrows

        allocate(AxxEntryNum(xmomSize))
        allocate(AxpEntryNum(xmomSize))
        allocate(AyyEntryNum(ymomSize))
        allocate(AypEntryNum(ymomSize))
        allocate(AcxEntryNum(conSize))
        allocate(AcyEntryNum(conSize))
        allocate(ACfEntryNum(CfSize))
        allocate(ASwEntryNum(SwSize))
        allocate(ATemEntryNum(TemSize))

        AxxEntryNum(:) = 5
        AxpEntryNum(:) = 2
        AyyEntryNum(:) = 5
        AypEntryNum(:) = 2
        AcxEntryNum(:) = 2
        AcyEntryNum(:) = 2
        ACfEntryNum(:) = 9
        ASwEntryNum(:) = 5
        AtemEntryNum(:) = 5

        if(pcol /= pncols) then
            indexr = localncols
        else
            indexr = localncols + 1
        end if

        if(pcol == 1) then
            do n = 0, localnrows-1
                AxxEntryNum(1+indexr*n) = AxxEntryNum(1+indexr*n) - 1
            end do
        end if
        if(pcol == pncols) then
            do n = 1, localnrows
                AxxEntryNum(indexr*n) = AxxEntryNum(indexr*n) - 1
            end do
        end if
        if(prow == 1) then
            do n = 1, indexr
                AxxEntryNum(n) = AxxEntryNum(n) - 1
            end do
        end if
        if(prow == pnrows) then
            do n = 1, indexr
                AxxEntryNum(indexr*(localnrows-1)+n) = AxxEntryNum(indexr*(localnrows-1)+n) - 1
            end do
        end if
        do j = 1, localnrows
            if((isDiriX0_p(ylower+j-1)==0).and.(pcol==1)) then
                if(((j==1).and.(prow==1)).or.((j==localnrows).and.(prow==pnrows))) then
                    AxxEntryNum(1+indexr*(j-1)) = AxxEntryNum(1+indexr*(j-1)) - 2
                else
                    AxxEntryNum(1+indexr*(j-1)) = AxxEntryNum(1+indexr*(j-1)) - 3
                end if
            end if
            if((isDiriX1_p(ylower+j-1)==0).and.(pcol==pncols)) then
                if(((j==1).and.(prow==1)).or.((j==localnrows).and.(prow==pnrows))) then
                    AxxEntryNum(indexr*j) = AxxEntryNum(indexr*j) - 2
                else
                    AxxEntryNum(indexr*j) = AxxEntryNum(indexr*j) - 3
                end if
            end if
        end do

        if(pcol == 1) then
            do n = 0, localnrows-1
                AxpEntryNum(1+indexr*n) = AxpEntryNum(1+indexr*n) - 1
            end do
        end if
        if(pcol == pncols) then
            do n = 1, localnrows
                AxpEntryNum(indexr*n) = AxpEntryNum(indexr*n) - 1
            end do
        end if
        do j = 1, localnrows
            if((pcol==1).and.(isDiriX0_p(ylower+j-1)==0)) then
                AxpEntryNum(1+indexr*(j-1)) = AxpEntryNum(1+indexr*(j-1)) - 1
            end if
            if((pcol==pncols).and.(isDiriX1_p(ylower+j-1)==0)) then
                AxpEntryNum(indexr*j) = AxpEntryNum(indexr*j) - 1
            end if
        end do

        if(prow /= pnrows) then
            indexu = localnrows
        else
            indexu = localnrows + 1
        end if

        if(prow == 1) then
            do n = 1, localncols
                AyyEntryNum(n) = AyyEntryNum(n) - 1
            end do
        end if
        if(prow == pnrows) then
            do n = 1, localncols
                AyyEntryNum(localncols*localnrows+n) = AyyEntryNum(localncols*localnrows+n) - 1
            end do
        end if
        if(pcol == 1) then
            do n = 1, indexu
                AyyEntryNum(1+(n-1)*localncols) = AyyEntryNum(1+(n-1)*localncols) - 1
            end do
        end if
        if(pcol == pncols) then
            do n = 1, indexu
                AyyEntryNum(n*localncols) = AyyEntryNum(n*localncols) - 1
            end do
        end if
        do i = 1, localncols
            if((isDiriY0_p(xlower+i-1)==0).and.(prow==1)) then
                if(((i==1).and.(pcol==1)).or.((i==localncols).and.(pcol==pncols))) then
                    AyyEntryNum(i) = AyyEntryNum(i) - 2
                else
                    AyyEntryNum(i) = AyyEntryNum(i) - 3
                end if
            end if
            if((isDiriY1_p(xlower+i-1)==0).and.(prow==pnrows)) then
                if(((i==1).and.(pcol==1)).or.((i==localncols).and.(pcol==pncols))) then
                    AyyEntryNum(localncols*localnrows+i) = AyyEntryNum(localncols*localnrows+i) - 2
                else
                    AyyEntryNum(localncols*localnrows+i) = AyyEntryNum(localncols*localnrows+i) - 3
                end if
            end if
        end do

        if(prow == 1) then
            do n = 1, localncols
                AypEntryNum(n) = AypEntryNum(n) - 1
            end do
        end if
        if(prow == pnrows) then
            do n = 1, localncols
                AypEntryNum(localncols*localnrows+n) = AypEntryNum(localncols*localnrows+n) - 1
            end do
        end if
        do i = 1, localncols
            if((prow==1).and.(isDiriY0_p(xlower+i-1)==0)) then
                AypEntryNum(i) = AypEntryNum(i) - 1
            end if
            if((prow==pnrows).and.(isDiriY1_p(xlower+i-1)==0)) then
                AypEntryNum(localncols*localnrows+i) = AypEntryNum(localncols*localnrows+i) - 1
            end if
        end do

        if(pcol == 1) then
            do n = 0, localnrows-1
                ACfEntryNum(1+localncols*n) = ACfEntryNum(1+localncols*n) - 3
            end do
        end if
        if(pcol == pncols) then
            do n = 1, localnrows
                ACfEntryNum(localncols*n) = ACfEntryNum(localncols*n) - 3
            end do
        end if
        if(prow == 1) then
            do n = 1, localncols
                ACfEntryNum(n) = ACfEntryNum(n) - 3
            end do
        end if
        if(prow == pnrows) then
            do n = 1, localncols
                ACfEntryNum(localncols*(localnrows-1)+n) = ACfEntryNum(localncols*(localnrows-1)+n) - 3
            end do
        end if
        if((pcol==1).and.(prow==1)) then
            ACfEntryNum(1) = ACfEntryNum(1) + 1
        end if
        if((pcol==1).and.(prow==pnrows)) then
            ACfEntryNum(1+localncols*(localnrows-1)) = ACfEntryNum(1+localncols*(localnrows-1)) + 1
        end if
        if((pcol==pncols).and.(prow==1)) then
            ACfEntryNum(localncols) = ACfEntryNum(localncols) + 1
        end if
        if((pcol==pncols).and.(prow==pnrows)) then
            ACfEntryNum(localncols*localnrows) = ACfEntryNum(localncols*localnrows) + 1
        end if

        if(pcol == 1) then
            do n = 0, localnrows-1
                ASwEntryNum(1+localncols*n) = ASwEntryNum(1+localncols*n) - 1
            end do
        end if
        if(pcol == pncols) then
            do n = 1, localnrows
                ASwEntryNum(localncols*n) = ASwEntryNum(localncols*n) - 1
            end do
        end if
        if(prow == 1) then
            do n = 1, localncols
                ASwEntryNum(n) = ASwEntryNum(n) - 1
            end do
        end if
        if(prow == pnrows) then
            do n = 1, localncols
                ASwEntryNum(localncols*(localnrows-1)+n) = ASwEntryNum(localncols*(localnrows-1)+n) - 1
            end do
        end if

        if(pcol == 1) then
            do n = 0, localnrows-1
                ATemEntryNum(1+localncols*n) = ATemEntryNum(1+localncols*n) - 1
            end do
        end if
        if(pcol == pncols) then
            do n = 1, localnrows
                ATemEntryNum(localncols*n) = ATemEntryNum(localncols*n) - 1
            end do
        end if
        if(prow == 1) then
            do n = 1, localncols
                ATemEntryNum(n) = ATemEntryNum(n) - 1
            end do
        end if
        if(prow == pnrows) then
            do n = 1, localncols
                ATemEntryNum(localncols*(localnrows-1)+n) = ATemEntryNum(localncols*(localnrows-1)+n) - 1
            end do
        end if

        allocate(AxxEntryBase(xmomSize))
        allocate(AxpEntryBase(xmomSize))
        allocate(AyyEntryBase(ymomSize))
        allocate(AypEntryBase(ymomSize))
        allocate(AcxEntryBase(conSize))
        allocate(AcyEntryBase(conSize))
        allocate(ACfEntryBase(CfSize))
        allocate(ASwEntryBase(SwSize))
        allocate(ATemEntryBase(TemSize))

        AxxEntryBase(1) = 1
        AxpEntryBase(1) = 1
        do n = 2, xmomSize
            AxxEntryBase(n) = AxxEntryBase(n-1) + AxxEntryNum(n-1)
            AxpEntryBase(n) = AxpEntryBase(n-1) + AxpEntryNum(n-1)
        end do

        AyyEntryBase(1) = 1
        AypEntryBase(1) = 1
        do n = 2, ymomSize
            AyyEntryBase(n) = AyyEntryBase(n-1) + AyyEntryNum(n-1)
            AypEntryBase(n) = AypEntryBase(n-1) + AypEntryNum(n-1)
        end do

        AcxEntryBase(1) = 1
        AcyEntryBase(1) = 1
        do n = 2, conSize
            AcxEntryBase(n) = AcxEntryBase(n-1) + AcxEntryNum(n-1)
            AcyEntryBase(n) = AcyEntryBase(n-1) + AcyEntryNum(n-1)
        end do

        ACfEntryBase(1) = 1
        do n = 2, CfSize
            ACfEntryBase(n) = ACfEntryBase(n-1) + ACfEntryNum(n-1)
        end do

        ASwEntryBase(1) = 1
        do n = 2, SwSize
            ASwEntryBase(n) = ASwEntryBase(n-1) + ASwEntryNum(n-1)
        end do
        
        ATemEntryBase(1) = 1
        do n = 2, TemSize
            ATemEntryBase(n) = ATemEntryBase(n-1) + ATemEntryNum(n-1)
        end do

        ! compute matrix size
        AxxSize = 0
        do n = 1, xmomSize
            AxxSize = AxxSize + AxxEntryNum(n)
        end do
        AxpSize = 0
        do n = 1, xmomSize
            AxpSize = AxpSize + AxpEntryNum(n)
        end do
        AyySize = 0
        do n = 1, ymomSize
            AyySize = AyySize + AyyEntryNum(n)
        end do
        AypSize = 0
        do n = 1, ymomSize
            AypSize = AypSize + AypEntryNum(n)
        end do
        AcxSize = 2*localncols*localnrows
        AcySize = 2*localncols*localnrows
        ACfSize = 0
        do n = 1, CfSize
            ACfSize = ACfSize + ACfEntryNum(n)
        end do
        ASwSize = 0
        do n = 1, SwSize
            ASwSize = ASwSize + ASwEntryNum(n)
        end do
        ATemSize = 0
        do n = 1, TemSize
            ATemSize = ATemSize + ATemEntryNum(n)
        end do

        ! initialize matrix
        allocate(local_rhs_vp(local_x_size_vp))
        allocate(local_rhs_vp_static(local_x_size_vp))
        allocate(local_rhs_Cf(local_x_size_Cf))
        allocate(local_rhs_Sw(local_x_size_Sw))
        allocate(local_rhs_Tem(local_x_size_Tem))

        allocate(AxxwCols(AxxSize))
        allocate(AxxwRows(AxxSize))
        allocate(AxxwStaticValues(AxxSize))
        allocate(AxxwDynValues(AxxSize))
        allocate(AxxwValues(AxxSize))
        allocate(AxpwCols(AxpSize))
        allocate(AxpwRows(AxpSize))
        allocate(AxpwValues(AxpSize))
        allocate(AyywCols(AyySize))
        allocate(AyywRows(AyySize))
        allocate(AyywStaticValues(AyySize))
        allocate(AyywDynValues(AyySize))
        allocate(AyywValues(AyySize))
        allocate(AypwCols(AypSize))
        allocate(AypwRows(AypSize))
        allocate(AypwValues(AypSize))

        allocate(AxxnCols(AxxSize))
        allocate(AxxnRows(AxxSize))
        allocate(AxxnStaticValues(AxxSize))
        allocate(AxxnDynValues(AxxSize))
        allocate(AxxnValues(AxxSize))
        allocate(AxpnCols(AxpSize))
        allocate(AxpnRows(AxpSize))
        allocate(AxpnValues(AxpSize))
        allocate(AyynCols(AyySize))
        allocate(AyynRows(AyySize))
        allocate(AyynStaticValues(AyySize))
        allocate(AyynDynValues(AyySize))
        allocate(AyynValues(AyySize))
        allocate(AypnCols(AypSize))
        allocate(AypnRows(AypSize))
        allocate(AypnValues(AypSize))

        allocate(AcxwCols(AcxSize))
        allocate(AcxwRows(AcxSize))
        allocate(AcxwValues(AcxSize))
        allocate(AcywCols(AcySize))
        allocate(AcywRows(AcySize))
        allocate(AcywValues(AcySize))

        allocate(AcxnCols(AcxSize))
        allocate(AcxnRows(AcxSize))
        allocate(AcxnValues(AcxSize))
        allocate(AcynCols(AcySize))
        allocate(AcynRows(AcySize))
        allocate(AcynValues(AcySize))
       
        allocate(ACfCols(ACfSize))
        allocate(ACfRows(ACfSize))
        allocate(ACfValues(ACfSize))

        allocate(ASwCols(ASwSize))
        allocate(ASwRows(ASwSize))
        allocate(ASwValues(ASwSize))

        allocate(ATemCols(ATemSize))
        allocate(ATemRows(ATemSize))
        allocate(ATemValues(ATemSize))

        AxxwCols(:) = 0
        AxxwRows(:) = 0
        AxxwStaticValues(:) = 0.D0
        AxxwDynValues(:) = 0.D0
        AxxwValues(:) = 0.D0
        AxpwCols(:) = 0
        AxpwRows(:) = 0
        AxpwValues(:) = 0.D0
        AyywCols(:) = 0
        AyywRows(:) = 0
        AyywStaticValues(:) = 0.D0
        AyywDynValues(:) = 0.D0
        AyywValues(:) = 0.D0
        AypwCols(:) = 0
        AypwRows(:) = 0
        AypwValues(:) = 0.D0

        AxxnCols(:) = 0
        AxxnRows(:) = 0
        AxxnStaticValues(:) = 0.D0
        AxxnDynValues(:) = 0.D0
        AxxnValues(:) = 0.D0
        AxpnCols(:) = 0
        AxpnRows(:) = 0
        AxpnValues(:) = 0.D0
        AyynCols(:) = 0
        AyynRows(:) = 0
        AyynStaticValues(:) = 0.D0
        AyynDynValues(:) = 0.D0
        AyynValues(:) = 0.D0
        AypnCols(:) = 0
        AypnRows(:) = 0
        AypnValues(:) = 0.D0

        AcxwCols(:) = 0
        AcxwRows(:) = 0
        AcxwValues(:) = 0.D0
        AcywCols(:) = 0
        AcywRows(:) = 0
        AcywValues(:) = 0.D0

        AcxnCols(:) = 0
        AcxnRows(:) = 0
        AcxnValues(:) = 0.D0
        AcynCols(:) = 0
        AcynRows(:) = 0
        AcynValues(:) = 0.D0
       
        ACfCols(:) = 0
        ACfRows(:) = 0
        ACfValues(:) = 0.D0

        ASwCols(:) = 0
        ASwRows(:) = 0
        ASwValues(:) = 0.D0

        ATemCols(:) = 0
        ATemRows(:) = 0
        ATemValues(:) = 0.D0

        presDropInit = 0.D0
        isFindPresDropInit = .false.

        ! open the data output files
        if(myid == 0) then
            inquire(file = trim(adjustl(soludoc)), exist = alive)
            if(.not.alive) then
                call system("mkdir "//trim(adjustl(soludoc)))
            end if
            open(unit=40, file=trim(adjustl(soludoc))//'/his_poro_avg.txt', status='replace', iostat=ierr)
            if(ierr /= 0) then
                print *, 'open file ', trim(adjustl(soludoc))//'/his_poro_avg.txt', ' error. ', ierr
                stop
            end if
            open(unit=41, file=trim(adjustl(soludoc))//'/his_Kxx_avg.txt', status='replace', iostat=ierr)
            if(ierr /= 0) then
                print *, 'open file ', trim(adjustl(soludoc))//'/his_Kxx_avg.txt', ' error. ', ierr
                stop
            end if
            open(unit=42, file=trim(adjustl(soludoc))//'/his_av_avg.txt', status='replace', iostat=ierr)
            if(ierr /= 0) then
                print *, 'open file ', trim(adjustl(soludoc))//'/his_av_avg.txt', ' error. ', ierr
                stop
            end if
            open(unit=43, file=trim(adjustl(soludoc))//'/his_p_avg.txt', status='replace', iostat=ierr)
            if(ierr /= 0) then
                print *, 'open file ', trim(adjustl(soludoc))//'/his_p_avg.txt', ' error. ', ierr
                stop
            end if
            open(unit=44, file=trim(adjustl(soludoc))//'/his_Cf_avg.txt', status='replace', iostat=ierr)
            if(ierr /= 0) then
                print *, 'open file ', trim(adjustl(soludoc))//'/his_Cf_avg.txt', ' error. ', ierr
                stop
            end if
            open(unit=45, file=trim(adjustl(soludoc))//'/his_Tem_avg.txt', status='replace', iostat=ierr)
            if(ierr /= 0) then
                print *, 'open file ', trim(adjustl(soludoc))//'/his_Tem_avg.txt', ' error. ', ierr
                stop
            end if
            open(unit=46, file=trim(adjustl(soludoc))//'/his_q_avg.txt', status='replace', iostat=ierr)
            if(ierr /= 0) then
                print *, 'open file ', trim(adjustl(soludoc))//'/his_q_avg.txt', ' error. ', ierr
                stop
            end if
            open(unit=47, file=trim(adjustl(soludoc))//'/his_lp_avg.txt', status='replace', iostat=ierr)
            if(ierr /= 0) then
                print *, 'open file ', trim(adjustl(soludoc))//'/his_lp_avg.txt', ' error. ', ierr
                stop
            end if
        end if

#ifdef LAPACK

        allocate(A_lapack_vp(local_x_size_vp,local_x_size_vp))
        allocate(b_lapack_vp(local_x_size_vp))
        allocate(IPIV_vp(local_x_size_vp))

        allocate(A_lapack_Cf(local_x_size_Cf,local_x_size_Cf))
        allocate(b_lapack_Cf(local_x_size_Cf))
        allocate(IPIV_Cf(local_x_size_Cf))

        allocate(A_lapack_Sw(local_x_size_Sw,local_x_size_Sw))
        allocate(b_lapack_Sw(local_x_size_Sw))
        allocate(IPIV_Sw(local_x_size_Sw))

        allocate(A_lapack_Tem(local_x_size_Tem,local_x_size_Tem))
        allocate(b_lapack_Tem(local_x_size_Tem))
        allocate(IPIV_Tem(local_x_size_Tem))

#elif defined(UMFPACK)

        allocate(Ap_vp(1:local_x_size_vp+1))
        allocate(Ai_vp(1:8*local_x_size_vp))
        allocate(Ax_vp(1:8*local_x_size_vp))

        allocate(Ap_Cf(1:local_x_size_Cf+1))
        allocate(Ai_Cf(1:9*local_x_size_Cf))
        allocate(Ax_Cf(1:9*local_x_size_Cf))

        allocate(Ap_Sw(1:local_x_size_Sw+1))
        allocate(Ai_Sw(1:5*local_x_size_Sw))
        allocate(Ax_Sw(1:5*local_x_size_Sw))

        allocate(Ap_Tem(1:local_x_size_Tem+1))
        allocate(Ai_Tem(1:5*local_x_size_Tem))
        allocate(Ax_Tem(1:5*local_x_size_Tem))

#elif defined(MUMPS)

        mumps_par_vp%COMM = MPI_COMM_WORLD
        mumps_par_vp%SYM = 0
        mumps_par_vp%PAR = 1
        mumps_par_vp%JOB = -1
        call DMUMPS(mumps_par_vp)

        mumps_par_vp%ICNTL(4) = 0
        mumps_par_vp%ICNTL(7) = 5
        mumps_par_vp%ICNTL(18) = 3
        if(mumps_par_vp%MYID == 0) then
            mumps_par_vp%N = (nx+1)*ny*2+nx*(ny+1)*2+nx*ny
        end if
        allocate(mumps_par_vp%RHS(mumps_par_vp%N))

        allocate(mumps_IRN_loc_vp(local_x_size_vp*8))
        allocate(mumps_JCN_loc_vp(local_x_size_vp*8))
        allocate(mumps_A_loc_vp(local_x_size_vp*8))

        mumps_par_Cf%COMM = MPI_COMM_WORLD
        mumps_par_Cf%SYM = 0
        mumps_par_Cf%PAR = 1
        mumps_par_Cf%JOB = -1
        call DMUMPS(mumps_par_Cf)

        mumps_par_Cf%ICNTL(4) = 0
        mumps_par_Cf%ICNTL(7) = 5
        mumps_par_Cf%ICNTL(18) = 3
        if(mumps_par_Cf%MYID == 0) then
            mumps_par_Cf%N = nx*ny
        end if
        allocate(mumps_par_Cf%RHS(mumps_par_Cf%N))

        allocate(mumps_IRN_loc_Cf(local_x_size_Cf*9))
        allocate(mumps_JCN_loc_Cf(local_x_size_Cf*9))
        allocate(mumps_A_loc_Cf(local_x_size_Cf*9))

        mumps_par_Sw%COMM = MPI_COMM_WORLD
        mumps_par_Sw%SYM = 0
        mumps_par_Sw%PAR = 1
        mumps_par_Sw%JOB = -1
        call DMUMPS(mumps_par_Sw)

        mumps_par_Sw%ICNTL(4) = 0
        mumps_par_Sw%ICNTL(7) = 5
        mumps_par_Sw%ICNTL(18) = 3
        if(mumps_par_Sw%MYID == 0) then
            mumps_par_Sw%N = nx*ny
        end if
        allocate(mumps_par_Sw%RHS(mumps_par_Sw%N))

        allocate(mumps_IRN_loc_Sw(local_x_size_Sw*5))
        allocate(mumps_JCN_loc_Sw(local_x_size_Sw*5))
        allocate(mumps_A_loc_Sw(local_x_size_Sw*5))

        mumps_par_Tem%COMM = MPI_COMM_WORLD
        mumps_par_Tem%SYM = 0
        mumps_par_Tem%PAR = 1
        mumps_par_Tem%JOB = -1
        call DMUMPS(mumps_par_Tem)

        mumps_par_Tem%ICNTL(4) = 0
        mumps_par_Tem%ICNTL(7) = 5
        mumps_par_Tem%ICNTL(18) = 3
        if(mumps_par_Tem%MYID == 0) then
            mumps_par_Tem%N = nx*ny
        end if
        allocate(mumps_par_Tem%RHS(mumps_par_Tem%N))

        allocate(mumps_IRN_loc_Tem(local_x_size_Tem*5))
        allocate(mumps_JCN_loc_Tem(local_x_size_Tem*5))
        allocate(mumps_A_loc_Tem(local_x_size_Tem*5))

#elif defined(HYPRE)

        if(prow /= 1) then
            call index_convert_local_global(myid-pncols, 11, 1, localnrows, jlower_vp)
        elseif(pcol /= 1) then
            call index_convert_local_global(myid-1, 11, localncols, 1, jlower_vp)
        else
            jlower_vp = ilower_vp
        end if
        if(prow /= pnrows) then
            call index_convert_local_global(myid+pncols, 3, localncols, 1, jupper_vp)
        elseif(pcol /= pncols) then
            call index_convert_local_global(myid+1, 3, 1, localnrows, jupper_vp)
        else
            jupper_vp = iupper_vp
        end if

        call HYPRE_IJMatrixCreate(MPI_COMM_WORLD, ilower_vp, iupper_vp, jlower_vp, jupper_vp, A_vp, ierr)
        call HYPRE_IJMatrixSetObjectType(A_vp, HYPRE_PARCSR, ierr)
        call HYPRE_IJMatrixInitialize(A_vp, ierr)

        call HYPRE_IJVectorCreate(MPI_COMM_WORLD, ilower_vp, iupper_vp, b_vp, ierr)
        call HYPRE_IJVectorSetObjectType(b_vp, HYPRE_PARCSR, ierr)
        call HYPRE_IJVectorInitialize(b_vp, ierr)

        call HYPRE_IJVectorCreate(MPI_COMM_WORLD, ilower_vp, iupper_vp, x_vp, ierr)
        call HYPRE_IJVectorSetObjectType(x_vp, HYPRE_PARCSR, ierr)
        call HYPRE_IJVectorInitialize(x_vp, ierr)

        if((pcol /= 1).and.(prow /= 1)) then
            call index_convert_local_global(myid-1-pncols, 4, localncols, localnrows, jlower_Cf)
        elseif((pcol /= 1).and.(prow == 1)) then
            call index_convert_local_global(myid-1, 4, localncols, 1, jlower_Cf)
        elseif((pcol == 1).and.(prow /= 1)) then
            call index_convert_local_global(myid-pncols, 4, 1, localnrows, jlower_Cf)
        else
            jlower_Cf = ilower_Cf
        end if
        if((pcol /= pncols).and.(prow /= pnrows)) then
            call index_convert_local_global(myid+1+pncols, 4, 1, 1, jupper_Cf)
        elseif((pcol /= pncols).and.(prow == pnrows)) then
            call index_convert_local_global(myid+1, 4, 1, localnrows, jupper_Cf)
        elseif((pcol == pncols).and.(prow /= pnrows)) then
            call index_convert_local_global(myid+pncols, 4, localncols, 1, jupper_Cf)
        else
            jupper_Cf = iupper_Cf
        end if

        call HYPRE_IJMatrixCreate(MPI_COMM_WORLD, ilower_Cf, iupper_Cf, jlower_Cf, jupper_Cf, A_Cf, ierr)
        call HYPRE_IJMatrixSetObjectType(A_Cf, HYPRE_PARCSR, ierr)
        call HYPRE_IJMatrixInitialize(A_Cf, ierr)

        call HYPRE_IJVectorCreate(MPI_COMM_WORLD, ilower_Cf, iupper_Cf, b_Cf, ierr)
        call HYPRE_IJVectorSetObjectType(b_Cf, HYPRE_PARCSR, ierr)
        call HYPRE_IJVectorInitialize(b_Cf, ierr)

        call HYPRE_IJVectorCreate(MPI_COMM_WORLD, ilower_Cf, iupper_Cf, x_Cf, ierr)
        call HYPRE_IJVectorSetObjectType(x_Cf, HYPRE_PARCSR, ierr)
        call HYPRE_IJVectorInitialize(x_Cf, ierr)

        if(prow /= 1) then
            call index_convert_local_global(myid-pncols, 5, 1, localnrows, jlower_Sw)
        elseif(pcol /= 1) then
            call index_convert_local_global(myid-1, 5, localncols, 1, jlower_Sw)
        else
            jlower_Sw = ilower_Sw
        end if
        if(prow /= pnrows) then
            call index_convert_local_global(myid+pncols, 5, localncols, 1, jupper_Sw)
        elseif(pcol /= pncols) then
            call index_convert_local_global(myid+1, 5, 1, localnrows, jupper_Sw)
        else
            jupper_Sw = iupper_Sw
        end if

        call HYPRE_IJMatrixCreate(MPI_COMM_WORLD, ilower_Sw, iupper_Sw, jlower_Sw, jupper_Sw, A_Sw, ierr)
        call HYPRE_IJMatrixSetObjectType(A_Sw, HYPRE_PARCSR, ierr)
        call HYPRE_IJMatrixInitialize(A_Sw, ierr)

        call HYPRE_IJVectorCreate(MPI_COMM_WORLD, ilower_Sw, iupper_Sw, b_Sw, ierr)
        call HYPRE_IJVectorSetObjectType(b_Sw, HYPRE_PARCSR, ierr)
        call HYPRE_IJVectorInitialize(b_Sw, ierr)

        call HYPRE_IJVectorCreate(MPI_COMM_WORLD, ilower_Sw, iupper_Sw, x_Sw, ierr)
        call HYPRE_IJVectorSetObjectType(x_Sw, HYPRE_PARCSR, ierr)
        call HYPRE_IJVectorInitialize(x_Sw, ierr)

        if(prow /= 1) then
            call index_convert_local_global(myid-pncols, 6, 1, localnrows, jlower_Tem)
        elseif(pcol /= 1) then
            call index_convert_local_global(myid-1, 6, localncols, 1, jlower_Tem)
        else
            jlower_Tem = ilower_Tem
        end if
        if(prow /= pnrows) then
            call index_convert_local_global(myid+pncols, 6, localncols, 1, jupper_Tem)
        elseif(pcol /= pncols) then
            call index_convert_local_global(myid+1, 6, 1, localnrows, jupper_Tem)
        else
            jupper_Tem = iupper_Tem
        end if

        call HYPRE_IJMatrixCreate(MPI_COMM_WORLD, ilower_Tem, iupper_Tem, jlower_Tem, jupper_Tem, A_Tem, ierr)
        call HYPRE_IJMatrixSetObjectType(A_Tem, HYPRE_PARCSR, ierr)
        call HYPRE_IJMatrixInitialize(A_Tem, ierr)

        call HYPRE_IJVectorCreate(MPI_COMM_WORLD, ilower_Tem, iupper_Tem, b_Tem, ierr)
        call HYPRE_IJVectorSetObjectType(b_Tem, HYPRE_PARCSR, ierr)
        call HYPRE_IJVectorInitialize(b_Tem, ierr)

        call HYPRE_IJVectorCreate(MPI_COMM_WORLD, ilower_Tem, iupper_Tem, x_Tem, ierr)
        call HYPRE_IJVectorSetObjectType(x_Tem, HYPRE_PARCSR, ierr)
        call HYPRE_IJVectorInitialize(x_Tem, ierr)

        call HYPRE_ParCSRPilutCreate(MPI_COMM_WORLD, precond_vp, ierr)
        call HYPRE_ParCSRGMRESCreate(MPI_COMM_WORLD, solver_vp, ierr)
        call HYPRE_ParCSRGMRESSetTol(solver_vp, 1.D-7, ierr)
        call HYPRE_ParCSRGMRESSetPrintLevel(solver_vp, 2, ierr)
        call HYPRE_ParCSRGMRESSetLogging(solver_vp, 1, ierr)
        ! 1 means the DS preconditioner, 4 means the ParaSails preconditioner
        call HYPRE_ParCSRGMRESSetPrecond(solver_vp, 3, precond_vp, ierr)

        call HYPRE_ParCSRPilutCreate(MPI_COMM_WORLD, precond_Cf, ierr)
        call HYPRE_ParCSRGMRESCreate(MPI_COMM_WORLD, solver_Cf, ierr)
        call HYPRE_ParCSRGMRESSetTol(solver_Cf, 1.D-7, ierr)
        call HYPRE_ParCSRGMRESSetPrintLevel(solver_Cf, 2, ierr)
        call HYPRE_ParCSRGMRESSetLogging(solver_Cf, 1, ierr)
        ! 1 means the DS preconditioner, 4 means the ParaSails preconditioner
        call HYPRE_ParCSRGMRESSetPrecond(solver_Cf, 1, precond_Cf, ierr)

        call HYPRE_EuclidCreate(MPI_COMM_WORLD, precond_Sw, ierr)
        call HYPRE_EuclidSetLevel(precond_Sw, 1, ierr)
        call HYPRE_EuclidSetBJ(precond_Sw, 1, ierr)
        call HYPRE_EuclidSetRowScale(precond_Sw, 1, ierr)
        call HYPRE_ParCSRGMRESCreate(MPI_COMM_WORLD, solver_Sw, ierr)
        call HYPRE_ParCSRGMRESSetTol(solver_Sw, 1.D-7, ierr)
        call HYPRE_ParCSRGMRESSetPrintLevel(solver_Sw, 3, ierr)
        call HYPRE_ParCSRGMRESSetLogging(solver_Sw, 3, ierr)
        ! 1 means the DS preconditioner, 4 means the ParaSails preconditioner
        call HYPRE_ParCSRGMRESSetPrecond(solver_Sw, 0, precond_Sw, ierr)

        call HYPRE_ParCSRGMRESCreate(MPI_COMM_WORLD, solver_Tem, ierr)
        call HYPRE_ParCSRGMRESSetPrintLevel(solver_Tem, 2, ierr)
        call HYPRE_ParCSRGMRESSetLogging(solver_Tem, 1, ierr)

        call HYPRE_ParaSailsCreate(MPI_COMM_WORLD, precond_Tem,ierr)
        call HYPRE_ParaSailsSetParams(precond_Tem, -9D-1, 2, ierr)
        call HYPRE_ParaSailsSetFilter(precond_Tem, -9D-1, ierr)
        ! Because the matrix A is nonsymmetric and indefinite, you must choose the parameter as 0.
        call HYPRE_ParaSailsSetSym(precond_Tem, 0)
        call HYPRE_ParaSailsSetLogging(precond_Tem, 1, ierr)

        ! 1 means the DS preconditioner, 4 means the ParaSails preconditioner
        call HYPRE_ParCSRGMRESSetPrecond(solver_Tem, 1, precond_Tem, ierr)

        allocate(rows_vp(local_x_size_vp))
        do n = 1, local_x_size_vp
            rows_vp(n) = ilower_vp + n - 1
        end do

        allocate(rows_Cf(local_x_size_Cf))
        do n = 1, local_x_size_Cf
            rows_Cf(n) = ilower_Cf + n - 1
        end do

        allocate(rows_Sw(local_x_size_Sw))
        do n = 1, local_x_size_Sw
            rows_Sw(n) = ilower_Sw + n - 1
        end do

        allocate(rows_Tem(local_x_size_Tem))
        do n = 1, local_x_size_Tem
            rows_Tem(n) = ilower_Tem + n - 1
        end do

        ! the initial guess of x in the solver iteration
        allocate(initial_x_guess_vp(local_x_size_vp))
        initial_x_guess_vp(:) = 0.D0
        c = local_x_size_vp - localncols*localnrows
        do j = 1, localnrows
            do i = 1, localncols
                c = c + 1
                initial_x_guess_vp(c) = pInit(xlower+i-1,ylower+j-1)
            end do
        end do

        allocate(initial_x_guess_Cf(local_x_size_Cf))
        c = 0
        do j = 1, localnrows
            do i = 1, localncols
                c = c + 1
                initial_x_guess_Cf(c) = CfInit(xlower+i-1,ylower+j-1)
            end do
        end do

        allocate(initial_x_guess_Sw(local_x_size_Sw))
        c = 0
        do j = 1, localnrows
            do i = 1, localncols
                c = c + 1
                initial_x_guess_Sw(c) = SwInit(xlower+i-1,ylower+j-1)
            end do
        end do

        allocate(initial_x_guess_Tem(local_x_size_Tem))
        c = 0
        do j = 1, localnrows
            do i = 1, localncols
                c = c + 1
                initial_x_guess_Tem(c) = TemInit(xlower+i-1,ylower+j-1)
            end do
        end do

#endif

        if(mod(nt, NUMFRAME) /= 0) then
            if(myid == 0) then
                print *, 'The number of time steps must be divided by the number of frames.'
            end if
            stop
        end if

    end subroutine initialize

    ! Generate the values on the right-hand side and the coefficients of the matrix A,
    ! and the subroutine will generate the values that will not change with the time iteration.
    subroutine genStaticPara_vp()

        integer :: findexl, findexr, findexd, findexu ! field index
        integer :: eindexl, eindexr, eindexd, eindexu ! equation index
        integer :: global_ind, phase
        integer, dimension(:,:), pointer :: velx, vely, pres
        logical :: isField
        real(kind=8), dimension(:,:), pointer :: rhs_velx_bw, rhs_velx_bn, rhs_dpdx, &!
            rhs_vely_bw, rhs_vely_bn, rhs_dpdy
        real(kind=8), dimension(:,:), pointer :: resiAxx_b, resiAxp, resiAyy_b, resiAyp, resitemp
        integer :: i, j, n, c

        if(pcol /= 1) then
            findexl = 0
        else
            findexl = 1
        end if
        findexr = localncols + 1
        if(prow /= 1) then
            findexd = 0
        else
            findexd = 1
        end if
        if(prow /= pnrows) then
            findexu = localnrows + 1
        else
            findexu = localnrows
        end if
        allocate(velx(findexl:findexr,findexd:findexu))

        eindexl = 1
        if(pcol /= pncols) then
            eindexr = localncols
        else
            eindexr = localncols + 1
        end if
        eindexd = 1
        eindexu = localnrows
        allocate(rhs_velx_bw(eindexl:eindexr,eindexd:eindexu))
        allocate(rhs_velx_bn(eindexl:eindexr,eindexd:eindexu))
        allocate(resiAxx_b(eindexl:eindexr,eindexd:eindexu))
        allocate(resitemp(eindexl:eindexr,eindexd:eindexu))
        do phase = 1, 2
            velx(:,:) = 0
            if(phase == 1) then
                call Resi_xmom_velx_b(velx, rhs_velx_bw, phase)
            else
                call Resi_xmom_velx_b(velx, rhs_velx_bn, phase)
            end if
            do j = findexd, findexd+2
                do i = findexl, findexl+2
                    call genExpField(i, j, findexr, findexu, velx, isField)
                    if(isField) then
                        call Resi_xmom_velx_b(velx, resiAxx_b, phase)
                        if(phase == 1) then
                            resitemp = resiAxx_b - rhs_velx_bw
                            call constructAxx(velx, resitemp, 1111)
                        else
                            resitemp = resiAxx_b - rhs_velx_bn
                            call constructAxx(velx, resitemp, 1121)
                        end if
                    end if
                end do
            end do
        end do
        deallocate(velx)
        deallocate(resiAxx_b)
        deallocate(resitemp)

        do j = 1, localnrows
            if((pcol==1).and.(isDiriX0_p(ylower+j-1)==0)) then
                call index_convert_local_global(myid, 11, 1, j, global_ind)
                do n = 1, AxxSize
                    if(AxxwRows(n) == global_ind) then
                        AxxwValues(n) = AxxwStaticValues(n)
                    end if
                end do
                call index_convert_local_global(myid, 12, 1, j, global_ind)
                do n = 1, AxxSize
                    if(AxxnRows(n) == global_ind) then
                        AxxnValues(n) = AxxnStaticValues(n)
                    end if
                end do
            end if
            if((pcol==pncols).and.(isDiriX1_p(ylower+j-1)==0)) then
                call index_convert_local_global(myid, 11, localncols+1, j, global_ind)
                do n = 1, AxxSize
                    if(AxxwRows(n)==global_ind) then
                        AxxwValues(n) = AxxwStaticValues(n)
                    end if
                end do
                call index_convert_local_global(myid, 12, localncols+1, j, global_ind)
                do n = 1, AxxSize
                    if(AxxnRows(n)==global_ind) then
                        AxxnValues(n) = AxxnStaticValues(n)
                    end if
                end do
            end if
        end do

        if(pcol /= 1) then
            findexl = 0
        else
            findexl = 1
        end if
        findexr = localncols
        findexd = 1
        findexu = localnrows
        allocate(pres(findexl:findexr,findexd:findexu))
        allocate(rhs_dpdx(eindexl:eindexr,eindexd:eindexu))
        allocate(resiAxp(eindexl:eindexr,eindexd:eindexu))
        allocate(resitemp(eindexl:eindexr,eindexd:eindexu))
        pres(:,:) = 0
        call Resi_xmom_pres(pres, rhs_dpdx)
        do j = findexd, findexd+2
            do i = findexl, findexl+2
                call genExpField(i, j, findexr, findexu, pres, isField)
                if(isField) then
                    call Resi_xmom_pres(pres, resiAxp)
                    resitemp = resiAxp - rhs_dpdx
                    call constructAxp(pres, resitemp, 121)
                    call constructAxp(pres, resitemp, 122)
                end if
            end do
        end do
        deallocate(pres)
        deallocate(resiAxp)
        deallocate(resitemp)

        c = 0
        do j = eindexd, eindexu
            do i = eindexl, eindexr
                c = c + 1
                local_rhs_vp_static(c) = -(rhs_velx_bw(i,j) + rhs_dpdx(i,j) + rhofw*gravX)
                if(((pcol==1).and.(i==1).and.(isDiriX0_p(ylower+j-1)==0)).or.((pcol==pncols).and. &!
                    (i==localncols+1).and.(isDiriX1_p(ylower+j-1)==0))) then
                    local_rhs_vp(c) = local_rhs_vp_static(c) + rhofw*gravX
                end if
            end do
        end do
        deallocate(rhs_velx_bw)

        do j = eindexd, eindexu
            do i = eindexl, eindexr
                c = c + 1
                local_rhs_vp_static(c) = -(rhs_velx_bn(i,j) + rhs_dpdx(i,j) + rhofn*gravX)
                if(((pcol==1).and.(i==1).and.(isDiriX0_p(ylower+j-1)==0)).or.((pcol==pncols).and. &!
                    (i==localncols+1).and.(isDiriX1_p(ylower+j-1)==0))) then
                    local_rhs_vp(c) = local_rhs_vp_static(c) + rhofn*gravX
                end if
            end do
        end do
        deallocate(rhs_velx_bn)

        deallocate(rhs_dpdx)

        if(pcol /= 1) then
            findexl = 0
        else
            findexl = 1
        end if
        if(pcol /= pncols) then
            findexr = localncols + 1
        else
            findexr = localncols
        end if
        if(prow /= 1) then
            findexd = 0
        else
            findexd = 1
        end if
        findexu = localnrows + 1
        allocate(vely(findexl:findexr,findexd:findexu))
        eindexl = 1
        eindexr = localncols
        eindexd = 1
        if(prow /= pnrows) then
            eindexu = localnrows
        else
            eindexu = localnrows + 1
        end if
        allocate(rhs_vely_bw(eindexl:eindexr,eindexd:eindexu))
        allocate(rhs_vely_bn(eindexl:eindexr,eindexd:eindexu))
        allocate(resiAyy_b(eindexl:eindexr,eindexd:eindexu))
        allocate(resitemp(eindexl:eindexr,eindexd:eindexu))
        do phase = 1, 2
            vely(:,:) = 0
            if(phase == 1) then
                call Resi_ymom_vely_b(vely, rhs_vely_bw, phase)
            else
                call Resi_ymom_vely_b(vely, rhs_vely_bn, phase)
            end if
            do j = findexd, findexd+2
                do i = findexl, findexl+2
                    call genExpField(i, j, findexr, findexu, vely, isField)
                    if(isField) then
                        call Resi_ymom_vely_b(vely, resiAyy_b, phase)
                        if(phase == 1) then
                            resitemp = resiAyy_b - rhs_vely_bw
                            call constructAyy(vely, resitemp, 2111)
                        else
                            resitemp = resiAyy_b - rhs_vely_bn
                            call constructAyy(vely, resitemp, 2121)
                        end if
                    end if
                end do
            end do
        end do
        deallocate(vely)
        deallocate(resiAyy_b)
        deallocate(resitemp)

        do i = 1, localncols
            if((prow==1).and.(isDiriY0_p(xlower+i-1)==0)) then
                call index_convert_local_global(myid, 21, i, 1, global_ind)
                do n = 1, AyySize
                    if(AyywRows(n)==global_ind) then
                        AyywValues(n) = AyywStaticValues(n)
                    end if
                end do
                call index_convert_local_global(myid, 22, i, 1, global_ind)
                do n = 1, AyySize
                    if(AyynRows(n)==global_ind) then
                        AyynValues(n) = AyynStaticValues(n)
                    end if
                end do
            end if
            if((prow==pnrows).and.(isDiriY1_p(xlower+i-1)==0)) then
                call index_convert_local_global(myid, 21, i, localnrows+1, global_ind)
                do n = 1, AyySize
                    if(AyywRows(n)==global_ind) then
                        AyywValues(n) = AyywStaticValues(n)
                    end if
                end do
                call index_convert_local_global(myid, 22, i, localnrows+1, global_ind)
                do n = 1, AyySize
                    if(AyynRows(n)==global_ind) then
                        AyynValues(n) = AyynStaticValues(n)
                    end if
                end do
            end if
        end do

        findexl = 1
        findexr = localncols
        if(prow /= 1) then
            findexd = 0
        else
            findexd = 1
        end if
        findexu = localnrows
        allocate(pres(findexl:findexr,findexd:findexu))
        allocate(rhs_dpdy(eindexl:eindexr,eindexd:eindexu))
        allocate(resiAyp(eindexl:eindexr,eindexd:eindexu))
        allocate(resitemp(eindexl:eindexr,eindexd:eindexu))
        pres(:,:) = 0
        call Resi_ymom_pres(pres, rhs_dpdy)
        do j = findexd, findexd+2
            do i = findexl, findexl+2
                call genExpField(i, j, findexr, findexu, pres, isField)
                if(isField) then
                    call Resi_ymom_pres(pres, resiAyp)
                    resitemp = resiAyp - rhs_dpdy
                    call constructAyp(pres, resitemp, 221)
                    call constructAyp(pres, resitemp, 222)
                end if
            end do
        end do
        deallocate(pres)
        deallocate(resiAyp)
        deallocate(resitemp)

        do j = eindexd, eindexu
            do i = eindexl, eindexr
                c = c + 1
                local_rhs_vp_static(c) = -(rhs_vely_bw(i,j) + rhs_dpdy(i,j) + rhofw*gravY)
                if(((prow==1).and.(j==1).and.(isDiriY0_p(xlower+i-1)==0)).or.((prow==pnrows).and.(j==localnrows+1) &!
                    .and.(isDiriY1_p(xlower+i-1)==0))) then
                    local_rhs_vp(c) = local_rhs_vp_static(c) + rhofw*gravY
                end if
            end do
        end do
        deallocate(rhs_vely_bw)

        do j = eindexd, eindexu
            do i = eindexl, eindexr
                c = c + 1
                local_rhs_vp_static(c) = -(rhs_vely_bn(i,j) + rhs_dpdy(i,j) + rhofn*gravY)
                if(((prow==1).and.(j==1).and.(isDiriY0_p(xlower+i-1)==0)).or.((prow==pnrows).and.(j==localnrows+1) &!
                    .and.(isDiriY1_p(xlower+i-1)==0))) then
                    local_rhs_vp(c) = local_rhs_vp_static(c) + rhofn*gravY
                end if
            end do
        end do
        deallocate(rhs_vely_bn)

        deallocate(rhs_dpdy)

    end subroutine genStaticPara_vp

    subroutine computedm()

        integer :: status(MPI_STATUS_SIZE)
        integer :: requestl, requestr, requestd, requestu
        real(kind=8), dimension(:), allocatable :: sent, recv
        integer :: sentSize, recvSize
        integer :: i, j, c, ierr

        do j = 1, localnrows
            do i = 1, localncols

                dm(i,j) = dmRef * exp(Eg/Rg*(1.D0/TemRef_dm-1.D0/Tem(i,j)))

            end do
        end do

        ! send
        ! the 4 edges
        if(pcol /= 1) then
            sentSize = localnrows
            allocate(sent(sentSize))
            c = 0
            do j = 1, localnrows
                c = c + 1
                sent(c) = dm(1,j)
            end do
            call MPI_IBSEND(sent, sentSize, MPI_DOUBLE_PRECISION, myid-1, myid, MPI_COMM_WORLD, requestl, ierr)
            deallocate(sent)
        end if

        if(pcol /= pncols) then
            sentSize = localnrows
            allocate(sent(sentSize))
            c = 0
            do j = 1, localnrows
                c = c + 1
                sent(c) = dm(localncols,j)
            end do
            call MPI_IBSEND(sent, sentSize, MPI_DOUBLE_PRECISION, myid+1, myid, MPI_COMM_WORLD, requestr, ierr)
            deallocate(sent)
        end if

        if(prow /= 1) then
            sentSize = localncols
            allocate(sent(sentSize))
            c = 0
            do i = 1, localncols
                c = c + 1
                sent(c) = dm(i,1)
            end do
            call MPI_IBSEND(sent, sentSize, MPI_DOUBLE_PRECISION, myid-pncols, myid, MPI_COMM_WORLD, requestd, ierr)
            deallocate(sent)
        end if

        if(prow /= pnrows) then
            sentSize = localncols
            allocate(sent(sentSize))
            c = 0
            do i = 1, localncols
                c = c + 1
                sent(c) = dm(i,localnrows)
            end do
            call MPI_IBSEND(sent, sentSize, MPI_DOUBLE_PRECISION, myid+pncols, myid, MPI_COMM_WORLD, requestu, ierr)
            deallocate(sent)
        end if

        ! receive
        ! the 4 edges
        if(pcol /= pncols) then
            recvSize = localnrows
            allocate(recv(recvSize))
            call MPI_RECV(recv, recvSize, MPI_DOUBLE_PRECISION, myid+1, myid+1, MPI_COMM_WORLD, status, ierr)
            c = 0
            do j = 1, localnrows
                c = c + 1
                dm(localncols+1,j) = recv(c)
            end do
            deallocate(recv)
        end if

        if(pcol /= 1) then
            recvSize = localnrows
            allocate(recv(recvSize))
            call MPI_RECV(recv, recvSize, MPI_DOUBLE_PRECISION, myid-1, myid-1, MPI_COMM_WORLD, status, ierr)
            c = 0
            do j = 1, localnrows
                c = c + 1
                dm(0,j) = recv(c)
            end do
            deallocate(recv)
        end if

        if(prow /= pnrows) then
            recvSize = localncols
            allocate(recv(recvSize))
            call MPI_RECV(recv, recvSize, MPI_DOUBLE_PRECISION, myid+pncols, myid+pncols, MPI_COMM_WORLD, status, ierr)
            c = 0
            do i = 1, localncols
                c = c + 1
                dm(i,localnrows+1) = recv(c)
            end do
            deallocate(recv)
        end if

        if(prow /= 1) then
            recvSize = localncols
            allocate(recv(recvSize))
            call MPI_RECV(recv, recvSize, MPI_DOUBLE_PRECISION, myid-pncols, myid-pncols, MPI_COMM_WORLD, status, ierr)
            c = 0
            do i = 1, localncols
                c = c + 1
                dm(i,0) = recv(c)
            end do
            deallocate(recv)
        end if

        if(pcol /= 1) then
            call MPI_WAIT(requestl, status, ierr)
        end if
        if(pcol /= pncols) then
            call MPI_WAIT(requestr, status, ierr)
        end if
        if(prow /= 1) then
            call MPI_WAIT(requestd, status, ierr)
        end if
        if(prow /= pnrows) then
            call MPI_WAIT(requestu, status, ierr)
        end if

    end subroutine computedm

    subroutine computeks()

        integer :: i, j

        do j = 1, localnrows
            do i = 1, localncols

                ks(i,j) = ksRef * exp(Eg/Rg*(1.D0/TemRef_ks-1.D0/Tem(i,j)))

            end do
        end do

    end subroutine computeks

    subroutine computekc()

        real(kind=8) :: Sc, vwmodulus, radius, Rep
        integer :: i, j

        do j = 1, localnrows
            do i = 1, localncols

                vwmodulus = dsqrt(((vxw(i+1,j)-vxw(i,j))**2.D0+(vyw(i,j+1)-vyw(i,j))**2.D0))

                radius = radiusInit*(1.D0-poroInit(xlower+i-1,ylower+j-1))/ &!
                    (poroInit(xlower+i-1,ylower+j-1)*(1.D0-poro(i,j)))*poro(i,j)

                Rep = 2.D0*vwmodulus*radius/(viscw/rhofw)

                Sc = viscw/rhofw*dm(i,j)

                kc(i,j) = (ShInfinity+7.D-1*Rep**(1.D0/2.D0)*Sc**(1.D0/3.D0))*dm(i,j)/2.D0/radius

            end do
        end do

    end subroutine computekc

    subroutine computePoro()

        integer :: status(MPI_STATUS_SIZE)
        integer :: requestl, requestr, requestd, requestu, requestlu, requestrd, requestru
        real(kind=8), dimension(:), allocatable :: sent, recv
        integer :: sentSize, recvSize
        integer :: i, j, c, ierr

        do j = 1, localnrows
            do i = 1, localncols

                poro_old(i,j) = poro(i,j)
                poro(i,j) = poro_old(i,j) + Sw(i,j)*av(i,j)*al*Cf(i,j)*kc(i,j)*ks(i,j)*(ts(t)-ts(t-1))/(rhos*(kc(i,j)+ks(i,j)))

            end do
        end do

        if(pcol /= 1) then
            sentSize = 2*localnrows
            allocate(sent(sentSize))
            c = 0
            do j = 1, localnrows
                c = c + 1
                sent(c) = poro(1,j)
                c = c + 1
                sent(c) = poro_old(1,j)
            end do
            call MPI_IBSEND(sent, sentSize, MPI_DOUBLE_PRECISION, myid-1, myid, MPI_COMM_WORLD, requestl, ierr)
            deallocate(sent)
        end if

        if(pcol /= pncols) then
            sentSize = 4*localnrows
            allocate(sent(sentSize))
            c = 0
            do j = 1, localnrows
                do i = localncols-1, localncols
                    c = c + 1
                    sent(c) = poro(i,j)
                    c = c + 1
                    sent(c) = poro_old(i,j)
                end do
            end do
            call MPI_IBSEND(sent, sentSize, MPI_DOUBLE_PRECISION, myid+1, myid, MPI_COMM_WORLD, requestr, ierr)
            deallocate(sent)
        end if

        if(prow /= 1) then
            sentSize = 2*localncols
            allocate(sent(sentSize))
            c = 0
            do i = 1, localncols
                c = c + 1
                sent(c) = poro(i,1)
                c = c + 1
                sent(c) = poro_old(i,1)
            end do
            call MPI_IBSEND(sent, sentSize, MPI_DOUBLE_PRECISION, myid-pncols, myid, MPI_COMM_WORLD, requestd, ierr)
            deallocate(sent)
        end if

        if(prow /= pnrows) then
            sentSize = localncols*4
            allocate(sent(sentSize))
            c = 0
            do j = localnrows-1, localnrows
                do i = 1, localncols
                    c = c + 1
                    sent(c) = poro(i,j)
                    c = c + 1
                    sent(c) = poro_old(i,j)
                end do
            end do
            call MPI_IBSEND(sent, sentSize, MPI_DOUBLE_PRECISION, myid+pncols, myid, MPI_COMM_WORLD, requestu, ierr)
            deallocate(sent)
        end if

        if((pcol/=1).and.(prow/=pnrows)) then
            allocate(sent(2))
            sent(1) = poro(1,localnrows)
            sent(2) = poro_old(1,localnrows)
            call MPI_IBSEND(sent, 2, MPI_DOUBLE_PRECISION, myid-1+pncols, myid, MPI_COMM_WORLD, requestlu, ierr)
            deallocate(sent)
        end if

        if((pcol/=pncols).and.(prow/=1)) then
            allocate(sent(2))
            sent(1) = poro(localncols,1)
            sent(2) = poro_old(localncols,1)
            call MPI_IBSEND(sent, 2, MPI_DOUBLE_PRECISION, myid+1-pncols, myid, MPI_COMM_WORLD, requestrd, ierr)
            deallocate(sent)
        end if

        if((pcol/=pncols).and.(prow/=pnrows)) then
            allocate(sent(2))
            sent(1) = poro(localncols,localnrows)
            sent(2) = poro_old(localncols,localnrows)
            call MPI_IBSEND(sent, 2, MPI_DOUBLE_PRECISION, myid+1+pncols, myid, MPI_COMM_WORLD, requestru, ierr)
            deallocate(sent)
        end if

        if(pcol /= pncols) then
            recvSize = 2*localnrows
            allocate(recv(recvSize))
            call MPI_RECV(recv, recvSize, MPI_DOUBLE_PRECISION, myid+1, myid+1, MPI_COMM_WORLD, status, ierr)
            c = 0
            do j = 1, localnrows
                c = c + 1
                poro(localncols+1,j) = recv(c)
                c = c + 1
                poro_old(localncols+1,j) = recv(c)
            end do
            deallocate(recv)
        end if

        if(pcol /= 1) then
            recvSize = 4*localnrows
            allocate(recv(recvSize))
            call MPI_RECV(recv, recvSize, MPI_DOUBLE_PRECISION, myid-1, myid-1, MPI_COMM_WORLD, status, ierr)
            c = 0
            do j = 1, localnrows
                do i = -1, 0
                    c = c + 1
                    poro(i,j) = recv(c)
                    c = c + 1
                    poro_old(i,j) = recv(c)
                end do
            end do
            deallocate(recv)
        end if

        if(prow /= pnrows) then
            recvSize = 2*localncols
            allocate(recv(recvSize))
            call MPI_RECV(recv, recvSize, MPI_DOUBLE_PRECISION, myid+pncols, myid+pncols, MPI_COMM_WORLD, status, ierr)
            c = 0
            do i = 1, localncols
                c = c + 1
                poro(i,localnrows+1) = recv(c)
                c = c + 1
                poro_old(i,localnrows+1) = recv(c)
            end do
            deallocate(recv)
        end if

        if(prow /= 1) then
            recvSize = localncols*4
            allocate(recv(recvSize))
            call MPI_RECV(recv, recvSize, MPI_DOUBLE_PRECISION, myid-pncols, myid-pncols, MPI_COMM_WORLD, status, ierr)
            c = 0
            do j = -1, 0
                do i = 1, localncols
                    c = c + 1
                    poro(i,j) = recv(c)
                    c = c + 1
                    poro_old(i,j) = recv(c)
                end do
            end do
            deallocate(recv)
        end if

        if((pcol/=pncols).and.(prow/=1)) then
            allocate(recv(2))
            call MPI_RECV(recv, 2, MPI_DOUBLE_PRECISION, myid+1-pncols, myid+1-pncols, MPI_COMM_WORLD, status, ierr)
            poro(localncols+1,0) = recv(1)
            poro_old(localncols+1,0) = recv(2)
            deallocate(recv)
        end if

        if((pcol/=1).and.(prow/=pnrows)) then
            allocate(recv(2))
            call MPI_RECV(recv, 2, MPI_DOUBLE_PRECISION, myid-1+pncols, myid-1+pncols, MPI_COMM_WORLD, status, ierr)
            poro(0,localnrows+1) = recv(1)
            poro_old(0,localnrows+1) = recv(2)
            deallocate(recv)
        end if

        if((pcol/=1).and.(prow/=1)) then
            allocate(recv(2))
            call MPI_RECV(recv, 2, MPI_DOUBLE_PRECISION, myid-1-pncols, myid-1-pncols, MPI_COMM_WORLD, status, ierr)
            poro(0,0) = recv(1)
            poro_old(0,0) = recv(2)
            deallocate(recv)
        end if

        if(pcol /= 1) then
            call MPI_WAIT(requestl, status, ierr)
        end if
        if(pcol /= pncols) then
            call MPI_WAIT(requestr, status, ierr)
        end if
        if(prow /= 1) then
            call MPI_WAIT(requestd, status, ierr)
        end if
        if(prow /= pnrows) then
            call MPI_WAIT(requestu, status, ierr)
        end if
        if((pcol/=1).and.(prow/=pnrows)) then
            call MPI_WAIT(requestlu, status, ierr)
        end if
        if((pcol/=pncols).and.(prow/=1)) then
            call MPI_WAIT(requestrd, status, ierr)
        end if
        if((pcol/=pncols).and.(prow/=pnrows)) then
            call MPI_WAIT(requestru, status, ierr)
        end if

    end subroutine computePoro

    subroutine computePoroEdge(var)

        integer, intent(in) :: var

        real(kind=8), dimension(:,:), pointer :: EdgeX, EdgeY
        integer :: indexl, indexr, indexd, indexu
        integer :: i, j

        if(var == 1) then
            EdgeX => poroEdgeX_old
            EdgeY => poroEdgeY_old
        elseif(var == 2) then
            EdgeX => poroEdgeX
            EdgeY => poroEdgeY
        else
            print *, 'The var value in computePoroEdge is wrong!'
            stop
        end if

        if(pcol /= 1) then
            indexl = 0
        else
            indexl = 1
        end if
        indexr = localncols + 1
        if(prow /= 1) then
            indexd = 0
        else
            indexd = 1
        end if
        if(prow /= pnrows) then
            indexu = localnrows + 1
        else
            indexu = localnrows
        end if

        if(pcol == 1) then
            EdgeX(1,indexd:indexu) = poro(1,indexd:indexu)
            indexl = indexl + 1
        end if

        if(pcol == pncols) then
            EdgeX(localncols+1,indexd:indexu) = poro(localncols,indexd:indexu)
            indexr = indexr - 1
        end if

        do j = indexd, indexu
            do i = indexl, indexr
                EdgeX(i,j) = (hx(i-1)+hx(i)) / (hx(i-1)/poro(i-1,j)+hx(i)/poro(i,j))
            end do
        end do

        if(pcol /= 1) then
            indexl = 0
        else
            indexl = 1
        end if
        if(pcol /= pncols) then
            indexr = localncols + 1
        else
            indexr = localncols
        end if
        if(prow /= 1) then
            indexd = 0
        else
            indexd = 1
        end if
        indexu = localnrows + 1

        if(prow == 1) then
            EdgeY(indexl:indexr,1) = poro(indexl:indexr,1)
            indexd = indexd + 1
        end if

        if(prow == pnrows) then
            EdgeY(indexl:indexr,localnrows+1) = poro(indexl:indexr,localnrows)
            indexu = indexu - 1
        end if

        do j = indexd, indexu
            do i = indexl, indexr
                EdgeY(i,j) = (hy(j-1)+hy(j)) / (hy(j-1)/poro(i,j-1)+hy(j)/poro(i,j))
            end do
        end do

        if((t == 2).and.(var == 1)) then
            poroEdgeXInit = EdgeX
            poroEdgeYInit = EdgeY
        end if

    end subroutine computePoroEdge

    subroutine genDynPara_Sw()

        integer :: findexl, findexr, findexd, findexu
        integer :: eindexl, eindexr, eindexd, eindexu
        integer, dimension(:,:), pointer :: satw
        logical :: isField
        real(kind=8), dimension(:,:), pointer :: rhs_Sw
        real(kind=8), dimension(:,:), pointer :: resiASw, resitemp
        integer :: i, j, c

        if(pcol /= 1) then
            findexl = 0
        else
            findexl = 1
        end if
        if(pcol /= pncols) then
            findexr = localncols + 1
        else
            findexr = localncols
        end if
        if(prow /= 1) then
            findexd = 0
        else
            findexd = 1
        end if
        if(prow /= pnrows) then
            findexu = localnrows + 1
        else
            findexu = localnrows
        end if
        allocate(satw(findexl:findexr,findexd:findexu))
        satw(:,:) = 0

        eindexl = 1
        eindexr = localncols
        eindexd = 1
        eindexu = localnrows
        allocate(rhs_Sw(eindexl:eindexr,eindexd:eindexu))
        rhs_Sw(:,:) = 0.D0

        call Resi_con_Sw(satw, rhs_Sw)
        allocate(resiASw(eindexl:eindexr,eindexd:eindexu))
        allocate(resitemp(eindexl:eindexr,eindexd:eindexu))
        do j = findexd, findexd+2
            do i = findexl, findexl+2
                call genExpField(i, j, findexr, findexu, satw, isField)
                if(isField) then
                    call Resi_con_Sw(satw, resiASw)
                    resitemp = resiASw - rhs_Sw
                    call constructASw(satw, resitemp)
                end if
            end do
        end do

        c = 0
        do j = 1, localnrows
            do i = 1, localncols
                c = c + 1
                local_rhs_Sw(c) = -rhs_Sw(i,j)
            end do
        end do

        deallocate(satw)
        deallocate(rhs_Sw)
        deallocate(resiASw)
        deallocate(resitemp)

    end subroutine genDynPara_Sw

    subroutine computeSw()

        integer, dimension(:), allocatable :: cols
        real(kind=8), dimension(:), allocatable :: values
        real(kind=8), dimension(:), pointer :: local_x
        integer :: ASwBeInd
        integer :: status(MPI_STATUS_SIZE)
        integer :: request
        integer, dimension(:), allocatable :: requestarray
        real(kind=8), dimension(:), allocatable :: slave_data
        real(kind=8) :: solvertimestart, solvertimefinish

        integer :: requestl, requestr, requestd, requestu, requestlu, requestrd, requestru
        real(kind=8), dimension(:), allocatable :: sent, recv
        integer :: sentSize, recvSize

        integer :: i, j, l, n, c, num_iter, ierr

        allocate(cols(5))
        allocate(values(5)) ! Each row of A has at most 5 nonzero values.
        allocate(local_x(local_x_size_Sw))

#ifdef LAPACK

        A_lapack_Sw(:,:) = 0.D0
        b_lapack_Sw(:) = local_rhs_Sw(:)
        IPIV_Sw(:) = 0.D0

#elif defined(UMFPACK)

        Ap_Sw(1) = 0
        Ai_Sw(:) = 0
        Ax_Sw(:) = 0.D0

#elif defined(MUMPS)

        mumps_NNZ_loc_Sw = 0
        mumps_IRN_loc_Sw(:) = 0
        mumps_JCN_loc_Sw(:) = 0
        mumps_A_loc_Sw(:) = 0.D0

        if(mumps_par_Sw%MYID /= 0) then
            call MPI_IBSEND(local_rhs_Sw, local_x_size_Sw, MPI_DOUBLE_PRECISION, 0, myid, &!
                MPI_COMM_WORLD, request, ierr)
        end if

        if(mumps_par_Sw%MYID == 0) then
            mumps_par_Sw%RHS(1:local_x_size_Sw) = local_rhs_Sw(1:local_x_size_Sw)
            c = local_x_size_Sw + 1
            do i = 1, nProcs-1
                allocate(slave_data(local_x_size_Sw))
                call MPI_RECV(slave_data, local_x_size_Sw, MPI_DOUBLE_PRECISION, i, i, MPI_COMM_WORLD, status, ierr)
                mumps_par_Sw%RHS(c:c+local_x_size_Sw-1) = slave_data(:)
                c = c + local_x_size_Sw
                deallocate(slave_data)
            end do
        end if

        if(mumps_par_Sw%MYID /= 0) then
            call MPI_WAIT(request, status, ierr)
        end if

#elif defined(HYPRE)

        call HYPRE_IJVectorSetValues(b_Sw, local_x_size_Sw, rows_Sw, local_rhs_Sw, ierr)
        call HYPRE_IJVectorSetValues(x_Sw, local_x_size_Sw, rows_Sw, initial_x_guess_Sw, ierr)

        call HYPRE_IJVectorAssemble(b_Sw, ierr)
        call HYPRE_IJVectorAssemble(x_Sw, ierr)

        call HYPRE_IJVectorGetObject(b_Sw, par_b_Sw, ierr)
        call HYPRE_IJVectorGetObject(x_Sw, par_x_Sw, ierr)

#endif

        ASwBeInd = 1
        do n = ilower_Sw, iupper_Sw

            cols(:) = 0
            values(:) = 0.D0

            c = 0
            do l = ASwBeInd, ASwSize
                if(ASwRows(l) == n) then
                    c = c + 1
                    cols(c) = ASwCols(l)
                    values(c) = ASwValues(l)
                else
                    ASwBeInd = l
                    exit
                end if
            end do

#ifdef LAPACK
            do l = 1, c
                A_lapack_Sw(n,cols(l)) = values(l)
            end do
#elif defined(UMFPACK)
            Ap_Sw(n+1) = Ap_Sw(n) + c
            Ai_Sw(Ap_Sw(n)+1:Ap_Sw(n+1)) = cols(1:c) - 1
            Ax_Sw(Ap_Sw(n)+1:Ap_Sw(n+1)) = values(1:c)
#elif defined(MUMPS)
            mumps_IRN_loc_Sw(mumps_NNZ_loc_Sw+1:mumps_NNZ_loc_Sw+c) = n
            mumps_JCN_loc_Sw(mumps_NNZ_loc_Sw+1:mumps_NNZ_loc_Sw+c) = cols(1:c)
            mumps_A_loc_Sw(mumps_NNZ_loc_Sw+1:mumps_NNZ_loc_Sw+c) = values(1:c)
            mumps_NNZ_loc_Sw = mumps_NNZ_loc_Sw + c
#elif defined(HYPRE)
            call HYPRE_IJMatrixSetValues(A_Sw, 1, c, n, cols, values, ierr)
#endif

        end do

#ifdef LAPACK

        ! you have to make sure that the number of processors is set to 1 when using such method.
        solvertimestart = MPI_Wtime()
        call dgesv(local_x_size_Sw, 1, A_lapack_Sw, local_x_size_Sw, IPIV_Sw, b_lapack_Sw, local_x_size_Sw, LAPACKINFO)
        solvertimefinish = MPI_Wtime()
        solvertime = solvertime + solvertimefinish - solvertimestart
        local_x(:) = b_lapack_Sw(:)
        if(LAPACKINFO /= 0) then
            print *, 'LAPACK solver error when solving Sw. INFO = ', LAPACKINFO
            stop
        end if

#elif defined(UMFPACK)

        solvertimestart = MPI_Wtime()
        call umf4def(control)
        call umf4sym(local_x_size_Sw, local_x_size_Sw, Ap_Sw, Ai_Sw, Ax_Sw, symbolic, control, umfinfo)
        call umf4num(Ap_Sw, Ai_Sw, Ax_Sw, symbolic, numeric, control, umfinfo)
        call umf4fsym(symbolic)
        call umf4sol(1, local_x, local_rhs_Sw, numeric, control, umfinfo) ! 1 means A'x=b
        if(umfinfo(1) < 0) then
            print *, 'UMFPACK solver error when solving Sw. Info: ', umfinfo(1)
            stop
        end if
        call umf4fnum(numeric)
        solvertimefinish = MPI_Wtime()
        solvertime = solvertime + solvertimefinish - solvertimestart

#elif defined(MUMPS)

        mumps_par_Sw%NNZ_loc = mumps_NNZ_loc_Sw
        allocate(mumps_par_Sw%IRN_loc(mumps_par_Sw%NNZ_loc))
        allocate(mumps_par_Sw%JCN_loc(mumps_par_Sw%NNZ_loc))
        allocate(mumps_par_Sw%A_loc(mumps_par_Sw%NNZ_loc))
        mumps_par_Sw%IRN_loc(1:mumps_par_Sw%NNZ_loc) = mumps_IRN_loc_Sw(1:mumps_NNZ_loc_Sw)
        mumps_par_Sw%JCN_loc(1:mumps_par_Sw%NNZ_loc) = mumps_JCN_loc_Sw(1:mumps_NNZ_loc_Sw)
        mumps_par_Sw%A_loc(1:mumps_par_Sw%NNZ_loc) = mumps_A_loc_Sw(1:mumps_NNZ_loc_Sw)

        mumps_par_Sw%JOB = 6
        solvertimestart = MPI_Wtime()
        call DMUMPS(mumps_par_Sw)
        solvertimefinish = MPI_Wtime()
        solvertime = solvertime + solvertimefinish - solvertimestart
        if(mumps_par_Sw%INFOG(1) < 0) then
            stop
        end if

        if(mumps_par_Sw%MYID == 0) then
            local_x(:) = mumps_par_Sw%RHS(1:local_x_size_Sw)
            allocate(requestarray(nProcs-1))
            c = local_x_size_Sw + 1
            do i = 1, nProcs-1
                allocate(slave_data(local_x_size_Sw))
                slave_data(:) = mumps_par_Sw%RHS(c:c+local_x_size_Sw-1)
                call MPI_IBSEND(slave_data, local_x_size_Sw, MPI_DOUBLE_PRECISION, i, myid, &!
                    MPI_COMM_WORLD, requestarray(i), ierr)
                c = c + local_x_size_Sw
                deallocate(slave_data)
            end do
        end if

        if(mumps_par_Sw%MYID /= 0) then
            allocate(slave_data(local_x_size_Sw))
            call MPI_RECV(slave_data, local_x_size_Sw, MPI_DOUBLE_PRECISION, 0, 0, MPI_COMM_WORLD, status, ierr)
            local_x(:) = slave_data(:)
            deallocate(slave_data)
        end if

        if(mumps_par_Sw%MYID == 0) then
            do i = 1, nProcs-1
                call MPI_WAIT(requestarray(i), status, ierr)
            end do
            deallocate(requestarray)
        end if

        deallocate(mumps_par_Sw%IRN_loc)
        deallocate(mumps_par_Sw%JCN_loc)
        deallocate(mumps_par_Sw%A_loc)
    
#elif defined(HYPRE)

        call HYPRE_IJMatrixAssemble(A_Sw, ierr)
        call HYPRE_IJMatrixGetObject(A_Sw, parcsr_A_Sw, ierr)
        solvertimestart = MPI_Wtime()
        call HYPRE_ParCSRGMRESSetup(solver_Sw, parcsr_A_Sw, par_b_Sw, par_x_Sw, ierr)
        call HYPRE_ParCSRGMRESSolve(solver_Sw, parcsr_A_Sw, par_b_Sw, par_x_Sw, ierr)
        solvertimefinish = MPI_Wtime()
        solvertime = solvertime + solvertimefinish - solvertimestart
        call HYPRE_ParCSRGMRESGetNumIteratio(solver_Sw, num_iter, ierr)
        if(ierr /= 0) then
            if(myid == 0) then
                print *, 'HYPRE solver error when solving Sw. ierr = ', ierr
                stop
            end if
        end if

        call HYPRE_IJVectorGetValues(x_Sw, local_x_size_Sw, rows_Sw, local_x, ierr)

        ! let the solution of this time step be the initial x guess in the next time step.
        ! by this way, the number of solver iteration steps can be reduced greatly.
        initial_x_guess_Sw(:) = local_x(:)

#endif

        c = 0
        do j = 1, localnrows
            do i = 1, localncols
                c = c + 1
                Sw_old(i,j) = Sw(i,j)
                Sw(i,j) = local_x(c)
            end do
        end do

        deallocate(cols)
        deallocate(values)
        deallocate(local_x)

        if(pcol /= 1) then
            sentSize = 2*localnrows
            allocate(sent(sentSize))
            c = 0
            do j = 1, localnrows
                c = c + 1
                sent(c) = Sw(1,j)
                c = c + 1
                sent(c) = Sw_old(1,j)
            end do
            call MPI_IBSEND(sent, sentSize, MPI_DOUBLE_PRECISION, myid-1, myid, MPI_COMM_WORLD, requestl, ierr)
            deallocate(sent)
        end if

        if(pcol /= pncols) then
            sentSize = 4*localnrows
            allocate(sent(sentSize))
            c = 0
            do j = 1, localnrows
                do i = localncols-1, localncols
                    c = c + 1
                    sent(c) = Sw(i,j)
                    c = c + 1
                    sent(c) = Sw_old(i,j)
                end do
            end do
            call MPI_IBSEND(sent, sentSize, MPI_DOUBLE_PRECISION, myid+1, myid, MPI_COMM_WORLD, requestr, ierr)
            deallocate(sent)
        end if

        if(prow /= 1) then
            sentSize = 2*localncols
            allocate(sent(sentSize))
            c = 0
            do i = 1, localncols
                c = c + 1
                sent(c) = Sw(i,1)
                c = c + 1
                sent(c) = Sw_old(i,1)
            end do
            call MPI_IBSEND(sent, sentSize, MPI_DOUBLE_PRECISION, myid-pncols, myid, MPI_COMM_WORLD, requestd, ierr)
            deallocate(sent)
        end if

        if(prow /= pnrows) then
            sentSize = localncols*4
            allocate(sent(sentSize))
            c = 0
            do j = localnrows-1, localnrows
                do i = 1, localncols
                    c = c + 1
                    sent(c) = Sw(i,j)
                    c = c + 1
                    sent(c) = Sw_old(i,j)
                end do
            end do
            call MPI_IBSEND(sent, sentSize, MPI_DOUBLE_PRECISION, myid+pncols, myid, MPI_COMM_WORLD, requestu, ierr)
            deallocate(sent)
        end if

        if((pcol/=1).and.(prow/=pnrows)) then
            allocate(sent(2))
            sent(1) = Sw(1,localnrows)
            sent(2) = Sw_old(1,localnrows)
            call MPI_IBSEND(sent, 2, MPI_DOUBLE_PRECISION, myid-1+pncols, myid, MPI_COMM_WORLD, requestlu, ierr)
            deallocate(sent)
        end if

        if((pcol/=pncols).and.(prow/=1)) then
            allocate(sent(2))
            sent(1) = Sw(localncols,1)
            sent(2) = Sw_old(localncols,1)
            call MPI_IBSEND(sent, 2, MPI_DOUBLE_PRECISION, myid+1-pncols, myid, MPI_COMM_WORLD, requestrd, ierr)
            deallocate(sent)
        end if

        if((pcol/=pncols).and.(prow/=pnrows)) then
            allocate(sent(2))
            sent(1) = Sw(localncols,localnrows)
            sent(2) = Sw_old(localncols,localnrows)
            call MPI_IBSEND(sent, 2, MPI_DOUBLE_PRECISION, myid+1+pncols, myid, MPI_COMM_WORLD, requestru, ierr)
            deallocate(sent)
        end if

        if(pcol /= pncols) then
            recvSize = 2*localnrows
            allocate(recv(recvSize))
            call MPI_RECV(recv, recvSize, MPI_DOUBLE_PRECISION, myid+1, myid+1, MPI_COMM_WORLD, status, ierr)
            c = 0
            do j = 1, localnrows
                c = c + 1
                Sw(localncols+1,j) = recv(c)
                c = c + 1
                Sw_old(localncols+1,j) = recv(c)
            end do
            deallocate(recv)
        end if

        if(pcol /= 1) then
            recvSize = 4*localnrows
            allocate(recv(recvSize))
            call MPI_RECV(recv, recvSize, MPI_DOUBLE_PRECISION, myid-1, myid-1, MPI_COMM_WORLD, status, ierr)
            c = 0
            do j = 1, localnrows
                do i = -1, 0
                    c = c + 1
                    Sw(i,j) = recv(c)
                    c = c + 1
                    Sw_old(i,j) = recv(c)
                end do
            end do
            deallocate(recv)
        end if

        if(prow /= pnrows) then
            recvSize = 2*localncols
            allocate(recv(recvSize))
            call MPI_RECV(recv, recvSize, MPI_DOUBLE_PRECISION, myid+pncols, myid+pncols, MPI_COMM_WORLD, status, ierr)
            c = 0
            do i = 1, localncols
                c = c + 1
                Sw(i,localnrows+1) = recv(c)
                c = c + 1
                Sw_old(i,localnrows+1) = recv(c)
            end do
            deallocate(recv)
        end if

        if(prow /= 1) then
            recvSize = localncols*4
            allocate(recv(recvSize))
            call MPI_RECV(recv, recvSize, MPI_DOUBLE_PRECISION, myid-pncols, myid-pncols, MPI_COMM_WORLD, status, ierr)
            c = 0
            do j = -1, 0
                do i = 1, localncols
                    c = c + 1
                    Sw(i,j) = recv(c)
                    c = c + 1
                    Sw_old(i,j) = recv(c)
                end do
            end do
            deallocate(recv)
        end if

        if((pcol/=pncols).and.(prow/=1)) then
            allocate(recv(2))
            call MPI_RECV(recv, 2, MPI_DOUBLE_PRECISION, myid+1-pncols, myid+1-pncols, MPI_COMM_WORLD, status, ierr)
            Sw(localncols+1,0) = recv(1)
            Sw_old(localncols+1,0) = recv(2)
            deallocate(recv)
        end if

        if((pcol/=1).and.(prow/=pnrows)) then
            allocate(recv(2))
            call MPI_RECV(recv, 2, MPI_DOUBLE_PRECISION, myid-1+pncols, myid-1+pncols, MPI_COMM_WORLD, status, ierr)
            Sw(0,localnrows+1) = recv(1)
            Sw_old(0,localnrows+1) = recv(2)
            deallocate(recv)
        end if

        if((pcol/=1).and.(prow/=1)) then
            allocate(recv(2))
            call MPI_RECV(recv, 2, MPI_DOUBLE_PRECISION, myid-1-pncols, myid-1-pncols, MPI_COMM_WORLD, status, ierr)
            Sw(0,0) = recv(1)
            Sw_old(0,0) = recv(2)
            deallocate(recv)
        end if

        if(pcol /= 1) then
            call MPI_WAIT(requestl, status, ierr)
        end if
        if(pcol /= pncols) then
            call MPI_WAIT(requestr, status, ierr)
        end if
        if(prow /= 1) then
            call MPI_WAIT(requestd, status, ierr)
        end if
        if(prow /= pnrows) then
            call MPI_WAIT(requestu, status, ierr)
        end if
        if((pcol/=1).and.(prow/=pnrows)) then
            call MPI_WAIT(requestlu, status, ierr)
        end if
        if((pcol/=pncols).and.(prow/=1)) then
            call MPI_WAIT(requestrd, status, ierr)
        end if
        if((pcol/=pncols).and.(prow/=pnrows)) then
            call MPI_WAIT(requestru, status, ierr)
        end if

    end subroutine computeSw

    subroutine computeSwEdge(var)

        integer, intent(in) :: var

        real(kind=8), dimension(:,:), pointer :: EdgeX, EdgeY
        integer :: indexl, indexr, indexd, indexu
        integer :: i, j

        if(var == 1) then
            EdgeX => SwEdgeX_old
            EdgeY => SwEdgeY_old
        elseif(var == 2) then
            EdgeX => SwEdgeX
            EdgeY => SwEdgeY
        else
            print *, 'The var value in computeSwEdge is wrong!'
            stop
        end if

        if(prow /= 1) then
            indexd = 0
        else
            indexd = 1
        end if
        if(prow /= pnrows) then
            indexu = localnrows + 1
        else
            indexu = localnrows
        end if

        do j = indexd, indexu
            if(pcol /= 1) then
                indexl = 0
            else
                indexl = 1
            end if
            indexr = localncols + 1

            if((pcol == 1).and.(isDiriX0_Sw(ylower+j-1) == 0)) then
                EdgeX(1,j) = Sw(1,j)
                indexl = indexl + 1
            elseif((pcol == 1).and.(isDiriX0_Sw(ylower+j-1) == 1).and.(vxw(1,j) > 0.D0)) then
                EdgeX(1,j) = SwBdryX0(ylower+j-1)
                indexl = indexl + 1
            end if

            if((pcol == pncols).and.(isDiriX1_Sw(ylower+j-1) == 0)) then
                EdgeX(localncols+1,j) = Sw(localncols,j)
                indexr = indexr - 1
            elseif((pcol == pncols).and.(isDiriX1_Sw(ylower+j-1) == 1).and. &!
                (vxw(localncols+1,j) < 0.D0)) then
                EdgeX(localncols+1,j) = SwBdryX1(ylower+j-1)
                indexr = indexr - 1
            end if

            do i = indexl, indexr
                if(vxw(i,j) > 0.D0) then
                    EdgeX(i,j) = Sw(i-1,j)
                else
                    EdgeX(i,j) = Sw(i,j)
                end if
            end do
        end do

        if(pcol /= 1) then
            indexl = 0
        else
            indexl = 1
        end if
        if(pcol /= pncols) then
            indexr = localncols + 1
        else
            indexr = localncols
        end if

        do i = indexl, indexr
            if(prow /= 1) then
                indexd = 0
            else
                indexd = 1
            end if
            indexu = localnrows + 1

            if((prow == 1).and.(isDiriY0_Sw(xlower+i-1) == 0)) then
                EdgeY(i,1) = Sw(i,1)
                indexd = indexd + 1
            elseif((prow == 1).and.(isDiriY0_Sw(xlower+i-1) == 1).and.(vyw(i,1) > 0.D0)) then
                EdgeY(i,1) = SwBdryY0(xlower+i-1)
                indexd = indexd + 1
            end if

            if((prow == pnrows).and.(isDiriY1_Sw(xlower+i-1) == 0)) then
                EdgeY(i,localnrows+1) = Sw(i,localnrows)
                indexu = indexu - 1
            elseif((prow == pnrows).and.(isDiriY1_Sw(xlower+i-1) == 1).and. &!
                (vyw(i,localnrows+1) < 0.D0)) then
                EdgeY(i,localnrows+1) = SwBdryY1(xlower+i-1)
                indexu = indexu - 1
            end if

            do j = indexd, indexu
                if(vyw(i,j) > 0.D0) then
                    EdgeY(i,j) = Sw(i,j-1)
                else
                    EdgeY(i,j) = Sw(i,j)
                end if
            end do
        end do

    end subroutine computeSwEdge

    subroutine computeK()

        integer :: status(MPI_STATUS_SIZE)
        integer :: requestr, requestu
        real(kind=8), dimension(:), allocatable :: sent, recv
        integer :: sentSize, recvSize
        integer :: i, j, c, ierr

        do j = 1, localnrows
            do i = 1, localncols
                Kxx(i,j) = poro(i,j)/poroInit(xlower+i-1,ylower+j-1)*(poro(i,j)*(1.D0-poroInit(xlower+i-1,ylower+j-1)) &!
                    /poroInit(xlower+i-1,ylower+j-1)/(1.D0-poro(i,j)))**2 * KxxInit(xlower+i-1,ylower+j-1)
                Kyy(i,j) = poro(i,j)/poroInit(xlower+i-1,ylower+j-1)*(poro(i,j)*(1.D0-poroInit(xlower+i-1,ylower+j-1)) &!
                    /poroInit(xlower+i-1,ylower+j-1)/(1.D0-poro(i,j)))**2 * KyyInit(xlower+i-1,ylower+j-1)
            end do
        end do

        if(pcol /= pncols) then
            sentSize = 2*localnrows
            allocate(sent(sentSize))
            c = 0
            do j = 1, localnrows
                c = c + 1
                sent(c) = Kxx(localncols,j)
            end do
            do j = 1, localnrows
                c = c + 1
                sent(c) = Kyy(localncols,j)
            end do
            call MPI_IBSEND(sent, sentSize, MPI_DOUBLE_PRECISION, myid+1, myid, MPI_COMM_WORLD, requestr, ierr)
            deallocate(sent)
        end if

        if(prow /= pnrows) then
            sentSize = 2*localncols
            allocate(sent(sentSize))
            c = 0
            do i = 1, localncols
                c = c + 1
                sent(c) = Kxx(i,localnrows)
            end do
            do i = 1, localncols
                c = c + 1
                sent(c) = Kyy(i,localnrows)
            end do
            call MPI_IBSEND(sent, sentSize, MPI_DOUBLE_PRECISION, myid+pncols, myid, MPI_COMM_WORLD, requestu, ierr)
            deallocate(sent)
        end if

        if(pcol /= 1) then
            recvSize = 2*localnrows
            allocate(recv(recvSize))
            call MPI_RECV(recv, recvSize, MPI_DOUBLE_PRECISION, myid-1, myid-1, MPI_COMM_WORLD, status, ierr)
            c = 0
            do j = 1, localnrows
                c = c + 1
                Kxx(0,j) = recv(c)
            end do
            do j = 1, localnrows
                c = c + 1
                Kyy(0,j) = recv(c)
            end do
            deallocate(recv)
        end if

        if(prow /= 1) then
            recvSize = 2*localncols
            allocate(recv(recvSize))
            call MPI_RECV(recv, recvSize, MPI_DOUBLE_PRECISION, myid-pncols, myid-pncols, MPI_COMM_WORLD, status, ierr)
            c = 0
            do i = 1, localncols
                c = c + 1
                Kxx(i,0) = recv(c)
            end do
            do i = 1, localncols
                c = c + 1
                Kyy(i,0) = recv(c)
            end do
            deallocate(recv)
        end if

        if(pcol /= pncols) then
            call MPI_WAIT(requestr, status, ierr)
        end if
        if(prow /= pnrows) then
            call MPI_WAIT(requestu, status, ierr)
        end if

    end subroutine computeK

    subroutine computeKEdge()

        integer :: indexl, indexr, indexd, indexu
        integer :: i, j

        indexd = 1
        indexu = localnrows
        if(pncols == 1) then
            KxxEdge(1,indexd:indexu) = Kxx(1,indexd:indexu)
            do j = indexd, indexu
                do i = 2, localncols
                    KxxEdge(i,j) = (hx(i-1)+hx(i)) / (hx(i-1)/Kxx(i-1,j)+hx(i)/Kxx(i,j))
                end do
            end do
            KxxEdge(localncols+1,indexd:indexu) = Kxx(localncols,indexd:indexu)
        elseif(pcol == 1) then
            KxxEdge(1,indexd:indexu) = Kxx(1,indexd:indexu)
            do j = indexd, indexu
                do i = 2, localncols
                    KxxEdge(i,j) = (hx(i-1)+hx(i)) / (hx(i-1)/Kxx(i-1,j)+hx(i)/Kxx(i,j))
                end do
            end do
        elseif(pcol == pncols) then
            do j = indexd, indexu
                do i = 1, localncols
                    KxxEdge(i,j) = (hx(i-1)+hx(i)) / (hx(i-1)/Kxx(i-1,j)+hx(i)/Kxx(i,j))
                end do
            end do
            KxxEdge(localncols+1,indexd:indexu) = Kxx(localncols,indexd:indexu)
        else
            do j = indexd, indexu
                do i = 1, localncols
                    KxxEdge(i,j) = (hx(i-1)+hx(i)) / (hx(i-1)/Kxx(i-1,j)+hx(i)/Kxx(i,j))
                end do
            end do
        end if

        indexl = 1
        indexr = localncols
        if(pnrows == 1) then
            KyyEdge(indexl:indexr,1) = Kyy(indexl:indexr,1)
            do j = 2, localnrows
                do i = indexl, indexr
                    KyyEdge(i,j) = (hy(j-1)+hy(j)) / (hy(j-1)/Kyy(i,j-1)+hy(j)/Kyy(i,j))
                end do
            end do
            KyyEdge(indexl:indexr,localnrows+1) = Kyy(indexl:indexr,localnrows)
        elseif(prow == 1) then
            KyyEdge(indexl:indexr,1) = Kyy(indexl:indexr,1)
            do j = 2, localnrows
                do i = indexl, indexr
                    KyyEdge(i,j) = (hy(j-1)+hy(j)) / (hy(j-1)/Kyy(i,j-1)+hy(j)/Kyy(i,j))
                end do
            end do
        elseif(prow == pnrows) then
            do j = 1, localnrows
                do i = indexl, indexr
                    KyyEdge(i,j) = (hy(j-1)+hy(j)) / (hy(j-1)/Kyy(i,j-1)+hy(j)/Kyy(i,j))
                end do
            end do
            KyyEdge(indexl:indexr,localnrows+1) = Kyy(indexl:indexr,localnrows)
        else
            do j = 1, localnrows
                do i = indexl, indexr
                    KyyEdge(i,j) = (hy(j-1)+hy(j)) / (hy(j-1)/Kyy(i,j-1)+hy(j)/Kyy(i,j))
                end do
            end do
        end if

    end subroutine computeKEdge

    subroutine computeav()

        integer :: i, j

        ! suppose kxx=kyy. if you find the anisotropic permeability equation to compute av, you can change the
        ! computation here.
        do j = 1, localnrows
            do i = 1, localncols
                av(i,j) = avInit(xlower+i-1,ylower+j-1)*poro(i,j)/poroInit(xlower+i-1,ylower+j-1)* &!
                    dsqrt(1.D0*KxxInit(xlower+i-1,ylower+j-1)*poro(i,j)/Kxx(i,j)/poroInit(xlower+i-1,ylower+j-1))
            end do
        end do

    end subroutine computeav

    ! Generate the values on the right-hand side and the coefficients of the matrix A,
    ! and the subroutine will generate the values that will change with the time iteration.
    subroutine genDynPara_vp()
       
        integer :: findexl, findexr, findexd, findexu
        integer :: eindexl, eindexr, eindexd, eindexu
        integer :: global_ind, phase
        integer, dimension(:,:), pointer :: velx, vely, pres
        logical :: isField
        real(kind=8) :: miuw_l, miuw_r, miuw_d, miuw_u, miun_l, miun_r, miun_d, miun_u
        real(kind=8), dimension(:,:), pointer :: dmiuwdx, dmiuwdy, dmiundx, dmiundy
        real(kind=8), dimension(:,:), pointer :: rhs_velxw, rhs_velxn, rhs_velyw, rhs_velyn
        real(kind=8), dimension(:,:), pointer :: resiAxx, resiAyy, resiAcp, resitemp
        integer :: AxxwBeInd, AxxnBeInd, AyywBeInd, AyynBeInd
        integer :: i, j, n, c
        
        if(pcol /= 1) then
            findexl = 0
        else
            findexl = 1
        end if
        findexr = localncols + 1
        if(prow /= 1) then
            findexd = 0
        else
            findexd = 1
        end if
        if(prow /= pnrows) then
            findexu = localnrows + 1
        else
            findexu = localnrows
        end if
        allocate(velx(findexl:findexr,findexd:findexu))

        eindexl = 1
        if(pcol /= pncols) then
            eindexr = localncols
        else
            eindexr = localncols + 1
        end if
        eindexd = 1
        eindexu = localnrows
        allocate(rhs_velxw(eindexl:eindexr,eindexd:eindexu))
        allocate(rhs_velxn(eindexl:eindexr,eindexd:eindexu))
        allocate(resiAxx(eindexl:eindexr,eindexd:eindexu))
        allocate(resitemp(eindexl:eindexr,eindexd:eindexu))
        do phase = 1, 2
            velx(:,:) = 0
            if(phase == 1) then
                call Resi_xmom_velx(velx, rhs_velxw, phase)
            else
                call Resi_xmom_velx(velx, rhs_velxn, phase)
            end if
            do j = findexd, findexd+2
                do i = findexl, findexl+2
                    call genExpField(i, j, findexr, findexu, velx, isField)
                    if(isField) then
                        call Resi_xmom_velx(velx, resiAxx, phase)
                        ! generate the dynamic coefficients of Axx. The scheme is just like the DNA copy according
                        ! to a template, here the AxxCols and the AxxRows are just the template.
                        if(phase == 1) then
                            resitemp = resiAxx - rhs_velxw
                            call constructAxx(velx, resitemp, 1112)
                        else
                            resitemp = resiAxx - rhs_velxn
                            call constructAxx(velx, resitemp, 1122)
                        end if
                    end if
                end do
            end do
        end do
        deallocate(velx)
        deallocate(resiAxx)
        deallocate(resitemp)

        AxxwBeInd = 1
        AxxnBeInd = 1
        do j = eindexd, eindexu
            do i = eindexl, eindexr
                if((.not.((pcol==1).and.(i==1).and.(isDiriX0_p(ylower+j-1)==0))).and.(.not.((pcol==pncols) &!
                    .and.(i==localncols+1).and.(isDiriX1_p(ylower+j-1)==0)))) then
                    call index_convert_local_global(myid, 11, i, j, global_ind)
                    do n = AxxwBeInd, AxxSize
                        if(AxxwRows(n)==global_ind) then
                            AxxwValues(n) = AxxwStaticValues(n) + AxxwDynValues(n)
                        else
                            AxxwBeInd = n
                            exit
                        end if
                    end do
                    call index_convert_local_global(myid, 12, i, j, global_ind)
                    do n = AxxnBeInd, AxxSize
                        if(AxxnRows(n)==global_ind) then
                            AxxnValues(n) = AxxnStaticValues(n) + AxxnDynValues(n)
                        else
                            AxxnBeInd = n
                            exit
                        end if
                    end do
                elseif(((pcol==1).and.(i==1).and.(isDiriX0_p(ylower+j-1)==0)).or.((pcol==pncols).and. &!
                    (i==localncols+1).and.(isDiriX1_p(ylower+j-1)==0))) then
                    AxxwBeInd = AxxwBeInd + 1
                    AxxnBeInd = AxxnBeInd + 1
                end if
            end do
        end do

        allocate(dmiuwdx(eindexl:eindexr,eindexd:eindexu))
        allocate(dmiundx(eindexl:eindexr,eindexd:eindexu))
        do j = eindexd, eindexu
            do i = eindexl, eindexr
                if((pcol==1).and.(i==1).and.(isDiriX0_Sw(ylower+j-1)==1)) then
                    miuw_l = gammaw*DLOG(max(SwBdryX0(ylower+j-1),1.D-7)) + gammawn*(1.D0-SwBdryX0(ylower+j-1))
                    miun_l = gamman*DLOG(max(1.D0-SwBdryX0(ylower+j-1),1.D-7)) + gammawn*SwBdryX0(ylower+j-1)
                elseif(.not.((pcol==1).and.(i==1).and.(isDiriX0_Sw(ylower+j-1)==0))) then
                    miuw_l = gammaw*DLOG(max(Sw(i-1,j),1.D-7)) + gammawn*(1.D0-Sw(i-1,j))
                    miun_l = gamman*DLOG(max(1.D0-Sw(i-1,j),1.D-7)) + gammawn*Sw(i-1,j)
                end if

                if((pcol==pncols).and.(i==localncols+1).and.(isDiriX1_Sw(ylower+j-1)==1)) then
                    miuw_r = gammaw*DLOG(max(SwBdryX1(ylower+j-1),1.D-7)) + gammawn*(1.D0-SwBdryX1(ylower+j-1))
                    miun_r = gamman*DLOG(max(1.D0-SwBdryX1(ylower+j-1),1.D-7)) + gammawn*SwBdryX1(ylower+j-1)
                elseif(.not.((pcol==pncols).and.(i==localncols+1).and.(isDiriX1_Sw(ylower+j-1)==0))) then
                    miuw_r = gammaw*DLOG(max(Sw(i,j),1.D-7)) + gammawn*(1.D0-Sw(i,j))
                    miun_r = gamman*DLOG(max(1.D0-Sw(i,j),1.D-7)) + gammawn*Sw(i,j)
                end if

                if((pcol==1).and.(i==1).and.(isDiriX0_Sw(ylower+j-1)==0)) then
                    dmiuwdx(i,j) = 0.D0
                    dmiundx(i,j) = 0.D0
                elseif((pcol==1).and.(i==1).and.(isDiriX0_Sw(ylower+j-1)==1)) then
                    dmiuwdx(i,j) = (miuw_r-miuw_l)/hx(i)
                    dmiundx(i,j) = (miun_r-miun_l)/hx(i)
                elseif((pcol==pncols).and.(i==localncols+1).and.(isDiriX1_Sw(ylower+j-1)==0)) then
                    dmiuwdx(i,j) = 0.D0
                    dmiundx(i,j) = 0.D0
                elseif((pcol==pncols).and.(i==localncols+1).and.(isDiriX1_Sw(ylower+j-1)==1)) then
                    dmiuwdx(i,j) = (miuw_r-miuw_l)/hx(i-1)
                    dmiundx(i,j) = (miun_r-miun_l)/hx(i-1)
                else
                    dmiuwdx(i,j) = (miuw_r-miuw_l)/((hx(i-1)+hx(i))/2.D0)
                    dmiundx(i,j) = (miun_r-miun_l)/((hx(i-1)+hx(i))/2.D0)
                end if
            end do
        end do

        c = 0
        do j = eindexd, eindexu
            do i = eindexl, eindexr
                c = c + 1
                if(.not.(((pcol==1).and.(i==1).and.(isDiriX0_p(ylower+j-1)==0)).or.((pcol==pncols).and. &!
                    (i==localncols+1).and.(isDiriX1_p(ylower+j-1)==0)))) then
                    local_rhs_vp(c) = local_rhs_vp_static(c) - rhs_velxw(i,j) + dmiuwdx(i,j)
                end if
            end do
        end do

        do j = eindexd, eindexu
            do i = eindexl, eindexr
                c = c + 1
                if(.not.(((pcol==1).and.(i==1).and.(isDiriX0_p(ylower+j-1)==0)).or.((pcol==pncols).and. &!
                    (i==localncols+1).and.(isDiriX1_p(ylower+j-1)==0)))) then
                    local_rhs_vp(c) = local_rhs_vp_static(c) - rhs_velxn(i,j) + dmiundx(i,j)
                end if
            end do
        end do

        deallocate(dmiuwdx)
        deallocate(dmiundx)
        deallocate(rhs_velxw)
        deallocate(rhs_velxn)

        if(pcol /= 1) then
            findexl = 0
        else
            findexl = 1
        end if
        if(pcol /= pncols) then
            findexr = localncols + 1
        else
            findexr = localncols
        end if
        if(prow /= 1) then
            findexd = 0
        else
            findexd = 1
        end if
        findexu = localnrows + 1
        allocate(vely(findexl:findexr,findexd:findexu))
        eindexl = 1
        eindexr = localncols
        eindexd = 1
        if(prow /= pnrows) then
            eindexu = localnrows
        else
            eindexu = localnrows + 1
        end if
        allocate(rhs_velyw(eindexl:eindexr,eindexd:eindexu))
        allocate(rhs_velyn(eindexl:eindexr,eindexd:eindexu))
        allocate(resiAyy(eindexl:eindexr,eindexd:eindexu))
        allocate(resitemp(eindexl:eindexr,eindexd:eindexu))
        do phase = 1, 2
            vely(:,:) = 0
            if(phase == 1) then
                call Resi_ymom_vely(vely, rhs_velyw, phase)
            else
                call Resi_ymom_vely(vely, rhs_velyn, phase)
            end if
            do j = findexd, findexd+2
                do i = findexl, findexl+2
                    call genExpField(i, j, findexr, findexu, vely, isField)
                    if(isField) then
                        call Resi_ymom_vely(vely, resiAyy, phase)
                        if(phase == 1) then
                            resitemp = resiAyy - rhs_velyw
                            call constructAyy(vely, resitemp, 2112)
                        else
                            resitemp = resiAyy - rhs_velyn
                            call constructAyy(vely, resitemp, 2122)
                        end if
                    end if
                end do
            end do
        end do
        deallocate(vely)
        deallocate(resiAyy)
        deallocate(resitemp)

        AyywBeInd = 1
        AyynBeInd = 1
        do j = eindexd, eindexu
            do i = eindexl, eindexr
                if((.not.((prow==1).and.(j==1).and.(isDiriY0_p(xlower+i-1)==0))).and.(.not.((prow==pnrows).and. &!
                    (j==localnrows+1).and.(isDiriY1_p(xlower+i-1)==0)))) then
                    call index_convert_local_global(myid, 21, i, j, global_ind)
                    do n = AyywBeInd, AyySize
                        if(AyywRows(n)==global_ind) then
                            AyywValues(n) = AyywStaticValues(n) + AyywDynValues(n)
                        else
                            AyywBeInd = n
                            exit
                        end if
                    end do
                    call index_convert_local_global(myid, 22, i, j, global_ind)
                    do n = AyynBeInd, AyySize
                        if(AyynRows(n)==global_ind) then
                            AyynValues(n) = AyynStaticValues(n) + AyynDynValues(n)
                        else
                            AyynBeInd = n
                            exit
                        end if
                    end do
                elseif(((prow==1).and.(j==1).and.(isDiriY0_p(xlower+i-1)==0)).or.((prow==pnrows).and. &!
                    (j==localnrows+1).and.(isDiriY1_p(xlower+i-1)==0))) then
                    AyywBeInd = AyywBeInd + 1
                    AyynBeInd = AyynBeInd + 1
                end if
            end do
        end do

        allocate(dmiuwdy(eindexl:eindexr,eindexd:eindexu))
        allocate(dmiundy(eindexl:eindexr,eindexd:eindexu))
        do j = eindexd, eindexu
            do i = eindexl, eindexr
                if((prow==1).and.(j==1).and.(isDiriY0_Sw(xlower+i-1)==1)) then
                    miuw_d = gammaw*DLOG(max(SwBdryY0(xlower+i-1),1.D-7)) + gammawn*(1.D0-SwBdryY0(xlower+i-1))
                    miun_d = gamman*DLOG(max(1.D0-SwBdryY0(xlower+i-1),1.D-7)) + gammawn*SwBdryY0(xlower+i-1)
                elseif(.not.((prow==1).and.(j==1).and.(isDiriY0_Sw(xlower+i-1)==0))) then
                    miuw_d = gammaw*DLOG(max(Sw(i,j-1),1.D-7)) + gammawn*(1.D0-Sw(i,j-1))
                    miun_d = gamman*DLOG(max(1.D0-Sw(i,j-1),1.D-7)) + gammawn*Sw(i,j-1)
                end if

                if((prow==pnrows).and.(j==localnrows+1).and.(isDiriY1_Sw(xlower+i-1)==1)) then
                    miuw_u = gammaw*DLOG(max(SwBdryY1(xlower+i-1),1.D-7)) + gammawn*(1.D0-SwBdryY1(xlower+i-1))
                    miun_u = gamman*DLOG(max(1.D0-SwBdryY1(xlower+i-1),1.D-7)) + gammawn*SwBdryY1(xlower+i-1)
                elseif(.not.((prow==pnrows).and.(j==localnrows+1).and.(isDiriY1_Sw(xlower+i-1)==0))) then
                    miuw_u = gammaw*DLOG(max(Sw(i,j),1.D-7)) + gammawn*(1.D0-Sw(i,j))
                    miun_u = gamman*DLOG(max(1.D0-Sw(i,j),1.D-7)) + gammawn*Sw(i,j)
                end if

                if((prow==1).and.(j==1).and.(isDiriY0_Sw(xlower+i-1)==0)) then
                    dmiuwdy(i,j) = 0.D0
                    dmiundy(i,j) = 0.D0
                elseif((prow==1).and.(j==1).and.(isDiriY0_Sw(xlower+i-1)==1)) then
                    dmiuwdy(i,j) = (miuw_u-miuw_d)/hy(j)
                    dmiundy(i,j) = (miun_u-miun_d)/hy(j)
                elseif((prow==pnrows).and.(j==localnrows+1).and.(isDiriY1_Sw(xlower+i-1)==0)) then
                    dmiuwdy(i,j) = 0.D0
                    dmiundy(i,j) = 0.D0
                elseif((prow==pnrows).and.(j==localnrows+1).and.(isDiriY1_Sw(xlower+i-1)==1)) then
                    dmiuwdy(i,j) = (miuw_u-miuw_d)/hy(j-1)
                    dmiundy(i,j) = (miun_u-miun_d)/hy(j-1)
                else
                    dmiuwdy(i,j) = (miuw_u-miuw_d)/((hy(j-1)+hy(j))/2.D0)
                    dmiundy(i,j) = (miun_u-miun_d)/((hy(j-1)+hy(j))/2.D0)
                end if
            end do
        end do

        do j = eindexd, eindexu
            do i = eindexl, eindexr
                c = c + 1
                if(.not.(((prow==1).and.(j==1).and.(isDiriY0_p(xlower+i-1)==0)).or.((prow==pnrows).and. &!
                    (j==localnrows+1).and.(isDiriY1_p(xlower+i-1)==0)))) then
                    local_rhs_vp(c) = local_rhs_vp_static(c) - rhs_velyw(i,j) + dmiuwdy(i,j)
                end if
            end do
        end do
        
        do j = eindexd, eindexu
            do i = eindexl, eindexr
                c = c + 1
                if(.not.(((prow==1).and.(j==1).and.(isDiriY0_p(xlower+i-1)==0)).or.((prow==pnrows).and. &!
                    (j==localnrows+1).and.(isDiriY1_p(xlower+i-1)==0)))) then
                    local_rhs_vp(c) = local_rhs_vp_static(c) - rhs_velyn(i,j) + dmiundy(i,j)
                end if
            end do
        end do

        deallocate(dmiuwdy)
        deallocate(dmiundy)
        deallocate(rhs_velyw)
        deallocate(rhs_velyn)

        allocate(velx(1:localncols+1,1:localnrows))
        allocate(vely(1:localncols,1:localnrows+1))
        allocate(rhs_velxw(1:localncols,1:localnrows))
        allocate(rhs_velxn(1:localncols,1:localnrows))
        allocate(rhs_velyw(1:localncols,1:localnrows))
        allocate(rhs_velyn(1:localncols,1:localnrows))
        allocate(resiAxx(1:localncols,1:localnrows))
        allocate(resiAyy(1:localncols,1:localnrows))
        allocate(resitemp(1:localncols,1:localnrows))

        do phase = 1, 2
            velx(:,:) = 0
            vely(:,:) = 0
            if(phase == 1) then
                call Resi_consum_velx(velx, rhs_velxw, phase)
                call Resi_consum_vely(vely, rhs_velyw, phase)
            else
                call Resi_consum_velx(velx, rhs_velxn, phase)
                call Resi_consum_vely(vely, rhs_velyn, phase)
            end if
            do j = 1, 3
                do i = 1, 3
                    call genExpField(i, j, localncols+1, localnrows, velx, isField)
                    if(isField) then
                        call Resi_consum_velx(velx, resiAxx, phase)
                        if(phase == 1) then
                            resitemp = resiAxx - rhs_velxw
                            call constructAcx(velx, resitemp, 311)
                        else
                            resitemp = resiAxx - rhs_velxn
                            call constructAcx(velx, resitemp, 312)
                        end if
                    end if
                end do
            end do
            do j = 1, 3
                do i = 1, 3
                    call genExpField(i, j, localncols, localnrows+1, vely, isField)
                    if(isField) then
                        call Resi_consum_vely(vely, resiAyy, phase)
                        if(phase == 1) then
                            resitemp = resiAyy - rhs_velyw
                            call constructAcy(vely, resitemp, 321)
                        else
                            resitemp = resiAyy - rhs_velyn
                            call constructAcy(vely, resitemp, 322)
                        end if
                    end if
                end do
            end do
        end do

        deallocate(velx)
        deallocate(vely)
        deallocate(resiAxx)
        deallocate(resiAyy)
        deallocate(resitemp)

        do j = 1, localnrows
            do i = 1, localncols
                c = c + 1
                local_rhs_vp(c) = - rhs_velxw(i,j) - rhs_velxn(i,j) - rhs_velyw(i,j) - rhs_velyn(i,j) - &!
                    (poro(i,j)-poro_old(i,j))/(timeEnd/nt)
            end do
        end do

        deallocate(rhs_velxw)
        deallocate(rhs_velxn)
        deallocate(rhs_velyw)
        deallocate(rhs_velyn)

    end subroutine genDynPara_vp

    subroutine computevp()

        integer :: indexr, indexu
        integer :: nVelx, nVely, nPres
        integer :: AxxwBeInd, AxpwBeInd, AyywBeInd, AypwBeInd, AcxwBeInd, AcywBeInd, &!
            AxxnBeInd, AxpnBeInd, AyynBeInd, AypnBeInd, AcxnBeInd, AcynBeInd
        integer, dimension(:), allocatable :: cols
        real(kind=8), dimension(:), allocatable :: values
        real(kind=8), dimension(:), pointer :: local_x
        integer :: status(MPI_STATUS_SIZE)
        integer :: request, requestl, requestr, requestd, requestu, requestld, requestlu, requestrd
        integer, dimension(:), allocatable :: requestarray
        real(kind=8), dimension(:), allocatable :: slave_data
        real(kind=8), dimension(:), allocatable :: sent, recv
        integer :: sentSize, recvSize
        real(kind=8) :: solvertimestart, solvertimefinish
        integer :: i, j, l, n, c, num_iter, ierr

        allocate(cols(8))
        allocate(values(8)) ! Each row of A has at most 8 nonzero values.
        allocate(local_x(local_x_size_vp))

#ifdef LAPACK

        A_lapack_vp(:,:) = 0.D0
        b_lapack_vp(:) = local_rhs_vp(:)
        IPIV_vp(:) = 0.D0

#elif defined(UMFPACK)

        Ap_vp(1) = 0
        Ai_vp(:) = 0
        Ax_vp(:) = 0.D0

#elif defined(MUMPS)

        mumps_NNZ_loc_vp = 0
        mumps_IRN_loc_vp(:) = 0
        mumps_JCN_loc_vp(:) = 0
        mumps_A_loc_vp(:) = 0.D0

        if(mumps_par_vp%MYID /= 0) then
            call MPI_IBSEND(local_rhs_vp, local_x_size_vp, MPI_DOUBLE_PRECISION, 0, myid, &!
                MPI_COMM_WORLD, request, ierr)
        end if

        if(mumps_par_vp%MYID == 0) then
            mumps_par_vp%RHS(1:local_x_size_vp) = local_rhs_vp(1:local_x_size_vp)
            c = local_x_size_vp + 1
            do i = 1, nProcs-1
                allocate(slave_data(slave_vp_data_size(i)))
                call MPI_RECV(slave_data, slave_vp_data_size(i), MPI_DOUBLE_PRECISION, i, i, MPI_COMM_WORLD, status, ierr)
                mumps_par_vp%RHS(c:c+slave_vp_data_size(i)-1) = slave_data(:)
                c = c + slave_vp_data_size(i)
                deallocate(slave_data)
            end do
        end if

        if(mumps_par_vp%MYID /= 0) then
            call MPI_WAIT(request, status, ierr)
        end if

#elif defined(HYPRE)

        call HYPRE_IJVectorSetValues(b_vp, local_x_size_vp, rows_vp, local_rhs_vp, ierr)
        call HYPRE_IJVectorSetValues(x_vp, local_x_size_vp, rows_vp, initial_x_guess_vp, ierr)

        call HYPRE_IJVectorAssemble(b_vp, ierr)
        call HYPRE_IJVectorAssemble(x_vp, ierr)

        call HYPRE_IJVectorGetObject(b_vp, par_b_vp, ierr)
        call HYPRE_IJVectorGetObject(x_vp, par_x_vp, ierr)

#endif

        if(pcol /= pncols) then
            nVelx = localncols*localnrows
        else
            nVelx = (localncols+1)*localnrows
        end if
        if(prow /= pnrows) then
            nVely = localncols*localnrows
        else
            nVely = localncols*(localnrows+1)
        end if
        nPres = localncols*localnrows

        AxxwBeInd = 1
        AxpwBeInd = 1
        AyywBeInd = 1
        AypwBeInd = 1
        AcxwBeInd = 1
        AcywBeInd = 1
        AxxnBeInd = 1
        AxpnBeInd = 1
        AyynBeInd = 1
        AypnBeInd = 1
        AcxnBeInd = 1
        AcynBeInd = 1

        do n = ilower_vp, iupper_vp

            cols(:) = 0
            values(:) = 0
            c = 0

            ! the line is in the water-phase x-momentum part
            if(n <= ilower_vp+nVelx-1) then

                do l = AxxwBeInd, AxxSize
                    if(AxxwRows(l) == n) then
                        c = c + 1
                        cols(c) = AxxwCols(l)
                        values(c) = AxxwValues(l)
                    else
                        AxxwBeInd = l
                        exit
                    end if
                end do
                do l = AxpwBeInd, AxpSize
                    if(AxpwRows(l) == n) then
                        c = c + 1
                        cols(c) = AxpwCols(l)
                        values(c) = AxpwValues(l)
                    else
                        AxpwBeInd = l
                        exit
                    end if
                end do

            ! the line is in the oil-phase x-momentum part
            elseif((n >= ilower_vp+nVelx).and.(n <= ilower_vp+nVelx*2-1)) then

                do l = AxxnBeInd, AxxSize
                    if(AxxnRows(l) == n) then
                        c = c + 1
                        cols(c) = AxxnCols(l)
                        values(c) = AxxnValues(l)
                    else
                        AxxnBeInd = l
                        exit
                    end if
                end do
                do l = AxpnBeInd, AxpSize
                    if(AxpnRows(l) == n) then
                        c = c + 1
                        cols(c) = AxpnCols(l)
                        values(c) = AxpnValues(l)
                    else
                        AxpnBeInd = l
                        exit
                    end if
                end do

            ! the line is in the water-phase y-momentum part
            elseif((n >= ilower_vp+nVelx*2).and.(n <= ilower_vp+nVelx*2+nVely-1)) then

                do l = AyywBeInd, AyySize
                    if(AyywRows(l) == n) then
                        c = c + 1
                        cols(c) = AyywCols(l)
                        values(c) = AyywValues(l)
                    else
                        AyywBeInd = l
                        exit
                    end if
                end do
                do l = AypwBeInd, AypSize
                    if(AypwRows(l) == n) then
                        c = c + 1
                        cols(c) = AypwCols(l)
                        values(c) = AypwValues(l)
                    else
                        AypwBeInd = l
                        exit
                    end if
                end do

            ! the line is in the oil-phase y-momentum part
            elseif((n >= ilower_vp+nVelx*2+nVely).and.(n <= ilower_vp+nVelx*2+nVely*2-1)) then

                do l = AyynBeInd, AyySize
                    if(AyynRows(l) == n) then
                        c = c + 1
                        cols(c) = AyynCols(l)
                        values(c) = AyynValues(l)
                    else
                        AyynBeInd = l
                        exit
                    end if
                end do
                do l = AypnBeInd, AypSize
                    if(AypnRows(l) == n) then
                        c = c + 1
                        cols(c) = AypnCols(l)
                        values(c) = AypnValues(l)
                    else
                        AypnBeInd = l
                        exit
                    end if
                end do

            ! the line is in the continuity part
            elseif(n >= ilower_vp+nVelx*2+nVely*2) then

                do l = AcxwBeInd, AcxSize
                    if(AcxwRows(l) == n) then
                        c = c + 1
                        cols(c) = AcxwCols(l)
                        values(c) = AcxwValues(l)
                    else
                        AcxwBeInd = l
                        exit
                    end if
                end do
                do l = AcxnBeInd, AcxSize
                    if(AcxnRows(l) == n) then
                        c = c + 1
                        cols(c) = AcxnCols(l)
                        values(c) = AcxnValues(l)
                    else
                        AcxnBeInd = l
                        exit
                    end if
                end do
                do l = AcywBeInd, AcySize
                    if(AcywRows(l) == n) then
                        c = c + 1
                        cols(c) = AcywCols(l)
                        values(c) = AcywValues(l)
                    else
                        AcywBeInd = l
                        exit
                    end if
                end do
                do l = AcynBeInd, AcySize
                    if(AcynRows(l) == n) then
                        c = c + 1
                        cols(c) = AcynCols(l)
                        values(c) = AcynValues(l)
                    else
                        AcynBeInd = l
                        exit
                    end if
                end do

            end if

#ifdef LAPACK
            do l = 1, c
                A_lapack_vp(n,cols(l)) = values(l)
            end do
#elif defined(UMFPACK)
            Ap_vp(n+1) = Ap_vp(n) + c
            Ai_vp(Ap_vp(n)+1:Ap_vp(n+1)) = cols(1:c) - 1
            Ax_vp(Ap_vp(n)+1:Ap_vp(n+1)) = values(1:c)
#elif defined(MUMPS)
            mumps_IRN_loc_vp(mumps_NNZ_loc_vp+1:mumps_NNZ_loc_vp+c) = n
            mumps_JCN_loc_vp(mumps_NNZ_loc_vp+1:mumps_NNZ_loc_vp+c) = cols(1:c)
            mumps_A_loc_vp(mumps_NNZ_loc_vp+1:mumps_NNZ_loc_vp+c) = values(1:c)
            mumps_NNZ_loc_vp = mumps_NNZ_loc_vp + c
#elif defined(HYPRE)
            call HYPRE_IJMatrixSetValues(A_vp, 1, c, n, cols, values, ierr)
#endif
        end do

#ifdef LAPACK

        ! you have to make sure that the number of processors is set to 1 when using such method.
        solvertimestart = MPI_Wtime()
        call dgesv(local_x_size_vp, 1, A_lapack_vp, local_x_size_vp, IPIV_vp, b_lapack_vp, local_x_size_vp, LAPACKINFO)
        solvertimefinish = MPI_Wtime()
        solvertime = solvertime + solvertimefinish - solvertimestart
        local_x(:) = b_lapack_vp(:)

        if(LAPACKINFO /= 0) then
            print *, 'LAPACK solver error. INFO = ', LAPACKINFO
            stop
        end if

#elif defined(UMFPACK)

        solvertimestart = MPI_Wtime()
        call umf4def(control)
        call umf4sym(local_x_size_vp, local_x_size_vp, Ap_vp, Ai_vp, Ax_vp, symbolic, control, umfinfo)
        call umf4num(Ap_vp, Ai_vp, Ax_vp, symbolic, numeric, control, umfinfo)
        call umf4fsym(symbolic)
        call umf4sol(1, local_x, local_rhs_vp, numeric, control, umfinfo) ! 1 means A'x=b
        if(umfinfo(1) < 0) then
            print *, 'UMFPACK solver error. Info: ', umfinfo(1)
            stop
        end if
        call umf4fnum(numeric)
        solvertimefinish = MPI_Wtime()
        solvertime = solvertime + solvertimefinish - solvertimestart

#elif defined(MUMPS)

        mumps_par_vp%NNZ_loc = mumps_NNZ_loc_vp
        allocate(mumps_par_vp%IRN_loc(mumps_par_vp%NNZ_loc))
        allocate(mumps_par_vp%JCN_loc(mumps_par_vp%NNZ_loc))
        allocate(mumps_par_vp%A_loc(mumps_par_vp%NNZ_loc))
        mumps_par_vp%IRN_loc(1:mumps_par_vp%NNZ_loc) = mumps_IRN_loc_vp(1:mumps_NNZ_loc_vp)
        mumps_par_vp%JCN_loc(1:mumps_par_vp%NNZ_loc) = mumps_JCN_loc_vp(1:mumps_NNZ_loc_vp)
        mumps_par_vp%A_loc(1:mumps_par_vp%NNZ_loc) = mumps_A_loc_vp(1:mumps_NNZ_loc_vp)
        mumps_par_vp%JOB = 6
        solvertimestart = MPI_Wtime()
        call DMUMPS(mumps_par_vp)
        solvertimefinish = MPI_Wtime()
        solvertime = solvertime + solvertimefinish - solvertimestart
        if(mumps_par_vp%INFOG(1) < 0) then
            stop
        end if

        if(mumps_par_vp%MYID == 0) then
            local_x(:) = mumps_par_vp%RHS(1:local_x_size_vp)
            allocate(requestarray(nProcs-1))
            c = local_x_size_vp + 1
            do i = 1, nProcs-1
                allocate(slave_data(slave_vp_data_size(i)))
                slave_data(:) = mumps_par_vp%RHS(c:c+slave_vp_data_size(i)-1)
                call MPI_IBSEND(slave_data, slave_vp_data_size(i), MPI_DOUBLE_PRECISION, i, myid, &!
                    MPI_COMM_WORLD, requestarray(i), ierr)
                c = c + slave_vp_data_size(i)
                deallocate(slave_data)
            end do
        end if

        if(mumps_par_vp%MYID /= 0) then
            allocate(slave_data(local_x_size_vp))
            call MPI_RECV(slave_data, local_x_size_vp, MPI_DOUBLE_PRECISION, 0, 0, MPI_COMM_WORLD, status, ierr)
            local_x(:) = slave_data(:)
            deallocate(slave_data)
        end if

        if(mumps_par_vp%MYID == 0) then
            do i = 1, nProcs-1
                call MPI_WAIT(requestarray(i), status, ierr)
            end do
            deallocate(requestarray)
        end if

        deallocate(mumps_par_vp%IRN_loc)
        deallocate(mumps_par_vp%JCN_loc)
        deallocate(mumps_par_vp%A_loc)

#elif defined(HYPRE)

        call HYPRE_IJMatrixAssemble(A_vp, ierr)
        call HYPRE_IJMatrixGetObject(A_vp, parcsr_A_vp, ierr)

        solvertimestart = MPI_Wtime()
        call HYPRE_ParCSRGMRESSetup(solver_vp, parcsr_A_vp, par_b_vp, par_x_vp, ierr)
        call HYPRE_ParCSRGMRESSolve(solver_vp, parcsr_A_vp, par_b_vp, par_x_vp, ierr)
        solvertimefinish = MPI_Wtime()
        solvertime = solvertime + solvertimefinish - solvertimestart
        call HYPRE_ParCSRGMRESGetNumIteratio(solver_vp, num_iter, ierr)
        if(ierr /= 0) then
            if(myid == 0) then
                print *, 'HYPRE solver error. ierr = ', ierr
                stop
            end if
        end if

        call HYPRE_IJVectorGetValues(x_vp, local_x_size_vp, rows_vp, local_x, ierr)

        ! let the solution of this time step be the initial x guess in the next time step.
        ! by this way, the number of solver iteration steps can be reduced greatly.
        initial_x_guess_vp(:) = local_x(:)

#endif

        c = 0
        indexr = localncols
        if(pcol == pncols) then
            indexr = localncols + 1
        end if
        do j = 1, localnrows
            do i = 1, indexr
                c = c + 1
                vxw(i,j) = local_x(c)
            end do
        end do
        do j = 1, localnrows
            do i = 1, indexr
                c = c + 1
                vxn(i,j) = local_x(c)
            end do
        end do

        indexu = localnrows
        if(prow == pnrows) then
            indexu = localnrows + 1
        end if
        do j = 1, indexu
            do i = 1, localncols
                c = c + 1
                vyw(i,j) = local_x(c)
            end do
        end do
        do j = 1, indexu
            do i = 1, localncols
                c = c + 1
                vyn(i,j) = local_x(c)
            end do
        end do

        do j = 1, localnrows
            do i = 1, localncols
                c = c + 1
                p(i,j) = local_x(c)
            end do
        end do

        if(pcol /= 1) then
            if(prow /= pnrows) then
                indexu = localnrows
            else
                indexu = localnrows + 1
            end if
            sentSize = localnrows*2 + indexu*2
            allocate(sent(sentSize))
            c = 0
            do j = 1, localnrows
                c = c + 1
                sent(c) = vxw(1,j)
            end do
            do j = 1, localnrows
                c = c + 1
                sent(c) = vxn(1,j)
            end do
            do j = 1, indexu
                c = c + 1
                sent(c) = vyw(1,j)
            end do
            do j = 1, indexu
                c = c + 1
                sent(c) = vyn(1,j)
            end do
            call MPI_IBSEND(sent, sentSize, MPI_DOUBLE_PRECISION, myid-1, myid, MPI_COMM_WORLD, requestl, ierr)
            deallocate(sent)
        end if

        if(pcol /= pncols) then
            if(prow /= pnrows) then
                indexu = localnrows
            else
                indexu = localnrows + 1
            end if
            sentSize = indexu*2
            allocate(sent(sentSize))
            c = 0
            do j = 1, indexu
                c = c + 1
                sent(c) = vyw(localncols,j)
            end do
            do j = 1, indexu
                c = c + 1
                sent(c) = vyn(localncols,j)
            end do
            call MPI_IBSEND(sent, sentSize, MPI_DOUBLE_PRECISION, myid+1, myid, MPI_COMM_WORLD, requestr, ierr)
            deallocate(sent)
        end if

        if(prow /= 1) then
            if(pcol /= pncols) then
                indexr = localncols
            else
                indexr = localncols + 1
            end if
            sentSize = indexr*2 + localncols*2
            allocate(sent(sentSize))
            c = 0
            do i = 1, indexr
                c = c + 1
                sent(c) = vxw(i,1)
            end do
            do i = 1, indexr
                c = c + 1
                sent(c) = vxn(i,1)
            end do
            do i = 1, localncols
                c = c + 1
                sent(c) = vyw(i,1)
            end do
            do i = 1, localncols
                c = c + 1
                sent(c) = vyn(i,1)
            end do
            call MPI_IBSEND(sent, sentSize, MPI_DOUBLE_PRECISION, myid-pncols, myid, MPI_COMM_WORLD, requestd, ierr)
            deallocate(sent)
        end if

        if(prow /= pnrows) then
            if(pcol /= pncols) then
                indexr = localncols
            else
                indexr = localncols + 1
            end if
            sentSize = indexr*2
            allocate(sent(sentSize))
            c = 0
            do i = 1, indexr
                c = c + 1
                sent(c) = vxw(i,localnrows)
            end do
            do i = 1, indexr
                c = c + 1
                sent(c) = vxn(i,localnrows)
            end do
            call MPI_IBSEND(sent, sentSize, MPI_DOUBLE_PRECISION, myid+pncols, myid, MPI_COMM_WORLD, requestu, ierr)
            deallocate(sent)
        end if

        if((pcol/=1).and.(prow/=1)) then
            allocate(sent(4))
            sent(1) = vxw(1,1)
            sent(2) = vxn(1,1)
            sent(3) = vyw(1,1)
            sent(4) = vyn(1,1)
            call MPI_IBSEND(sent, 4, MPI_DOUBLE_PRECISION, myid-1-pncols, myid, MPI_COMM_WORLD, requestld, ierr)
            deallocate(sent)
        end if

        if((pcol/=1).and.(prow/=pnrows)) then
            allocate(sent(2))
            sent(1) = vxw(1,localnrows)
            sent(2) = vxn(1,localnrows)
            call MPI_IBSEND(sent, 2, MPI_DOUBLE_PRECISION, myid-1+pncols, myid, MPI_COMM_WORLD, requestlu, ierr)
            deallocate(sent)
        end if

        if((pcol/=pncols).and.(prow/=1)) then
            allocate(sent(2))
            sent(1) = vyw(localncols,1)
            sent(2) = vyn(localncols,1)
            call MPI_IBSEND(sent, 2, MPI_DOUBLE_PRECISION, myid+1-pncols, myid, MPI_COMM_WORLD, requestrd, ierr)
            deallocate(sent)
        end if

        if(pcol /= pncols) then
            if(prow /= pnrows) then
                indexu = localnrows
            else
                indexu = localnrows + 1
            end if
            recvSize = localnrows*2 + indexu*2
            allocate(recv(recvSize))
            call MPI_RECV(recv, recvSize, MPI_DOUBLE_PRECISION, myid+1, myid+1, MPI_COMM_WORLD, status, ierr)
            c = 0
            do j = 1, localnrows
                c = c + 1
                vxw(localncols+1,j) = recv(c)
            end do
            do j = 1, localnrows
                c = c + 1
                vxn(localncols+1,j) = recv(c)
            end do
            do j = 1, indexu
                c = c + 1
                vyw(localncols+1,j) = recv(c)
            end do
            do j = 1, indexu
                c = c + 1
                vyn(localncols+1,j) = recv(c)
            end do
            deallocate(recv)
        end if

        if(pcol /= 1) then
            if(prow /= pnrows) then
                indexu = localnrows
            else
                indexu = localnrows + 1
            end if
            recvSize = indexu*2
            allocate(recv(recvSize))
            call MPI_RECV(recv, recvSize, MPI_DOUBLE_PRECISION, myid-1, myid-1, MPI_COMM_WORLD, status, ierr)
            c = 0
            do j = 1, indexu
                c = c + 1
                vyw(0,j) = recv(c)
            end do
            do j = 1, indexu
                c = c + 1
                vyn(0,j) = recv(c)
            end do
            deallocate(recv)
        end if

        if(prow /= pnrows) then
            if(pcol /= pncols) then
                indexr = localncols
            else
                indexr = localncols + 1
            end if
            recvSize = indexr*2 + localncols*2
            allocate(recv(recvSize))
            call MPI_RECV(recv, recvSize, MPI_DOUBLE_PRECISION, myid+pncols, myid+pncols, MPI_COMM_WORLD, status, ierr)
            c = 0
            do i = 1, indexr
                c = c + 1
                vxw(i,localnrows+1) = recv(c)
            end do
            do i = 1, indexr
                c = c + 1
                vxn(i,localnrows+1) = recv(c)
            end do
            do i = 1, localncols
                c = c + 1
                vyw(i,localnrows+1) = recv(c)
            end do
            do i = 1, localncols
                c = c + 1
                vyn(i,localnrows+1) = recv(c)
            end do
            deallocate(recv)
        end if

        if(prow /= 1) then
            if(pcol /= pncols) then
                indexr = localncols
            else
                indexr = localncols + 1
            end if
            recvSize = indexr*2
            allocate(recv(recvSize))
            call MPI_RECV(recv, recvSize, MPI_DOUBLE_PRECISION, myid-pncols, myid-pncols, MPI_COMM_WORLD, status, ierr)
            c = 0
            do i = 1, indexr
                c = c + 1
                vxw(i,0) = recv(c)
            end do
            do i = 1, indexr
                c = c + 1
                vxn(i,0) = recv(c)
            end do
            deallocate(recv)
        end if

        if((pcol/=pncols).and.(prow/=pnrows)) then
            allocate(recv(4))
            call MPI_RECV(recv, 4, MPI_DOUBLE_PRECISION, myid+1+pncols, myid+1+pncols, MPI_COMM_WORLD, status, ierr)
            vxw(localncols+1,localnrows+1) = recv(1)
            vxn(localncols+1,localnrows+1) = recv(2)
            vyw(localncols+1,localnrows+1) = recv(3)
            vyn(localncols+1,localnrows+1) = recv(4)
            deallocate(recv)
        end if

        if((pcol/=pncols).and.(prow/=1)) then
            allocate(recv(2))
            call MPI_RECV(recv, 2, MPI_DOUBLE_PRECISION, myid+1-pncols, myid+1-pncols, MPI_COMM_WORLD, status, ierr)
            vxw(localncols+1,0) = recv(1)
            vxn(localncols+1,0) = recv(2)
            deallocate(recv)
        end if

        if((pcol/=1).and.(prow/=pnrows)) then
            allocate(recv(2))
            call MPI_RECV(recv, 2, MPI_DOUBLE_PRECISION, myid-1+pncols, myid-1+pncols, MPI_COMM_WORLD, status, ierr)
            vyw(0,localnrows+1) = recv(1)
            vyn(0,localnrows+1) = recv(2)
            deallocate(recv)
        end if

        if(pcol /= 1) then
            call MPI_WAIT(requestl, status, ierr)
        end if
        if(pcol /= pncols) then
            call MPI_WAIT(requestr, status, ierr)
        end if
        if(prow /= 1) then
            call MPI_WAIT(requestd, status, ierr)
        end if
        if(prow /= pnrows) then
            call MPI_WAIT(requestu, status, ierr)
        end if
        if((pcol/=1).and.(prow/=1)) then
            call MPI_WAIT(requestld, status, ierr)
        end if
        if((pcol/=1).and.(prow/=pnrows)) then
            call MPI_WAIT(requestlu, status, ierr)
        end if
        if((pcol/=pncols).and.(prow/=1)) then
            call MPI_WAIT(requestrd, status, ierr)
        end if

        deallocate(cols)
        deallocate(values)
        deallocate(local_x)

    end subroutine computevp

    subroutine genDynPara_Cf()

        integer :: findexl, findexr, findexd, findexu
        integer :: eindexl, eindexr, eindexd, eindexu
        integer, dimension(:,:), pointer :: conc
        logical :: isField
        real(kind=8), dimension(:,:), pointer :: rhs_Cf
        real(kind=8), dimension(:,:), pointer :: resiACf, resitemp
        integer :: i, j, c

        findexl = 0
        findexr = localncols + 1
        findexd = 0
        findexu = localnrows + 1
        allocate(conc(findexl:findexr,findexd:findexu))
        conc(:,:) = 0

        eindexl = 1
        eindexr = localncols
        eindexd = 1
        eindexu = localnrows
        allocate(rhs_Cf(eindexl:eindexr,eindexd:eindexu))
        rhs_Cf(:,:) = 0.D0

        call Resi_tran_Cf(conc, rhs_Cf)
        allocate(resiACf(eindexl:eindexr,eindexd:eindexu))
        allocate(resitemp(eindexl:eindexr,eindexd:eindexu))
        do j = findexd, findexd+2
            do i = findexl, findexl+2
                call genExpField(i, j, findexr, findexu, conc, isField)
                if(isField) then
                    call Resi_tran_Cf(conc, resiACf)
                    resitemp = resiACf - rhs_Cf
                    call constructACf(conc, resitemp)
                end if
            end do
        end do

        c = 0
        do j = 1, localnrows
            do i = 1, localncols
                c = c + 1
                local_rhs_Cf(c) = -rhs_Cf(i,j)
            end do
        end do

        deallocate(conc)
        deallocate(rhs_Cf)
        deallocate(resiACf)
        deallocate(resitemp)

    end subroutine genDynPara_Cf

    subroutine computeCf()

        integer, dimension(:), allocatable :: cols
        real(kind=8), dimension(:), allocatable :: values
        real(kind=8), dimension(:), pointer :: local_x
        integer :: ACfBeInd
        integer :: status(MPI_STATUS_SIZE)
        integer :: request
        integer, dimension(:), allocatable :: requestarray
        real(kind=8), dimension(:), allocatable :: slave_data
        real(kind=8) :: solvertimestart, solvertimefinish
        integer :: i, j, l, n, c, num_iter, ierr

        allocate(cols(9))
        allocate(values(9)) ! Each row of A has at most 9 nonzero values.
        allocate(local_x(local_x_size_Cf))

#ifdef LAPACK

        A_lapack_Cf(:,:) = 0.D0
        b_lapack_Cf(:) = local_rhs_Cf(:)
        IPIV_Cf(:) = 0.D0

#elif defined(UMFPACK)

        Ap_Cf(1) = 0
        Ai_Cf(:) = 0
        Ax_Cf(:) = 0.D0

#elif defined(MUMPS)

        mumps_NNZ_loc_Cf = 0
        mumps_IRN_loc_Cf(:) = 0
        mumps_JCN_loc_Cf(:) = 0
        mumps_A_loc_Cf(:) = 0.D0

        if(mumps_par_Cf%MYID /= 0) then
            call MPI_IBSEND(local_rhs_Cf, local_x_size_Cf, MPI_DOUBLE_PRECISION, 0, myid, &!
                MPI_COMM_WORLD, request, ierr)
        end if

        if(mumps_par_Cf%MYID == 0) then
            mumps_par_Cf%RHS(1:local_x_size_Cf) = local_rhs_Cf(1:local_x_size_Cf)
            c = local_x_size_Cf + 1
            do i = 1, nProcs-1
                allocate(slave_data(local_x_size_Cf))
                call MPI_RECV(slave_data, local_x_size_Cf, MPI_DOUBLE_PRECISION, i, i, MPI_COMM_WORLD, status, ierr)
                mumps_par_Cf%RHS(c:c+local_x_size_Cf-1) = slave_data(:)
                c = c + local_x_size_Cf
                deallocate(slave_data)
            end do
        end if

        if(mumps_par_Cf%MYID /= 0) then
            call MPI_WAIT(request, status, ierr)
        end if

#elif defined(HYPRE)

        call HYPRE_IJVectorSetValues(b_Cf, local_x_size_Cf, rows_Cf, local_rhs_Cf, ierr)
        call HYPRE_IJVectorSetValues(x_Cf, local_x_size_Cf, rows_Cf, initial_x_guess_Cf, ierr)

        call HYPRE_IJVectorAssemble(b_Cf, ierr)
        call HYPRE_IJVectorAssemble(x_Cf, ierr)

        call HYPRE_IJVectorGetObject(b_Cf, par_b_Cf, ierr)
        call HYPRE_IJVectorGetObject(x_Cf, par_x_Cf, ierr)

#endif

        ACfBeInd = 1
        do n = ilower_Cf, iupper_Cf

            cols(:) = 0
            values(:) = 0.D0

            c = 0
            do l = ACfBeInd, ACfSize
                if(ACfRows(l) == n) then
                    c = c + 1
                    cols(c) = ACfCols(l)
                    values(c) = ACfValues(l)
                else
                    ACfBeInd = l
                    exit
                end if
            end do

#ifdef LAPACK
            do l = 1, c
                A_lapack_Cf(n,cols(l)) = values(l)
            end do
#elif defined(UMFPACK)
            Ap_Cf(n+1) = Ap_Cf(n) + c
            Ai_Cf(Ap_Cf(n)+1:Ap_Cf(n+1)) = cols(1:c) - 1
            Ax_Cf(Ap_Cf(n)+1:Ap_Cf(n+1)) = values(1:c)
#elif defined(MUMPS)
            mumps_IRN_loc_Cf(mumps_NNZ_loc_Cf+1:mumps_NNZ_loc_Cf+c) = n
            mumps_JCN_loc_Cf(mumps_NNZ_loc_Cf+1:mumps_NNZ_loc_Cf+c) = cols(1:c)
            mumps_A_loc_Cf(mumps_NNZ_loc_Cf+1:mumps_NNZ_loc_Cf+c) = values(1:c)
            mumps_NNZ_loc_Cf = mumps_NNZ_loc_Cf + c
#elif defined(HYPRE)
            call HYPRE_IJMatrixSetValues(A_Cf, 1, c, n, cols, values, ierr)
#endif

        end do

#ifdef LAPACK

        ! you have to make sure that the number of processors is set to 1 when using such method.
        solvertimestart = MPI_Wtime()
        call dgesv(local_x_size_Cf, 1, A_lapack_Cf, local_x_size_Cf, IPIV_Cf, b_lapack_Cf, local_x_size_Cf, LAPACKINFO)
        solvertimefinish = MPI_Wtime()
        solvertime = solvertime + solvertimefinish - solvertimestart
        local_x(:) = b_lapack_Cf(:)
        if(LAPACKINFO /= 0) then
            print *, 'LAPACK solver error when solving Cf. INFO = ', LAPACKINFO
            stop
        end if

#elif defined(UMFPACK)

        solvertimestart = MPI_Wtime()
        call umf4def(control)
        call umf4sym(local_x_size_Cf, local_x_size_Cf, Ap_Cf, Ai_Cf, Ax_Cf, symbolic, control, umfinfo)
        call umf4num(Ap_Cf, Ai_Cf, Ax_Cf, symbolic, numeric, control, umfinfo)
        call umf4fsym(symbolic)
        call umf4sol(1, local_x, local_rhs_Cf, numeric, control, umfinfo) ! 1 means A'x=b
        if(umfinfo(1) < 0) then
            print *, 'UMFPACK solver error when solving Cf. Info: ', umfinfo(1)
            stop
        end if
        call umf4fnum(numeric)
        solvertimefinish = MPI_Wtime()
        solvertime = solvertime + solvertimefinish - solvertimestart

#elif defined(MUMPS)

        mumps_par_Cf%NNZ_loc = mumps_NNZ_loc_Cf
        allocate(mumps_par_Cf%IRN_loc(mumps_par_Cf%NNZ_loc))
        allocate(mumps_par_Cf%JCN_loc(mumps_par_Cf%NNZ_loc))
        allocate(mumps_par_Cf%A_loc(mumps_par_Cf%NNZ_loc))
        mumps_par_Cf%IRN_loc(1:mumps_par_Cf%NNZ_loc) = mumps_IRN_loc_Cf(1:mumps_NNZ_loc_Cf)
        mumps_par_Cf%JCN_loc(1:mumps_par_Cf%NNZ_loc) = mumps_JCN_loc_Cf(1:mumps_NNZ_loc_Cf)
        mumps_par_Cf%A_loc(1:mumps_par_Cf%NNZ_loc) = mumps_A_loc_Cf(1:mumps_NNZ_loc_Cf)

        mumps_par_Cf%JOB = 6
        solvertimestart = MPI_Wtime()
        call DMUMPS(mumps_par_Cf)
        solvertimefinish = MPI_Wtime()
        solvertime = solvertime + solvertimefinish - solvertimestart
        if(mumps_par_Cf%INFOG(1) < 0) then
            stop
        end if

        if(mumps_par_Cf%MYID == 0) then
            local_x(:) = mumps_par_Cf%RHS(1:local_x_size_Cf)
            allocate(requestarray(nProcs-1))
            c = local_x_size_Cf + 1
            do i = 1, nProcs-1
                allocate(slave_data(local_x_size_Cf))
                slave_data(:) = mumps_par_Cf%RHS(c:c+local_x_size_Cf-1)
                call MPI_IBSEND(slave_data, local_x_size_Cf, MPI_DOUBLE_PRECISION, i, myid, &!
                    MPI_COMM_WORLD, requestarray(i), ierr)
                c = c + local_x_size_Cf
                deallocate(slave_data)
            end do
        end if

        if(mumps_par_Cf%MYID /= 0) then
            allocate(slave_data(local_x_size_Cf))
            call MPI_RECV(slave_data, local_x_size_Cf, MPI_DOUBLE_PRECISION, 0, 0, MPI_COMM_WORLD, status, ierr)
            local_x(:) = slave_data(:)
            deallocate(slave_data)
        end if

        if(mumps_par_Cf%MYID == 0) then
            do i = 1, nProcs-1
                call MPI_WAIT(requestarray(i), status, ierr)
            end do
            deallocate(requestarray)
        end if

        deallocate(mumps_par_Cf%IRN_loc)
        deallocate(mumps_par_Cf%JCN_loc)
        deallocate(mumps_par_Cf%A_loc)

#elif defined(HYPRE)

        call HYPRE_IJMatrixAssemble(A_Cf, ierr)
        call HYPRE_IJMatrixGetObject(A_Cf, parcsr_A_Cf, ierr)

        solvertimestart = MPI_Wtime()
        call HYPRE_ParCSRGMRESSetup(solver_Cf, parcsr_A_Cf, par_b_Cf, par_x_Cf, ierr)
        call HYPRE_ParCSRGMRESSolve(solver_Cf, parcsr_A_Cf, par_b_Cf, par_x_Cf, ierr)
        solvertimefinish = MPI_Wtime()
        solvertime = solvertime + solvertimefinish - solvertimestart
        call HYPRE_ParCSRGMRESGetNumIteratio(solver_Cf, num_iter, ierr)
        if(ierr /= 0) then
            if(myid == 0) then
                print *, 'HYPRE solver error when solving Cf. ierr = ', ierr
                stop
            end if
        end if

        call HYPRE_IJVectorGetValues(x_Cf, local_x_size_Cf, rows_Cf, local_x, ierr)

        ! let the solution of this time step be the initial x guess in the next time step.
        ! by this way, the number of solver iteration steps can be reduced greatly.
        initial_x_guess_Cf(:) = local_x(:)

#endif

        c = 0
        do j = 1, localnrows
            do i = 1, localncols
                c = c + 1
                Cf(i,j) = local_x(c)
            end do
        end do

        deallocate(cols)
        deallocate(values) 
        deallocate(local_x)

    end subroutine computeCf

    subroutine genDynPara_Tem()

        integer :: findexl, findexr, findexd, findexu
        integer :: eindexl, eindexr, eindexd, eindexu
        integer, dimension(:,:), pointer :: tempe
        logical :: isField
        real(kind=8), dimension(:,:), pointer :: rhs_Tem
        real(kind=8), dimension(:,:), pointer :: resiATem, resitemp
        integer :: i, j, c

        findexl = 0
        findexr = localncols + 1
        findexd = 0
        findexu = localnrows + 1
        allocate(tempe(findexl:findexr,findexd:findexu))
        tempe(:,:) = 0

        eindexl = 1
        eindexr = localncols
        eindexd = 1
        eindexu = localnrows
        allocate(rhs_Tem(eindexl:eindexr,eindexd:eindexu))
        rhs_Tem(:,:) = 0.D0

        call Resi_ener_Tem(tempe, rhs_Tem)

        allocate(resiATem(eindexl:eindexr,eindexd:eindexu))
        allocate(resitemp(eindexl:eindexr,eindexd:eindexu))
        do j = findexd, findexd+2
            do i = findexl, findexl+2
                call genExpField(i, j, findexr, findexu, tempe, isField)
                if(isField) then
                    call Resi_ener_Tem(tempe, resiATem)
                    resitemp = resiATem - rhs_Tem
                    call constructATem(tempe, resitemp)
                end if
            end do
        end do

        c = 0
        do j = 1, localnrows
            do i = 1, localncols
                c = c + 1
                local_rhs_Tem(c) = -rhs_Tem(i,j)
            end do
        end do

        deallocate(tempe)
        deallocate(rhs_Tem)
        deallocate(resiATem)
        deallocate(resitemp)

    end subroutine genDynPara_Tem

    subroutine computeTem()

        integer, dimension(:), allocatable :: cols
        real(kind=8), dimension(:), allocatable :: values
        real(kind=8), dimension(:), pointer :: local_x
        integer :: ATemBeInd
        integer :: status(MPI_STATUS_SIZE)
        integer :: request
        integer, dimension(:), allocatable :: requestarray
        real(kind=8), dimension(:), allocatable :: slave_data
        real(kind=8) :: solvertimestart, solvertimefinish
        integer :: i, j, l, n, c, num_iter, ierr

        allocate(cols(5))
        allocate(values(5)) ! Each row of A has at most 5 nonzero values.
        allocate(local_x(local_x_size_Tem))

#ifdef LAPACK

        A_lapack_Tem(:,:) = 0.D0
        b_lapack_Tem(:) = local_rhs_Tem(:)
        IPIV_Tem(:) = 0.D0

#elif defined(UMFPACK)

        Ap_Tem(1) = 0
        Ai_Tem(:) = 0
        Ax_Tem(:) = 0.D0

#elif defined(MUMPS)

        mumps_NNZ_loc_Tem= 0
        mumps_IRN_loc_Tem(:) = 0
        mumps_JCN_loc_Tem(:) = 0
        mumps_A_loc_Tem(:) = 0.D0

        if(mumps_par_Tem%MYID /= 0) then
            call MPI_IBSEND(local_rhs_Tem, local_x_size_Tem, MPI_DOUBLE_PRECISION, 0, myid, &!
                MPI_COMM_WORLD, request, ierr)
        end if

        if(mumps_par_Tem%MYID == 0) then
            mumps_par_Tem%RHS(1:local_x_size_Tem) = local_rhs_Tem(1:local_x_size_Tem)
            c = local_x_size_Tem + 1
            do i = 1, nProcs-1
                allocate(slave_data(local_x_size_Tem))
                call MPI_RECV(slave_data, local_x_size_Tem, MPI_DOUBLE_PRECISION, i, i, MPI_COMM_WORLD, status, ierr)
                    mumps_par_Tem%RHS(c:c+local_x_size_Tem-1) = slave_data(:)
                c = c + local_x_size_Tem
                deallocate(slave_data)
            end do
        end if

        if(mumps_par_Tem%MYID /= 0) then
            call MPI_WAIT(request, status, ierr)
        end if

#elif defined(HYPRE)

        call HYPRE_IJVectorSetValues(b_Tem, local_x_size_Tem, rows_Tem, local_rhs_Tem, ierr)
        call HYPRE_IJVectorSetValues(x_Tem, local_x_size_Tem, rows_Tem, initial_x_guess_Tem, ierr)

        call HYPRE_IJVectorAssemble(b_Tem, ierr)
        call HYPRE_IJVectorAssemble(x_Tem, ierr)

        call HYPRE_IJVectorGetObject(b_Tem, par_b_Tem, ierr)
        call HYPRE_IJVectorGetObject(x_Tem, par_x_Tem, ierr)

#endif

        ATemBeInd = 1
        do n = ilower_Tem, iupper_Tem

            cols(:) = 0
            values(:) = 0.D0

            c = 0
            do l = ATemBeInd, ATemSize
                if(ATemRows(l) == n) then
                    c = c + 1
                    cols(c) = ATemCols(l)
                    values(c) = ATemValues(l)
                else
                    ATemBeInd = l
                    exit
                end if
            end do

#ifdef LAPACK
            do l = 1, c
                A_lapack_Tem(n,cols(l)) = values(l)
            end do
#elif defined(UMFPACK)
            Ap_Tem(n+1) = Ap_Tem(n) + c
            Ai_Tem(Ap_Tem(n)+1:Ap_Tem(n+1)) = cols(1:c) - 1
            Ax_Tem(Ap_Tem(n)+1:Ap_Tem(n+1)) = values(1:c)
#elif defined(MUMPS)
            mumps_IRN_loc_Tem(mumps_NNZ_loc_Tem+1:mumps_NNZ_loc_Tem+c) = n
            mumps_JCN_loc_Tem(mumps_NNZ_loc_Tem+1:mumps_NNZ_loc_Tem+c) = cols(1:c)
            mumps_A_loc_Tem(mumps_NNZ_loc_Tem+1:mumps_NNZ_loc_Tem+c) = values(1:c)
            mumps_NNZ_loc_Tem = mumps_NNZ_loc_Tem + c
#elif defined(HYPRE)
            call HYPRE_IJMatrixSetValues(A_Tem, 1, c, n, cols, values, ierr)
#endif

        end do

#ifdef LAPACK

        ! you have to make sure that the number of processors is set to 1 when using such method.
        solvertimestart = MPI_Wtime()
        call dgesv(local_x_size_Tem, 1, A_lapack_Tem, local_x_size_Tem, IPIV_Tem, b_lapack_Tem, local_x_size_Tem, LAPACKINFO)
        solvertimefinish = MPI_Wtime()
        solvertime = solvertime + solvertimefinish - solvertimestart
        local_x(:) = b_lapack_Tem(:)
        if(LAPACKINFO /= 0) then
            print *, 'LAPACK solver error when solving Tem. INFO = ', LAPACKINFO
            stop
        end if

#elif defined(UMFPACK)

        solvertimestart = MPI_Wtime()
        call umf4def(control)
        call umf4sym(local_x_size_Tem, local_x_size_Tem, Ap_Tem, Ai_Tem, Ax_Tem, symbolic, control, umfinfo)
        call umf4num(Ap_Tem, Ai_Tem, Ax_Tem, symbolic, numeric, control, umfinfo)
        call umf4fsym(symbolic)
        call umf4sol(1, local_x, local_rhs_Tem, numeric, control, umfinfo) ! 1 means A'x=b
        if(umfinfo(1) < 0) then
            print *, 'UMFPACK solver error when solving Tem. Info: ', umfinfo(1)
            stop
        end if
        call umf4fnum(numeric)
        solvertimefinish = MPI_Wtime()
        solvertime = solvertime + solvertimefinish - solvertimestart

#elif defined(MUMPS)

        mumps_par_Tem%NNZ_loc = mumps_NNZ_loc_Tem
        allocate(mumps_par_Tem%IRN_loc(mumps_par_Tem%NNZ_loc))
        allocate(mumps_par_Tem%JCN_loc(mumps_par_Tem%NNZ_loc))
        allocate(mumps_par_Tem%A_loc(mumps_par_Tem%NNZ_loc))
        mumps_par_Tem%IRN_loc(1:mumps_par_Tem%NNZ_loc) = mumps_IRN_loc_Tem(1:mumps_NNZ_loc_Tem)
        mumps_par_Tem%JCN_loc(1:mumps_par_Tem%NNZ_loc) = mumps_JCN_loc_Tem(1:mumps_NNZ_loc_Tem)
        mumps_par_Tem%A_loc(1:mumps_par_Tem%NNZ_loc) = mumps_A_loc_Tem(1:mumps_NNZ_loc_Tem)

        mumps_par_Tem%JOB = 6
        solvertimestart = MPI_Wtime()
        call DMUMPS(mumps_par_Tem)
        solvertimefinish = MPI_Wtime()
        solvertime = solvertime + solvertimefinish - solvertimestart
        if(mumps_par_Tem%INFOG(1) < 0) then
            stop
        end if

        if(mumps_par_Tem%MYID == 0) then
            local_x(:) = mumps_par_Tem%RHS(1:local_x_size_Tem)
            allocate(requestarray(nProcs-1))
            c = local_x_size_Tem + 1
            do i = 1, nProcs-1
                allocate(slave_data(local_x_size_Tem))
                slave_data(:) = mumps_par_Tem%RHS(c:c+local_x_size_Tem-1)
                call MPI_IBSEND(slave_data, local_x_size_Tem, MPI_DOUBLE_PRECISION, i, myid, &!
                    MPI_COMM_WORLD, requestarray(i), ierr)
                c = c + local_x_size_Tem
                deallocate(slave_data)
            end do
        end if

        if(mumps_par_Tem%MYID /= 0) then
            allocate(slave_data(local_x_size_Tem))
            call MPI_RECV(slave_data, local_x_size_Tem, MPI_DOUBLE_PRECISION, 0, 0, MPI_COMM_WORLD, status, ierr)
            local_x(:) = slave_data(:)
            deallocate(slave_data)
        end if

        if(mumps_par_Tem%MYID == 0) then
            do i = 1, nProcs-1
                call MPI_WAIT(requestarray(i), status, ierr)
            end do
            deallocate(requestarray)
        end if

        deallocate(mumps_par_Tem%IRN_loc)
        deallocate(mumps_par_Tem%JCN_loc)
        deallocate(mumps_par_Tem%A_loc)

#elif defined(HYPRE)

        call HYPRE_IJMatrixAssemble(A_Tem, ierr)
        call HYPRE_IJMatrixGetObject(A_Tem, parcsr_A_Tem, ierr)

        solvertimestart = MPI_Wtime()
        call HYPRE_ParCSRGMRESSetup(solver_Tem, parcsr_A_Tem, par_b_Tem, par_x_Tem, ierr)
        call HYPRE_ParCSRGMRESSolve(solver_Tem, parcsr_A_Tem, par_b_Tem, par_x_Tem, ierr)
        solvertimefinish = MPI_Wtime()
        solvertime = solvertime + solvertimefinish - solvertimestart
        call HYPRE_ParCSRGMRESGetNumIteratio(solver_Tem, num_iter, ierr)
        if(ierr /= 0) then
            if(myid == 0) then
                print *, 'HYPRE solver error when solving Tem. ierr = ', ierr
                stop
            end if
        end if

        call HYPRE_IJVectorGetValues(x_Tem, local_x_size_Tem, rows_Tem, local_x, ierr)

        ! let the solution of this time step be the initial x guess in the next time step.
        ! by this way, the number of solver iteration steps can be reduced greatly.
        initial_x_guess_Tem(:) = local_x(:)

#endif

        c = 0
        do j = 1, localnrows
            do i = 1, localncols
                c = c + 1
                Tem(i,j) = local_x(c)
            end do
        end do

        deallocate(cols)
        deallocate(values)
        deallocate(local_x)

    end subroutine computeTem

    subroutine outputHisData()

        real(kind=8), dimension(:), allocatable :: local_data, recv
        real(kind=8), dimension(:), pointer :: global_data
        integer :: local_data_size
        integer :: status(MPI_STATUS_SIZE)
        integer :: request
        real(kind=8) :: poroavg, Kxxavg, avavg, qavg, pavg, Cfavg, Temavg
        real(kind=8) :: localqsum, localpsum
        integer :: i, j, c, ierr

        local_data_size = 6
        allocate(local_data(local_data_size))
        local_data(:) = 0.D0
        do j = 1, localnrows
            do i = 1, localncols
                local_data(1) = local_data(1) + poro(i,j)
                local_data(2) = local_data(2) + Kxx(i,j)
                local_data(3) = local_data(3) + av(i,j)
                local_data(4) = local_data(4) + p(i,j)
                local_data(5) = local_data(5) + Cf(i,j)
                local_data(6) = local_data(6) + Tem(i,j)
            end do
        end do

        ! output the average values of poro, Kxx, av, p, Cf, Tem
        if(myid /= 0) then
            call MPI_IBSEND(local_data, local_data_size, MPI_DOUBLE_PRECISION, 0, myid, MPI_COMM_WORLD, request, ierr)
        end if

        if(myid == 0) then

            allocate(global_data(local_data_size*nProcs))
            allocate(recv(local_data_size))

            global_data(1:local_data_size) = local_data(1:local_data_size)
            c = local_data_size + 1
            do i = 1, nProcs-1
                call MPI_RECV(recv, local_data_size, MPI_DOUBLE_PRECISION, i, i, MPI_COMM_WORLD, status, ierr)
                global_data(c:c+local_data_size-1) = recv(:)
                c = c + local_data_size
            end do
            
            poroavg = 0.D0
            do c = 1, local_data_size*nProcs, local_data_size
                poroavg = poroavg + global_data(c)
            end do
            poroavg = poroavg/(nx*ny)

            Kxxavg = 0.D0
            do c = 2, local_data_size*nProcs, local_data_size
                Kxxavg = Kxxavg + global_data(c)
            end do
            Kxxavg = Kxxavg/(nx*ny)

            avavg = 0.D0
            do c = 3, local_data_size*nProcs, local_data_size
                avavg = avavg + global_data(c)
            end do
            avavg = avavg/(nx*ny)

            pavg = 0.D0
            do c = 4, local_data_size*nProcs, local_data_size
                pavg = pavg + global_data(c)
            end do
            pavg = pavg/(nx*ny)

            Cfavg = 0.D0
            do c = 5, local_data_size*nProcs, local_data_size
                Cfavg = Cfavg + global_data(c)
            end do
            Cfavg = Cfavg/(nx*ny)

            Temavg = 0.D0
            do c = 6, local_data_size*nProcs, local_data_size
                Temavg = Temavg + global_data(c)
            end do
            Temavg = Temavg/(nx*ny)

            write(40, fmt='(es24.16)', iostat=ierr) poroavg
            if(ierr /= 0) then
                print *, 'write file error. ', ierr
                stop
            end if
            write(41, fmt='(es24.16)', iostat=ierr) Kxxavg
            if(ierr /= 0) then
                print *, 'write file error. ', ierr
                stop
            end if
            write(42, fmt='(es24.16)', iostat=ierr) avavg
            if(ierr /= 0) then
                print *, 'write file error. ', ierr
                stop
            end if
            write(43, fmt='(es24.16)', iostat=ierr) pavg
            if(ierr /= 0) then
                print *, 'write file error. ', ierr
                stop
            end if
            write(44, fmt='(es24.16)', iostat=ierr) Cfavg
            if(ierr /= 0) then
                print *, 'write file error. ', ierr
                stop
            end if
            write(45, fmt='(es24.16)', iostat=ierr) Temavg
            if(ierr /= 0) then
                print *, 'write file error. ', ierr
                stop
            end if

            deallocate(global_data)
            deallocate(recv)

        end if

        if(myid /= 0) then
            call MPI_WAIT(request, status, ierr)
        end if

        deallocate(local_data)

        ! output Q at the exit
        if(pcol == pncols) then

            localqsum = 0.D0
            do j = 1, localnrows
                localqsum = localqsum + vxw(localncols+1,j)
            end do

            if(myid /= 0) then
                call MPI_IBSEND(localqsum, 1, MPI_DOUBLE_PRECISION, 0, myid, MPI_COMM_WORLD, request, ierr)
            end if

        end if

        if(myid == 0) then

            allocate(global_data(pnrows))
            allocate(recv(1))

            if(pcol == pncols) then
                global_data(1) = localqsum
                c = 1
                do i = pncols-1, nProcs-1, pncols
                    if(i > 0) then
                        call MPI_RECV(recv, 1, MPI_DOUBLE_PRECISION, i, i, MPI_COMM_WORLD, status, ierr)
                        c = c + 1
                        global_data(c) = recv(1)
                    end if
                end do
            else
                c = 0
                do i = pncols-1, nProcs-1, pncols
                    call MPI_RECV(recv, 1, MPI_DOUBLE_PRECISION, i, i, MPI_COMM_WORLD, status, ierr)
                    c = c + 1
                    global_data(c) = recv(1)
                end do
            end if

            qavg = 0.D0
            do c = 1, pnrows
                qavg = qavg + global_data(c)
            end do
            qavg = qavg/ny

            write(46, fmt='(es24.16)', iostat=ierr) qavg
            if(ierr /= 0) then
                print *, 'write file error. ', ierr
                stop
            end if

            deallocate(global_data)
            deallocate(recv)

        end if

        if((pcol==pncols).and.(myid/=0)) then
            call MPI_WAIT(request, status, ierr)
        end if

        ! output p at the entry
        if(pcol == 1) then

            localpsum = 0.D0
            do j = 1, localnrows
                localpsum = localpsum + p(1,j)
            end do

            if(myid /= 0) then
                call MPI_IBSEND(localpsum, 1, MPI_DOUBLE_PRECISION, 0, myid, MPI_COMM_WORLD, request, ierr)
            end if

        end if

        if(myid == 0) then

            allocate(global_data(pnrows))
            allocate(recv(1))

            global_data(1) = localpsum
            c = 1
            do i = 0, nProcs-1, pncols
                if(i > 0) then
                    call MPI_RECV(recv, 1, MPI_DOUBLE_PRECISION, i, i, MPI_COMM_WORLD, status, ierr)
                    c = c + 1
                    global_data(c) = recv(1)
                end if
            end do

            pavg = 0.D0
            do c = 1, pnrows
                pavg = pavg + global_data(c)
            end do
            pavg = pavg/ny

            write(47, fmt='(es24.16)', iostat=ierr) pavg
            if(ierr /= 0) then
                print *, 'write file error. ', ierr
                stop
            end if

            deallocate(global_data)
            deallocate(recv)

        end if

        if((pcol==1).and.(myid/=0)) then
            call MPI_WAIT(request, status, ierr)
        end if

    end subroutine outputHisData

    ! output the raw results
    subroutine outputRawData()

        real(kind=8), dimension(:), allocatable :: local_data, recv
        real(kind=8), dimension(:), pointer :: global_data
        integer :: local_data_size, slave_data_size
        integer :: status(MPI_STATUS_SIZE)
        integer :: request
        integer :: indexr, indexu
        integer :: i, j, c, ierr

        local_data_size = local_x_size_vp+5*localncols*localnrows
        allocate(local_data(local_data_size))

        c = 0

        do j = 1, localnrows
            do i = 1, localncols
                c = c + 1
                local_data(c) = poro(i,j)
            end do
        end do

        do j = 1, localnrows
            do i = 1, localncols
                c = c + 1
                local_data(c) = Sw(i,j)
            end do
        end do

        do j = 1, localnrows
            do i = 1, localncols
                c = c + 1
                local_data(c) = Kxx(i,j)
            end do
        end do

        if(pcol /= pncols) then
            indexr = localncols
        else
            indexr = localncols + 1
        end if
        do j = 1, localnrows
            do i = 1, indexr
                c = c + 1
                local_data(c) = vxw(i,j)
            end do
        end do

        do j = 1, localnrows
            do i = 1, indexr
                c = c + 1
                local_data(c) = vxn(i,j)
            end do
        end do

        if(prow /= pnrows) then
            indexu = localnrows
        else
            indexu = localnrows + 1
        end if
        do j = 1, indexu
            do i = 1, localncols
                c = c + 1
                local_data(c) = vyw(i,j)
            end do
        end do

        do j = 1, indexu
            do i = 1, localncols
                c = c + 1
                local_data(c) = vyn(i,j)
            end do
        end do

        do j = 1, localnrows
            do i = 1, localncols
                c = c + 1
                local_data(c) = p(i,j)
            end do
        end do

        do j = 1, localnrows
            do i = 1, localncols
                c = c + 1
                local_data(c) = Cf(i,j)
            end do
        end do

        do j = 1, localnrows
            do i = 1, localncols
                c = c + 1
                local_data(c) = Tem(i,j)
            end do
        end do

        if(myid /= 0) then
            call MPI_IBSEND(local_data, local_data_size, MPI_DOUBLE_PRECISION, 0, myid, &!
                MPI_COMM_WORLD, request, ierr)
        end if

        if(myid == 0) then

            allocate(global_data((nx+1)*ny*2+nx*(ny+1)*2+6*nx*ny))

            global_data(1:local_data_size) = local_data(1:local_data_size)
            c = local_data_size + 1
            do i = 1, nProcs-1
                slave_data_size = slave_vp_data_size(i)+5*localncols*localnrows
                allocate(recv(slave_data_size))
                call MPI_RECV(recv, slave_data_size, MPI_DOUBLE_PRECISION, i, i, MPI_COMM_WORLD, status, ierr)
                global_data(c:c+slave_data_size-1) = recv(:)
                c = c + slave_data_size
                deallocate(recv)
            end do

            call exportResults(global_data)

            deallocate(global_data)

        end if

        if(myid /= 0) then
            call MPI_WAIT(request, status, ierr)
        end if

        deallocate(local_data)

    end subroutine outputRawData

    subroutine computePresDrop(isBT)

        logical, intent(out) :: isBT
        real(kind=8) :: inpavg, outpavg, localpsum, presDrop
        real(kind=8), dimension(:), allocatable :: global_data, recv
        integer :: status(MPI_STATUS_SIZE)
        integer :: request
        integer :: i, j, c, ierr

        if(pcol == 1) then

            localpsum = 0.D0
            do j = 1, localnrows
                localpsum = localpsum + p(1,j)
            end do

            if(myid /= 0) then
                call MPI_IBSEND(localpsum, 1, MPI_DOUBLE_PRECISION, 0, myid, MPI_COMM_WORLD, request, ierr)
            end if

        end if

        if(myid == 0) then

            allocate(global_data(pnrows))
            allocate(recv(1))

            global_data(1) = localpsum
            c = 1
            do i = 0, nProcs-1, pncols
                if(i > 0) then
                    call MPI_RECV(recv, 1, MPI_DOUBLE_PRECISION, i, i, MPI_COMM_WORLD, status, ierr)
                    c = c + 1
                    global_data(c) = recv(1)
                end if
            end do

            inpavg = 0.D0
            do c = 1, pnrows
                inpavg = inpavg + global_data(c)
            end do
            inpavg = inpavg/ny

            deallocate(global_data)
            deallocate(recv)

        end if

        if((pcol==1).and.(myid/=0)) then
            call MPI_WAIT(request, status, ierr)
        end if

        if(pcol == pncols) then

            localpsum = 0.D0
            do j = 1, localnrows
                localpsum = localpsum + p(localncols,j)
            end do

            if(myid /= 0) then
                call MPI_IBSEND(localpsum, 1, MPI_DOUBLE_PRECISION, 0, myid, MPI_COMM_WORLD, request, ierr)
            end if

        end if

        if(myid == 0) then

            allocate(global_data(pnrows))
            allocate(recv(1))

            if(pcol == pncols) then
                global_data(1) = localpsum
                c = 1
                do i = pncols-1, nProcs-1, pncols
                    if(i > 0) then
                        call MPI_RECV(recv, 1, MPI_DOUBLE_PRECISION, i, i, MPI_COMM_WORLD, status, ierr)
                        c = c + 1
                        global_data(c) = recv(1)
                    end if
                end do
            else
                c = 0
                do i = pncols-1, nProcs-1, pncols
                    call MPI_RECV(recv, 1, MPI_DOUBLE_PRECISION, i, i, MPI_COMM_WORLD, status, ierr)
                    c = c + 1
                    global_data(c) = recv(1)
                end do
            end if

            outpavg = 0.D0
            do c = 1, pnrows
                outpavg = outpavg + global_data(c)
            end do
            outpavg = outpavg/ny

            deallocate(global_data)
            deallocate(recv)

        end if

        if((pcol==pncols).and.(myid/=0)) then
            call MPI_WAIT(request, status, ierr)
        end if

        if(myid == 0) then
            isBT = .false.
            presDrop = abs(outpavg-inpavg)
            if(.not.isFindPresDropInit) then
                presDropInit = presDrop
                isFindPresDropInit = .true.
            end if
            if(presDrop/presDropInit < 1.D-2) then
                isBT = .true.
                print *, 'Breakthrough has been achieved! Program stops now.'
                print *, 'The breakthrough time is ', (t-1)*timeEnd/nt, ' seconds.'
                print *, 'The pore volume to breakthrough is ', (t-1)*timeEnd/nt*abs(vxwBdryX0(1))
            elseif(mod(t,100) == 0) then
                print *, 'The normalized pressure drop = ', presDrop/presDropInit*1.D2, '%'
            end if
        end if
        if(nProcs > 1) then
            call MPI_BCAST(isBT, 1, MPI_LOGICAL, 0, MPI_COMM_WORLD, ierr)
        end if

    end subroutine computePresDrop

    subroutine finalize()

        real(kind=8) :: timefinish
        integer :: ierr

        deallocate(xs)
        deallocate(ys)
        deallocate(ts)
        deallocate(src)
        deallocate(poroInit)
        deallocate(SwInit)
        deallocate(KxxInit)
        deallocate(KyyInit)
        deallocate(avInit)
        deallocate(vxwBdryX0)
        deallocate(vxwBdryX1)
        deallocate(vxwBdryY0)
        deallocate(vxwBdryY1)
        deallocate(vxnBdryX0)
        deallocate(vxnBdryX1)
        deallocate(vxnBdryY0)
        deallocate(vxnBdryY1)
        deallocate(vywBdryX0)
        deallocate(vywBdryX1)
        deallocate(vywBdryY0)
        deallocate(vywBdryY1)
        deallocate(vynBdryX0)
        deallocate(vynBdryX1)
        deallocate(vynBdryY0)
        deallocate(vynBdryY1)
        deallocate(isDiriX0_p)
        deallocate(isDiriX1_p)
        deallocate(isDiriY0_p)
        deallocate(isDiriY1_p)
        deallocate(pBdryX0)
        deallocate(pBdryX1)
        deallocate(pBdryY0)
        deallocate(pBdryY1)
        deallocate(pInit)
        deallocate(isDiriX0_Cf)
        deallocate(isDiriX1_Cf)
        deallocate(isDiriY0_Cf)
        deallocate(isDiriY1_Cf)
        deallocate(CfBdryX0)
        deallocate(CfBdryX1)
        deallocate(CfBdryY0)
        deallocate(CfBdryY1)
        deallocate(CfInit)
        deallocate(isDiriX0_Tem)
        deallocate(isDiriX1_Tem)
        deallocate(isDiriY0_Tem)
        deallocate(isDiriY1_Tem)
        deallocate(TemBdryX0)
        deallocate(TemBdryX1)
        deallocate(TemBdryY0)
        deallocate(TemBdryY1)
        deallocate(TemInit)

        deallocate(dm)
        deallocate(kc)
        deallocate(ks)
        deallocate(hx)
        deallocate(hy)
        deallocate(poro)
        deallocate(poro_old)
        deallocate(poroEdgeX)
        deallocate(poroEdgeX_old)
        deallocate(poroEdgeXInit)
        deallocate(poroEdgeY)
        deallocate(poroEdgeY_old)
        deallocate(poroEdgeYInit)
        deallocate(Sw)
        deallocate(Sw_old)
        deallocate(SwEdgeX)
        deallocate(SwEdgeX_old)
        deallocate(SwEdgeY)
        deallocate(SwEdgeY_old)
        deallocate(Kxx)
        deallocate(Kyy)
        deallocate(KxxEdge)
        deallocate(KyyEdge)
        deallocate(av)
        deallocate(vxw)
        deallocate(vxn)
        deallocate(vyw)
        deallocate(vyn)
        deallocate(p)
        deallocate(Cf)
        deallocate(Tem)
        
        deallocate(local_rhs_vp)
        deallocate(local_rhs_vp_static)
        deallocate(local_rhs_Cf)
        deallocate(local_rhs_Sw)
        deallocate(local_rhs_Tem)

        deallocate(AxxwCols)
        deallocate(AxxwRows)
        deallocate(AxxwStaticValues)
        deallocate(AxxwDynValues)
        deallocate(AxxwValues)
        deallocate(AxpwCols)
        deallocate(AxpwRows)
        deallocate(AxpwValues)
        deallocate(AyywCols)
        deallocate(AyywRows)
        deallocate(AyywStaticValues)
        deallocate(AyywDynValues)
        deallocate(AyywValues)
        deallocate(AypwCols)
        deallocate(AypwRows)
        deallocate(AypwValues)
        deallocate(AxxnCols)
        deallocate(AxxnRows)
        deallocate(AxxnStaticValues)
        deallocate(AxxnDynValues)
        deallocate(AxxnValues)
        deallocate(AxpnCols)
        deallocate(AxpnRows)
        deallocate(AxpnValues)
        deallocate(AyynCols)
        deallocate(AyynRows)
        deallocate(AyynStaticValues)
        deallocate(AyynDynValues)
        deallocate(AyynValues)
        deallocate(AypnCols)
        deallocate(AypnRows)
        deallocate(AypnValues)
        deallocate(AcxwCols)
        deallocate(AcxwRows)
        deallocate(AcxwValues)
        deallocate(AcywCols)
        deallocate(AcywRows)
        deallocate(AcywValues)
        deallocate(AcxnCols)
        deallocate(AcxnRows)
        deallocate(AcxnValues)
        deallocate(AcynCols)
        deallocate(AcynRows)
        deallocate(AcynValues)
        deallocate(ACfCols)
        deallocate(ACfRows)
        deallocate(ACfValues)
        deallocate(ASwCols)
        deallocate(ASwRows)
        deallocate(ASwValues)
        deallocate(ATemCols)
        deallocate(ATemRows)
        deallocate(ATemValues)

        deallocate(AxxEntryNum)
        deallocate(AxpEntryNum)
        deallocate(AyyEntryNum)
        deallocate(AypEntryNum)
        deallocate(AcxEntryNum)
        deallocate(AcyEntryNum)
        deallocate(ACfEntryNum)
        deallocate(ASwEntryNum)
        deallocate(ATemEntryNum)

        deallocate(AxxEntryBase)
        deallocate(AxpEntryBase)
        deallocate(AyyEntryBase)
        deallocate(AypEntryBase)
        deallocate(AcxEntryBase)
        deallocate(AcyEntryBase)
        deallocate(ACfEntryBase)
        deallocate(ASwEntryBase)
        deallocate(ATemEntryBase)

        if((nProcs>1).and.(myid==0)) then
            deallocate(slave_vp_data_size)
        end if

        if(myid == 0) then
            close(40)
            close(41)
            close(42)
            close(43)
            close(44)
            close(45)
            close(46)
            close(47)
        end if

#ifdef LAPACK

        deallocate(A_lapack_vp)
        deallocate(b_lapack_vp)
        deallocate(IPIV_vp)
        deallocate(A_lapack_Cf)
        deallocate(b_lapack_Cf)
        deallocate(IPIV_Cf)
        deallocate(A_lapack_Sw)
        deallocate(b_lapack_Sw)
        deallocate(IPIV_Sw)
        deallocate(A_lapack_Tem)
        deallocate(b_lapack_Tem)
        deallocate(IPIV_Tem)

#elif defined(UMFPACK)

        deallocate(Ap_vp)
        deallocate(Ai_vp)
        deallocate(Ax_vp)
        deallocate(Ap_Cf)
        deallocate(Ai_Cf)
        deallocate(Ax_Cf)
        deallocate(Ap_Sw)
        deallocate(Ai_Sw)
        deallocate(Ax_Sw)
        deallocate(Ap_Tem)
        deallocate(Ai_Tem)
        deallocate(Ax_Tem)

#elif defined(MUMPS)

        deallocate(mumps_par_vp%RHS)
        deallocate(mumps_par_Cf%RHS)
        deallocate(mumps_par_Sw%RHS)
        deallocate(mumps_par_Tem%RHS)

        mumps_par_vp%JOB = -2
        call DMUMPS(mumps_par_vp)
        mumps_par_Cf%JOB = -2
        call DMUMPS(mumps_par_Cf)
        mumps_par_Sw%JOB = -2
        call DMUMPS(mumps_par_Sw)
        mumps_par_Tem%JOB = -2
        call DMUMPS(mumps_par_Tem)

        deallocate(mumps_IRN_loc_vp)
        deallocate(mumps_JCN_loc_vp)
        deallocate(mumps_A_loc_vp)
        deallocate(mumps_IRN_loc_Cf)
        deallocate(mumps_JCN_loc_Cf)
        deallocate(mumps_A_loc_Cf)
        deallocate(mumps_IRN_loc_Sw)
        deallocate(mumps_JCN_loc_Sw)
        deallocate(mumps_A_loc_Sw)
        deallocate(mumps_IRN_loc_Tem)
        deallocate(mumps_JCN_loc_Tem)
        deallocate(mumps_A_loc_Tem)

#elif defined(HYPRE)

        call HYPRE_ParCSRPilutDestroy(precond_vp, ierr)
        call HYPRE_ParCSRPilutDestroy(precond_Cf, ierr)
        call HYPRE_ParCSRPilutDestroy(precond_Sw, ierr)
        call HYPRE_ParaSailsDestroy(precond_Tem, ierr)
        call HYPRE_ParCSRGMRESDestroy(solver_vp, ierr)
        call HYPRE_ParCSRGMRESDestroy(solver_Cf, ierr)
        call HYPRE_ParCSRGMRESDestroy(solver_Sw, ierr)
        call HYPRE_ParCSRGMRESDestroy(solver_Tem, ierr)

        call HYPRE_IJMatrixDestroy(A_vp, ierr)
        call HYPRE_IJVectorDestroy(b_vp, ierr)
        call HYPRE_IJVectorDestroy(x_vp, ierr)
        call HYPRE_IJMatrixDestroy(A_Cf, ierr)
        call HYPRE_IJVectorDestroy(b_Cf, ierr)
        call HYPRE_IJVectorDestroy(x_Cf, ierr)
        call HYPRE_IJMatrixDestroy(A_Sw, ierr)
        call HYPRE_IJVectorDestroy(b_Sw, ierr)
        call HYPRE_IJVectorDestroy(x_Sw, ierr)
        call HYPRE_IJMatrixDestroy(A_Tem, ierr)
        call HYPRE_IJVectorDestroy(b_Tem, ierr)
        call HYPRE_IJVectorDestroy(x_Tem, ierr)

        deallocate(rows_vp)
        deallocate(initial_x_guess_vp)
        deallocate(rows_Cf)
        deallocate(initial_x_guess_Cf)
        deallocate(rows_Sw)
        deallocate(initial_x_guess_Sw)
        deallocate(rows_Tem)
        deallocate(initial_x_guess_Tem)

#endif

#if defined(MUMPS) || defined(HYPRE)
        call MPI_BUFFER_DETACH(buffer, buffer_size, ierr)
#endif

        if(myid == 0) then
            print *, 'Solver time = ', solvertime, ' seconds.'
        end if
        call MPI_Barrier(MPI_COMM_WORLD, ierr)
        timefinish = MPI_Wtime()
        if(myid == 0) then
            print *, 'Elapsed time = ', timefinish-timestart, ' seconds.'
        end if 

        call MPI_Finalize(ierr)

    end subroutine finalize

    subroutine driver()

        logical :: isBT

        call initialize()

        call computedm()
        call computekc()

        ! time iteration
        do t = 2, nt+1

            call computeks()
            call computePoroEdge(1) ! old
            call computePoro()
            call computePoroEdge(2) ! new

            call computeSwEdge(1) ! old
            call genDynPara_Sw()
            call computeSw()
            call computeSwEdge(2) ! new

            call computeK()
            call computeKEdge()

            call computeav()

            call computekc()
            if(t == 2) then
                call genStaticPara_vp()
            end if
            call genDynPara_vp()
            call computevp()

            call computekc()
            call genDynPara_Cf()
            call computeCf()

            call genDynPara_Tem()
            call computeTem()
            call computedm()

if(Tem(1,1)/=Tem(1,1)) then
PRINT *, t
pause
end if

            call outputHisData()
            call computePresDrop(isBT)
            if(isBT) then
                call outputRawData()
                exit
            elseif(t == nt+1) then
                if(myid == 0) then
                    print *, 'Program stops now, but breakthrough has NOT been achieved.'
                    print *, 'More simulation time is needed to achieve breakthrough!'
                end if
            end if
            if((t==2).or.(mod(t-1,nt/NUMFRAME)==0)) then
                call outputRawData()
            end if

        end do

        call finalize()

    end subroutine driver

end module DBF_driver

