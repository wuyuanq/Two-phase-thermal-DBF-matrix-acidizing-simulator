
!!$ Author:
!!$   Yuanqing Wu, DGUT, P.R.China
!!$
!!$ History:
!!$   2025-5-9 by Yuanqing Wu
!!$
!!$ Support:
!!$   wuyuanq@gmail.com

module DBF_constructMat

    use DBF_globalData
    implicit none

Contains

    subroutine constructAxx(velx, resi, m_kind)

        integer, dimension(:,:), pointer, intent(in) :: velx
        real(kind=8), dimension(:,:), pointer, intent(in) :: resi
        integer, intent(in) :: m_kind

        integer :: lg_kind, fieldInd, equInd
        integer :: indexl, indexr, indexd, indexu, indextemp
        integer :: i, j

        if((m_kind == 1111).or.(m_kind == 1112)) then ! AxxwStatic, AxxwDyn
            lg_kind = 11
        elseif((m_kind == 1121).or.(m_kind == 1122)) then ! AxxnStatic, AxxnDyn
            lg_kind = 12
        else
            print *, 'm_kind in constructAxx is wrong!'
            stop
        end if

        ! the field index
        if(pcol /= 1) then
            indexl = 0
        else
            indexl = 1
        end if
        indexr = localncols+1
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

        if(pcol /= pncols) then
            indextemp = localncols
        else
            indextemp = localncols + 1
        end if

        do j = indexd, indexu
            do i = indexl, indexr

                if(velx(i,j) == 1) then

                    if((i==0).and.(j/=0).and.(j/=localnrows+1)) then
                        call index_convert_local_global(myid-1, lg_kind, localncols, j, fieldInd)
                        call index_convert_local_global(myid, lg_kind, 1, j, equInd)
                        call setMatValue(fieldInd, equInd, resi(i+1,j), m_kind, 1, j)
                    elseif((pcol/=pncols).and.(i==localncols+1).and.(j/=0).and.(j/=localnrows+1)) then
                        call index_convert_local_global(myid+1, lg_kind, 1, j, fieldInd)
                        call index_convert_local_global(myid, lg_kind, localncols, j, equInd)
                        call setMatValue(fieldInd, equInd, resi(i-1,j), m_kind, localncols, j)
                    elseif((j==0).and.(i/=0).and.(i/=indextemp+1).and.(.not.((pcol==1).and.(i==1).and. &!
                        (isDiriX0_p(ylower+j)==0))).and.(.not.((pcol==pncols).and.(i==indextemp).and. &!
                        (isDiriX1_p(ylower+j)==0)))) then
                        call index_convert_local_global(myid-pncols, lg_kind, i, localnrows, fieldInd)
                        call index_convert_local_global(myid, lg_kind, i, 1, equInd)
                        call setMatValue(fieldInd, equInd, resi(i,j+1), m_kind, i, 1)
                    elseif((j==localnrows+1).and.(i/=0).and.(i/=indextemp+1).and.(.not.((pcol==1).and. &!
                        (i==1).and.(isDiriX0_p(ylower+j-2)==0))).and.(.not.((pcol==pncols).and.(i==indextemp).and. &!
                        (isDiriX1_p(ylower+j-2)==0)))) then
                        call index_convert_local_global(myid+pncols, lg_kind, i, 1, fieldInd)
                        call index_convert_local_global(myid, lg_kind, i, localnrows, equInd)
                        call setMatValue(fieldInd, equInd, resi(i,j-1), m_kind, i, localnrows)
                    elseif((i>=1).and.(i<=indextemp).and.(j>=1).and.(j<=localnrows)) then
                        call index_convert_local_global(myid, lg_kind, i, j, fieldInd)
                        if(((j/=1).and.(i>=2).and.(i<=indextemp-1)).or.((j/=1).and.(i==1).and.(pcol/=1)).or. &!
                            ((j/=1).and.(i==1).and.(pcol==1).and.(isDiriX0_p(ylower+j-2)/=0)).or.((j/=1).and.(i==indextemp) &!
                            .and.(pcol/=pncols)).or.((j/=1).and.(i==indextemp).and.(pcol==pncols).and. &!
                            (isDiriX1_p(ylower+j-2)/=0))) then
                            call index_convert_local_global(myid, lg_kind, i, j-1, equInd)
                            call setMatValue(fieldInd, equInd, resi(i,j-1), m_kind, i, j-1)
                        end if
                        if((i>=3).or.((i==2).and.(pcol/=1)).or.((i==2).and.(pcol==1).and.(isDiriX0_p(ylower+j-1)/=0))) then
                            call index_convert_local_global(myid, lg_kind, i-1, j, equInd)
                            call setMatValue(fieldInd, equInd, resi(i-1,j), m_kind, i-1, j)
                        end if
                        call setMatValue(fieldInd, fieldInd, resi(i,j), m_kind, i, j)
                        if((i<=indextemp-2).or.((i==indextemp-1).and.(pcol/=pncols)).or.((i==indextemp-1).and. &!
                            (pcol==pncols).and.(isDiriX1_p(ylower+j-1)/=0))) then
                            call index_convert_local_global(myid, lg_kind, i+1, j, equInd)
                            call setMatValue(fieldInd, equInd, resi(i+1,j), m_kind, i+1, j)
                        end if
                        if(((j/=localnrows).and.(i>=2).and.(i<=indextemp-1)).or.((j/=localnrows).and.(i==1).and. &!
                            (pcol/=1)).or.((j/=localnrows).and.(i==1).and.(pcol==1).and.(isDiriX0_p(ylower+j)/=0)).or. &!
                            ((j/=localnrows).and.(i==indextemp).and.(pcol/=pncols)).or.((j/=localnrows).and. &!
                            (i==indextemp).and.(pcol==pncols).and.(isDiriX1_p(ylower+j)/=0))) then
                            call index_convert_local_global(myid, lg_kind, i, j+1, equInd)
                            call setMatValue(fieldInd, equInd, resi(i,j+1), m_kind, i, j+1)
                        end if
                    end if

                end if

            end do
        end do

    end subroutine constructAxx

    subroutine constructAxp(pres, resi, m_kind)

        ! notice that the dimensions of field are different from the dimensions of resi
        ! in the function
        integer, dimension(:,:), pointer, intent(in) :: pres
        real(kind=8), dimension(:,:), pointer, intent(in) :: resi
        integer, intent(in) :: m_kind

        integer :: lg_kind, fieldInd, equInd
        integer :: indexl, indexr, indexd, indexu
        integer :: i, j

        if(m_kind == 121) then ! Axpw
            lg_kind = 11
        elseif(m_kind == 122) then ! Axpn
            lg_kind = 12
        else
            print *, 'm_kind in constructAxp is wrong!'
            stop
        end if

        ! the field index
        if(pcol /= 1) then
            indexl = 0
        else
            indexl = 1
        end if
        indexr = localncols
        indexd = 1
        indexu = localnrows

        do j = indexd, indexu
            do i = indexl, indexr

                if(pres(i,j) == 1) then
                    if(i == 0) then
                        call index_convert_local_global(myid-1, 3, localncols, j, fieldInd)
                        call index_convert_local_global(myid, lg_kind, 1, j, equInd)
                        call setMatValue(fieldInd, equInd, resi(i+1,j), m_kind, 1, j)
                    elseif((i==1).and.(pcol==1)) then
                        call index_convert_local_global(myid, 3, i, j, fieldInd)
                        if(isDiriX0_p(ylower+j-1) /= 0) then
                            call index_convert_local_global(myid, lg_kind, i, j, equInd)
                            call setMatValue(fieldInd, equInd, resi(i,j), m_kind, i, j)
                        end if
                        call index_convert_local_global(myid, lg_kind, i+1, j, equInd)
                        call setMatValue(fieldInd, equInd, resi(i+1,j), m_kind, i+1, j)
                    elseif((i==localncols).and.(pcol/=pncols)) then
                        call index_convert_local_global(myid, 3, localncols, j, fieldInd)
                        call index_convert_local_global(myid, lg_kind, localncols, j, equInd)
                        call setMatValue(fieldInd, equInd, resi(i,j), m_kind, localncols, j)
                    elseif((i==localncols).and.(pcol==pncols)) then
                        call index_convert_local_global(myid, 3, i, j, fieldInd)
                        call index_convert_local_global(myid, lg_kind, i, j, equInd)
                        call setMatValue(fieldInd, equInd, resi(i,j), m_kind, i, j)
                        if(isDiriX1_p(ylower+j-1) /= 0) then
                            call index_convert_local_global(myid, lg_kind, i+1, j, equInd)
                            call setMatValue(fieldInd, equInd, resi(i+1,j), m_kind, i+1, j)
                        end if
                    else
                        call index_convert_local_global(myid, 3, i, j, fieldInd)
                        call index_convert_local_global(myid, lg_kind, i, j, equInd)
                        call setMatValue(fieldInd, equInd, resi(i,j), m_kind, i, j)
                        call index_convert_local_global(myid, lg_kind, i+1, j, equInd)
                        call setMatValue(fieldInd, equInd, resi(i+1,j), m_kind, i+1, j)
                    end if
                end if

            end do
        end do

    end subroutine constructAxp

    subroutine constructAyy(vely, resi, m_kind)

        integer, dimension(:,:), pointer, intent(in) :: vely
        real(kind=8), dimension(:,:), pointer, intent(in) :: resi
        integer, intent(in) :: m_kind

        integer :: lg_kind, fieldInd, equInd
        integer :: indexl, indexr, indexd, indexu, indextemp
        integer :: i, j

        if((m_kind == 2111).or.(m_kind == 2112)) then ! AyywStatic, AyywDyn
            lg_kind = 21
        elseif((m_kind == 2121).or.(m_kind == 2122)) then ! AyynStatic, AyynDyn
            lg_kind = 22
        else
            print *, 'm_kind in constructAyy is wrong!'
            stop
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
        if(prow /= 1) then
            indexd = 0
        else
            indexd = 1
        end if
        indexu = localnrows + 1

        if(prow /= pnrows) then
            indextemp = localnrows
        else
            indextemp = localnrows + 1
        end if

        do j = indexd, indexu
            do i = indexl, indexr

                if(vely(i,j) == 1) then

                    if((j==0).and.(i/=0).and.(i/=localncols+1)) then
                        call index_convert_local_global(myid-pncols, lg_kind, i, localnrows, fieldInd)
                        call index_convert_local_global(myid, lg_kind, i, 1, equInd)
                        call setMatValue(fieldInd, equInd, resi(i,j+1), m_kind, i, 1)
                    elseif((prow/=pnrows).and.(j==localnrows+1).and.(i/=0).and.(i/=localncols+1)) then
                        call index_convert_local_global(myid+pncols, lg_kind, i, 1, fieldInd)
                        call index_convert_local_global(myid, lg_kind, i, localnrows, equInd)
                        call setMatValue(fieldInd, equInd, resi(i,j-1), m_kind, i, localnrows)
                    elseif((i==0).and.(j/=0).and.(j/=indextemp+1).and.(.not.((prow==1).and.(j==1).and. &!
                        (isDiriY0_p(xlower+i)==0))).and.(.not.((prow==pnrows).and.(j==indextemp).and. &!
                        (isDiriY1_p(xlower+i)==0)))) then
                        call index_convert_local_global(myid-1, lg_kind, localncols, j, fieldInd)
                        call index_convert_local_global(myid, lg_kind, 1, j, equInd)
                        call setMatValue(fieldInd, equInd, resi(i+1,j), m_kind, 1, j)
                    elseif((i==localncols+1).and.(j/=0).and.(j/=indextemp+1).and.(.not.((prow==1).and.(j==1).and. &!
                        (isDiriY0_p(xlower+i-2)==0))).and.(.not.((prow==pnrows).and.(j==indextemp).and. &!
                        (isDiriY1_p(xlower+i-2)==0)))) then
                        call index_convert_local_global(myid+1, lg_kind, 1, j, fieldInd)
                        call index_convert_local_global(myid, lg_kind, localncols, j, equInd)
                        call setMatValue(fieldInd, equInd, resi(i-1,j), m_kind, localncols, j)
                    elseif((j>=1).and.(j<=indextemp).and.(i>=1).and.(i<=localncols)) then
                        call index_convert_local_global(myid, lg_kind, i, j, fieldInd)
                        if((j>=3).or.((j==2).and.(prow/=1)).or.((j==2).and.(prow==1).and.(isDiriY0_p(xlower+i-1)/=0))) then
                            call index_convert_local_global(myid, lg_kind, i, j-1, equInd)
                            call setMatValue(fieldInd, equInd, resi(i,j-1), m_kind, i, j-1)
                        end if
                        if(((i/=1).and.(j>=2).and.(j<=indextemp-1)).or.((i/=1).and.(j==1).and.(prow/=1)).or.((i/=1).and. &!
                            (j==1).and.(prow==1).and.(isDiriY0_p(xlower+i-2)/=0)).or.((i/=1).and.(j==indextemp).and. &!
                            (prow/=pnrows)).or.((i/=1).and.(j==indextemp).and.(prow==pnrows).and. &!
                            (isDiriY1_p(xlower+i-2)/=0))) then
                            call index_convert_local_global(myid, lg_kind, i-1, j, equInd)
                            call setMatValue(fieldInd, equInd, resi(i-1,j), m_kind, i-1, j)
                        end if
                        call setMatValue(fieldInd, fieldInd, resi(i,j), m_kind, i, j)
                        if(((i/=localncols).and.(j>=2).and.(j<=indextemp-1)).or.((i/=localncols).and.(j==1).and.(prow/=1)).or. &!
                            ((i/=localncols).and.(j==1).and.(prow==1).and.(isDiriY0_p(xlower+i)/=0)).or.((i/=localncols).and. &!
                            (j==indextemp).and.(prow/=pnrows)).or.((i/=localncols).and.(j==indextemp).and.(prow==pnrows) &!
                            .and.(isDiriY1_p(xlower+i)/=0))) then
                            call index_convert_local_global(myid, lg_kind, i+1, j, equInd)
                            call setMatValue(fieldInd, equInd, resi(i+1,j), m_kind, i+1, j)
                        end if
                        if((j<=indextemp-2).or.((j==indextemp-1).and.(prow/=pnrows)).or.((j==indextemp-1).and. &!
                            (prow==pnrows).and.(isDiriY1_p(xlower+i-1)/=0))) then
                            call index_convert_local_global(myid, lg_kind, i, j+1, equInd)
                            call setMatValue(fieldInd, equInd, resi(i,j+1), m_kind, i, j+1)
                        end if
                    end if

                end if

            end do
        end do

    end subroutine constructAyy

    subroutine constructAyp(pres, resi, m_kind)

        integer, dimension(:,:), pointer, intent(in) :: pres
        real(kind=8), dimension(:,:), pointer, intent(in) :: resi
        integer, intent(in) :: m_kind

        integer :: lg_kind, fieldInd, equInd
        integer :: indexl, indexr, indexd, indexu
        integer :: i, j

        if(m_kind == 221) then ! Aypw
            lg_kind = 21
        elseif(m_kind == 222) then ! Aypn
            lg_kind = 22
        else
            print *, 'm_kind in constructAyp is wrong!'
            stop
        end if

        ! the field index
        indexl = 1
        indexr = localncols
        if(prow /= 1) then
            indexd = 0
        else
            indexd = 1
        end if
        indexu = localnrows

        do j = indexd, indexu
            do i = indexl, indexr

                if(pres(i,j) == 1) then
                    if(j == 0) then
                        call index_convert_local_global(myid-pncols, 3, i, localnrows, fieldInd)
                        call index_convert_local_global(myid, lg_kind, i, 1, equInd)
                        call setMatValue(fieldInd, equInd, resi(i,j+1), m_kind, i, 1)
                    elseif((j==1).and.(prow==1)) then
                        call index_convert_local_global(myid, 3, i, j, fieldInd)
                        if(isDiriY0_p(xlower+i-1) /= 0) then
                            call index_convert_local_global(myid, lg_kind, i, j, equInd)
                            call setMatValue(fieldInd, equInd, resi(i,j), m_kind, i, j)
                        end if
                        call index_convert_local_global(myid, lg_kind, i, j+1, equInd)
                        call setMatValue(fieldInd, equInd, resi(i,j+1), m_kind, i, j+1)
                    elseif((j==localnrows).and.(prow/=pnrows)) then
                        call index_convert_local_global(myid, 3, i, localnrows, fieldInd)
                        call index_convert_local_global(myid, lg_kind, i, localnrows, equInd)
                        call setMatValue(fieldInd, equInd, resi(i,j), m_kind, i, localnrows)
                    elseif((j==localnrows).and.(prow==pnrows)) then
                        call index_convert_local_global(myid, 3, i, localnrows, fieldInd)
                        call index_convert_local_global(myid, lg_kind, i, localnrows, equInd)
                        call setMatValue(fieldInd, equInd, resi(i,j), m_kind, i, localnrows)
                        if(isDiriY1_p(xlower+i-1) /= 0) then
                            call index_convert_local_global(myid, lg_kind, i, localnrows+1, equInd)
                            call setMatValue(fieldInd, equInd, resi(i,j+1), m_kind, i, localnrows+1)
                        end if
                    else
                        call index_convert_local_global(myid, 3, i, j, fieldInd)
                        call index_convert_local_global(myid, lg_kind, i, j, equInd)
                        call setMatValue(fieldInd, equInd, resi(i,j), m_kind, i, j)
                        call index_convert_local_global(myid, lg_kind, i, j+1, equInd)
                        call setMatValue(fieldInd, equInd, resi(i,j+1), m_kind, i, j+1)
                    end if
                end if

            end do
        end do

    end subroutine constructAyp

    subroutine constructAcx(velx, resi, m_kind)

        integer, dimension(:,:), pointer, intent(in) :: velx
        real(kind=8), dimension(:,:), pointer, intent(in) :: resi
        integer, intent(in) :: m_kind
        
        integer :: lg_kind, fieldInd, equInd
        integer :: indexl, indexr, indexd, indexu
        integer :: i, j

        if(m_kind == 311) then ! Acxw
            lg_kind = 11
        elseif(m_kind == 312) then ! Acxn
            lg_kind = 12
        else
            print *, 'm_kind in constructAcx is wrong!'
            stop
        end if

        ! the field index
        indexl = 1
        indexr = localncols + 1
        indexd = 1
        indexu = localnrows

        do j = indexd, indexu
            do i = indexl, indexr

                if(velx(i,j) == 1) then
                    if(i == 1) then
                        call index_convert_local_global(myid, lg_kind, 1, j, fieldInd)
                        call index_convert_local_global(myid, 3, 1, j, equInd)
                        call setMatValue(fieldInd, equInd, resi(i,j), m_kind, 1, j)
                    elseif((i==localncols+1).and.(pcol/=pncols)) then
                        call index_convert_local_global(myid+1, lg_kind, 1, j, fieldInd)
                        call index_convert_local_global(myid, 3, localncols, j, equInd)
                        call setMatValue(fieldInd, equInd, resi(i-1,j), m_kind, localncols, j)
                    elseif((i==localncols+1).and.(pcol==pncols)) then
                        call index_convert_local_global(myid, lg_kind, localncols+1, j, fieldInd)
                        call index_convert_local_global(myid, 3, localncols, j, equInd)
                        call setMatValue(fieldInd, equInd, resi(i-1,j), m_kind, localncols, j)
                    else
                        call index_convert_local_global(myid, lg_kind, i, j, fieldInd)
                        call index_convert_local_global(myid, 3, i-1, j, equInd)
                        call setMatValue(fieldInd, equInd, resi(i-1,j), m_kind, i-1, j)
                        call index_convert_local_global(myid, 3, i, j, equInd)
                        call setMatValue(fieldInd, equInd, resi(i,j), m_kind, i, j)
                    end if
                end if

            end do
        end do

    end subroutine constructAcx

    subroutine constructAcy(vely, resi, m_kind)

        integer, dimension(:,:), pointer, intent(in) :: vely
        real(kind=8), dimension(:,:), pointer, intent(in) :: resi
        integer, intent(in) :: m_kind
        
        integer :: lg_kind, fieldInd, equInd
        integer :: indexl, indexr, indexd, indexu
        integer :: i, j

        if(m_kind == 321) then ! Acyw
            lg_kind = 21
        elseif(m_kind == 322) then ! Acyn
            lg_kind = 22
        else
            print *, 'm_kind in constructAcy is wrong!'
            stop
        end if

        ! the field index
        indexl = 1
        indexr = localncols
        indexd = 1
        indexu = localnrows + 1

        do j = indexd, indexu
            do i = indexl, indexr

                if(vely(i,j) == 1) then
                    if(j == 1) then
                        call index_convert_local_global(myid, lg_kind, i, 1, fieldInd)
                        call index_convert_local_global(myid, 3, i, 1, equInd)
                        call setMatValue(fieldInd, equInd, resi(i,j), m_kind, i, 1)
                    elseif((j==localnrows+1).and.(prow/=pnrows)) then
                        call index_convert_local_global(myid+pncols, lg_kind, i, 1, fieldInd)
                        call index_convert_local_global(myid, 3, i, localnrows, equInd)
                        call setMatValue(fieldInd, equInd, resi(i,j-1), m_kind, i, localnrows)
                    elseif((j==localnrows+1).and.(prow==pnrows)) then
                        call index_convert_local_global(myid, lg_kind, i, localnrows+1, fieldInd)
                        call index_convert_local_global(myid, 3, i, localnrows, equInd)
                        call setMatValue(fieldInd, equInd, resi(i,j-1), m_kind, i, localnrows)
                    else
                        call index_convert_local_global(myid, lg_kind, i, j, fieldInd)
                        call index_convert_local_global(myid, 3, i, j-1, equInd)
                        call setMatValue(fieldInd, equInd, resi(i,j-1), m_kind, i, j-1)
                        call index_convert_local_global(myid, 3, i, j, equInd)
                        call setMatValue(fieldInd, equInd, resi(i,j), m_kind, i, j)
                    end if
                end if

            end do
        end do

    end subroutine constructAcy

    subroutine constructACf(conc, resi)

        integer, dimension(:,:), pointer, intent(in) :: conc
        real(kind=8), dimension(:,:), pointer, intent(in) :: resi
        
        integer :: fieldInd, equInd
        integer :: indexl, indexr, indexd, indexu
        integer :: i, j

        ! the field index
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
        if(prow /= pnrows) then
            indexu = localnrows + 1
        else
            indexu = localnrows
        end if

        do j = indexd, indexu
            do i = indexl, indexr

                if(conc(i,j) == 1) then

                    if((i==0).and.(j/=0).and.(j/=localnrows+1)) then
                        call index_convert_local_global(myid-1, 4, localncols, j, fieldInd)
                        call index_convert_local_global(myid, 4, i+1, j, equInd)
                        call setMatValue(fieldInd, equInd, resi(i+1,j), 41, i+1, j)
                        if(j-1 >= 1) then
                            call index_convert_local_global(myid, 4, i+1, j-1, equInd)
                            call setMatValue(fieldInd, equInd, resi(i+1,j-1), 41, i+1, j-1)
                        end if
                        if(j+1 <= localnrows) then
                            call index_convert_local_global(myid, 4, i+1, j+1, equInd)
                            call setMatValue(fieldInd, equInd, resi(i+1,j+1), 41, i+1, j+1)
                        end if
                    elseif((i==localncols+1).and.(j/=0).and.(j/=localnrows+1)) then
                        call index_convert_local_global(myid+1, 4, 1, j, fieldInd)
                        call index_convert_local_global(myid, 4, i-1, j, equInd)
                        call setMatValue(fieldInd, equInd, resi(i-1,j), 41, i-1, j)
                        if(j-1 >= 1) then
                            call index_convert_local_global(myid, 4, i-1, j-1, equInd)
                            call setMatValue(fieldInd, equInd, resi(i-1,j-1), 41, i-1, j-1)
                        end if
                        if(j+1 <= localnrows) then
                            call index_convert_local_global(myid, 4, i-1, j+1, equInd)
                            call setMatValue(fieldInd, equInd, resi(i-1,j+1), 41, i-1, j+1)
                        end if
                    elseif((j==0).and.(i/=0).and.(i/=localncols+1)) then
                        call index_convert_local_global(myid-pncols, 4, i, localnrows, fieldInd)
                        call index_convert_local_global(myid, 4, i, j+1, equInd)
                        call setMatValue(fieldInd, equInd, resi(i,j+1), 41, i, j+1)
                        if(i-1 >= 1) then
                            call index_convert_local_global(myid, 4, i-1, j+1, equInd)
                            call setMatValue(fieldInd, equInd, resi(i-1,j+1), 41, i-1, j+1)
                        end if
                        if(i+1 <= localncols) then
                            call index_convert_local_global(myid, 4, i+1, j+1, equInd)
                            call setMatValue(fieldInd, equInd, resi(i+1,j+1), 41, i+1, j+1)
                        end if
                    elseif((j==localnrows+1).and.(i/=0).and.(i/=localncols+1)) then
                        call index_convert_local_global(myid+pncols, 4, i, 1, fieldInd)
                        call index_convert_local_global(myid, 4, i, j-1, equInd)
                        call setMatValue(fieldInd, equInd, resi(i,j-1), 41, i, j-1)
                        if(i-1 >= 1) then
                            call index_convert_local_global(myid, 4, i-1, j-1, equInd)
                            call setMatValue(fieldInd, equInd, resi(i-1,j-1), 41, i-1, j-1)
                        end if
                        if(i+1 <= localncols) then
                            call index_convert_local_global(myid, 4, i+1, j-1, equInd)
                            call setMatValue(fieldInd, equInd, resi(i+1,j-1), 41, i+1, j-1)
                        end if
                    elseif((i==0).and.(j==0)) then
                        call index_convert_local_global(myid-1-pncols, 4, localncols, localnrows, fieldInd)
                        call index_convert_local_global(myid, 4, 1, 1, equInd)
                        call setMatValue(fieldInd, equInd, resi(1,1), 41, 1, 1)
                    elseif((i==0).and.(j==localnrows+1)) then
                        call index_convert_local_global(myid-1+pncols, 4, localncols, 1, fieldInd)
                        call index_convert_local_global(myid, 4, 1, localnrows, equInd)
                        call setMatValue(fieldInd, equInd, resi(1,localnrows), 41, 1, localnrows)
                    elseif((i==localncols+1).and.(j==0)) then
                        call index_convert_local_global(myid+1-pncols, 4, 1, localnrows, fieldInd)
                        call index_convert_local_global(myid, 4, localncols, 1, equInd)
                        call setMatValue(fieldInd, equInd, resi(localncols,1), 41, localncols, 1)
                    elseif((i==localncols+1).and.(j==localnrows+1)) then
                        call index_convert_local_global(myid+1+pncols, 4, 1, 1, fieldInd)
                        call index_convert_local_global(myid, 4, localncols, localnrows, equInd)
                        call setMatValue(fieldInd, equInd, resi(localncols,localnrows), 41, localncols, localnrows)
                    else
                        call index_convert_local_global(myid, 4, i, j, fieldInd)
                        call index_convert_local_global(myid, 4, i, j, equInd)
                        call setMatValue(fieldInd, equInd, resi(i,j), 41, i, j)
                        if((i-1>=1).and.(j-1>=1)) then
                            call index_convert_local_global(myid, 4, i-1, j-1, equInd)
                            call setMatValue(fieldInd, equInd, resi(i-1,j-1), 41, i-1, j-1)
                        end if
                        if(j-1 >= 1) then
                            call index_convert_local_global(myid, 4, i, j-1, equInd)
                            call setMatValue(fieldInd, equInd, resi(i,j-1), 41, i, j-1)
                        end if
                        if((i+1<=localncols).and.(j-1>=1)) then
                            call index_convert_local_global(myid, 4, i+1, j-1, equInd)
                            call setMatValue(fieldInd, equInd, resi(i+1,j-1), 41, i+1, j-1)
                        end if
                        if(i-1 >= 1) then
                            call index_convert_local_global(myid, 4, i-1, j, equInd)
                            call setMatValue(fieldInd, equInd, resi(i-1,j), 41, i-1, j)
                        end if
                        if(i+1 <= localncols) then
                            call index_convert_local_global(myid, 4, i+1, j, equInd)
                            call setMatValue(fieldInd, equInd, resi(i+1,j), 41, i+1, j)
                        end if
                        if((i-1>=1).and.(j+1<=localnrows)) then
                            call index_convert_local_global(myid, 4, i-1, j+1, equInd)
                            call setMatValue(fieldInd, equInd, resi(i-1,j+1), 41, i-1, j+1)
                        end if
                        if(j+1 <= localnrows) then
                            call index_convert_local_global(myid, 4, i, j+1, equInd)
                            call setMatValue(fieldInd, equInd, resi(i,j+1), 41, i, j+1)
                        end if
                        if((i+1<=localncols).and.(j+1<=localnrows)) then
                            call index_convert_local_global(myid, 4, i+1, j+1, equInd)
                            call setMatValue(fieldInd, equInd, resi(i+1,j+1), 41, i+1, j+1)
                        end if
                    end if

                end if

            end do
        end do

    end subroutine constructACf

    subroutine constructASw(satw, resi)

        integer, dimension(:,:), pointer, intent(in) :: satw
        real(kind=8), dimension(:,:), pointer, intent(in) :: resi
    
        integer :: fieldInd, equInd
        integer :: indexl, indexr, indexd, indexu
        integer :: i, j

        ! the field index
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
        if(prow /= pnrows) then
            indexu = localnrows + 1
        else
            indexu = localnrows
        end if

        do j = indexd, indexu
            do i = indexl, indexr

                if(satw(i,j) == 1) then

                    if((i==0).and.(j/=0).and.(j/=localnrows+1)) then
                        call index_convert_local_global(myid-1, 5, localncols, j, fieldInd)
                        call index_convert_local_global(myid, 5, i+1, j, equInd)
                        call setMatValue(fieldInd, equInd, resi(i+1,j), 51, i+1, j)
                    elseif((i==localncols+1).and.(j/=0).and.(j/=localnrows+1)) then
                        call index_convert_local_global(myid+1, 5, 1, j, fieldInd)
                        call index_convert_local_global(myid, 5, i-1, j, equInd)
                        call setMatValue(fieldInd, equInd, resi(i-1,j), 51, i-1, j)
                    elseif((j==0).and.(i/=0).and.(i/=localncols+1)) then
                        call index_convert_local_global(myid-pncols, 5, i, localnrows, fieldInd)
                        call index_convert_local_global(myid, 5, i, j+1, equInd)
                        call setMatValue(fieldInd, equInd, resi(i,j+1), 51, i, j+1)
                    elseif((j==localnrows+1).and.(i/=0).and.(i/=localncols+1)) then
                        call index_convert_local_global(myid+pncols, 5, i, 1, fieldInd)
                        call index_convert_local_global(myid, 5, i, j-1, equInd)
                        call setMatValue(fieldInd, equInd, resi(i,j-1), 51, i, j-1)
                    elseif((i>=1).and.(i<=localncols).and.(j>=1).and.(j<=localnrows)) then
                        call index_convert_local_global(myid, 5, i, j, fieldInd)
                        call index_convert_local_global(myid, 5, i, j, equInd)
                        call setMatValue(fieldInd, equInd, resi(i,j), 51, i, j)
                        if(i-1 >= 1) then
                            call index_convert_local_global(myid, 5, i-1, j, equInd)
                            call setMatValue(fieldInd, equInd, resi(i-1,j), 51, i-1, j)
                        end if
                        if(j-1 >= 1) then
                            call index_convert_local_global(myid, 5, i, j-1, equInd)
                            call setMatValue(fieldInd, equInd, resi(i,j-1), 51, i, j-1)
                        end if
                        if(i+1 <= localncols) then
                            call index_convert_local_global(myid, 5, i+1, j, equInd)
                            call setMatValue(fieldInd, equInd, resi(i+1,j), 51, i+1, j)
                        end if
                        if(j+1 <= localnrows) then
                            call index_convert_local_global(myid, 5, i, j+1, equInd)
                            call setMatValue(fieldInd, equInd, resi(i,j+1), 51, i, j+1)
                        end if
                    end if

                end if

            end do
        end do

    end subroutine constructASw

    subroutine constructATem(tempe, resi)

        integer, dimension(:,:), pointer, intent(in) :: tempe
        real(kind=8), dimension(:,:), pointer, intent(in) :: resi

        integer :: fieldInd, equInd
        integer :: indexl, indexr, indexd, indexu
        integer :: i, j

        ! the field index
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
        if(prow /= pnrows) then
            indexu = localnrows + 1
        else
            indexu = localnrows
        end if

        do j = indexd, indexu
            do i = indexl, indexr

                if(tempe(i,j) == 1) then

                    if((i==0).and.(j/=0).and.(j/=localnrows+1)) then
                        call index_convert_local_global(myid-1, 6, localncols, j, fieldInd)
                        call index_convert_local_global(myid, 6, 1, j, equInd)
                        call setMatValue(fieldInd, equInd, resi(i+1,j), 61, 1, j)
                    elseif((i==localncols+1).and.(j/=0).and.(j/=localnrows+1)) then
                        call index_convert_local_global(myid+1, 6, 1, j, fieldInd)
                        call index_convert_local_global(myid, 6, localncols, j, equInd)
                        call setMatValue(fieldInd, equInd, resi(i-1,j), 61, localncols, j)
                    elseif((j==0).and.(i/=0).and.(i/=localncols+1)) then
                        call index_convert_local_global(myid-pncols, 6, i, localnrows, fieldInd)
                        call index_convert_local_global(myid, 6, i, 1, equInd)
                        call setMatValue(fieldInd, equInd, resi(i,j+1), 61, i, 1)
                    elseif((j==localnrows+1).and.(i/=0).and.(i/=localncols+1)) then
                        call index_convert_local_global(myid+pncols, 6, i, 1, fieldInd)
                        call index_convert_local_global(myid, 6, i, localnrows, equInd)
                        call setMatValue(fieldInd, equInd, resi(i,j-1), 61, i, localnrows)
                    elseif((i>=1).and.(i<=localncols).and.(j>=1).and.(j<=localnrows)) then
                        call index_convert_local_global(myid, 6, i, j, fieldInd)
                        call setMatValue(fieldInd, fieldInd, resi(i,j), 61, i, j)
                        if(i /= 1) then
                            call index_convert_local_global(myid, 6, i-1, j, equInd)
                            call setMatValue(fieldInd, equInd, resi(i-1,j), 61, i-1, j)
                        end if
                        if(i /= localncols) then
                            call index_convert_local_global(myid, 6, i+1, j, equInd)
                            call setMatValue(fieldInd, equInd, resi(i+1,j), 61, i+1, j)
                        end if
                        if(j /= 1) then
                            call index_convert_local_global(myid, 6, i, j-1, equInd)
                            call setMatValue(fieldInd, equInd, resi(i,j-1), 61, i, j-1)
                        end if
                        if(j /= localnrows) then
                            call index_convert_local_global(myid, 6, i, j+1, equInd)
                            call setMatValue(fieldInd, equInd, resi(i,j+1), 61, i, j+1)
                        end if
                    end if

                end if

            end do
        end do

    end subroutine constructATem

    subroutine setMatValue(col, row, value, m_kind, eq_i, eq_j)

        integer, intent(in) :: col ! global column index
        integer, intent(in) :: row ! global row index
        real(kind=8), intent(in) :: value
        integer, intent(in) :: m_kind ! matrix kind
        integer, intent(in) :: eq_i ! equation x-direction coordinate
        integer, intent(in) :: eq_j ! equation y-direction coordinate

        integer, dimension(:), pointer :: Acols
        integer, dimension(:), pointer :: Arows
        real(kind=8), dimension(:), pointer :: Avalues
        integer, dimension(:), pointer :: AEntryBase
        integer, dimension(:), pointer :: AEntryNum

        integer :: indexr, pos, base, tail, shend, m ,n
        integer :: left, right, mid

        ! matrix kind
        if(m_kind == 1111) then ! AxxwStatic
            Acols => AxxwCols
            Arows => AxxwRows
            Avalues => AxxwStaticValues
            AEntryBase => AxxEntryBase
            AEntryNum => AxxEntryNum
        elseif(m_kind == 1112) then ! AxxwDyn
            Acols => AxxwCols
            Arows => AxxwRows
            Avalues => AxxwDynValues
            AEntryBase => AxxEntryBase
            AEntryNum => AxxEntryNum
        elseif(m_kind == 1121) then ! AxxnStatic
            Acols => AxxnCols
            Arows => AxxnRows
            Avalues => AxxnStaticValues
            AEntryBase => AxxEntryBase
            AEntryNum => AxxEntryNum
        elseif(m_kind == 1122) then ! AxxnDyn
            Acols => AxxnCols
            Arows => AxxnRows
            Avalues => AxxnDynValues
            AEntryBase => AxxEntryBase
            AEntryNum => AxxEntryNum
        elseif(m_kind == 121) then ! Axpw
            Acols => AxpwCols
            Arows => AxpwRows
            Avalues => AxpwValues
            AEntryBase => AxpEntryBase
            AEntryNum => AxpEntryNum
        elseif(m_kind == 122) then ! Axpn
            Acols => AxpnCols
            Arows => AxpnRows
            Avalues => AxpnValues
            AEntryBase => AxpEntryBase
            AEntryNum => AxpEntryNum
        elseif(m_kind == 2111) then ! AyywStatic
            Acols => AyywCols
            Arows => AyywRows
            Avalues => AyywStaticValues
            AEntryBase => AyyEntryBase
            AEntryNum => AyyEntryNum
        elseif(m_kind == 2112) then ! AyywDyn
            Acols => AyywCols
            Arows => AyywRows
            Avalues => AyywDynValues
            AEntryBase => AyyEntryBase
            AEntryNum => AyyEntryNum
        elseif(m_kind == 2121) then ! AyynStatic
            Acols => AyynCols
            Arows => AyynRows
            Avalues => AyynStaticValues
            AEntryBase => AyyEntryBase
            AEntryNum => AyyEntryNum
        elseif(m_kind == 2122) then ! AyynDyn
            Acols => AyynCols
            Arows => AyynRows
            Avalues => AyynDynValues
            AEntryBase => AyyEntryBase
            AEntryNum => AyyEntryNum
        elseif(m_kind == 221) then ! Aypw
            Acols => AypwCols
            Arows => AypwRows
            Avalues => AypwValues
            AEntryBase => AypEntryBase
            AEntryNum => AypEntryNum
        elseif(m_kind == 222) then ! Aypn
            Acols => AypnCols
            Arows => AypnRows
            Avalues => AypnValues
            AEntryBase => AypEntryBase
            AEntryNum => AypEntryNum
        elseif(m_kind == 311) then ! Acxw
            Acols => AcxwCols
            Arows => AcxwRows
            Avalues => AcxwValues
            AEntryBase => AcxEntryBase
            AEntryNum => AcxEntryNum
        elseif(m_kind == 312) then ! Acxn
            Acols => AcxnCols
            Arows => AcxnRows
            Avalues => AcxnValues
            AEntryBase => AcxEntryBase
            AEntryNum => AcxEntryNum
        elseif(m_kind == 321) then ! Acyw
            Acols => AcywCols
            Arows => AcywRows
            Avalues => AcywValues
            AEntryBase => AcyEntryBase
            AEntryNum => AcyEntryNum
        elseif(m_kind == 322) then ! Acyn
            Acols => AcynCols
            Arows => AcynRows
            Avalues => AcynValues
            AEntryBase => AcyEntryBase
            AEntryNum => AcyEntryNum
        elseif(m_kind == 41) then ! ACf
            Acols => ACfCols
            Arows => ACfRows
            Avalues => ACfValues
            AEntryBase => ACfEntryBase
            AEntryNum => ACfEntryNum
        elseif(m_kind == 51) then ! ASw
            Acols => ASwCols
            Arows => ASwRows
            Avalues => ASwValues
            AEntryBase => ASwEntryBase
            AEntryNum => ASwEntryNum
        elseif(m_kind == 61) then ! ATem
            Acols => ATemCols
            Arows => ATemRows
            Avalues => ATemValues
            AEntryBase => ATemEntryBase
            AEntryNum => ATemEntryNum
        else
            print *, 'm_kind in setMatValue is wrong!'
            stop
        end if

        ! equation kind
        if((m_kind/1000==1).or.(m_kind/100==1)) then
            if(pcol /= pncols) then
                indexr = localncols
            else
                indexr = localncols + 1
            end if
            pos = indexr*(eq_j-1) + eq_i
        else
            pos = localncols*(eq_j-1) + eq_i
        end if

        base = AEntryBase(pos)
        tail = base + AEntryNum(pos) - 1

        if(t == 2) then

            do n = base, tail
                if(Acols(n) > col) then
                    do m = n+1, tail
                        if(Acols(m) == 0) then
                            shend = m - 1
                            exit
                        end if
                    end do
                    do m = shend, n, -1
                        Acols(m+1) = Acols(m)
                        Arows(m+1) = Arows(m)
                        Avalues(m+1) = Avalues(m)
                    end do
                    Acols(n) = col
                    Arows(n) = row
                    Avalues(n) = value
                    exit
                elseif(Acols(n) == col) then
                    Avalues(n) = value
                    exit
                elseif(Acols(n) == 0) then
                    Acols(n) = col
                    Arows(n) = row
                    Avalues(n) = value
                    exit
                end if
            end do

        else

            left = base
            right = tail
            mid = (left+right)/2
            do n = 1, AEntryNum(pos)
                if(Acols(mid) == col) then
                    Avalues(mid) = value
                    exit
                elseif(Acols(mid) < col) then
                    left = mid
                    mid = (left+right)/2
                    if((right-left) == 1) then
                        if(Acols(left) == col) then
                            Avalues(left) = value
                            exit
                        elseif(Acols(right) == col) then
                            Avalues(right) = value
                            exit
                        end if
                    end if
                elseif(Acols(mid) > col) then
                    right = mid
                    mid = (left+right)/2
                    if((right-left) == 1) then
                        if(Acols(left) == col) then
                            Avalues(left) = value
                            exit
                        elseif(Acols(right) == col) then
                            Avalues(right) = value
                            exit
                        end if
                    end if
                end if
            end do

        end if

    end subroutine setMatValue

    ! change the index of the unknowns from the local index to the global index
    subroutine index_convert_local_global(pid, lg_kind, local_i, local_j, global_ind)

        integer, intent(in) :: pid
        integer, intent(in) :: lg_kind
        integer, intent(in) :: local_i
        integer, intent(in) :: local_j
        integer, intent(out) :: global_ind

        integer :: p_pcol, p_prow
        integer :: base, ubase, vbase, CSTbase

        p_pcol = mod(pid,pncols)+1
        p_prow = pid/pncols+1

        base = (p_prow-1)*((nx+1)*localnrows*2+nx*localnrows*2+nx*localnrows)
        base = base + (p_pcol-1)*(localncols*localnrows*5)
        if(p_prow == pnrows) then
            base = base + (p_pcol-1)*localncols*2
        end if

        if(p_pcol == pncols) then
            ubase = localnrows*(localncols+1)*2
        else
            ubase = localnrows*localncols*2
        end if

        if(p_prow == pnrows) then
            vbase = (localnrows+1)*localncols*2
        else
            vbase = localnrows*localncols*2
        end if

        CSTbase = (p_prow-1)*nx*localnrows + (p_pcol-1)*localncols*localnrows

        ! the 'vx' index
        if(lg_kind == 11) then ! vxw
            if(p_pcol == pncols) then
                global_ind = base + (local_j-1)*(localncols+1) + local_i
            else
                global_ind = base + (local_j-1)*localncols + local_i
            end if

        elseif(lg_kind == 12) then ! vxn
            if(p_pcol == pncols) then
                global_ind = base + (localncols+1)*localnrows + (local_j-1)*(localncols+1) + local_i
            else
                global_ind = base + localncols*localnrows + (local_j-1)*localncols + local_i
            end if

        ! the 'vy' index
        elseif(lg_kind == 21) then ! vyw
            global_ind = base + ubase + (local_j-1)*localncols + local_i

        elseif(lg_kind == 22) then ! vyn
            if(p_prow == pnrows) then
                global_ind = base + ubase + localncols*(localnrows+1) + (local_j-1)*localncols + local_i
            else
                global_ind = base + ubase + localncols*localnrows + (local_j-1)*localncols + local_i
            end if

        ! the 'p' index
        elseif(lg_kind == 3) then ! p
            global_ind = base + ubase + vbase + (local_j-1)*localncols+local_i

        ! the 'Cf' index
        elseif(lg_kind == 4) then
            global_ind = CSTbase + (local_j-1)*localncols + local_i

        ! the 'Sw' index
        elseif(lg_kind == 5) then
            global_ind = CSTbase + (local_j-1)*localncols + local_i

        ! the 'Tem' index
        elseif(lg_kind == 6) then
            global_ind = CSTbase + (local_j-1)*localncols + local_i

        else
            print *, 'lg_kind in index_convert_local_global is wrong!', lg_kind
            stop
        end if

    end subroutine index_convert_local_global

    subroutine genExpField(bx, by, local_nx, local_ny, field, isField)

        integer, intent(in) :: bx, by, local_nx, local_ny
        integer, dimension(:,:), pointer, intent(in out) :: field
        logical, intent(out) :: isField

        field(:,:) = 0

        if((bx>local_nx).or.(by>local_ny)) then
            isField = .false.
        else
            field(bx:local_nx:3, by:local_ny:3) = 1
            isField = .true.
        end if

    end subroutine genExpField

end module DBF_constructMat
