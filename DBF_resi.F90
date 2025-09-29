
!!$ Author:
!!$   Yuanqing Wu, DGUT, P.R.China
!!$
!!$ History:
!!$   2025-5-9 by Yuanqing Wu
!!$
!!$ Support:
!!$   wuyuanq@gmail.com

module DBF_resi

    use DBF_globalData
    implicit none

Contains

    ! Calc Resi x-momentum
    ! The derivatives are evaluated on edges

    subroutine Resi_xmom_velx_b(velx, resi, phase)

        integer, dimension(:,:), pointer, intent(in) :: velx
        real(kind=8), dimension(:,:), pointer, intent(in out) :: resi
        integer, intent(in) :: phase

        real(kind=8), dimension(:), pointer :: vxBdryX0, vxBdryX1
        integer :: j

        if(phase == 1) then
            vxBdryX0 => vxwBdryX0
            vxBdryX1 => vxwBdryX1
        elseif(phase == 2) then
            vxBdryX0 => vxnBdryX0
            vxBdryX1 => vxnBdryX1
        else
            print *, 'The phase value in Resi_xmom_velx_b is wrong!'
            stop
        end if

        resi(:,:) = 0.D0

        ! Applying Neumann BC
        do j = 1, localnrows
            if((pcol==1).and.(isDiriX0_p(ylower+j-1)==0)) then
                resi(1,j) = velx(1,j) - vxBdryX0(ylower+j-1)
            end if
            if((pcol==pncols).and.(isDiriX1_p(ylower+j-1)==0)) then
                resi(localncols+1,j) = velx(localncols+1,j) - vxBdryX1(ylower+j-1)
            end if
        end do

    end subroutine Resi_xmom_velx_b

    subroutine Resi_xmom_velx(velx, resi, phase)

        integer, dimension(:,:), pointer, intent(in) :: velx
        real(kind=8), dimension(:,:), pointer, intent(in out) :: resi
        integer, intent(in) :: phase
        
        real(kind=8) :: velxLaplace, Forchh
        real(kind=8) :: vyAv, vAbs, DVxDx, DVxDy
        real(kind=8) :: poro_d, poro_u, S_d, S_u

        real(kind=8), dimension(:,:), pointer :: vx, vy
        real(kind=8), dimension(:), pointer :: vxBdryX0, vxBdryX1, vxBdryY0, vxBdryY1, &!
            vyBdryX0, vyBdryX1, vyBdryY0, vyBdryY1
        real(kind=8), dimension(:,:), pointer :: S, SEdgeX, SEdgeX_old
        real(kind=8), dimension(:), pointer :: SBdryX0, SBdryX1, SBdryY0, SBdryY1
        real(kind=8), dimension(:,:), pointer :: temp, temp_old, temp_Sw
        real(kind=8), dimension(:), pointer :: temp_SBdryX0, temp_SBdryX1, temp_SBdryY0, temp_SBdryY1
        real(kind=8), pointer :: rhof, visc

        integer :: indexl, indexr, indexd, indexu
        integer :: i, j, c

        if(phase == 1) then
            vx => vxw
            vy => vyw
            vxBdryX0 => vxwBdryX0
            vxBdryX1 => vxwBdryX1
            vxBdryY0 => vxwBdryY0
            vxBdryY1 => vxwBdryY1
            vyBdryX0 => vywBdryX0
            vyBdryX1 => vywBdryX1
            vyBdryY0 => vywBdryY0
            vyBdryY1 => vywBdryY1
            S => Sw
            SEdgeX => SwEdgeX
            SEdgeX_old => SwEdgeX_old
            SBdryX0 => SwBdryX0
            SBdryX1 => SwBdryX1
            SBdryY0 => SwBdryY0
            SBdryY1 => SwBdryY1
            rhof => rhofw
            visc => viscw
        elseif(phase == 2) then
            vx => vxn
            vy => vyn
            vxBdryX0 => vxnBdryX0
            vxBdryX1 => vxnBdryX1
            vxBdryY0 => vxnBdryY0
            vxBdryY1 => vxnBdryY1
            vyBdryX0 => vynBdryX0
            vyBdryX1 => vynBdryX1
            vyBdryY0 => vynBdryY0
            vyBdryY1 => vynBdryY1
            allocate(temp_Sw(lbound(Sw,dim=1):ubound(Sw,dim=1),lbound(Sw,dim=2):ubound(Sw,dim=2)))
            allocate(temp(lbound(SwEdgeX,dim=1):ubound(SwEdgeX,dim=1), &!
                lbound(SwEdgeX,dim=2):ubound(SwEdgeX,dim=2)))
            allocate(temp_old(lbound(SwEdgeX_old,dim=1):ubound(SwEdgeX_old,dim=1), &!
                lbound(SwEdgeX_old,dim=2):ubound(SwEdgeX_old,dim=2)))
            temp_Sw = 1.D0 - Sw
            temp = 1.D0 - SwEdgeX
            temp_old = 1.D0 - SwEdgeX_old
            S => temp_Sw
            SEdgeX => temp
            SEdgeX_old => temp_old
            allocate(temp_SBdryX0(size(SwBdryX0, dim=1)))
            allocate(temp_SBdryX1(size(SwBdryX1, dim=1)))
            allocate(temp_SBdryY0(size(SwBdryY0, dim=1)))
            allocate(temp_SBdryY1(size(SwBdryY1, dim=1)))
            temp_SBdryX0 = 1.D0 - SwBdryX0
            temp_SBdryX1 = 1.D0 - SwBdryX1
            temp_SBdryY0 = 1.D0 - SwBdryY0
            temp_SBdryY1 = 1.D0 - SwBdryY1
            SBdryX0 => temp_SBdryX0
            SBdryX1 => temp_SBdryX1
            SBdryY0 => temp_SBdryY0
            SBdryY1 => temp_SBdryY1
            rhof => rhofn
            visc => viscn
        else
            print *, 'The phase value in Resi_xmom_velx is wrong!'
            stop
        end if

        indexl = 1
        if(pcol /= pncols) then
            indexr = localncols
        else
            indexr = localncols + 1
        end if
        indexd = 1
        indexu = localnrows

        do j = indexd, indexu
            do i = indexl, indexr

                ! Determine the advection term
                if((pcol==1).and.(i==1)) then
                    vyAv = (vyBdryX0(ylower+j-1)+vyBdryX0(ylower+j)+vy(i,j)+vy(i,j+1)) / 4.D0
                elseif((pcol==pncols).and.(i == localncols+1)) then
                    vyAv = (vy(i-1,j)+vy(i-1,j+1)+vyBdryX1(ylower+j-1)+vyBdryX1(ylower+j)) / 4.D0
                else
                    vyAv = (vy(i-1,j)+vy(i-1,j+1)+vy(i,j)+vy(i,j+1)) / 4.D0
                end if

                if((pcol==1).and.(i==1).and.(isDiriX0_p(ylower+j-1)==0)) then
                    DVxDx = (velx(i+1,j)*poroEdgeX(i+1,j)*SEdgeX(i+1,j) - &!
                        vxBdryX0(ylower+j-1)*poroEdgeX(i,j)*SEdgeX(i,j))/hx(1)
                elseif((pcol==1).and.(i==1).and.(isDiriX0_p(ylower+j-1)==1)) then
                    DVxDx = (velx(i+1,j)*poroEdgeX(i+1,j)*SEdgeX(i+1,j) - &!
                        velx(i,j)*poroEdgeX(i,j)*SEdgeX(i,j))/hx(1)
                elseif((pcol==pncols).and.(i==localncols+1).and.(isDiriX1_p(ylower+j-1)==0)) then
                    DVxDx = (vxBdryX1(ylower+j-1)*poroEdgeX(i,j)*SEdgeX(i,j) - &!
                        velx(i-1,j)*poroEdgeX(i-1,j)*SEdgeX(i-1,j))/hx(localncols)
                elseif((pcol==pncols).and.(i==localncols+1).and.(isDiriX1_p(ylower+j-1)==1)) then
                    DVxDx = (velx(i,j)*poroEdgeX(i,j)*SEdgeX(i,j) - &!
                        velx(i-1,j)*poroEdgeX(i-1,j)*SEdgeX(i-1,j))/hx(localncols)
                else
                    if(vx(i,j)>0.D0) then
                        DVxDx = (velx(i,j)*poroEdgeX(i,j)*SEdgeX(i,j) - &!
                            velx(i-1,j)*poroEdgeX(i-1,j)*SEdgeX(i-1,j))/hx(i-1)
                    else
                        DVxDx = (velx(i+1,j)*poroEdgeX(i+1,j)*SEdgeX(i+1,j) - &!
                            velx(i,j)*poroEdgeX(i,j)*SEdgeX(i,j))/hx(i)
                    end if
                end if

                if((prow==1).and.(j==1)) then
                    DVxDy = (velx(i,j)*poroEdgeX(i,j)*SEdgeX(i,j) - &!
                        vxBdryY0(xlower+i-1)*poroEdgeX(i,j)*SEdgeX(i,j))/hy(j)
                elseif((prow==pnrows).and.(j==localnrows)) then
                    DVxDy = (vxBdryY1(xlower+i-1)*poroEdgeX(i,j)*SEdgeX(i,j) - &!
                        velx(i,j)*poroEdgeX(i,j)*SEdgeX(i,j))/hy(j)
                else
                    if(vyAv>0.D0) then
                        DVxDy = (velx(i,j)*poroEdgeX(i,j)*SEdgeX(i,j) - &!
                            velx(i,j-1)*poroEdgeX(i,j-1)*SEdgeX(i,j-1))/((hy(j-1)+hy(j))/2.D0)
                    else
                        DVxDy = (velx(i,j+1)*poroEdgeX(i,j+1)*SEdgeX(i,j+1) - &!
                            velx(i,j)*poroEdgeX(i,j)*SEdgeX(i,j))/((hy(j)+hy(j+1))/2.D0)
                    end if
                end if

                resi(i,j) = - rhof * (vx(i,j)*DVxDx+vyAv*DVxDy)

                ! Determine the Darcian term
                if(isDarcy) then
                    if((.not.((pcol==1).and.(i==1).and.(isDiriX0_p(ylower+j-1)==0))).and. &!
                        (.not.((pcol==pncols).and.(i==localncols+1).and.(isDiriX1_p(ylower+j-1)==0)))) then
                        resi(i,j) = resi(i,j) - &!
                            visc*poroEdgeX(i,j)*SEdgeX(i,j)/(SEdgeX(i,j)**3.D0*KxxEdge(i,j))*velx(i,j)
                    end if
                end if

                ! Determine the Brinkman term
                if(isBrinkman) then
                    if((pcol==1).and.(i==1).and.(isDiriX0_p(ylower+j-1)==0)) then
                        velxLaplace = (velx(i+1,j)-vxBdryX0(ylower+j-1))/(hx(1)**2.D0)*poro(i,j)*S(i,j)
                    elseif((pcol==1).and.(i==1).and.(isDiriX0_p(ylower+j-1)==1)) then
                        velxLaplace = (velx(i+1,j)-velx(i,j))/(hx(1)**2.D0)*poro(i,j)*S(i,j)
                    elseif((pcol==pncols).and.(i == localncols+1).and.(isDiriX1_p(ylower+j-1)==0)) then
                        velxLaplace = (vxBdryX1(ylower+j-1)-velx(i-1,j))/(hx(localncols)**2.D0)*poro(i-1,j)*S(i-1,j)
                    elseif((pcol==pncols).and.(i == localncols+1).and.(isDiriX1_p(ylower+j-1)==1)) then
                        velxLaplace = (velx(i,j)-velx(i-1,j))/(hx(localncols)**2.D0)*poro(i-1,j)*S(i-1,j)
                    else
                        velxLaplace = ((velx(i+1,j)-velx(i,j))/hx(i)*poro(i,j)*S(i,j) - &!
                            (velx(i,j)-velx(i-1,j))/hx(i-1)*poro(i-1,j)*S(i-1,j)) / ((hx(i-1)+hx(i))/2.D0)
                    end if

                    poro_d = 0.D0
                    c = 0
                    if(.not.((pcol==1).and.(i==1))) then
                        poro_d = poro_d + poro(i-1,j)
                        c = c + 1
                        if(.not.((prow==1).and.(j==1))) then
                            poro_d = poro_d + poro(i-1,j-1)
                            c = c + 1
                        end if
                    end if
                    if(.not.((pcol==pncols).and.(i==localncols+1))) then
                        poro_d = poro_d + poro(i,j)
                        c = c + 1
                        if(.not.((prow==1).and.(j==1))) then
                            poro_d = poro_d + poro(i,j-1)
                            c = c + 1
                        end if
                    end if
                    poro_d = poro_d / c

                    poro_u = 0.D0
                    c = 0
                    if(.not.((pcol==1).and.(i==1))) then
                        poro_u = poro_u + poro(i-1,j)
                        c = c + 1
                        if(.not.((prow==pnrows).and.(j==localnrows))) then
                            poro_u = poro_u + poro(i-1,j+1)
                            c = c + 1
                        end if
                    end if
                    if(.not.((pcol==pncols).and.(i==localncols+1))) then
                        poro_u = poro_u + poro(i,j)
                        c = c + 1
                        if(.not.((prow==pnrows).and.(j==localnrows))) then
                            poro_u = poro_u + poro(i,j+1)
                            c = c + 1
                        end if
                    end if
                    poro_u = poro_u / c

                    S_d = 0.D0
                    c = 0
                    if((pcol==1).and.(i==1).and.(isDiriX0_Sw(ylower+j-1)==1)) then
                        S_d = S_d + SBdryX0(ylower+j-1)
                        c = c + 1
                    elseif(.not.((pcol==1).and.(i==1))) then
                        S_d = S_d + S(i-1,j)
                        c = c + 1
                    end if
                    if(.not.((pcol==1).and.(i==1).and.(prow==1).and.(j==1))) then
                        if((pcol==1).and.(i==1).and.(j>=2).and.(isDiriX0_Sw(ylower+j-2)==1)) then
                            S_d = S_d + SBdryX0(ylower+j-2)
                            c = c + 1
                        elseif((prow==1).and.(j==1).and.(i>=2).and.(isDiriY0_Sw(xlower+i-2)==1)) then
                            S_d = S_d + SBdryY0(xlower+i-2)
                            c = c + 1
                        elseif((i>=2).and.(j>=2)) then
                            S_d = S_d + S(i-1,j-1)
                            c = c + 1
                        end if
                    end if
                    if((pcol==pncols).and.(i==localncols+1).and.(isDiriX1_Sw(ylower+j-1)==1)) then
                        S_d = S_d + SBdryX1(ylower+j-1)
                        c = c + 1
                    elseif(.not.((pcol==pncols).and.(i==localncols+1))) then
                        S_d = S_d + S(i,j)
                        c = c + 1
                    end if
                    if(.not.((pcol==pncols).and.(i==localncols+1).and.(prow==1).and.(j==1))) then
                        if((pcol==pncols).and.(i==localncols+1).and.(j>=2).and.(isDiriX1_Sw(ylower+j-2)==1)) then
                            S_d = S_d + SBdryX1(ylower+j-2)
                            c = c + 1
                        elseif((prow==1).and.(j==1).and.(i<=localncols).and.(isDiriY0_Sw(xlower+i-1)==1)) then
                            S_d = S_d + SBdryY0(xlower+i-1)
                            c = c + 1
                        elseif((i<=localncols).and.(j>=2)) then
                            S_d = S_d + S(i,j-1)
                            c = c + 1
                        end if
                    end if
                    S_d = S_d / c

                    S_u = 0.D0
                    c = 0
                    if((pcol==1).and.(i==1).and.(isDiriX0_Sw(ylower+j-1)==1)) then
                        S_u = S_u + SBdryX0(ylower+j-1)
                        c = c + 1
                    elseif(.not.((pcol==1).and.(i==1))) then
                        S_u = S_u + S(i-1,j)
                        c = c + 1
                    end if
                    if(.not.((pcol==1).and.(i==1).and.(prow==pnrows).and.(j==localnrows))) then
                        if((pcol==1).and.(i==1).and.(j<=localnrows-1).and.(isDiriX0_Sw(ylower+j)==1)) then
                            S_u = S_u + SBdryX0(ylower+j)
                            c = c + 1
                        elseif((prow==pnrows).and.(j==localnrows).and.(i>=2).and.(isDiriY1_Sw(xlower+i-2)==1)) then
                            S_u = S_u + SBdryY1(xlower+i-2)
                            c = c + 1
                        elseif((i>=2).and.(j<=localnrows-1)) then
                            S_u = S_u + S(i-1,j+1)
                            c = c + 1
                        end if
                    end if
                    if((pcol==pncols).and.(i==localncols+1).and.(isDiriX1_Sw(ylower+j-1)==1)) then
                        S_u = S_u + SBdryX1(ylower+j-1)
                        c = c + 1
                    elseif(.not.((pcol==pncols).and.(i==localncols+1))) then
                        S_u = S_u + S(i,j)
                        c = c + 1
                    end if
                    if(.not.((pcol==pncols).and.(i==localncols+1).and.(prow==pnrows).and.(j==localnrows))) then
                        if((pcol==pncols).and.(i==localncols+1).and.(j<=localnrows-1).and.(isDiriX1_Sw(ylower+j)==1)) then
                            S_u = S_u + SBdryX1(ylower+j)
                            c = c + 1
                        elseif((prow==pnrows).and.(j==localnrows).and.(i<=localncols).and.(isDiriY1_Sw(xlower+i-1)==1)) then
                            S_u = S_u + SBdryY1(xlower+i-1)
                            c = c + 1
                        elseif((i<=localncols).and.(j<=localnrows-1)) then
                            S_u = S_u + S(i,j+1)
                            c = c + 1
                        end if
                    end if
                    S_u = S_u / c

                    if((prow==1).and.(j==1)) then
                        velxLaplace = velxLaplace + ((velx(i,j+1)-velx(i,j))/ &!
                            ((hy(j+1)+hy(j))/2.D0)*poro_u*S_u - &!
                            (velx(i,j)-vxBdryY0(xlower+i-1))/hy(1)*poro_d*S_d) / hy(j)
                    elseif((prow==pnrows).and.(j==localnrows)) then
                        velxLaplace = velxLaplace + ((vxBdryY1(xlower+i-1)-velx(i,j))/hy(localnrows)*poro_u*S_u - &!
                            (velx(i,j)-velx(i,j-1))/((hy(j)+hy(j-1))/2.D0)*poro_d*S_d) / hy(j)
                    else
                        velxLaplace = velxLaplace + ((velx(i,j+1)-velx(i,j))/((hy(j+1)+hy(j))/2.D0)*poro_u*S_u - &!
                            (velx(i,j)-velx(i,j-1))/((hy(j)+hy(j-1))/2.D0)*poro_d*S_d) / hy(j)
                    end if

                    resi(i,j) = resi(i,j) + visc*velxLaplace

                end if

                ! Determine the Forchhimer term
                if(isForchheimer) then
                    if((.not.((pcol==1).and.(i==1).and.(isDiriX0_p(ylower+j-1)==0))).and. &!
                        (.not.((pcol==pncols).and.(i==localncols+1).and.(isDiriX1_p(ylower+j-1)==0)))) then
                        vAbs = dsqrt(vx(i,j)**2.D0 + vyAv**2.D0)
                        Forchh = 1.75D0/dsqrt(1.5D2*(poroEdgeX(i,j)*SEdgeX(i,j))**3.D0)
                        resi(i,j) = resi(i,j) - poroEdgeX(i,j)*SEdgeX(i,j)*rhof*Forchh*vAbs*velx(i,j)
                    end if
                end if

                ! Determine the time term
                resi(i,j) = resi(i,j) - rhof*(velx(i,j)*poroEdgeX(i,j)*SEdgeX(i,j) - &!
                    vx(i,j)*poroEdgeX_old(i,j)*SEdgeX_old(i,j))/(timeEnd/nt)

            end do
        end do

        if(phase == 2) then
            deallocate(temp_Sw)
            deallocate(temp)
            deallocate(temp_old)
            deallocate(temp_SBdryX0)
            deallocate(temp_SBdryX1)
            deallocate(temp_SBdryY0)
            deallocate(temp_SBdryY1)
        end if

    end subroutine Resi_xmom_velx

    subroutine Resi_xmom_pres(pres, resi)

        integer, dimension(:,:), pointer, intent(in) :: pres
        real(kind=8), dimension(:,:), pointer, intent(in out) :: resi

        real(kind=8), dimension(:,:), allocatable :: pPad
        real(kind=8), dimension(:), allocatable :: xc
        integer :: indexl, indexr, indexd, indexu
        integer :: i, j

        ! the index for equation
        indexl = 1
        if(pcol /= pncols) then
            indexr = localncols
        else
            indexr = localncols + 1
        end if
        indexd = 1
        indexu = localnrows

        allocate(pPad(0:localncols+1,1:localnrows))
        allocate(xc(0:localncols+1))

        xc(1:localncols) = (xs(xlower:xlower+localncols-1)+xs(xlower+1:xlower+localncols))/2.D0
        if(pcol /= 1) then
            xc(0) = (xs(xlower-1)+xs(xlower))/2.D0
        else
            xc(0) = xs(xlower)
        end if
        if(pcol /= pncols) then
            xc(localncols+1) = (xs(xlower+localncols)+xs(xlower+localncols+1))/2.D0
        else
            xc(localncols+1) = xs(xlower+localncols)
        end if

        pPad(1:localncols,1:localnrows) = pres(1:localncols,1:localnrows)
        ! insert BC
        if(pcol /= 1) then
            pPad(0,1:localnrows) = pres(0,1:localnrows)
        else
            pPad(0,1:localnrows) = pBdryX0(ylower:ylower+localnrows-1) * isDiriX0_p(ylower:ylower+localnrows-1)
        end if
        if(pcol /= pncols) then
            pPad(localncols+1,1:localnrows) = pres(localncols+1,1:localnrows)
        else
            pPad(localncols+1,1:localnrows) = pBdryX1(ylower:ylower+localnrows-1) * isDiriX1_p(ylower:ylower+localnrows-1)
        end if

        do j = indexd, indexu
            do i = indexl, indexr

                ! Applying Neumann BC
                if(((pcol==1).and.(i==1).and.(isDiriX0_p(ylower+j-1)==0)).or.((pcol==pncols).and.(i==localncols+1) &!
                    .and.(isDiriX1_p(ylower+j-1)==0))) then
                    resi(i,j) = 0.D0
                else
                    resi(i,j) = -(pPad(i,j) - pPad(i-1,j))/(xc(i) - xc(i-1))
                end if

            end do
        end do

        deallocate(pPad)
        deallocate(xc)

    end subroutine Resi_xmom_pres

    ! End Calc Resi x-momentum

    ! Calc Resi y-momentum
    ! The derivatives are evaluated on edges

    subroutine Resi_ymom_vely_b(vely, resi, phase)

        integer, dimension(:,:), pointer, intent(in) :: vely
        real(kind=8), dimension(:,:), pointer, intent(in out) :: resi
        integer, intent(in) :: phase

        real(kind=8), dimension(:), pointer :: vyBdryY0, vyBdryY1
        integer :: i

        if(phase == 1) then
            vyBdryY0 => vywBdryY0
            vyBdryY1 => vywBdryY1
        elseif(phase == 2) then
            vyBdryY0 => vynBdryY0
            vyBdryY1 => vynBdryY1
        else
            print *, 'The phase value in Resi_ymom_vely_b is wrong!'
            stop
        end if

        resi(:,:) = 0.D0

        ! Applying Neumann BC
        do i = 1, localncols
            if((prow==1).and.(isDiriY0_p(xlower+i-1)==0)) then
                resi(i,1) = vely(i,1) - vyBdryY0(xlower+i-1)
            end if
            if((prow==pnrows).and.(isDiriY1_p(xlower+i-1)==0)) then
                resi(i,localnrows+1) = vely(i,localnrows+1) - vyBdryY1(xlower+i-1)
            end if
        end do
        
    end subroutine Resi_ymom_vely_b

    subroutine Resi_ymom_vely(vely, resi, phase)

        integer, dimension(:,:), pointer, intent(in) :: vely
        real(kind=8), dimension(:,:), pointer, intent(in out) :: resi
        integer, intent(in) :: phase
        
        real(kind=8) :: velyLaplace, Forchh
        real(kind=8) :: vxAv, vAbs, DVyDx, DVyDy
        real(kind=8) :: poro_l, poro_r, S_l, S_r

        real(kind=8), dimension(:,:), pointer :: vx, vy
        real(kind=8), dimension(:), pointer :: vxBdryX0, vxBdryX1, vxBdryY0, vxBdryY1, &!
            vyBdryX0, vyBdryX1, vyBdryY0, vyBdryY1
        real(kind=8), dimension(:,:), pointer :: S, SEdgeY, SEdgeY_old
        real(kind=8), dimension(:), pointer :: SBdryX0, SBdryX1, SBdryY0, SBdryY1
        real(kind=8), dimension(:,:), pointer :: temp, temp_old, temp_Sw
        real(kind=8), dimension(:), pointer :: temp_SBdryX0, temp_SBdryX1, temp_SBdryY0, temp_SBdryY1
        real(kind=8), pointer :: rhof, visc

        integer :: indexl, indexr, indexd, indexu
        integer :: i, j, c 

        if(phase == 1) then
            vx => vxw
            vy => vyw
            vxBdryX0 => vxwBdryX0
            vxBdryX1 => vxwBdryX1
            vxBdryY0 => vxwBdryY0
            vxBdryY1 => vxwBdryY1
            vyBdryX0 => vywBdryX0
            vyBdryX1 => vywBdryX1
            vyBdryY0 => vywBdryY0
            vyBdryY1 => vywBdryY1
            S => Sw
            SEdgeY => SwEdgeY
            SEdgeY_old => SwEdgeY_old
            SBdryX0 => SwBdryX0
            SBdryX1 => SwBdryX1
            SBdryY0 => SwBdryY0
            SBdryY1 => SwBdryY1
            rhof => rhofw
            visc => viscw
        elseif(phase == 2) then
            vx => vxn
            vy => vyn
            vxBdryX0 => vxnBdryX0
            vxBdryX1 => vxnBdryX1
            vxBdryY0 => vxnBdryY0
            vxBdryY1 => vxnBdryY1
            vyBdryX0 => vynBdryX0
            vyBdryX1 => vynBdryX1
            vyBdryY0 => vynBdryY0
            vyBdryY1 => vynBdryY1
            allocate(temp_Sw(lbound(Sw,dim=1):ubound(Sw,dim=1),lbound(Sw,dim=2):ubound(Sw,dim=2)))
            allocate(temp(lbound(SwEdgeY,dim=1):ubound(SwEdgeY,dim=1), &!
                lbound(SwEdgeY,dim=2):ubound(SwEdgeY,dim=2)))
            allocate(temp_old(lbound(SwEdgeY_old,dim=1):ubound(SwEdgeY_old,dim=1), &!
                lbound(SwEdgeY_old,dim=2):ubound(SwEdgeY_old,dim=2)))
            temp_Sw = 1.D0 - Sw
            temp = 1.D0 - SwEdgeY
            temp_old = 1.D0 - SwEdgeY_old
            S => temp_Sw
            SEdgeY => temp
            SEdgeY_old => temp_old
            allocate(temp_SBdryX0(size(SwBdryX0, dim=1)))
            allocate(temp_SBdryX1(size(SwBdryX1, dim=1)))
            allocate(temp_SBdryY0(size(SwBdryY0, dim=1)))
            allocate(temp_SBdryY1(size(SwBdryY1, dim=1)))
            temp_SBdryX0 = 1.D0 - SwBdryX0
            temp_SBdryX1 = 1.D0 - SwBdryX1
            temp_SBdryY0 = 1.D0 - SwBdryY0
            temp_SBdryY1 = 1.D0 - SwBdryY1
            SBdryX0 => temp_SBdryX0
            SBdryX1 => temp_SBdryX1
            SBdryY0 => temp_SBdryY0
            SBdryY1 => temp_SBdryY1
            rhof => rhofn
            visc => viscn
        else
            print *, 'The phase value in Resi_ymom_vely is wrong!'
            stop
        end if

        indexl = 1
        indexr = localncols
        indexd = 1
        if(prow /= pnrows) then
            indexu = localnrows
        else
            indexu = localnrows + 1
        end if

        do j = indexd, indexu
            do i = indexl, indexr

                ! Determine the advection term
                if((prow==1).and.(j==1)) then
                    vxAv = (vxBdryY0(xlower+i-1)+vxBdryY0(xlower+i)+vx(i,j)+vx(i+1,j)) / 4.D0
                elseif((prow==pnrows).and.(j==localnrows+1)) then
                    vxAv = (vx(i,j-1)+vx(i+1,j-1)+vxBdryY1(xlower+i-1)+vxBdryY1(xlower+i)) / 4.D0
                else
                    vxAv = (vx(i,j-1)+vx(i+1,j-1)+vx(i,j)+vx(i+1,j)) / 4.D0
                end if

                if((prow==1).and.(j==1).and.(isDiriY0_p(xlower+i-1)==0)) then
                    DVyDy = (vely(i,j+1)*poroEdgeY(i,j+1)*SEdgeY(i,j+1) - &!
                        vyBdryY0(xlower+i-1)*poroEdgeY(i,j)*SEdgeY(i,j))/hy(j)
                elseif((prow==1).and.(j==1).and.(isDiriY0_p(xlower+i-1)==1)) then
                    DVyDy = (vely(i,j+1)*poroEdgeY(i,j+1)*SEdgeY(i,j+1) - &!
                        vely(i,j)*poroEdgeY(i,j)*SEdgeY(i,j))/hy(j)
                elseif((prow==pnrows).and.(j==localnrows+1).and.(isDiriY1_p(xlower+i-1)==0)) then
                    DVyDy = (vyBdryY1(xlower+i-1)*poroEdgeY(i,j)*SEdgeY(i,j) - &!
                        vely(i,j-1)*poroEdgeY(i,j-1)*SEdgeY(i,j-1))/hy(j-1)
                elseif((prow==pnrows).and.(j==localnrows+1).and.(isDiriY1_p(xlower+i-1)==1)) then
                    DVyDy = (vely(i,j)*poroEdgeY(i,j)*SEdgeY(i,j) - &!
                        vely(i,j-1)*poroEdgeY(i,j-1)*SEdgeY(i,j-1))/hy(j-1)
                else
                    if(vy(i,j)>0.D0) then
                        DVyDy = (vely(i,j)*poroEdgeY(i,j)*SEdgeY(i,j) - &!
                            vely(i,j-1)*poroEdgeY(i,j-1)*SEdgeY(i,j-1))/hy(j-1)
                    else
                        DVyDy = (vely(i,j+1)*poroEdgeY(i,j+1)*SEdgeY(i,j+1) - &!
                            vely(i,j)*poroEdgeY(i,j)*SEdgeY(i,j))/hy(j)
                    end if
                end if
                if((pcol==1).and.(i==1)) then
                    DVyDx = (vely(i,j)*poroEdgeY(i,j)*SEdgeY(i,j) - &!
                        vyBdryX0(ylower+j-1)*poroEdgeY(i,j)*SEdgeY(i,j))/hx(i)
                elseif((pcol==pncols).and.(i==localncols)) then
                    DVyDx = (vyBdryX1(ylower+j-1)*poroEdgeY(i,j)*SEdgeY(i,j) - &!
                        vely(i,j)*poroEdgeY(i,j)*SEdgeY(i,j))/hx(i)
                else
                    if(vxAv>0.D0) then
                        DVyDx = (vely(i,j)*poroEdgeY(i,j)*SEdgeY(i,j) - &!
                            vely(i-1,j)*poroEdgeY(i-1,j)*SEdgeY(i-1,j))/((hx(i-1)+hx(i))/2.D0)
                    else
                        DVyDx = (vely(i+1,j)*poroEdgeY(i+1,j)*SEdgeY(i+1,j) - &!
                            vely(i,j)*poroEdgeY(i,j)*SEdgeY(i,j))/((hx(i)+hx(i+1))/2.D0)
                    end if
                end if

                resi(i,j) = - rhof * (vy(i,j)*DVyDy+vxAv*DVyDx)

                ! Determine the Darcian term
                if(isDarcy) then
                    if((.not.((prow==1).and.(j==1).and.(isDiriY0_p(xlower+i-1)==0))).and. &!
                        (.not.((prow==pnrows).and.(j==localnrows+1).and.(isDiriY1_p(xlower+i-1)==0)))) then
                        resi(i,j) = resi(i,j) - &!
                            visc*poroEdgeY(i,j)*SEdgeY(i,j)/(SEdgeY(i,j)**3.D0*KyyEdge(i,j))*vely(i,j)
                    end if
                end if

                ! Determine the Brinkman term
                if(isBrinkman) then
                    if((prow==1).and.(j==1).and.(isDiriY0_p(xlower+i-1)==0)) then
                        velyLaplace = (vely(i,j+1)-vyBdryY0(xlower+i-1))/(hy(j)**2.D0)*poro(i,j)*S(i,j)
                    elseif((prow==1).and.(j==1).and.(isDiriY0_p(xlower+i-1)==1)) then
                        velyLaplace = (vely(i,j+1)-vely(i,j))/(hy(j)**2.D0)*poro(i,j)*S(i,j)
                    elseif((prow==pnrows).and.(j==localnrows+1).and.(isDiriY1_p(xlower+i-1)==0)) then
                        velyLaplace = (vyBdryY1(xlower+i-1)-vely(i,j-1))/(hy(j-1)**2.D0)*poro(i,j-1)*S(i,j-1)
                    elseif((prow==pnrows).and.(j==localnrows+1).and.(isDiriY1_p(xlower+i-1)==1)) then
                        velyLaplace = (vely(i,j)-vely(i,j-1))/(hy(j-1)**2.D0)*poro(i,j-1)*S(i,j-1)
                    else
                        velyLaplace = ((vely(i,j+1)-vely(i,j))/hy(j)*poro(i,j)*S(i,j) - &!
                            (vely(i,j)-vely(i,j-1))/hy(j-1)*poro(i,j-1)*S(i,j-1)) / ((hy(j-1)+hy(j))/2.D0)
                    end if

                    poro_l = 0.D0
                    c = 0
                    if(.not.((prow==1).and.(j==1))) then
                        poro_l = poro_l + poro(i,j-1)
                        c = c + 1
                        if(.not.((pcol==1).and.(i==1))) then
                            poro_l = poro_l + poro(i-1,j-1)
                            c = c + 1
                        end if
                    end if
                    if(.not.((prow==pnrows).and.(j==localnrows+1))) then
                        poro_l = poro_l + poro(i,j)
                        c = c + 1
                        if(.not.((pcol==1).and.(i==1))) then
                            poro_l = poro_l + poro(i-1,j)
                            c = c + 1
                        end if
                    end if
                    poro_l = poro_l / c

                    poro_r = 0.D0
                    c = 0
                    if(.not.((prow==1).and.(j==1))) then
                        poro_r = poro_r + poro(i,j-1)
                        c = c + 1
                        if(.not.((pcol==pncols).and.(i==localncols))) then
                            poro_r = poro_r + poro(i+1,j-1)
                            c = c + 1
                        end if
                    end if
                    if(.not.((prow==pnrows).and.(j==localnrows+1))) then
                        poro_r = poro_r + poro(i,j)
                        c = c + 1
                        if(.not.((pcol==pncols).and.(i==localncols))) then
                            poro_r = poro_r + poro(i+1,j)
                            c = c + 1
                        end if
                    end if
                    poro_r = poro_r / c

                    S_l = 0.D0
                    c = 0
                    if((prow==1).and.(j==1).and.(isDiriY0_Sw(xlower+i-1)==1)) then
                        S_l = S_l + SBdryY0(xlower+i-1)
                        c = c + 1
                    elseif(.not.((prow==1).and.(j==1))) then
                        S_l = S_l + S(i,j-1)
                        c = c + 1
                    end if
                    if(.not.((prow==1).and.(j==1).and.(pcol==1).and.(i==1))) then
                        if((prow==1).and.(j==1).and.(i>=2).and.(isDiriY0_Sw(xlower+i-2)==1)) then
                            S_l = S_l + SBdryY0(xlower+i-2)
                            c = c + 1
                        elseif((pcol==1).and.(i==1).and.(j>=2).and.(isDiriX0_Sw(ylower+j-2)==1)) then
                            S_l = S_l + SBdryX0(ylower+j-2)
                            c = c + 1
                        elseif((i>=2).and.(j>=2)) then
                            S_l = S_l + S(i-1,j-1)
                            c = c + 1
                        end if
                    end if
                    if((prow==pnrows).and.(j==localnrows+1).and.(isDiriY1_Sw(xlower+i-1)==1)) then
                        S_l = S_l + SBdryY1(xlower+i-1)
                        c = c + 1
                    elseif(.not.((prow==pnrows).and.(j==localnrows+1))) then
                        S_l = S_l + S(i,j)
                        c = c + 1
                    end if
                    if(.not.((prow==pnrows).and.(j==localnrows+1).and.(pcol==1).and.(i==1))) then
                        if((prow==pnrows).and.(j==localnrows+1).and.(i>=2).and.(isDiriY1_Sw(xlower+i-2)==1)) then
                            S_l = S_l + SBdryY1(xlower+i-2)
                            c = c + 1
                        elseif((pcol==1).and.(i==1).and.(j<=localnrows).and.(isDiriX0_Sw(ylower+j-1)==1)) then
                            S_l = S_l + SBdryX0(ylower+j-1)
                            c = c + 1
                        elseif((i>=2).and.(j<=localnrows)) then
                            S_l = S_l + S(i-1,j)
                            c = c + 1
                        end if
                    end if
                    S_l = S_l / c

                    S_r = 0.D0
                    c = 0
                    if((prow==1).and.(j==1).and.(isDiriY0_Sw(xlower+i-1)==1)) then
                        S_r = S_r + SBdryY0(xlower+i-1)
                        c = c + 1
                    elseif(.not.((prow==1).and.(j==1))) then
                        S_r = S_r + S(i,j-1)
                        c = c + 1
                    end if
                    if(.not.((prow==1).and.(j==1).and.(pcol==pncols).and.(i==localncols))) then
                        if((prow==1).and.(j==1).and.(i<=localncols-1).and.(isDiriY0_Sw(xlower+i)==1)) then
                            S_r = S_r + SBdryY0(xlower+i)
                            c = c + 1
                        elseif((pcol==pncols).and.(i==localncols).and.(j>=2).and.(isDiriX1_Sw(ylower+j-2)==1)) then
                            S_r = S_r + SBdryX1(ylower+j-2)
                            c = c + 1
                        elseif((i<=localncols-1).and.(j>=2)) then
                            S_r = S_r + S(i+1,j-1)
                            c = c + 1
                        end if
                    end if
                    if((prow==pnrows).and.(j==localnrows+1).and.(isDiriY1_Sw(xlower+i-1)==1)) then
                        S_r = S_r + SBdryY1(xlower+i-1)
                        c = c + 1
                    elseif(.not.((prow==pnrows).and.(j==localnrows+1))) then
                        S_r = S_r + S(i,j)
                        c = c + 1
                    end if
                    if(.not.((prow==pnrows).and.(j==localnrows+1).and.(pcol==pncols).and.(i==localncols))) then
                        if((prow==pnrows).and.(j==localnrows+1).and.(i<=localncols-1).and.(isDiriY1_Sw(xlower+i)==1)) then
                            S_r = S_r + SBdryY1(xlower+i)
                            c = c + 1
                        elseif((pcol==pncols).and.(i==localncols).and.(j<=localnrows).and.(isDiriX1_Sw(ylower+j-1)==1)) then
                            S_r = S_r + SBdryX1(ylower+j-1)
                            c = c + 1
                        elseif((i<=localncols-1).and.(j<=localnrows)) then
                            S_r = S_r + S(i+1,j)
                            c = c + 1
                        end if
                    end if
                    S_r = S_r / c

                    if((pcol==1).and.(i==1)) then
                        velyLaplace = velyLaplace + ((vely(i+1,j)-vely(i,j))/((hx(i+1)+hx(i))/2.D0)*poro_r*S_r - &!
                            (vely(i,j)-vyBdryX0(ylower+j-1))/hx(i)*poro_l*S_l) / hx(i)
                    elseif((pcol==pncols).and.(i==localncols)) then
                        velyLaplace = velyLaplace + ((vyBdryX1(ylower+j-1)-vely(i,j))/hx(i)*poro_r*S_r - &!
                            (vely(i,j)-vely(i-1,j))/((hx(i)+hx(i-1))/2.D0)*poro_l*S_l) / hx(i)
                    else
                        velyLaplace = velyLaplace + ((vely(i+1,j)-vely(i,j))/((hx(i+1)+hx(i))/2.D0)*poro_r*S_r - &!
                            (vely(i,j)-vely(i-1,j))/((hx(i)+hx(i-1))/2.D0)*poro_l*S_l) / hx(i)
                    end if

                    resi(i,j) = resi(i,j) + visc*velyLaplace
                end if

                ! Determine the Forchheimer term
                if(isForchheimer) then
                    if((.not.((prow==1).and.(j==1).and.(isDiriY0_p(xlower+i-1)==0))).and. &!
                        (.not.((prow==pnrows).and.(j==localnrows+1).and.(isDiriY1_p(xlower+i-1)==0)))) then
                        vAbs = dsqrt(vxAv**2.D0 + vy(i,j)**2.D0)
                        Forchh = 1.75D0/dsqrt(1.5D2*(poroEdgeY(i,j)*SEdgeY(i,j))**3.D0)
                        resi(i,j) = resi(i,j) -poroEdgeY(i,j)*SEdgeY(i,j)*rhof*Forchh*vAbs*vely(i,j)
                    end if
                end if

                ! Determine the time term
                resi(i,j) = resi(i,j) - rhof*(vely(i,j)*poroEdgeY(i,j)*SEdgeY(i,j) - &!
                    vy(i,j)*poroEdgeY_old(i,j)*SEdgeY_old(i,j))/(timeEnd/nt)

            end do
        end do

        if(phase == 2) then
            deallocate(temp_Sw)
            deallocate(temp)
            deallocate(temp_old)
            deallocate(temp_SBdryX0)
            deallocate(temp_SBdryX1)
            deallocate(temp_SBdryY0)
            deallocate(temp_SBdryY1)
        end if

    end subroutine Resi_ymom_vely

    subroutine Resi_ymom_pres(pres, resi)

        integer, dimension(:,:), pointer, intent(in) :: pres
        real(kind=8), dimension(:,:), pointer, intent(in out) :: resi
        
        real(kind=8), dimension(:,:), allocatable :: pPad
        real(kind=8), dimension(:), allocatable :: yc
        integer :: indexl, indexr, indexd, indexu
        integer :: i, j
 
        indexl = 1
        indexr = localncols
        indexd = 1
        if(prow /= pnrows) then
            indexu = localnrows 
        else
            indexu = localnrows + 1
        end if

        allocate(pPad(1:localncols,0:localnrows+1))
        allocate(yc(0:localnrows+1))

        yc(1:localnrows) = (ys(ylower:ylower+localnrows-1)+ys(ylower+1:ylower+localnrows))/2.D0
        if(prow /= 1) then
            yc(0) = (ys(ylower-1)+ys(ylower))/2.D0
        else
            yc(0) = ys(1)
        end if
        if(prow /= pnrows) then
            yc(localnrows+1) = (ys(ylower+localnrows)+ys(ylower+localnrows+1))/2.D0
        else
            yc(localnrows+1) = ys(ylower+localnrows)
        end if

        pPad(1:localncols,1:localnrows) = pres(1:localncols,1:localnrows)
        ! insert BC
        if(prow /= 1) then
            pPad(1:localncols,0) = pres(1:localncols,0)
        else
            pPad(1:localncols,0) = pBdryY0(xlower:xlower+localncols-1) * isDiriY0_p(xlower:xlower+localncols-1)
        end if
        if(prow /= pnrows) then
            pPad(1:localncols,localnrows+1) = pres(1:localncols,localnrows+1)
        else
            pPad(1:localncols,localnrows+1) = pBdryY1(xlower:xlower+localncols-1) * isDiriY1_p(xlower:xlower+localncols-1)
        end if

        do j = indexd, indexu
            do i = indexl, indexr

                ! Applying Neumann BC
                if(((prow==1).and.(j==1).and.(isDiriY0_p(xlower+i-1)==0)).or.((prow==pnrows).and.(j==localnrows+1) &!
                    .and.(isDiriY1_p(xlower+i-1)==0))) then
                    resi(i,j) = 0
                else
                    resi(i,j) = -(pPad(i,j) - pPad(i,j-1))/(yc(j) - yc(j-1))
                end if

            end do
        end do

        deallocate(pPad)
        deallocate(yc)

    end subroutine Resi_ymom_pres

    ! End Calc Resi y-momentum

    ! Calc Resi Continuity

    subroutine Resi_consum_velx(velx, resi, phase)

        integer, dimension(:,:), pointer, intent(in) :: velx
        real(kind=8), dimension(:,:), pointer, intent(in out) :: resi
        integer, intent(in) :: phase
        
        real(kind=8), dimension(:), pointer :: vxBdryX0, vxBdryX1
        real(kind=8), dimension(:,:), pointer :: SEdgeX, temp
        integer :: i, j

        if(phase == 1) then
            vxBdryX0 => vxwBdryX0
            vxBdryX1 => vxwBdryX1
            SEdgeX => SwEdgeX
        elseif(phase == 2) then
            vxBdryX0 => vxnBdryX0
            vxBdryX1 => vxnBdryX1
            allocate(temp(lbound(SwEdgeX,dim=1):ubound(SwEdgeX,dim=1), &!
                lbound(SwEdgeX,dim=2):ubound(SwEdgeX,dim=2)))
            temp = 1.D0 - SwEdgeX
            SEdgeX => temp
        else
            print *, 'The phase value in Resi_consum_velx is wrong!'
            stop
        end if

        do j = 1, localnrows
            do i = 1, localncols
                if((pcol==1).and.(i==1).and.(isDiriX0_p(ylower+j-1)==0)) then
                    resi(i,j) = (poroEdgeX(i+1,j)*SEdgeX(i+1,j)*velx(i+1,j) - &!
                        poroEdgeX(i,j)*SEdgeX(i,j)*vxBdryX0(ylower+j-1))/hx(i)
                elseif((pcol==pncols).and.(i==localncols).and.(isDiriX1_p(ylower+j-1)==0)) then
                    resi(i,j) = (poroEdgeX(i+1,j)*SEdgeX(i+1,j)*vxBdryX1(ylower+j-1) - &!
                        poroEdgeX(i,j)*SEdgeX(i,j)*velx(i,j))/hx(i)
                else
                    resi(i,j) = (poroEdgeX(i+1,j)*SEdgeX(i+1,j)*velx(i+1,j) - &!
                        poroEdgeX(i,j)*SEdgeX(i,j)*velx(i,j))/hx(i)
                end if
            end do
        end do

        if(phase == 2) then
            deallocate(temp)
        end if

    end subroutine Resi_consum_velx

    subroutine Resi_consum_vely(vely, resi, phase)

        integer, dimension(:,:), pointer, intent(in) :: vely
        real(kind=8), dimension(:,:), pointer, intent(in out) :: resi
        integer, intent(in) :: phase
    
        real(kind=8), dimension(:), pointer :: vyBdryY0, vyBdryY1
        real(kind=8), dimension(:,:), pointer :: SEdgeY, temp
        integer :: i, j

        if(phase == 1) then
            vyBdryY0 => vywBdryY0
            vyBdryY1 => vywBdryY1
            SEdgeY => SwEdgeY
        elseif(phase == 2) then
            vyBdryY0 => vynBdryY0
            vyBdryY1 => vynBdryY1
            allocate(temp(lbound(SwEdgeY,dim=1):ubound(SwEdgeY,dim=1), &!
                lbound(SwEdgeY,dim=2):ubound(SwEdgeY,dim=2)))
            temp = 1.D0 - SwEdgeY
            SEdgeY => temp
        else
            print *, 'The phase value in Resi_consum_vely is wrong!'
            stop
        end if
 
        do j = 1, localnrows
            do i = 1, localncols
                if((prow==1).and.(j==1).and.(isDiriY0_p(xlower+i-1)==0)) then
                    resi(i,j) = (poroEdgeY(i,j+1)*SEdgeY(i,j+1)*vely(i,j+1) - &!
                        poroEdgeY(i,j)*SEdgeY(i,j)*vyBdryY0(xlower+i-1))/hy(j)
                elseif((prow==pnrows).and.(j==localnrows).and.(isDiriY1_p(xlower+i-1)==0)) then
                    resi(i,j) = (poroEdgeY(i,j+1)*SEdgeY(i,j+1)*vyBdryY1(xlower+i-1) - &!
                        poroEdgeY(i,j)*SEdgeY(i,j)*vely(i,j))/hy(j)
                else
                    resi(i,j) = (poroEdgeY(i,j+1)*SEdgeY(i,j+1)*vely(i,j+1) - &!
                        poroEdgeY(i,j)*SEdgeY(i,j)*vely(i,j))/hy(j)
                end if
            end do
        end do

        if(phase == 2) then
            deallocate(temp)
        end if

    end subroutine Resi_consum_vely

    subroutine Resi_con_Sw(satw, resi)

        integer, dimension(:,:), pointer, intent(in) :: satw
        real(kind=8), dimension(:,:), pointer, intent(in out) :: resi

        real(kind=8) :: Sw_l, Sw_r, Sw_d, Sw_u
        integer :: i, j

        do j = 1, localnrows
            do i = 1, localncols

                resi(i,j) = (poro(i,j)*satw(i,j)-poro_old(i,j)*Sw(i,j))/(timeEnd/nt)

                if((pcol == 1).and.(i == 1).and.(isDiriX0_Sw(ylower+j-1)==0)) then
                    Sw_l = satw(i,j)
                elseif((pcol == 1).and.(i == 1).and.(isDiriX0_Sw(ylower+j-1)==1)) then
                    Sw_l = SwBdryX0(ylower+j-1)
                else
                    if(vxw(i,j) > 0.D0) then
                        Sw_l = satw(i-1,j)
                    else
                        Sw_l = satw(i,j)
                    end if
                end if

                if((pcol == pncols).and.(i == localncols).and.(isDiriX1_Sw(ylower+j-1)==0)) then
                    Sw_r = satw(i,j)
                elseif((pcol == pncols).and.(i == localncols).and.(isDiriX1_Sw(ylower+j-1)==1)) then
                    Sw_r = SwBdryX1(ylower+j-1)
                else
                    if(vxw(i+1,j) > 0.D0) then
                        Sw_r = satw(i,j)
                    else
                        Sw_r = satw(i+1,j)
                    end if
                end if

                if((prow == 1).and.(j == 1).and.(isDiriY0_Sw(xlower+i-1)==0)) then
                    Sw_d = satw(i,j)
                elseif((prow == 1).and.(j == 1).and.(isDiriY0_Sw(xlower+i-1)==1)) then
                    Sw_d = SwBdryY0(xlower+i-1)
                else
                    if(vyw(i,j) > 0.D0) then
                        Sw_d = satw(i,j-1)
                    else
                        Sw_d = satw(i,j)
                    end if
                end if

                if((prow == pnrows).and.(j == localnrows).and.(isDiriY1_Sw(xlower+i-1)==0)) then
                    Sw_u = satw(i,j)
                elseif((prow == pnrows).and.(j == localnrows).and.(isDiriY1_Sw(xlower+i-1)==1)) then
                    Sw_u = SwBdryY1(xlower+i-1)
                else
                    if(vyw(i,j+1) > 0.D0) then
                        Sw_u = satw(i,j)
                    else
                        Sw_u = satw(i,j+1)
                    end if
                end if

                resi(i,j) = resi(i,j) + (poroEdgeX(i+1,j)*vxw(i+1,j)*Sw_r - &!
                    poroEdgeX(i,j)*vxw(i,j)*Sw_l)/hx(i) + &!
                    (poroEdgeY(i,j+1)*vyw(i,j+1)*Sw_u - &!
                    poroEdgeY(i,j)*vyw(i,j)*Sw_d)/hy(j)

            end do
        end do

    end subroutine Resi_con_Sw

    ! End Calc Resi Continuity

    ! Calc Resi Concentration

    subroutine Resi_tran_Cf(conc, resi)

        integer, dimension(:,:), pointer, intent(in) :: conc
        real(kind=8), dimension(:,:), pointer, intent(in out) :: resi
        
        integer :: indexl, indexr, indexd, indexu
        real(kind=8) :: dl, dt
        real(kind=8), dimension(:,:,:), allocatable :: Ex, Ey, Dx, Dy
        real(kind=8), dimension(:,:), allocatable :: CfEdgeX, CfEdgeY
        real(kind=8) :: vxwAv, vywAv, vwmodulus
        real(kind=8) :: dmEdgeX, dmEdgeY
        real(kind=8), dimension(:,:), allocatable :: dCfdx, dCfdy
        real(kind=8) :: dCfdxleft, dCfdxright, dCfdxdown, dCfdxup
        real(kind=8) :: dCfdyleft, dCfdyright, dCfdydown, dCfdyup
        real(kind=8) :: div1, div2, reaction, timeTerm
        integer :: i, j

        indexl = 1
        indexr = localncols + 1
        indexd = 1
        indexu = localnrows

        allocate(Ex(indexl:indexr,indexd:indexu,2))
        allocate(Dx(indexl:indexr,indexd:indexu,2))
        allocate(CfEdgeX(indexl:indexr,indexd:indexu))

        do j = indexd, indexu
            do i = indexl, indexr

                if((pcol==1).and.(i==1)) then
                    vywAv = (vywBdryX0(ylower+j-1)+vywBdryX0(ylower+j)+vyw(i,j)+vyw(i,j+1)) / 4.D0
                elseif((pcol==pncols).and.(i == localncols+1)) then
                    vywAv = (vyw(i-1,j)+vyw(i-1,j+1)+vywBdryX1(ylower+j-1)+vywBdryX1(ylower+j)) / 4.D0
                else
                    vywAv = (vyw(i-1,j)+vyw(i-1,j+1)+vyw(i,j)+vyw(i,j+1)) / 4.D0
                end if

                vwmodulus = dsqrt((vxw(i,j)**2.D0+vywAv**2.D0)*1.D0)

                if(vwmodulus /= 0.D0) then
                    Ex(i,j,1) = vxw(i,j)**2.D0 / vwmodulus**2.D0
                    Ex(i,j,2) = vxw(i,j)*vywAv / vwmodulus**2.D0  ! the second column and the first row of the matrix
                else
                    Ex(i,j,1) = 0.D0
                    Ex(i,j,2) = 0.D0
                end if

                if((pcol==1).and.(i==1)) then
                    dmEdgeX = dm(i,j)
                elseif((pcol==pncols).and.(i==localncols+1)) then
                    dmEdgeX = dm(i-1,j)
                else
                    dmEdgeX = (hx(i-1)+hx(i)) / (hx(i-1)/dm(i-1,j)+hx(i)/dm(i,j))
                end if

                dl = alphaOS*dmEdgeX + 2.D0*lamdaX*vwmodulus*radiusInit*(1.D0-poroEdgeXInit(i,j))/ &!
                    (poroEdgeXInit(i,j)*(1.D0-poroEdgeX(i,j)))
                dt = alphaOS*dmEdgeX + 2.D0*lamdaT*vwmodulus*radiusInit*(1.D0-poroEdgeXInit(i,j))/ &!
                    (poroEdgeXInit(i,j)*(1.D0-poroEdgeX(i,j)))

                Dx(i,j,1) = alphaOS*dmEdgeX + vwmodulus/(poroEdgeX(i,j)*SwEdgeX(i,j))*((dl-dt)*Ex(i,j,1)+dt)
                Dx(i,j,2) = vwmodulus/(poroEdgeX(i,j)*SwEdgeX(i,j))*(dl-dt)*Ex(i,j,2)

            end do
        end do

        do j = indexd, indexu
            do i = indexl, indexr

                if((pcol==1).and.(i==1).and.(isDiriX0_Cf(ylower+j-1)==1)) then
                    CfEdgeX(i,j) = CfBdryX0(ylower+j-1)
                elseif((pcol==1).and.(i==1).and.(isDiriX0_Cf(ylower+j-1)==0)) then
                    CfEdgeX(i,j) = conc(i,j)
                elseif((pcol==pncols).and.(i==localncols+1).and.(isDiriX1_Cf(ylower+j-1)==1)) then
                    CfEdgeX(i,j) = CfBdryX1(ylower+j-1)
                elseif((pcol==pncols).and.(i==localncols+1).and.(isDiriX1_Cf(ylower+j-1)==0)) then
                    CfEdgeX(i,j) = conc(i-1,j)
                elseif(vxw(i,j) > 0.D0) then
                    CfEdgeX(i,j) = conc(i-1,j)
                else
                    CfEdgeX(i,j) = conc(i,j)
                end if

            end do
        end do

        indexl = 1
        indexr = localncols
        indexd = 1
        indexu = localnrows + 1

        allocate(Ey(indexl:indexr,indexd:indexu,2))
        allocate(Dy(indexl:indexr,indexd:indexu,2))
        allocate(CfEdgeY(indexl:indexr,indexd:indexu))

        do j = indexd, indexu
            do i = indexl, indexr

                if((prow==1).and.(j==1)) then
                    vxwAv = (vxwBdryY0(xlower+i-1)+vxwBdryY0(xlower+i)+vxw(i,j)+vxw(i+1,j)) / 4.D0
                elseif((prow==pnrows).and.(j==localnrows+1)) then
                    vxwAv = (vxw(i,j-1)+vxw(i+1,j-1)+vxwBdryY1(xlower+i-1)+vxwBdryY1(xlower+i)) / 4.D0
                else
                    vxwAv = (vxw(i,j-1)+vxw(i+1,j-1)+vxw(i,j)+vxw(i+1,j)) / 4.D0
                end if

                vwmodulus = dsqrt((vyw(i,j)**2.D0+vxwAv**2.D0)*1.D0)

                if(vwmodulus /= 0.D0) then
                    Ey(i,j,1) = vyw(i,j)*vxwAv / vwmodulus**2.D0
                    Ey(i,j,2) = vyw(i,j)**2.D0 / vwmodulus**2.D0
                else
                    Ey(i,j,1) = 0.D0
                    Ey(i,j,2) = 0.D0
                end if

                if((prow==1).and.(j==1)) then
                    dmEdgeY = dm(i,j)
                elseif((prow==pnrows).and.(j==localnrows+1)) then
                    dmEdgeY = dm(i,j-1)
                else
                    dmEdgeY = (hy(j-1)+hy(j)) / (hy(j-1)/dm(i,j-1)+hy(j)/dm(i,j))
                end if

                dl = alphaOS*dmEdgeY + 2.D0*lamdaX*vwmodulus*radiusInit*(1.D0-poroEdgeYInit(i,j))/ &!
                    (poroEdgeYInit(i,j)*(1.D0-poroEdgeY(i,j)))
                dt = alphaOS*dmEdgeY + 2.D0*lamdaT*vwmodulus*radiusInit*(1.D0-poroEdgeYInit(i,j))/ &!
                    (poroEdgeYInit(i,j)*(1.D0-poroEdgeY(i,j)))

                Dy(i,j,1) = vwmodulus/(poroEdgeY(i,j)*SwEdgeY(i,j))*(dl-dt)*Ey(i,j,1)
                Dy(i,j,2) = alphaOS*dmEdgeY + vwmodulus/(poroEdgeY(i,j)*SwEdgeY(i,j))*((dl-dt)*Ey(i,j,2)+dt)

            end do
        end do

        do j = indexd, indexu
            do i = indexl, indexr

                if((prow==1).and.(j==1).and.(isDiriY0_Cf(xlower+i-1)==1)) then
                    CfEdgeY(i,j) = CfBdryY0(xlower+i-1)
                elseif((prow==1).and.(j==1).and.(isDiriY0_Cf(xlower+i-1)==0)) then
                    CfEdgeY(i,j) = conc(i,j)
                elseif((prow==pnrows).and.(j==localnrows+1).and.(isDiriY1_Cf(xlower+i-1)==1)) then
                    CfEdgeY(i,j) = CfBdryY1(xlower+i-1)
                elseif((prow==pnrows).and.(j==localnrows+1).and.(isDiriY1_Cf(xlower+i-1)==0)) then
                    CfEdgeY(i,j) = conc(i,j-1)
                elseif(vyw(i,j) > 0.D0) then
                    CfEdgeY(i,j) = conc(i,j-1)
                else
                    CfEdgeY(i,j) = conc(i,j)
                end if

            end do
        end do

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
        allocate(dCfdx(1:localncols+1,indexd:indexu))

        do j = indexd, indexu
            do i = 1, localncols+1

                if((pcol==1).and.(i==1).and.(isDiriX0_Cf(ylower+j-1)==0)) then
                    dCfdx(i,j) = 0.D0
                elseif((pcol==1).and.(i==1).and.(isDiriX0_Cf(ylower+j-1)==1)) then
                    dCfdx(i,j) = (conc(i,j)-CfBdryX0(ylower+j-1))/hx(i)
                elseif((pcol==pncols).and.(i==localncols+1).and.(isDiriX1_Cf(ylower+j-1)==0)) then
                    dCfdx(i,j) = 0.D0
                elseif((pcol==pncols).and.(i==localncols+1).and.(isDiriX1_Cf(ylower+j-1)==1)) then
                    dCfdx(i,j) = (CfBdryX1(ylower+j-1)-conc(i-1,j))/hx(i-1)
                else
                    dCfdx(i,j) = (conc(i,j)-conc(i-1,j))/((hx(i)+hx(i-1))/2.D0)
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

        allocate(dCfdy(indexl:indexr,1:localnrows+1))

        do j = 1, localnrows+1
            do i = indexl, indexr

                if((prow==1).and.(j==1).and.(isDiriY0_Cf(xlower+i-1)==0)) then
                    dCfdy(i,j) = 0.D0
                elseif((prow==1).and.(j==1).and.(isDiriY0_Cf(xlower+i-1)==1)) then
                    dCfdy(i,j) = (conc(i,j)-CfBdryY0(xlower+i-1))/hy(j)
                elseif((prow==pnrows).and.(j==localnrows+1).and.(isDiriY1_Cf(xlower+i-1)==0)) then
                    dCfdy(i,j) = 0.D0
                elseif((prow==pnrows).and.(j==localnrows+1).and.(isDiriY1_Cf(xlower+i-1)==1)) then
                    dCfdy(i,j) = (CfBdryY1(xlower+i-1)-conc(i,j-1))/hy(j-1)
                else
                    dCfdy(i,j) = (conc(i,j)-conc(i,j-1))/((hy(j)+hy(j-1))/2.D0)
                end if

            end do
        end do

        do j = 1, localnrows
            do i = 1, localncols

                div1 = (vxw(i+1,j)*CfEdgeX(i+1,j)*poroEdgeX(i+1,j)*SwEdgeX(i+1,j) - &!
                    vxw(i,j)*CfEdgeX(i,j)*poroEdgeX(i,j)*SwEdgeX(i,j))/hx(i) + &!
                    (vyw(i,j+1)*CfEdgeY(i,j+1)*poroEdgeY(i,j+1)*SwEdgeY(i,j+1) - &!
                    vyw(i,j)*CfEdgeY(i,j)*poroEdgeY(i,j)*SwEdgeY(i,j))/hy(j)

                dCfdxleft = dCfdx(i,j)
                dCfdxright = dCfdx(i+1,j)
                if((prow==1).and.(j==1)) then
                    dCfdxdown = (dCfdx(i,j)+dCfdx(i+1,j)) / 2.D0
                else
                    dCfdxdown = (dCfdx(i,j)+dCfdx(i+1,j)+dCfdx(i,j-1)+dCfdx(i+1,j-1)) / 4.D0
                end if
                if((prow==pnrows).and.(j==localnrows)) then
                    dCfdxup = (dCfdx(i,j)+dCfdx(i+1,j)) / 2.D0
                else
                    dCfdxup = (dCfdx(i,j)+dCfdx(i+1,j)+dCfdx(i,j+1)+dCfdx(i+1,j+1)) / 4.D0
                end if

                if((pcol==1).and.(i==1)) then
                    dCfdyleft = (dCfdy(i,j)+dCfdy(i,j+1)) / 2.D0
                else
                    dCfdyleft = (dCfdy(i-1,j)+dCfdy(i-1,j+1)+dCfdy(i,j)+dCfdy(i,j+1)) / 4.D0
                end if
                if((pcol==pncols).and.(i==localncols)) then
                    dCfdyright = (dCfdy(i,j)+dCfdy(i,j+1)) / 2.D0
                else
                    dCfdyright = (dCfdy(i,j)+dCfdy(i,j+1)+dCfdy(i+1,j)+dCfdy(i+1,j+1)) / 4.D0
                end if
                dCfdydown = dCfdy(i,j)
                dCfdyup = dCfdy(i,j+1)

                div2 = (dCfdxright*poroEdgeX(i+1,j)*SwEdgeX(i+1,j)*Dx(i+1,j,1) - &!
                    dCfdxleft*poroEdgeX(i,j)*SwEdgeX(i,j)*Dx(i,j,1))/hx(i) + &!
                    (dCfdyright*poroEdgeX(i+1,j)*SwEdgeX(i+1,j)*Dx(i+1,j,2) - &!
                    dCfdyleft*poroEdgeX(i,j)*SwEdgeX(i,j)*Dx(i,j,2))/hx(i) + &!
                    (dCfdxup*poroEdgeY(i,j+1)*SwEdgeY(i,j+1)*Dy(i,j+1,1) - &!
                    dCfdxdown*poroEdgeY(i,j)*SwEdgeY(i,j)*Dy(i,j,1))/hy(j) + &!
                    (dCfdyup*poroEdgeY(i,j+1)*SwEdgeY(i,j+1)*Dy(i,j+1,2) - &!
                    dCfdydown*poroEdgeY(i,j)*SwEdgeY(i,j)*Dy(i,j,2))/hy(j)

                reaction = Sw(i,j)*av(i,j)*conc(i,j)*kc(i,j)*ks(i,j)/(kc(i,j)+ks(i,j))

                timeTerm = (conc(i,j)*poro(i,j)*Sw(i,j)-Cf(i,j)*poro_old(i,j)*Sw_old(i,j))/(timeEnd/nt)

                resi(i,j) = src(xlower+i-1,ylower+j-1) - reaction + div2 - div1 - timeTerm

            end do
        end do

        deallocate(Ex)
        deallocate(Ey)
        deallocate(Dx)
        deallocate(Dy)
        deallocate(CfEdgeX)
        deallocate(CfEdgeY)
        deallocate(dCfdx)
        deallocate(dCfdy)

    end subroutine Resi_tran_Cf

    ! End Calc Resi Concentration

    ! Calc Resi Temperature

    subroutine Resi_ener_Tem(tempe, resi)

        integer, dimension(:,:), pointer, intent(in) :: tempe
        real(kind=8), dimension(:,:), pointer, intent(in out) :: resi

        integer :: indexl, indexr, indexd, indexu
        real(kind=8), dimension(:,:), allocatable :: TemEdgeX, TemEdgeY
        real(kind=8), dimension(:,:), allocatable :: dTemdx, dTemdy
        real(kind=8) :: coe_l, coe_r, coe_d, coe_u
        real(kind=8) :: Hr, ke, vwabs2, vnabs2, Forchhw, Forchhn
        real(kind=8) :: dvdt, div1, div2, rhs1, rhs2, rhs3, rhs4
        integer :: i, j

        indexl = 1
        indexr = localncols + 1
        indexd = 1
        indexu = localnrows
        allocate(TemEdgeX(indexl:indexr,indexd:indexu))

        do j = indexd, indexu
            do i = indexl, indexr

                if((pcol==1).and.(i==1).and.(isDiriX0_Tem(ylower+j-1)==0)) then
                    TemEdgeX(i,j) = tempe(i,j)
                elseif((pcol==1).and.(i==1).and.(isDiriX0_Tem(ylower+j-1)==1)) then
                    TemEdgeX(i,j) = TemBdryX0(ylower+j-1)
                elseif((pcol==pncols).and.(i==localncols+1).and.(isDiriX1_Tem(ylower+j-1)==0)) then
                    TemEdgeX(i,j) = tempe(i-1,j)
                elseif((pcol==pncols).and.(i==localncols+1).and.(isDiriX1_Tem(ylower+j-1)==1)) then
                    TemEdgeX(i,j) = TemBdryX1(ylower+j-1)
                else
                    TemEdgeX(i,j) = (hx(i-1)+hx(i)) / (hx(i-1)/tempe(i-1,j)+hx(i)/tempe(i,j))
                end if

            end do
        end do

        indexl = 1
        indexr = localncols
        indexd = 1
        indexu = localnrows + 1
        allocate(TemEdgeY(indexl:indexr,indexd:indexu))

        do j = indexd, indexu
            do i = indexl, indexr

                if((prow==1).and.(j==1).and.(isDiriY0_Tem(xlower+i-1)==0)) then
                    TemEdgeY(i,j) = tempe(i,j)
                elseif((prow==1).and.(j==1).and.(isDiriY0_Tem(xlower+i-1)==1)) then
                    TemEdgeY(i,j) = TemBdryY0(xlower+i-1)
                elseif((prow==pnrows).and.(j==localnrows+1).and.(isDiriY1_Tem(xlower+i-1)==0)) then
                    TemEdgeY(i,j) = tempe(i,j-1)
                elseif((prow==pnrows).and.(j==localnrows+1).and.(isDiriY1_Tem(xlower+i-1)==1)) then
                    TemEdgeY(i,j) = TemBdryY1(xlower+i-1)
                else
                    TemEdgeY(i,j) = (hy(j-1)+hy(j)) / (hy(j-1)/tempe(i,j-1)+hy(j)/tempe(i,j))
                end if

            end do
        end do

        indexd = 1
        indexu = localnrows
        allocate(dTemdx(1:localncols+1,indexd:indexu))

        do j = indexd, indexu
            do i = 1, localncols+1

                if((pcol==1).and.(i==1).and.(isDiriX0_Tem(ylower+j-1)==0)) then
                    dTemdx(i,j) = 0.D0
                elseif((pcol==1).and.(i==1).and.(isDiriX0_Tem(ylower+j-1)==1)) then
                    dTemdx(i,j) = (tempe(i,j)-TemBdryX0(ylower+j-1))/hx(i)
                elseif((pcol==pncols).and.(i==localncols+1).and.(isDiriX1_Tem(ylower+j-1)==0)) then
                    dTemdx(i,j) = 0.D0
                elseif((pcol==pncols).and.(i==localncols+1).and.(isDiriX1_Tem(ylower+j-1)==1)) then
                    dTemdx(i,j) = (TemBdryX1(ylower+j-1)-tempe(i-1,j))/hx(i-1)
                else
                    dTemdx(i,j) = (tempe(i,j)-tempe(i-1,j))/((hx(i)+hx(i-1))/2.D0)
                end if

            end do
        end do

        indexl = 1
        indexr = localncols
        allocate(dTemdy(indexl:indexr,1:localnrows+1))

        do j = 1, localnrows+1
            do i = indexl, indexr

                if((prow==1).and.(j==1).and.(isDiriY0_Tem(xlower+i-1)==0)) then
                    dTemdy(i,j) = 0.D0
                elseif((prow==1).and.(j==1).and.(isDiriY0_Tem(xlower+i-1)==1)) then
                    dTemdy(i,j) = (tempe(i,j)-TemBdryY0(xlower+i-1))/hy(j)
                elseif((prow==pnrows).and.(j==localnrows+1).and.(isDiriY1_Tem(xlower+i-1)==0)) then
                    dTemdy(i,j) = 0.D0
                elseif((prow==pnrows).and.(j==localnrows+1).and.(isDiriY1_Tem(xlower+i-1)==1)) then
                    dTemdy(i,j) = (TemBdryY1(xlower+i-1)-tempe(i,j-1))/hy(j-1)
                else
                    dTemdy(i,j) = (tempe(i,j)-tempe(i,j-1))/((hy(j)+hy(j-1))/2.D0)
                end if

            end do
        end do

        do j = 1, localnrows
            do i = 1, localncols

                dvdt = rhofw*thetafw*(poro(i,j)*Sw(i,j)*tempe(i,j)-poro_old(i,j)*Sw_old(i,j)*Tem(i,j))/(timeEnd/nt) + &!
                    rhofn*thetafn*(poro(i,j)*(1.D0-Sw(i,j))*tempe(i,j)-poro_old(i,j)*(1.D0-Sw_old(i,j))*Tem(i,j))/(timeEnd/nt) + &!
                    rhos*thetas*((1.D0-poro(i,j))*tempe(i,j)-(1.D0-poro_old(i,j))*Tem(i,j))/(timeEnd/nt)

                div1 = rhofw*thetafw * ((poroEdgeX(i+1,j)*SwEdgeX(i+1,j)*vxw(i+1,j)*TemEdgeX(i+1,j) - &!
                    poroEdgeX(i,j)*SwEdgeX(i,j)*vxw(i,j)*TemEdgeX(i,j))/hx(i) + &!
                    (poroEdgeY(i,j+1)*SwEdgeY(i,j+1)*vyw(i,j+1)*TemEdgeY(i,j+1) - &!
                    poroEdgeY(i,j)*SwEdgeY(i,j)*vyw(i,j)*TemEdgeY(i,j))/hy(j)) + &!
                    rhofn*thetafn * ((poroEdgeX(i+1,j)*(1.D0-SwEdgeX(i+1,j))*vxn(i+1,j)*TemEdgeX(i+1,j) - &!
                    poroEdgeX(i,j)*(1.D0-SwEdgeX(i,j))*vxn(i,j)*TemEdgeX(i,j))/hx(i) + &!
                    (poroEdgeY(i,j+1)*(1.D0-SwEdgeY(i,j+1))*vyn(i,j+1)*TemEdgeY(i,j+1) - &!
                    poroEdgeY(i,j)*(1.D0-SwEdgeY(i,j))*vyn(i,j)*TemEdgeY(i,j))/hy(j))

                coe_l = poroEdgeX(i,j)*SwEdgeX(i,j)*Mw + poroEdgeX(i,j)*(1.D0-SwEdgeX(i,j))*Mn + (1.D0-poroEdgeX(i,j))*Ms
                coe_r = poroEdgeX(i+1,j)*SwEdgeX(i+1,j)*Mw + poroEdgeX(i+1,j)*(1.D0-SwEdgeX(i+1,j))*Mn + (1.D0-poroEdgeX(i+1,j))*Ms
                coe_d = poroEdgeY(i,j)*SwEdgeY(i,j)*Mw + poroEdgeY(i,j)*(1.D0-SwEdgeY(i,j))*Mn + (1.D0-poroEdgeY(i,j))*Ms
                coe_u = poroEdgeY(i,j+1)*SwEdgeY(i,j+1)*Mw + poroEdgeY(i,j+1)*(1.D0-SwEdgeY(i,j+1))*Mn + (1.D0-poroEdgeY(i,j+1))*Ms
                div2 = (dTemdx(i+1,j)*coe_r - dTemdx(i,j)*coe_l)/hx(i) + (dTemdy(i,j+1)*coe_u - dTemdy(i,j)*coe_d)/hy(j)

                Hr = -6.846D3 + 8.038*Tem(i,j) - 3.22D-3*Tem(i,j)**2.D0 - 8.703D2/Tem(i,j)

                ke = 2.D0*kc(i,j)*ks(i,j)/(kc(i,j)+2.D0*ks(i,j))
                rhs1 = -2.D0*Hr*av(i,j)*ke*Cf(i,j)

                vwabs2 = ((vxw(i,j)+vxw(i+1,j))/2.D0)**2.D0 + ((vyw(i,j)+vyw(i,j+1))/2.D0)**2.D0
                vnabs2 = ((vxn(i,j)+vxn(i+1,j))/2.D0)**2.D0 + ((vyn(i,j)+vyn(i,j+1))/2.D0)**2.D0
                rhs2 = viscw*poro(i,j)*Sw(i,j)/(Sw(i,j)**3.D0*Kxx(i,j))*vwabs2 + &!
                    viscn*poro(i,j)*(1.D0-Sw(i,j))/((1.D0-Sw(i,j))**3.D0*Kxx(i,j))*vnabs2

                rhs3 = viscw*poro(i,j)*Sw(i,j) * ((((vxw(i+1,j)-vxw(i,j))/hx(i))**2.D0 + ((vyw(i,j+1)-vyw(i,j))/hy(j)))**2.D0) + &!
                    viscn*poro(i,j)*(1.D0-Sw(i,j)) * ((((vxn(i+1,j)-vxn(i,j))/hx(i))**2.D0 + ((vyn(i,j+1)-vyn(i,j))/hy(j)))**2.D0)

                Forchhw = 1.75D0/dsqrt(1.5D2*(poro(i,j)*Sw(i,j))**3.D0)
                Forchhn = 1.75D0/dsqrt(1.5D2*(poro(i,j)*(1.D0-Sw(i,j)))**3.D0)
                rhs4 = rhofw*poro(i,j)*Sw(i,j)*Forchhw * (dsqrt(vwabs2))**3.D0 + &!
                    rhofn*poro(i,j)*(1.D0-Sw(i,j))*Forchhn * (dsqrt(vnabs2))**3.D0

                resi(i,j) = dvdt + div1 - div2 - rhs1 - rhs2 - rhs3 - rhs4

            end do
        end do

        deallocate(TemEdgeX)
        deallocate(TemEdgeY)
        deallocate(dTemdx)
        deallocate(dTemdy)

    end subroutine Resi_ener_Tem

    ! End Calc Resi Temperature

end module DBF_resi
