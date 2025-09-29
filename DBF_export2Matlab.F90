
!!$ Author:
!!$   Yuanqing Wu, DGUT, P.R.China
!!$
!!$ History:
!!$   2025-5-9 by Yuanqing Wu
!!$
!!$ Support:
!!$   wuyuanq@gmail.com

module DBF_export2Matlab
    
    use DBF_globalData
    implicit none

contains

    subroutine export2Matlab()

        character(len=30) :: fpmatlabm
        character(len=8) :: charLx, charLy
        character(len=5) :: charnx, charny
        character(len=11) :: chartimeEnd, charnt, chart, charnframe
        logical :: alive
        integer :: ierr

        fpmatlabm = trim(adjustl(soludoc))//"/matlabplot.m"

        open(unit=10, file=fpmatlabm, status='replace', iostat=ierr)
        if(ierr /= 0) then
            print *, 'open file ', fpmatlabm, ' error. ', ierr
            stop
        end if

        write(10, fmt="(a)") "path('/Users/yuanqingwu/research/MatrixAcidization/V_5.0/2D', path);"
        write(charLx,'(f8.4)') Lx
        write(10, fmt="(a)") "model.Lx = "//trim(adjustl(charLx))//";"
        write(charLy,'(f8.4)') Ly
        write(10, fmt="(a)") "model.Ly = "//trim(adjustl(charLy))//";"
        write(chartimeEnd,'(f11.1)') timeEnd
        write(10, fmt="(a)") "model.timeEnd = "//trim(adjustl(chartimeEnd))//";"
        write(charnx,'(i5)') nx
        write(10, fmt="(a)") "model.nx = "//trim(adjustl(charnx))//";"
        write(charny,'(i5)') ny
        write(10, fmt="(a)") "model.ny = "//trim(adjustl(charny))//";"
        write(charnt,'(i10)') nt
        write(10, fmt="(a)") "model.nt = "//trim(adjustl(charnt))//";"
        write(10, fmt="(a)") "model.xs = (0:model.nx)*model.Lx/model.nx;"
        write(10, fmt="(a)") "model.ys = (0:model.ny)*model.Ly/model.ny;"
        write(10, fmt="(a)") "model.ts = (0:model.nt)*model.timeEnd/model.nt;"
        write(charnframe,'(i10)') NUMFRAME
        write(10, fmt="(a)") "model.NUMFRAME = "//trim(adjustl(charnframe))//";"
        write(chart,'(i10)') t
        write(10, fmt="(a)") "model.fporotxt = 'soln_poro_raw_"//trim(adjustl(chart))//".txt';"
        write(10, fmt="(a)") "model.fSwtxt = 'soln_Sw_raw_"//trim(adjustl(chart))//".txt';"
        write(10, fmt="(a)") "model.fKxxtxt = 'soln_Kxx_raw_"//trim(adjustl(chart))//".txt';"
        write(10, fmt="(a)") "model.fvxwtxt = 'soln_vxw_raw_"//trim(adjustl(chart))//".txt';"
        write(10, fmt="(a)") "model.fvxntxt = 'soln_vxn_raw_"//trim(adjustl(chart))//".txt';"
        write(10, fmt="(a)") "model.fvywtxt = 'soln_vyw_raw_"//trim(adjustl(chart))//".txt';"
        write(10, fmt="(a)") "model.fvyntxt = 'soln_vyn_raw_"//trim(adjustl(chart))//".txt';"
        write(10, fmt="(a)") "model.fptxt = 'soln_p_raw_"//trim(adjustl(chart))//".txt';"
        write(10, fmt="(a)") "model.fCftxt = 'soln_Cf_raw_"//trim(adjustl(chart))//".txt';"
        write(10, fmt="(a)") "model.fTemtxt = 'soln_Tem_raw_"//trim(adjustl(chart))//".txt';"
        write(10, fmt="(a)") "model.fporohistxt = 'his_poro_avg.txt';"
        write(10, fmt="(a)") "model.fKxxhistxt = 'his_Kxx_avg.txt';"
        write(10, fmt="(a)") "model.favhistxt = 'his_av_avg.txt';"
        write(10, fmt="(a)") "model.fphistxt = 'his_p_avg.txt';"
        write(10, fmt="(a)") "model.fCfhistxt = 'his_Cf_avg.txt';"
        write(10, fmt="(a)") "model.fTemhistxt = 'his_Tem_avg.txt';"
        write(10, fmt="(a)") "model.fqhistxt = 'his_q_avg.txt';"
        write(10, fmt="(a)") "model.flphistxt = 'his_lp_avg.txt';"
        write(10, fmt="(a)") "model.soludoc = 'matlabplots';"
        inquire(file = trim(adjustl(soludoc))//'/matlabplots', exist = alive)
        if(.not.alive) then
            call system("mkdir "//trim(adjustl(soludoc))//"/matlabplots")
        end if
        write(10, fmt="(a)") "DBF_plot(model);"
        write(10, fmt="(a)") "DBF_his(model);"
        write(10, fmt="(a)") "rmpath('/Users/yuanqingwu/research/MatrixAcidization/V_5.0/2D');"

        close(10)

    end subroutine export2Matlab

end module DBF_export2Matlab
