Author: Yuanqing Wu, DGUT, P.R.China

History: 2025-9-29 by Yuanqing Wu

Support: wuyuanq@gmail.com

This is the parallel program to simulate the wormhole issue in 2D condition.

If you want to simulate with DBF framework, you must make sure that the three parameters isDarcy, isBrinkman and isForchheimer in the file "DBF_globalData.F90" are true. If you want to simulate with only Darcy framework, you must make sure the following things: isDarcy = .true., isBrinkman = .false. and isForchheimer = .false..

There are many solvers to choose. You can choose the solver in the Makefile. For example, if you want to use Hypre solve, you can uncomment the statement "SOLVER = HYPRE" in the Makefile.mac.

If you want to run the code in parallel, please set the number of the processors in the variable "NP" in the file Makefile.mac. You should do the same thing in the DBF.infile.F90. In the file, there is a statement "call proceAlloc(1, nx, ny, pncols, pnrows)", and please fill the number of the processors in the first parameter.

You can see the results in the document "case". In "case", there is a file called "matlabplot.m". You can run the file to generate the matlab figures in the document "matlabplots". The matlab figures only show the physical results at the end of the simulation.

You can also see the results in Tecplot. A series of .plt files have been generated in "case". The numbers at the tail of the name of the .plt files stand for the serial number of the time step. So in Tecplot, you can see the simulation history results and not only the results at the end of the simulation.

For a quick test of the code, please just open a terminal on Mac, enter the directory of the code, and input "make run -f Makefile.mac" in the terminal.
