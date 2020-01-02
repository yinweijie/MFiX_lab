@echo on

setlocal
set PYTHONPATH=F:/CFD_workdir/mfix_workdir/2019/fluid_bed_dem_2d;;%PYTHONPATH%
call "D:/Users/ywj123450/Anaconda3/envs/mfix-19.1/python.exe" -m mfixgui.pymfix %*
endlocal