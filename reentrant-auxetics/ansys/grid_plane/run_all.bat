@echo off
rem Runs every cell in inp\ with ANSYS Mechanical APDL in batch mode, one after another.
rem Each run appends one line to grid_plane_results.csv. Logs go to out\.
setlocal
set ANSYS="C:\Program Files\ANSYS Inc\v251\ansys\bin\winx64\ANSYS251.exe"
cd /d "%~dp0"
if not exist out mkdir out
if not exist work mkdir work
for %%f in (inp\*.inp) do (
  echo Running %%~nf ...
  %ANSYS% -b -np 2 -dir "%~dp0work" -j %%~nf -i "%~dp0%%f" -o "%~dp0out\%%~nf.out"
)
copy /Y "%~dp0work\grid_plane_results.csv" "%~dp0grid_plane_results.csv" >nul
echo.
echo Done. Results: grid_plane_results.csv
pause
