@echo off
cls

echo =================================

echo 设置Palabos根目录...
set PLBROOTDIR="D:\Palabos-V2.0R0"
echo PLB root = %PLBROOTDIR%

set MMGT_OPT=1
set MMGT_CLEAR=1 
set MMGT_REENTRANT="1" 

Set /p PARMNAME=在这里输入您想打开的参数文件名称:
If /I "%PARMNAME%"=="" set PARMNAME=Fibers-090-BFC-Bounceback-03

echo --------------------------------
echo 即将打开仿真文件Palabos.exe。
Set /p MODE=选择d -debug 或 r - release版本：
If /I "%MODE%"=="d" set MODE=Debug
if "%MODE%"=="r"  set MODE=Release
if "%MODE%"==""  set MODE=Release

echo 运行模式= %MODE%
pause

set COMMANDLINE="%PLBROOTDIR%\codeblocks\bin\%Mode%\Palabos.exe %~dp0%PARMNAME%.xml"

echo %COMMANDLINE%
pause

call %PLBROOTDIR%\codeblocks\bin\%Mode%\Palabos.exe %~dp0%PARMNAME%.xml
pause
