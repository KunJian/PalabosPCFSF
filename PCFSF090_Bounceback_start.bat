@echo off
cls

echo =================================

echo ����Palabos��Ŀ¼...
set PLBROOTDIR="D:\Palabos-V2.0R0"
echo PLB root = %PLBROOTDIR%

set MMGT_OPT=1
set MMGT_CLEAR=1 
set MMGT_REENTRANT="1" 

Set /p PARMNAME=��������������򿪵Ĳ����ļ�����:
If /I "%PARMNAME%"=="" set PARMNAME=Fibers-090-BFC-Bounceback-03

echo --------------------------------
echo �����򿪷����ļ�Palabos.exe��
Set /p MODE=ѡ��d -debug �� r - release�汾��
If /I "%MODE%"=="d" set MODE=Debug
if "%MODE%"=="r"  set MODE=Release
if "%MODE%"==""  set MODE=Release

echo ����ģʽ= %MODE%
pause

set COMMANDLINE="%PLBROOTDIR%\codeblocks\bin\%Mode%\Palabos.exe %~dp0%PARMNAME%.xml"

echo %COMMANDLINE%
pause

call %PLBROOTDIR%\codeblocks\bin\%Mode%\Palabos.exe %~dp0%PARMNAME%.xml
pause
