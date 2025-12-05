@echo off
:: This file detects MSVC preprocessor macros.
:: Synopsis: %0 <compiler> <language> <out_dir>
:: Params:
::   - <compiler> The MSVC compiler (cl.exe)
::   - <language> Language ('c' or 'cpp')
::   - <out_dir> The output directory
:: Created by TONGYI LINGMA

setlocal enabledelayedexpansion

if "%~3"=="" (
    echo Usage: %0 ^<compiler^> ^<language^> ^<output_dir^>
    exit /b 1
)

set CEU_DCPM_CC=%~1
set CEU_DCPM_LANG=%~2
set CEU_DCPM_OUT=%~3

:: Validate language
if /i not "%CEU_DCPM_LANG%"=="c" (
    if /i not "%CEU_DCPM_LANG%"=="cpp" (
        echo Language must be 'c' or 'cpp'
        exit /b 1
    )
)

:: Create output directory
mkdir "%CEU_DCPM_OUT%" 2>nul
if errorlevel 1 (
    echo Failed to create output directory: %CEU_DCPM_OUT%
    exit /b 1
)

:: Create temporary source file
echo int main(){return 0;} > "%CEU_DCPM_OUT%\temp_source.c"

:: Try different MSVC invocation methods
set COMMANDS[0]="%CEU_DCPM_CC%" /B /EP "%CEU_DCPM_OUT%\temp_source.c"
set COMMANDS[1]="%CEU_DCPM_CC%" /showIncludes "%CEU_DCPM_OUT%\temp_source.c"
set COMMANDS[2]="%CEU_DCPM_CC%" /P "%CEU_DCPM_OUT%\temp_source.c"

:: Execute commands
for %%i in (0,1,2) do (
    if defined COMMANDS[%%i] (
        echo Executing: !COMMANDS[%%i]!
        !COMMANDS[%%i]! >"%CEU_DCPM_OUT%\test_%CEU_DCPM_LANG%_cclog.log" 2>"%CEU_DCPM_OUT%\test_%CEU_DCPM_LANG%_cclog.err"
        if !errorlevel! equ 0 (
            echo Success with command %%i
            exit /b 0
        )
    )
)

echo All MSVC invocation attempts failed
exit /b 1
