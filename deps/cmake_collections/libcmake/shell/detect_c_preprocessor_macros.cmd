@echo off
:: This file detects MSVC preprocessor macros.
:: Synopsis: %0 <compiler> <language> <out_dir>
:: Params:
::   - <compiler> The MSVC compiler (cl.exe)
::   - <language> Language ('c' or 'c++')
::   - <out_dir> The output directory
:: Created by TONGYI LINGMA

setlocal enabledelayedexpansion

:: Check arguments
if "%~3"=="" (
    echo Usage: %0 ^<compiler^> ^<language^> ^<output_dir^>
    exit /b 1
)

set CEU_DCPM_CC=%~1
set CEU_DCPM_LANG=%~2
set CEU_DCPM_OUT=%~3

:: Validate language parameter
if /i not "%CEU_DCPM_LANG%"=="c" (
    if /i not "%CEU_DCPM_LANG%"=="c++" (
        echo Language must be 'c' or 'c++'
        exit /b 1
    )
)

:: Create output directory
mkdir "%CEU_DCPM_OUT%" 2>nul
if errorlevel 1 (
    echo Failed to create output directory: %CEU_DCPM_OUT%
    exit /b 1
)

:: Function to perform command execution
:perform_cmd
set CMD_EXEC=%~1
echo Executing: %CMD_EXEC%
%CMD_EXEC% >"%CEU_DCPM_OUT%\test_%CEU_DCPM_LANG%_cpp.h" 2>"%CEU_DCPM_OUT%\test_%CEU_DCPM_LANG%_cpp.err"
if %errorlevel% equ 0 (
    exit /b 0
)
exit /b 1

:: Try different MSVC invocation methods for preprocessing
set COMMANDS[0]="%CEU_DCPM_CC%" /nologo /EP /C /D _MSC_VER /showIncludes "%CEU_DCPM_OUT%\temp_dummy.c"
set COMMANDS[1]="%CEU_DCPM_CC%" /nologo /EP /C /D _MSC_VER NUL
set COMMANDS[2]="%CEU_DCPM_CC%" /nologo /E /C /D _MSC_VER NUL

:: Create dummy file for preprocessing
echo. > "%CEU_DCPM_OUT%\temp_dummy.c"

:: Execute commands
for %%i in (0,1,2) do (
    call :perform_cmd !COMMANDS[%%i]!
    if !errorlevel! equ 0 (
        del "%CEU_DCPM_OUT%\temp_dummy.c" 2>nul
        exit /b 0
    )
)

:: Cleanup
del "%CEU_DCPM_OUT%\temp_dummy.c" 2>nul

echo All MSVC invocation attempts failed
exit /b 1
