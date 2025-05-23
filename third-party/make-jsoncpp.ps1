$Root = $PSScriptRoot
$JsonCppVersion = "1.9.5"
$DownloadDir = "$Root/download"
$SourcesDir = "$DownloadDir/jsoncpp/$JsonCppVersion"

$Remote = "https://github.com/open-source-parsers/jsoncpp/archive/refs/tags/$JsonCppVersion.zip"
$Local = "$DownloadDir/jsoncpp-$JsonCppVersion.zip"

if (Test-Path "$Local" -PathType Leaf) {
    Write-Output "Skipping download"
}
else {
    New-Item -ItemType "directory" -Path "$DownloadDir"
    Invoke-WebRequest -URI "$Remote" -OutFile "$Local"
}

if (Test-Path "$SourcesDir" -PathType Container) {
    Remove-Item -LiteralPath "$SourcesDir" -Force -Recurse
}
New-Item -ItemType "directory" -Path "$SourcesDir"
Expand-Archive "$Local" -DestinationPath "$SourcesDir"


Get-ChildItem "$SourcesDir\jsoncpp-$JsonCppVersion" | ForEach-Object {
    Move-Item -Path "$_" -Destination "$SourcesDir" -Force
}
Remove-Item -LiteralPath "$SourcesDir\jsoncpp-$JsonCppVersion" -Force

$InstallPrefix = "$Root/jsoncpp"
if (Test-Path "$InstallPrefix" -PathType Container) {
    Remove-Item -LiteralPath "$InstallPrefix" -Force -Recurse
}
New-Item -ItemType "directory" -Path "$InstallPrefix"

$BuildDir = "$SourcesDir/build"
if (Test-Path "$BuildDir" -PathType Container) {
    Remove-Item -LiteralPath "$BuildDir" -Force -Recurse
}
New-Item -ItemType "directory" -Path "$BuildDir"

$Cmake = "C:\Program Files\Microsoft Visual Studio\2022\Community\Common7\IDE\CommonExtensions\Microsoft\CMake\CMake\bin\cmake.exe"
$MSBuild = "C:\Program Files\Microsoft Visual Studio\2022\Community\MSBuild\Current\Bin\MSBuild.exe"

Push-Location "$BuildDir"
& "$Cmake" @(
    "-G", "Visual Studio 17 2022", 
    "-DCMAKE_INSTALL_PREFIX='$InstallPrefix'", 
    "-DARCHIVE_INSTALL_DIR='$InstallPrefix'",
    "-DCMAKE_VS_INCLUDE_INSTALL_TO_DEFAULT_BUILD=ON", 
    "-DBUILD_STATIC_LIBS=ON",
    "-DBUILD_SHARED_LIBS=OFF",
    "-DCMAKE_BUILD_TYPE=debug", 
    "../"
)

& "$MSBuild" @(
    "/target:Clean"
    "/target:Build"
    "/property:Configuration=Debug"
    ".\jsoncpp.sln"
)

Pop-Location

