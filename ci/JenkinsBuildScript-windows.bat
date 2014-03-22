if not exist build (
    mkdir build
)

cd build

set CMAKE="C:\Program Files\CMake 2.8\bin\cmake.exe"

%CMAKE% -G"Visual Studio 11" -DUSE_MKL=OFF -DENABLE_PYTHON=ON -DBUILD_TESTING=ON -DENABLE_PACKAGING=ON -DMATLAB_TESTING=OFF ..
%CMAKE% --build . --config RelWithDebInfo
%CMAKE% --build . --config RelWithDebInfo --target update
%CMAKE% --build . --config RelWithDebInfo --target build

del build\JUnitTestResults.xml
pushd build
ctest.exe -T Test -C RelWithDebInfo --output-on-failure
verify >nul
popd
"C:\Python27\python.exe" ci/CTest2JUnit.py build ci/CTest2JUnit.xsl > build/JUnitTestResults.xml
