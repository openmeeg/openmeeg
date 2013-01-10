del build_32\JUnitTestResults.xml
pushd build_32\Tests
"C:\Program Files\CMake 2.8\bin\ctest.exe" -T Test -C RelWithDebInfo --output-on-failure
popd
verify >nul
C:\Python27\python.exe ci/CTest2JUnit.py build_32/Tests ci/CTest2JUnit.xsl > build_32/JUnitTestResults.xml
