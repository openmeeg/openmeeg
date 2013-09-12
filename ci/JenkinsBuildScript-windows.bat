del build\JUnitTestResults.xml
pushd build
ctest -T Test -C RelWithDebInfo --output-on-failure
verify >nul
popd
"C:\Python27\python.exe" ci/CTest2JUnit.py build ci/CTest2JUnit.xsl > build/JUnitTestResults.xml
