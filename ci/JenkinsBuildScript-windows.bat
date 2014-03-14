cd build

cmake -DUSE_ATLAS=ON -DUSE_MKL=OFF -DENABLE_PYTHON=ON -DBUILD_TESTING=ON -DENABLE_PACKAGING=ON -DMATLAB_TESTING=OFF ..
cmake --build . --config RelWithDebInfo
cmake --build . --config RelWithDebInfo --target update
cmake --build . --config RelWithDebInfo --target build

del build\JUnitTestResults.xml
pushd build
ctest.exe -T Test -C RelWithDebInfo --output-on-failure
verify >nul
popd
"C:\Python27\python.exe" ci/CTest2JUnit.py build ci/CTest2JUnit.xsl > build/JUnitTestResults.xml
