cd ..
tar -cjvf pack/multidirsolver.tar.bz2 mds.cpp multidirsolver.{cpp,h} CMakeLists.txt
cd pack

rm -rf tmp
mkdir -p tmp
cd tmp
tar -xjvf ../multidirsolver.tar.bz2
mkdir build
cd build
cmake ../
make
cd ../..
rm -rf tmp
