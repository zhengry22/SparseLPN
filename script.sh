mkdir -p build && cd build
rm -rf *
cmake -DCMAKE_BUILD_TYPE=Debug ..
make -j
# 运行时通过环境变量控制线程数（防止内存溢出）
export OMP_NUM_THREADS=4
#./testSparseMatrix
#./testEncScheme > debug.txt
./testEncScheme
# gdb ./testEncScheme