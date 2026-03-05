#include <NTL/ZZ.h>
#include <iostream>

int main() {
#ifdef NTL_THREADS
    std::cout << "NTL 线程安全已开启！" << std::endl;
#else
    std::cout << "警告：NTL 仍处于单线程模式。" << std::endl;
#endif
    return 0;
}