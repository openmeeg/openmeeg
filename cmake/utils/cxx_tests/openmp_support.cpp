#include <iostream>
#include <vector>

int
main() {
    std::vector<unsigned> V;
    for (unsigned i=0;i<100;++i)    
        V.push_back(i); 
    #pragma omp parallel for
    #if defined RANGEFOR
    for (const auto& val : V) {
    #elif defined ITERATOR
    for (std::vector<unsigned>::const_iterator v<V.begin();vi!=V.end();++vi) {
        const unsigned val = *vi;
    #else
    for (unsigned i=0;i<V.size();++i) {
        const unsigned val = V[i];
    #endif
        std::cerr << val << std::endl;
    }

    return 0;
}
