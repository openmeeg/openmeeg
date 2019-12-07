#include <memory>

int main() {
    std::shared_ptr<double[]> ptr(new double[16]);
    return 0;
}
