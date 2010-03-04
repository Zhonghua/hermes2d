#include <iostream>
#include <stdexcept>

#include "matrix.h"

#define ERROR_SUCCESS                               0
#define ERROR_FAILURE                              -1

void _assert(bool a)
{
    if (!a) throw std::runtime_error("Assertion failed.");
}

void test_matrix1()
{
    CooMatrix m(5);
    m.add(1, 3, 3.5);
    m.add(2, 3, 4.5);
    m.add(3, 4, 1.5);
    m.add(4, 2, 1.5);
    m.add(2, 3, 1);
    m.print();

    printf("----\n");
    // convert from COO
    CSRMatrix n1(&m);
    n1.print();
    CSCMatrix n2(&m);
    n2.print();
    // convert CSR <-> CSC
    CSRMatrix n3(&n2);
    n3.print();
    CSCMatrix n4(&n1);
    n4.print();
}

int main(int argc, char* argv[])
{
    try {
        test_matrix1();

        return ERROR_SUCCESS;
    } catch(std::exception const &ex) {
        std::cout << "Exception raised: " << ex.what() << "\n";
        return ERROR_FAILURE;
    } catch(...) {
        std::cout << "Exception raised." << "\n";
        return ERROR_FAILURE;
    }
}