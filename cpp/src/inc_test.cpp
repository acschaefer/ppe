#include <iostream>


void cpp_fun() {
  std::cout << "this is cpp" << std::endl;
}

extern "C" {

#include "inc_test.h"

void inc_test() {
  std::cout << "jello" << std::endl;
  cpp_fun();
}

}
