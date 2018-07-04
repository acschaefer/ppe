#include <iostream>

#include <ptxread.h>
#include <pcdraw.h>

int main (int argc, char** argv)
{
  ple::pointCloud* pc = ple::ptxread("../../data/cottage_01_head.ptx");
  ple::pcdraw(pc);
  return 0;
}
