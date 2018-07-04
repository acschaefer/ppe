#include <iostream>

#include <ptxread.h>
#include <pcextrpln.h>
#include <pcdraw.h>

int main (int argc, char** argv)
{
  ppe::pointCloud* pc = ppe::ptxread("../../data/cottage_01_head.ptx");
  pc = ppe::pcextrpln(pc);
  ppe::pcdraw(pc);
  delete pc;
  return 0;
}
