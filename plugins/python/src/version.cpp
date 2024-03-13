#include <Python.h>
#include <iostream>
int main() {
  std::string include("2.9.2");
  if (PY_MAJOR_VERSION >= 3 && PY_MINOR_VERSION > 5) include = "2.10.4"; 
  std::cout << include << "\n";
  return 0;
}
