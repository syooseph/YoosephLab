#include <iostream>
#include <list>
#include <string>
#include <boost/filesystem.hpp>

using namespace std;

std::string GetFileStem(const std::string& path)  {
  // going backward untill the '\/' character
  cout << path << endl;
  int i;
  for(i = path.length() - 1; i >= 0; -- i) {
    cout << i << endl;
    if(path[i] == '/')  {
      break;
    }
  }
  cout << i + 1 << "  " << path.length() - i - 1 << endl;
  return path.substr(i + 1, path.length() - i - 1);
}

int main()  {
  string path = "sim.1VDX.32.fa";
  string read_stem = GetFileStem(path);
  cout << read_stem << endl;
  return 0;
}
