#include <cstdlib>  // For system()
#include <iostream>  // For cout
using namespace std;
int main() {
    freopen("output.txt","w",stdout);
    for (int i = 0; i < 30; i++) {
        system("GEN-EPWOCD.exe");  // Runs x.exe 30 times
    }
    return 0;
}
 