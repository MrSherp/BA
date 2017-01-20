#include <iostream>
#include <fstream>
#include <sstream>
#include <string>

using namespace std;

template < typename RealType, typename VectorNd >
void loadSignal (VectorNd& Destination, const char *Filename ){
    
    string line;
    ifstream myfile;
    myfile.open(Filename);
    int N = 0;
    if (myfile.is_open()){
        while ( getline (myfile,line) ){
            ++N;
        }
    }
    else cout << "unable to open file" << endl;
    myfile.close();
    Destination.resize( N );
    int i = 0;    
    myfile.open(Filename);
    while ( getline (myfile,line) ){
        Destination(i) = std::stod (line);
        ++i;
    }
    myfile.close();
}