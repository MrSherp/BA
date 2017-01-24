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
    
    RealType max = Destination.maxCoeff();
    RealType min = Destination.minCoeff();
    RealType scale = 1. / (max - min);    
    for( int i = 0; i < N; ++i){
        Destination ( i ) -= min;
        Destination ( i ) *= scale;
    }
}



template < typename RealType, typename VectorNd >
void safeSignal ( VectorNd& Signal, const char *Filename ){
    
    ofstream myfile;
    myfile.open (Filename);
    int N = Signal.size();
    for (int i = 0; i < N; ++i){
        myfile << Signal(i) << std::endl;
    }
}



template < typename RealType, typename VectorNd >
void threshholding ( VectorNd Image, VectorNd& Signal, const RealType threshhold ){
    const unsigned int N = Signal.size();
    for ( int i = 0; i < N; ++i){
        if( threshhold > 1 || threshhold < 0 ){
            std::cout << "Threshhold not between 0 and 1!" << std::endl;
            break;
        }
        int j = 0;
        while( Image ( i + (N - 1 - j) * N ) >= threshhold ){            
            ++j;
        }
        Signal( i ) = static_cast < RealType > (j-1) / static_cast < RealType > (N-1);
    }    
}


