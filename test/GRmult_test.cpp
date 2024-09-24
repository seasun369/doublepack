#include "../PSS/pss.h"
#include "../RMFE/RMFE.h"


using namespace std;
using namespace packed_shamir;

int main()
{
    
    ZZ gr_prime = conv<ZZ>(2);  //p 
    long positive_integer = 16; //k
    long gr_degree = 4;         //r

    gr galoisring(gr_prime, positive_integer, gr_degree);

    ZZ_pX f = galoisring.Fps_poly;

    ZZ_pE a;
    ZZ_pE b;

    a = galoisring.set_T[5];
    b = galoisring.set_T[12];

    ZZ_pE cc = a*b;

    vector<vector<long>> matrix_b;

    generateMatrix(b, gr_degree, matrix_b);

    vector<long> vec_a;
    ZzpE2Veclong(a, vec_a, gr_degree);
    vector<long> vec_cc;
    ZzpE2Veclong(cc, vec_cc, gr_degree);

    vector<long> c_ = multiplyMatrixByVector(matrix_b, vec_a);

    vector<long> c = compression(f, c_);

    print(c);

    for(int i = 0; i<c.size(); i++){
        c[i] = fresh(c[i], positive_integer);
    }
    print(c);
    print(vec_cc);


}