#include "../PSS/pss.h"
#include "../RMFE/RMFE.h"

using namespace std;
using namespace packed_shamir;

int main()
{
    ZZ gr_prime = conv<ZZ>(2);  //p 
    long positive_integer = 4; //k
    long gr_degree = 4;         //r

    //cout << "Current modulus: " << ZZ_p::modulus() << endl;

    gr galoisring(gr_prime, positive_integer, gr_degree);

    cout << "Current modulus: " << ZZ_p::modulus() << endl;
    cout << "Current modulus polynomial: " << ZZ_pE::modulus() << endl;

    int members = 5;
	int packed_number = 2;
	int poly_degree = 3;
	int corrupt_parties = 2;

    packed_shamir::scheme Scheme(members,packed_number,poly_degree,corrupt_parties,galoisring);

    //std::cout << members << std::endl;

    vec_ZZ_pE test;
    test.SetLength(2);

    test[0] = random_ZZ_pE();
    test[1] = random_ZZ_pE();

    std::cout << test[0] << std::endl;
    std::cout << test[1] << std::endl;

    vec_ZZ_pE v = Scheme.create_shares(test);

    for(int i=0; i< members; i++)
    {
        std::cout << v[i] << std::endl;
    }

    vector<int> party = {1,2,3,4,5};

    vec_ZZ_pE u = Scheme.packed_reconstruct_shares(party, v);

    for(int i=0; i< packed_number; i++)
    {
        std::cout << u[i] << std::endl;
    }

}