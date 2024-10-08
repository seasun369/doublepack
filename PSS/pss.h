#ifndef PSS_GR
#define PSS_GR

#include <iostream>
#include <vector>
#include <random>
#include "../GR/gr.h"
using namespace std;

namespace packed_shamir
{
	typedef struct point
	{
		ZZ_pE x;
		ZZ_pE y;
	}point;

	typedef struct int_point
	{
     	int x;
		int y;
	}int_point;

	typedef std::vector<point> shares; 

	class scheme
	{
		int n; //members
		int m; //packed number
		int d; //poly degree
		//int t; //corrupt parties
		gr GR;
		//vec_ZZ_pE alpha_set;
		//vec_ZZ_pE beta_set;

		public:

			scheme(int members,int packed_number,int degree, gr galois_ring):
			        n(members), m(packed_number), d(degree), GR(galois_ring)
			{
				
				alpha_set.SetLength(n);
				beta_set.SetLength(m);

				for(int i=0; i<n; i++)
				{		
					alpha_set[i] = GR.set_T[i];
				}

				for(int i=0; i<m; i++)
				{
					beta_set[i] = GR.set_T[i+n];
				}
				
			}

			scheme() {}

			vec_ZZ_pE alpha_set;
			vec_ZZ_pE beta_set;


			vec_ZZ_pE create_shares(vec_ZZ_pE a);
			vector<vec_ZZ_pE> packed_create_shares(vec_ZZ_pE secret);
			vec_ZZ_pE packed_reconstruct_shares(vector<int> party, vec_ZZ_pE shares);
			vec_ZZ_pE create_one_shares(ZZ_pE a, long i);
			ZZ_pE reconstruct_one_shares(vec_ZZ_pE shares, long i);
			vec_ZZ_pE create_shares_with_points(vector<ZZ_pE> a, vector<ZZ_pE>b);
	};
}


#endif