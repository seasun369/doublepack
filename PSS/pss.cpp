#include "pss.h"
#include "../RMFE/utils.hpp"
using namespace std;
namespace packed_shamir
{

//#define MEMBERS 3//17//5//5//9//3//
//#define THRESHOLD 2 //people can reveal the secret

/*
ZZ gr_prime = conv<ZZ>(2);  //p 
long positive_integer = 64; //k
long gr_degree = 4;         //r
*/

typedef vector<point> shares;

/*
caculate the coefficients of product of two polynomials
*/
vec_ZZ_pE multiplyPolynomials(const vec_ZZ_pE& A, const vec_ZZ_pE& B) {
    long m = A.length();
    long n = B.length();
    vec_ZZ_pE C;
	C.SetLength(m+n-1);
    
    for (int i = 0; i < m; ++i) {
        for (int j = 0; j < n; ++j) {
            C[i + j] = C[i+j] + A[i] * B[j];
        } 
    }
    
    return C;
}

/*
poly = gh+f, f satisfies f(beta_i) = a_i, 
degree(g) = d-m, degree(h) = m.
*/
vec_ZZ_pE scheme::create_shares(vec_ZZ_pE a)
{
	vec_ZZ_pE value;
	value.SetLength(n);
	long _m = a.length();
	if(m != _m) Error("PSS: vector length mismatch");

	ZZ p = power(GR.p, GR.k);
	ZZ_p::init(p);
    ZZ_pX F = GR.Fps_poly;
    ZZ_pE::init(F);

	ZZ_pEX f;

	interpolate_for_GR(f, beta_set, a, GR.p, GR.k, GR.r);

	vec_ZZ_pE res;
	res.SetLength(d-m);

	while(res[d-m-1] == conv<ZZ_pE>(0))
	{
		for(int i=0; i<d-m; i++)
		{
			ZZ_pE temp = random_ZZ_pE();
			res[i] = temp;
		}
	}

	ZZ_pEX g;
	g.rep = res;

	vec_ZZ_pE ret;
    ret.SetLength(m);
	ret[0] = conv<ZZ_pE>(1);
	for (int i = 1; i < m; i++)
	{
		ret[i] = conv<ZZ_pE>(0);
	}
	
	for(int i=0; i<m; i++)
	{
		vec_ZZ_pE temp;
		temp.SetLength(2);
		temp[0] = -beta_set[i];
		temp[1] = conv<ZZ_pE>(1);
		ret = multiplyPolynomials(ret, temp);
	}

	ZZ_pEX h;
	h.rep = res;

	ZZ_pEX poly = g*h + f;

	for(int i=0; i<n; i++)
	{
		value[i] = eval(poly, alpha_set[i]);
	}

	return value;

}

vector<vec_ZZ_pE> scheme::packed_create_shares(vec_ZZ_pE secret)
{
	long number = secret.length();
	long newnum = (number/m)*m+m;
	/*
	Zero padding
	*/
	if(number % m != 0)
	{
		secret.SetLength(newnum);
		for(int i=number; i<newnum; i++)
		{
			secret[i] = ZZ_pE();
		}
	}

	long len = newnum / m;

	vector<vec_ZZ_pE> shares;
	for(int i=0; i<m; i++)
	{
		vec_ZZ_pE temp;
		temp.SetLength(m);
		for(int j=0; j<m; j++)
		{
			temp[j] = secret[i*m+j];
		}
		vec_ZZ_pE ret = create_shares(temp);
		shares.emplace_back(ret);
	}

	return shares;
}

vec_ZZ_pE scheme::packed_reconstruct_shares(vector<int> party, vec_ZZ_pE shares)
{
	long num = party.size();

	if(num != shares.length()) Error("PSS: reconstruct vector length mismatch");

	vec_ZZ_pE points;
	points.SetLength(num);

	long j = 0;
	for(auto i : party)
	{
		points[j] = alpha_set[i];
		j++;
	}

	ZZ p = power(GR.p, GR.k);
	ZZ_p::init(p);
    ZZ_pX F = GR.Fps_poly;
    ZZ_pE::init(F);

	ZZ_pEX f;

	interpolate_for_GR(f, points, shares, GR.p, GR.k, GR.r);

	vec_ZZ_pE value;
	value.SetLength(m);

	for(int i=0; i<m; i++)
	{
		value[i] = eval(f, beta_set[i]);
	}

	return value;
}

vec_ZZ_pE scheme::create_one_shares(ZZ_pE a, long i){
	
	vec_ZZ_pE shares;
	shares.SetLength(m);

	vector<ZZ_pE> coeff(d);

	for(long i=0;i<d;i++) 
	{
		coeff[i] = random_ZZ_pE();
	}

	ZZ_pE z(0);
	for(long j=0 ; j<d ; j++) 
	{
		z = z + (coeff[j] * power(beta_set[i],j));
	}

	coeff[0] = a-z+coeff[0];

	for(long k=0 ; k<n ; k++) 
	{
		ZZ_pE x = alpha_set[k];
		ZZ_pE y(0);
		for(long j=0 ; j<d ; j++) 
		{
			y = y + (coeff[j] * power(x,j));
		}
		shares[k] = y;
	}
	return shares;

}

}