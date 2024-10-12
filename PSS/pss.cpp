#include "pss.h"
#include "../RMFE/utils.hpp"
#include <iostream>
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
	std::cout << "Entering create_shares" << std::endl;
	std::cout << "Input vector a: ";
	for (long i = 0; i < a.length(); i++) {
		std::cout << a[i] << " ";
	}
	std::cout << std::endl;

	vec_ZZ_pE value;
	value.SetLength(n);
	long _m = a.length();
	if(m != _m) Error("PSS: vector length mismatch");

	ZZ p = power(GR.p, GR.k);
	ZZ_p::init(p);
    ZZ_pX F = GR.Fps_poly;
    ZZ_pE::init(F);

	ZZ_pEX f;

	
	std::cout << "Beta set used for interpolation:" << std::endl;
	for(int i=0; i<m; i++)
	{
		std::cout << "beta_set[" << i << "] = " << beta_set[i] << std::endl;
	}

	std::cout << "Interpolating in create_shares" << std::endl;
	interpolate_for_GR(f, beta_set, a, GR.p, GR.k, GR.r);

	// 添加这行来打印初始的多项式 f
	std::cout << "Initial polynomial f for secret sharing: " << f << std::endl;

	// 检查初始秘密是否匹配 eval(f, beta_set[i])
	std::cout << "Checking if initial secret matches eval(f, beta_set[i]):" << std::endl;
	for(int i=0; i<m; i++)
	{
		ZZ_pE result = eval(f, beta_set[i]);
		std::cout << "beta_set[" << i << "] = " << beta_set[i] << std::endl;
		std::cout << "Initial secret: " << a[i] << std::endl;
		std::cout << "eval(f, beta_set[" << i << "]) = " << result << std::endl;
		if (a[i] == result) {
			std::cout << "Match for index " << i << std::endl;
		} else {
			std::cout << "Mismatch for index " << i << std::endl;
		}
	}

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
	h.rep = ret;

	// 添加这段代码来检查 h(beta_i) = 0
	std::cout << "Checking if h(beta_i) = 0 for all i:" << std::endl;
	for(int i=0; i<m; i++)
	{
		ZZ_pE h_beta_i = eval(h, beta_set[i]);
		std::cout << "h(beta_" << i << ") = " << h_beta_i << std::endl;
		if (IsZero(h_beta_i)) {
			std::cout << "h(beta_" << i << ") is zero" << std::endl;
		} else {
			std::cout << "h(beta_" << i << ") is NOT zero" << std::endl;
		}
	}

	ZZ_pEX poly = g*h + f;

	// 添加这段代码来检查 poly(beta_i) = f(beta_i)
	std::cout << "Checking if poly(beta_i) = f(beta_i) for all i:" << std::endl;
	for(int i=0; i<m; i++)
	{
		ZZ_pE poly_beta_i = eval(poly, beta_set[i]);
		ZZ_pE f_beta_i = eval(f, beta_set[i]);
		std::cout << "poly(beta_" << i << ") = " << poly_beta_i << std::endl;
		std::cout << "f(beta_" << i << ") = " << f_beta_i << std::endl;
		if (poly_beta_i == f_beta_i) {
			std::cout << "poly(beta_" << i << ") = f(beta_" << i << ")" << std::endl;
		} else {
			std::cout << "poly(beta_" << i << ") != f(beta_" << i << ")" << std::endl;
		}
	}

	// 添加这行来打印最终的多项式 poly
	std::cout << "Final polynomial poly before share creation: " << poly << std::endl;

	std::cout << "Evaluating shares:" << std::endl;
	for(int i=0; i<n; i++)
	{
		ZZ_pE result = eval(poly, alpha_set[i]);
		std::cout << "Evaluating poly(" << alpha_set[i] << ") = " << result << std::endl;
		value[i] = result;
		std::cout << "Share " << i << ": " << value[i] << std::endl;
	}

	std::cout << "Polynomial f after interpolation: " << f << std::endl;

	std::cout << "Final polynomial poly: " << poly << std::endl;

	std::cout << "Exiting create_shares" << std::endl;
	return value;

}

vec_ZZ_pE scheme::create_shares(vector<ZZ_pE> b){
	std::cout << "Entering create_shares" << std::endl;
	std::cout << "Input vector a: ";

	vec_ZZ_pE a;
	a.SetLength(b.size());

	for (long i = 0; i < b.size(); i++) {
		a[i]=b[i];
		std::cout << a[i] << " ";
	}
	std::cout << std::endl;

	vec_ZZ_pE value;
	value.SetLength(n);
	long _m = a.length();
	if(m != _m) Error("PSS: vector length mismatch");

	ZZ p = power(GR.p, GR.k);
	ZZ_p::init(p);
    ZZ_pX F = GR.Fps_poly;
    ZZ_pE::init(F);

	ZZ_pEX f;

	
	std::cout << "Beta set used for interpolation:" << std::endl;
	for(int i=0; i<m; i++)
	{
		std::cout << "beta_set[" << i << "] = " << beta_set[i] << std::endl;
	}

	std::cout << "Interpolating in create_shares" << std::endl;
	interpolate_for_GR(f, beta_set, a, GR.p, GR.k, GR.r);

	// 添加这行来打印初始的多项式 f
	std::cout << "Initial polynomial f for secret sharing: " << f << std::endl;

	// 检查初始秘密是否匹配 eval(f, beta_set[i])
	std::cout << "Checking if initial secret matches eval(f, beta_set[i]):" << std::endl;
	for(int i=0; i<m; i++)
	{
		ZZ_pE result = eval(f, beta_set[i]);
		std::cout << "beta_set[" << i << "] = " << beta_set[i] << std::endl;
		std::cout << "Initial secret: " << a[i] << std::endl;
		std::cout << "eval(f, beta_set[" << i << "]) = " << result << std::endl;
		if (a[i] == result) {
			std::cout << "Match for index " << i << std::endl;
		} else {
			std::cout << "Mismatch for index " << i << std::endl;
		}
	}

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
	h.rep = ret;

	// 添加这段代码来检查 h(beta_i) = 0
	std::cout << "Checking if h(beta_i) = 0 for all i:" << std::endl;
	for(int i=0; i<m; i++)
	{
		ZZ_pE h_beta_i = eval(h, beta_set[i]);
		std::cout << "h(beta_" << i << ") = " << h_beta_i << std::endl;
		if (IsZero(h_beta_i)) {
			std::cout << "h(beta_" << i << ") is zero" << std::endl;
		} else {
			std::cout << "h(beta_" << i << ") is NOT zero" << std::endl;
		}
	}

	ZZ_pEX poly = g*h + f;

	// 添加这段代码来检查 poly(beta_i) = f(beta_i)
	std::cout << "Checking if poly(beta_i) = f(beta_i) for all i:" << std::endl;
	for(int i=0; i<m; i++)
	{
		ZZ_pE poly_beta_i = eval(poly, beta_set[i]);
		ZZ_pE f_beta_i = eval(f, beta_set[i]);
		std::cout << "poly(beta_" << i << ") = " << poly_beta_i << std::endl;
		std::cout << "f(beta_" << i << ") = " << f_beta_i << std::endl;
		if (poly_beta_i == f_beta_i) {
			std::cout << "poly(beta_" << i << ") = f(beta_" << i << ")" << std::endl;
		} else {
			std::cout << "poly(beta_" << i << ") != f(beta_" << i << ")" << std::endl;
		}
	}

	// 添加这行来打印最终的多项式 poly
	std::cout << "Final polynomial poly before share creation: " << poly << std::endl;

	std::cout << "Evaluating shares:" << std::endl;
	for(int i=0; i<n; i++)
	{
		ZZ_pE result = eval(poly, alpha_set[i]);
		std::cout << "Evaluating poly(" << alpha_set[i] << ") = " << result << std::endl;
		value[i] = result;
		std::cout << "Share " << i << ": " << value[i] << std::endl;
	}

	std::cout << "Polynomial f after interpolation: " << f << std::endl;

	std::cout << "Final polynomial poly: " << poly << std::endl;

	std::cout << "Exiting create_shares" << std::endl;
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
	std::cout << "Entering packed_reconstruct_shares" << std::endl;
	std::cout << "Party size: " << party.size() << std::endl;
	std::cout << "Shares size: " << shares.length() << std::endl;

	long num = party.size();

	if(num != shares.length()) {
		std::cout << "Error: size mismatch" << std::endl;
		Error("PSS: reconstruct vector length mismatch");
	}

	vec_ZZ_pE points;
	points.SetLength(num);
	std::cout << "Points vector created with size: " << points.length() << std::endl;

	long j = 0;
	for(auto i : party)
	{
		std::cout << "Processing party member: " << i << std::endl;
		if(i <= 0 || i > alpha_set.length()) {
			std::cout << "Error: Invalid party member index" << std::endl;
			Error("PSS: Invalid party member index");
		}
		points[j] = alpha_set[i-1];  // 使用 i-1 作为索引
		j++;
	}

	std::cout << "Points after assignment:" << std::endl;
	for (long i = 0; i < points.length(); i++) {
		std::cout << "Point " << i << ": " << points[i] << std::endl;
	}

	std::cout << "Shares:" << std::endl;
	for (long i = 0; i < shares.length(); i++) {
		std::cout << "Share " << i << ": " << shares[i] << std::endl;
	}

	ZZ p = power(GR.p, GR.k);
	ZZ_p::init(p);
    ZZ_pX F = GR.Fps_poly;
    ZZ_pE::init(F);

	ZZ_pEX f;

	std::cout << "Interpolating in packed_reconstruct_shares" << std::endl;
	interpolate_for_GR(f, points, shares, GR.p, GR.k, GR.r);

	std::cout << "Interpolated polynomial f: " << f << std::endl;

	vec_ZZ_pE value;
	value.SetLength(m);

	std::cout << "Beta set used for reconstruction:" << std::endl;
	for(int i=0; i<m; i++)
	{
		std::cout << "beta_set[" << i << "] = " << beta_set[i] << std::endl;
	}

	for(int i=0; i<m; i++)
	{
		ZZ_pE result = eval(f, beta_set[i]);
		std::cout << "Evaluating f(" << beta_set[i] << ") = " << result << std::endl;
		value[i] = result;
		std::cout << "Reconstructed value " << i << ": " << value[i] << std::endl;
	}

	std::cout << "Exiting packed_reconstruct_shares" << std::endl;
	return value;
}

vec_ZZ_pE scheme::create_one_shares(ZZ_pE a, long i){
    
    vec_ZZ_pE shares;
    shares.SetLength(n);

    vector<ZZ_pE> coeff(d);

    for(long j=0;j<d;j++) 
    {
        coeff[j] = random_ZZ_pE();
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

ZZ_pE scheme::reconstruct_one_shares(vec_ZZ_pE shares, long i){

	if(shares.length() != n) Error("PSS: reconstruct vector length mismatch");
	
	ZZ p = power(GR.p, GR.k);
	ZZ_p::init(p);
    ZZ_pX F = GR.Fps_poly;
    ZZ_pE::init(F);

	ZZ_pEX f;

	interpolate_for_GR(f, alpha_set, shares, GR.p, GR.k, GR.r);

	return eval(f, beta_set[i]);
}

vec_ZZ_pE scheme::create_shares_with_points(vector<ZZ_pE> a, vector<ZZ_pE>b){
    if (a.size() != b.size() || a.size() == 0) {
        throw std::invalid_argument("Input vectors must have the same non-zero size");
    }

    std::cout << "create_shares_with_points: Input sizes: " << a.size() << std::endl;

    ZZ p = power(GR.p, GR.k);
    ZZ_p::init(p);
    ZZ_pX F = GR.Fps_poly;
    ZZ_pE::init(F);

    ZZ_pEX f;
    vec_ZZ_pE value;
    value.SetLength(n);

    bool success = false;
    int max_attempts = 10;  // 最大尝试次数
    
    for (int attempt = 0; attempt < max_attempts && !success; ++attempt) {
        try {
            vec_ZZ_pE aa, bb;
            aa.SetLength(a.size());
            bb.SetLength(b.size());

            for (size_t i = 0; i < a.size(); i++) {
                aa[i] = (attempt == 0) ? a[i] : random_ZZ_pE();
                bb[i] = (attempt == 0) ? b[i] : random_ZZ_pE();
                std::cout << "Attempt " << attempt << ", Point " << i << ": (" << aa[i] << ", " << bb[i] << ")" << std::endl;
            }

            std::cout << "Calling interpolate_for_GR..." << std::endl;
            interpolate_for_GR(f, aa, bb, GR.p, GR.k, GR.r);
            std::cout << "Interpolation successful. Polynomial: " << f << std::endl;

            for(int i=0; i<n; i++) {
                value[i] = eval(f, alpha_set[i]);
                std::cout << "f(" << alpha_set[i] << ") = " << value[i] << std::endl;
            }

            success = true;
        } catch (const std::exception& e) {
            std::cerr << "Attempt " << attempt << " failed: " << e.what() << std::endl;
            if (attempt == max_attempts - 1) {
                throw std::runtime_error("Failed to find invertible points after max attempts");
            }
        }
    }

    return value;
}

}