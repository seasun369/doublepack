#include "Inverse.hpp"


void Inv_matrix_1(mat_ZZ_p& A, ZZ_pE& a, long d){


    ZZ_pX poly = rep(a);
    A.SetDims(d,2*d-1);
    clear(A);

    for(int i = 0 ;i < d;i++){
        
        int index = i;
        for(int j = 0; j< d;j++){
            A[i][j + index] = coeff(poly, j);
        }

    }
}


void Inv_matrix_2(mat_ZZ_p& A, long d){

    A.SetDims(2*d-1,d);
    clear(A);

    for(int i=0;i<d;i++){
        A[i][i] = 1;
    } 

    ZZ_pX H;
    SetCoeff(H,1,1);
    ZZ_pE b;
    conv(b,H);

    for(int i=d;i<2*d-1;i++){
        ZZ_pE tmp = power(b,i);
        ZZ_pX poly = rep(tmp);

        for(int j=0; j<d;j++){
            A[i][j] = coeff(poly,j);
        }
    }

}

/*
    output the inverse of a ; s is the degree of galois ring
*/
ZZ_pE Inv(ZZ_pE a, long s){

    try{
    mat_ZZ_p A1,A2;
    vec_ZZ x,b;

    x.SetLength(s);
    b.SetLength(s);
    clear(b);
    b[0] = 1;

    Inv_matrix_1(A1,a,s);
    Inv_matrix_2(A2,s);

    ZZ d;

    mat_ZZ AB;
    conv(AB,A1*A2);   


    solve(d,x,AB,b);

    vec_ZZ_p x_p;
    conv(x_p,x);

    ZZ_pX poly;

    for(int i=0;i<s;i++){
        SetCoeff(poly,i,x_p[i]);
    }

    ZZ_pE a_inverse;
    conv(a_inverse, poly);

    ZZ_pE m = a * a_inverse;
    ZZ_p r;
    ZZ_pX m_poly = rep(m);
    r = coeff(m_poly, 0);

    ZZ_p r_inverse = inv(r);

    ZZ_pE r_pE;
    conv(r_pE,r_inverse);

    return a_inverse * r_pE;

    } catch (const NTL::InvModErrorObject& e){
        ZZ_pE a;
        clear(a);
        return a;
    }

}



ZZ_pX ZZpXMod(ZZ_pX& a, ZZ n){

    long d = deg(a);
    ZZ_pX res;
    for(int i=0;i<=d;i++){
        ZZ m = rep(coeff(a,i));
        m = m % n;
        ZZ_p m_p;
        conv(m_p,m);
        SetCoeff(res,i,m_p);
    }
    return res;
}


ZZ_pX ZZpXDiv(ZZ_pX& a, ZZ n){

    long d = deg(a);
    ZZ_pX res;
    for(int i=0;i<=d;i++){
        ZZ m = rep(coeff(a,i));
        m = m / n;
        ZZ_p m_p;
        conv(m_p,m);
        SetCoeff(res,i,m_p);
    }
    return res;
}


ZZ_pX ZZpXMul(ZZ_pX& a, ZZ n){

    long d = deg(a);
    ZZ_pX res;
    for(int i=0;i<=d;i++){
        ZZ m = rep(coeff(a,i));
        m = m * n;
        ZZ_p m_p;
        conv(m_p,m);
        SetCoeff(res,i,m_p);
    }
    return res;
}


Vec<ZZ_pX> Inv_ai(ZZ_pX a, ZZ p, long s, long k){

    Vec<ZZ_pX> res, v1;

    for(int i = 1; i <= k; i++){
        ZZ mod = power(p,i);
        ZZ_pX tmp = ZZpXMod(a,mod);
        v1.append(tmp);
    }
    res.append(v1[0]);

    for(int i = 1; i < k; i++){
        ZZ_pX b = v1[i] - v1[i-1];
        long d = deg(b);

        ZZ n = power(p,i);
        ZZ_pX c = ZZpXDiv(b,n);
        res.append(c);
    }
    return res;
}


ZZ_pE Inv2(ZZ_pE a, ZZ_pX F, ZZ p, long s, long k){

    ZZ_pX a_x = rep(a);

    Vec<ZZ_pX> a_i = Inv_ai(a_x, p,s,k);

    //cout<<"a_i:"<<a_i<<endl;

    Vec<ZZ_pX> u_i, v_i;

    ZZ_pX d,u_0,v_0;
    ZZ_p::init(p);
    XGCD(d,u_0,v_0,a_i[0],F);
    ZZ_pX t = u_0 / F;
    u_0 = u_0 % F;
    v_0 = v_0 + t * a_i[0];
    u_i.append(u_0);
    v_i.append(v_0);
    //cout<<"u_0:"<<u_0<<"v_0:"<<v_0<<endl;

    for(int i = 2; i < k + 1;i++){
        ZZ mod = power(p,i);

        ZZ_p::init(mod);
        int n = i-1;
        ZZ_pX c;


        // for(int j = 0; j < n; j++){
        //     for(int l = 0; l < n; l++){
        //         ZZ_pX tmp = a_i[j] * u_i[l];
        //         ZZ_pX tmp2 = ZZpXMul(tmp,power(p,j+l));
        //         c = c + tmp2;
        //     }
        // }

        for(int y=0;y<i-1;y++){
            for(int j=0;j<=y;j++){
                ZZ_pX tmp = a_i[j] * u_i[y-j];
                ZZ_pX tmp2 = ZZpXMul(tmp, power(p,y));
                c = c + tmp2;
            }
        }


        for(int j = 0; j < n; j++){
            ZZ_pX tmp = F * v_i[j];
            ZZ_pX tmp2 = ZZpXMul(tmp,power(p,j));
            c = c + tmp2;
        }

        ZZ_pX e;
        SetCoeff(e,0,ZZ_p(1));
        c = c - e;
        ZZ_pX g = ZZpXMod(c, power(p,i));
        ZZ_pX h = ZZpXDiv(g,power(p,n));
        //cout<<"g:"<<g<<endl;
        //cout<<"h:"<<h<<endl;
        ZZ_pX o;


        ZZ_pX d,u,v;
        ZZ_p::init(p);
        for(int j=0;j<n;j++){
            ZZ_pX r = u_i[j] * a_i[n-j];
            o = o + r;
        }

        XGCD(d,u,v,a_i[0],F);
        u = u * (-1 * (h + o));
        ZZ_pX t = u / F;
        u = u % F;
        v = v * (-1 * (h + o));
        //cout<<"t"<<t * a_i[0]<<endl;
        
        //cout<<"o:"<<o<<"h+o:"<<(-1 * (h + o))<<endl;
        v = v + t * a_i[0];
        //cout<<"u:"<<u<<"v:"<<v<<endl;
        u_i.append(u);
        v_i.append(v);
    }

    ZZ_p::init(power(p,k));
    ZZ_pX res;
    //cout<<"u_i:"<<u_i<<"v_i:"<<v_i<<endl;
        for(int i=0;i<k;i++){
        ZZ_pX tmp = ZZpXMul(u_i[i], power(p,i));
        res = res + tmp;
    }
    ZZ_pE result;
    conv(result, res);

    return result;
}
