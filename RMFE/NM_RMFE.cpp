#include "NM_RMFE.hpp"


/*
    constructor of Non-malleable RMFE class 
    RMFE:  phase 1 : GR(p^k,s)^n1 -> GR(p^k,r)
           phase 2 : GR(p^k,r1)^n2 -> GR(p^k, d)
    
    Extended RMFE:  GR(p^k,s)^n2*GR(p^k,d) -> GR(p^k, n)

    !!!:user need to modify the coefficient as required in function NM_RMFE_GR_INIT2()
*/

NM_RMFE_GR::NM_RMFE_GR(ZZ p_, long s_, long k_, long D_, long n1_, long n2_,
            long r_, long d_, long n_)
        :p(p_), s(s_), k(k_), D(D_), n1(n1_), n2(n2_),
         r(r_), d(d_),n(n_)
{
    NM_RMFE_GR_INIT1();
    NM_EX_RMFE_GR_INIT2();
}


/*
    initialize the phase 1 of NM-RMFE
*/
void NM_RMFE_GR::NM_RMFE_GR_INIT1(){
    if(D*(n1-1)>r/s){
        throw invalid_argument("NM_RMFE_GR::NM_RMFE_GR_INIT: Note:D(n-1)<=r/s!");
    }
    if(power(p,s)<n1){
        throw invalid_argument("NM_RMFE_GR::NM_RMFE_GR_INIT: Note:(p,s)<n1 Lenstra Lemma!");
    }
    ZZ q1 = power(p,s);
    ZZ q2 = power(p,r);
    ZZ q3 = power(p,k);
    long q4;
    conv(q4,q1-ZZ(1));

    ZZ_p::init(p);
    ZZ_pX F;
    FindPrimitivePoly(F,p,s);

    ZZ_pX f;
    SetCoeff(f,q4,1);
    SetCoeff(f,0,-1);

    ZZ_pX g_;
    HenselLift(g_,f,F,p,k-1);

    fillIrred(g_,PrimitiveIrred1);
    ZZ_p::init(q3);
    ZZ_pE::init(g_);

    ZZ_pX H;
    SetCoeff(H,1,1);
    ZZ_pE a;
    conv(a,H);

    vec_ZZ_pE v1;
    v1.append(a);
    if(IsOne(a) && n1 == 2)
    {
        v1.append(ZZ_pE(0));
    }

    else{
        ZZ_pE a_tmp = a;
        for(int i=1;i<n1;i++){
            a_tmp = a_tmp * a;
            v1.append(a_tmp);
        }
    }
    fillInterpolation(v1,Interpolation1,s);
}


/*
    initialize the phase 2 of NM-RMFE ; user need to modify the coefficient as required
*/
void NM_RMFE_GR::NM_EX_RMFE_GR_INIT2(){
    if(D*(n2-1)>d/r){
        throw invalid_argument("NM_RMFE_GR::NM_RMFE_GR_INIT2: Note:D(n-1)<=r/s!");
    }
    if(power(p,r)<n2){
        throw invalid_argument("NM_RMFE_GR::NM_RMFE_GR_INIT:2 Note:(p,s)<n1 Lenstra Lemma!");
    }
    ZZ q1 = power(p,r);
    ZZ q2 = power(p,d);
    ZZ q3 = power(p,k);
    long q4;
    conv(q4,q1-ZZ(1));

    ZZ_p::init(q3);
    ZZ_pX g_;

    /////  user need to modify here
    SetCoeff(g_,12,1);
    SetCoeff(g_,7,1);
    SetCoeff(g_,4,1);
    SetCoeff(g_,3,1);
    SetCoeff(g_,0,1);
    /////

    fillIrred(g_,PrimitiveIrred2);
    ZZ_pE::init(g_);

    ZZ_pX H;
    SetCoeff(H,1,1);
    ZZ_pE a;
    conv(a,H);


    vec_ZZ_pE v1;
    v1.append(a);
    if(IsOne(a) && n2==2)
    {
        v1.append(ZZ_pE(0));
    }

    else{
        ZZ_pE a_tmp = a;
        for(int i=1;i<n2;i++){
            a_tmp = a_tmp * a;
            v1.append(a_tmp);
        }
    }
    fillInterpolation(v1,Interpolation2,r);

}


/*
    set up the input of NM-RMFE
*/
void NM_RMFE_GR::set_input(const vector<long>& Input){

    this->input.clear();
    
    long size = Input.size();
    if ( Input.size() != n1*n2*s )
        throw invalid_argument("NM_RMFE_GR::set_message: Input size error!");
    

    for(int i = 0;i<size;i++){
        this->input.push_back(Input[i]);

    }


}


/*
    the composite RMFE map phi
*/
void NM_RMFE_GR::NM_RMFE_GR_PHI(){

    vector<long> Input2; //D*n2*n1*s
    for(int i = 0; i < n2; i++){
        vector<long> Input1(input.begin()+i*n1*s,input.begin()+i*n1*s+n1*s);
        vector<long> result1 = NM_RMFE_GR_PHI1(Input1);
        vector<long> result1_pad;
        result1_pad = PadVectorToLength(result1,r);
        Input2.insert(Input2.end(),result1_pad.begin(),result1_pad.end());
    }
    NM_EX_RMFE_GR_PHI2(Input2);

}


/*
    the first RMFE map phi1
*/
vector<long> NM_RMFE_GR::NM_RMFE_GR_PHI1(vector<long>& input){

    if(input.size() != s * n1)
        throw invalid_argument("NM_RMFE_GR::NM_RMFE_GR_PHI1: input size error!");
    ZZ q3 = power(p,k);
    ZZ_p::init(q3);
    ZZ_pX F = long2ZZpX(PrimitiveIrred1);
    ZZ_pE::init(F);


    // interpolation points
    vec_ZZ_pE v1 ;
    v1.SetLength(n1);

    for(int i = 0; i < n1; i++){
        vector<long> coeff(Interpolation1.begin()+i*s, Interpolation1.begin()+i*s+s);
        v1[i] = long2ZZpE(coeff);
    }

    // input points

    vec_ZZ_pE v2;
    for (int i = 0; i < n1; i++){
        ZZ_pX f;
        ZZ_pE g;
        for (int j = 0; j < s; j++){
            SetCoeff(f,j,input[i*s+j]);
        }
        conv(g,f);
        v2.append(g);
    }
    ZZ_pEX h;
    interpolate_for_GR(h,v1,v2,p,k,s);


    vector<long> V;
    ZZpEX2long(h,V,s);
    return V;

}


/*
    the second RMFE map phi2
*/
void NM_RMFE_GR::NM_EX_RMFE_GR_PHI2(vector<long>& input){

    if(input.size() != r * n2)
        throw invalid_argument("NM_RMFE_GR::NM_RMFE_GR_PHI2: input size error!");
    

    ZZ q3 = power(p,k);
    ZZ_p::init(q3);
    ZZ_pX F = long2ZZpX(PrimitiveIrred2);
    ZZ_pE::init(F);


    //interpolation points
    vec_ZZ_pE v1 ;
    v1.SetLength(n2);
    for(int i = 0; i < n2; i++){
        vector<long> coeff(Interpolation2.begin()+i*r, Interpolation2.begin()+i*r+r);
        v1[i] = long2ZZpE(coeff);
    }

    // input points

    vec_ZZ_pE v2;
    for (int i = 0; i < n2; i++){
        ZZ_pX f;
        ZZ_pE g;
        for (int j = 0; j < r; j++){
            SetCoeff(f,j,input[i*r+j]);
        }
        conv(g,f);
        v2.append(g);
    }

    ZZ_pEX h;
    interpolate_for_GR(h,v1,v2,p,k,r);


    vector<long> h_long;
    ZZpEX2long(h,h_long,r);
    // ZZ_pX FF;

    // SetCoeff(FF,42,1);
    // SetCoeff(FF,7,1);
    // SetCoeff(FF,0,1);

    // ZZ_pE::init(FF);

    int nr = n/r;

    vector<long> res(n,0);

    for(int i=0;i<n2;i++){
        for(int j = 0; j<r;j++)
        {
            res[i + j * nr] = h_long[r * i + j];
        }
    }
    this->output = res;

}


/*
    the composite RMFE map psi
*/
vector<long> NM_RMFE_GR::NM_RMFE_GR_PSI(vector<long> Input){

    ZZ q3 = power(p,k);
    ZZ_p::init(q3);
    ZZ_pX F = long2ZZpX(PrimitiveIrred2);
    ZZ_pE::init(F);

    ZZ_pX H;
    SetCoeff(H,1,1);
    ZZ_pE beta;
    conv(beta,H);

    long nr = n/r;
    long dr = d/r;

    vector<long> input2;
    for(int i = 0 ; i < n / r ; i++){
        vector<long> tmp1(Input.begin()+i*nr,Input.begin()+i*nr+nr);
        vector<long> tmp2(tmp1.begin(),tmp1.begin() + dr);
        input2.insert(input2.end(),tmp2.begin(),tmp2.end());
    }

    
    vector<long> input3(input2.begin(), input2.begin() + dr * r);

    vector<vector<long>> tmp = splitVector(input3,dr);
    long numGroups = (input3.size() + dr - 1) / dr;

    //interpolation points
    vec_ZZ_pE v1 ;
    v1.SetLength(n2);
    for(int i = 0; i < n2; i++){
        vector<long> coeff(Interpolation2.begin()+i*r, Interpolation2.begin()+i*r+r);
        v1[i] = long2ZZpE(coeff);
    }

    vec_ZZ_pE v2;
    // for(int i = 0; i < n2; i++){
    //     v2.append(eval(V,v1[i]));
    // }

    for(int i=0;i<n2;i++){
        ZZ_pE as;
        clear(as);
        for(int j=0;j<dr;j++){
            ZZ_pE aa;
            clear(aa);
            for(int o=0;o<numGroups;o++){
                aa += tmp[o][j] * power(beta,o);

            }
            aa = aa *power(v1[i],j);
            as = as + aa;
        }
        v2.append(as);
    }

    //  phase 2
    q3 = power(p,k);
    ZZ_p::init(q3);
    F = long2ZZpX(PrimitiveIrred1);
    ZZ_pE::init(F);


    v1.SetLength(n1);

    for(int i = 0; i < n1; i++){
        vector<long> coeff(Interpolation1.begin()+i*s, Interpolation1.begin()+i*s+s);
        v1[i] = long2ZZpE(coeff);
    }

    vec_ZZ_pE v3 ;
    for(int i=0;i<n2;i++){
        ZZ_pE a = v2[i];
        ZZ_pX a_rep = rep(a);  // 获取 ZZ_pE 的内部表示 ZZ_pX

        ZZ_pEX result2;
        for(long j = 0; j <= deg(a_rep); j++){
            SetCoeff(result2,j,a_rep[j]);
        }

        for(int j=0;j<n1;j++){
            ZZ_pE c = eval(result2,v1[j]);
            v3.append(c);
        }

    }

    vector<long> m;
    VeczzpE2Veclong(v3, m, s);

    this->psi_output_1 = m;


    ////////////////////////////////////////////////

    F = long2ZZpX(PrimitiveIrred2);
    ZZ_pE::init(F);

    vector<vector<long>> tmp2;
    tmp2 = splitVector(Input,nr);
    numGroups = (Input.size() + nr - 1) / nr;

    //interpolation points


    v1.SetLength(n2);
    for(int i = 0; i < n2; i++){
        vector<long> coeff(Interpolation2.begin()+i*r, Interpolation2.begin()+i*r+r);
        v1[i] = long2ZZpE(coeff);
    }

    // for(int i = 0; i < n2; i++){
    //     v2.append(eval(V,v1[i]));
    // }

    vec_ZZ_pE v4;

    for(int i=0;i<n2;i++){
        ZZ_pE as;
        clear(as);
        for(int j=0;j<nr;j++){
            ZZ_pE aa;
            clear(aa);
            for(int o=0;o<numGroups;o++){
                aa += tmp2[o][j] * power(beta,o);

            }
            aa = aa *power(v1[i],j);
            as = as + aa;
        }
        v4.append(as);
    }
    //  phase 2
    q3 = power(p,k);
    ZZ_p::init(q3);
    F = long2ZZpX(PrimitiveIrred1);
    ZZ_pE::init(F);


    v1.SetLength(n1);

    for(int i = 0; i < n1; i++){
        vector<long> coeff(Interpolation1.begin()+i*s, Interpolation1.begin()+i*s+s);
        v1[i] = long2ZZpE(coeff);
    }

    vec_ZZ_pE v5;

    for(int i=0;i<n2;i++){
        ZZ_pE a = v4[i];
        ZZ_pX a_rep = rep(a);  // 获取 ZZ_pE 的内部表示 ZZ_pX

        ZZ_pEX result2;
        for(long j = 0; j <= deg(a_rep); j++){
            SetCoeff(result2,j,a_rep[j]);
        }

        for(int j=0;j<n1;j++){
            ZZ_pE c = eval(result2,v1[j]);
            v5.append(c);
        }

    }

    vector<long> m2;
    VeczzpE2Veclong(v5, m2, s);

    this->psi_output_2 = m2;

    bool flag = 1;

    int length = psi_output_1.size();

    if(length != psi_output_2.size())
        flag = 0;

    for(int i = 0 ; i < length;i++){
        if(psi_output_1[i]!= psi_output_2[i])
            flag = 0;
    }

    if(flag)
        return m;
    else {
        vector<long> a(1,0);
        return a;
    }

}

