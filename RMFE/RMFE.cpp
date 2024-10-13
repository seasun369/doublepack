#include "RMFE.hpp"


/*
    constructor of RMFE_GR class 
    RMFE:  phase 1 : GR(p^k,s)^n1 -> GR(p^k,r1)
           phase 2 : GR(p^k,r1)^n2 -> GR(p^k, r2)
    !!!:user need to modify the coefficient as required in function RMFE_GR_INIT2_cache()
*/
RMFE_GR::RMFE_GR(ZZ p_, long s_, long k_, long D_, long n1_, long n2_,
            long r1_, long r2_)
        :p(p_), s(s_), k(k_), D(D_), n1(n1_), n2(n2_),
         r1(r1_), r2(r2_)
{
    RMFE_GR_INIT1();
    RMFE_GR_INIT2_cache();
}


/*
    initialize the phase 1 of RMFE
*/
void RMFE_GR::RMFE_GR_INIT1(){
    if(D*(n1-1)>r1/s){
        throw invalid_argument("RMFE_GR::RMFE_GR_INIT1: Note:D(n-1)<=r/s!");
    }
    if(power(p,s) + 1<n1){
        throw invalid_argument("RMFE_GR::RMFE_GR_INIT1: Note:(p,s) + 1<n1 Lenstra Lemma!");
    }
    ZZ q1 = power(p,s);
    ZZ q2 = power(p,r1);
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
    if(k > 1)
        HenselLift(g_,f,F,p,k-1);
    else g_ = F;

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

    else if(IsZero(a) && n1 == 2){
        v1.append(ZZ_pE(1));
    }
    else if(power(p,s) + 1 == n1){
        v1.append(ZZ_pE(0));
    }
    else
        {
        ZZ_pE a_tmp = a;
        for(int i=1;i<n1;i++){
            a_tmp = a_tmp * a;
            v1.append(a_tmp);
        }
    }
    fillInterpolation(v1,Interpolation1,s);

}


/*
    initialize the phase 2 of RMFE
*/
void RMFE_GR::RMFE_GR_INIT2(){
    if(D*(n2-1)>r2/r1){
        throw invalid_argument("RMFE_GR::RMFE_GR_INIT2: Note:D(n-1)<=r/s!");
    }
    if(power(p,r1) + 1 < n2){
        throw invalid_argument("RMFE_GR::RMFE_GR_INIT2: Note:power(p,s)<n1 Lenstra Lemma!");
    }
    ZZ q1 = power(p,r1);
    ZZ q2 = power(p,r2);
    ZZ q3 = power(p,k);
    long q4;
    conv(q4,q1-ZZ(1));

    ZZ_p::init(p);
    ZZ_pX F;
    FindPrimitivePoly(F,p,r1);

    ZZ_pX f;
    SetCoeff(f,q4,1);
    SetCoeff(f,0,-1);

    ZZ_pX g_;
    HenselLift(g_,f,F,p,k-1);


    fillIrred(g_,PrimitiveIrred2);
    ZZ_p::init(q3);
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
    fillInterpolation(v1,Interpolation2,r1);
}


/*
    initialize the phase 2 of RMFE ; user need to modify the coefficient as required
*/
void RMFE_GR::RMFE_GR_INIT2_cache(){
    if(D*(n2-1)>r2/r1){
        throw invalid_argument("RMFE_GR::RMFE_GR_INIT2: Note:D(n-1)<=r/s!");
    }
    if(power(p,r1) + 1 < n2){
        throw invalid_argument("RMFE_GR::RMFE_GR_INIT2: Note:power(p,s)<n2 Lenstra Lemma!");
    }
    ZZ q1 = power(p,r1);
    ZZ q2 = power(p,r2);
    ZZ q3 = power(p,k);
    ZZ q5 = power(p,k-1);
    long q4;
    conv(q4,q1-ZZ(1));

    ZZ_p::init(q3);
    ZZ_pX g_;

    /////  user need to modify here
    SetCoeff(g_,5,1);
    SetCoeff(g_,2,1);
    SetCoeff(g_,0,1);
    /////

    fillIrred(g_,PrimitiveIrred2);
    ZZ_pE::init(g_);


    ZZ_pX H;
    SetCoeff(H,1,1);
    ZZ_pE aa;
    conv(aa,H);

    ZZ_pE a = power(aa,q5);

    vec_ZZ_pE v1;
    v1.append(a);
    if(IsOne(a) && n2==2)
    {
        v1.append(ZZ_pE(0));
    }

    else if(power(p,r1) + 1 == n2){
        ZZ_pE a_tmp = a;
        for(int i=1;i<n2-2;i++){
            a_tmp = a_tmp * a;
            v1.append(a_tmp);
        }
        v1.append(ZZ_pE(0));
    }
    else if(power(p,r1) == n2){
        ZZ_pE a_tmp = a;
        for(int i=1;i<n2-1;i++){
            a_tmp = a_tmp * a;
            v1.append(a_tmp);
        }
        v1.append(ZZ_pE(0));
    }
    else{
        ZZ_pE a_tmp = a;
        for(int i=1;i<n2;i++){
            a_tmp = a_tmp * a;
            v1.append(a_tmp);
        }
    }
    fillInterpolation(v1,Interpolation2,r1);

}


/*
    set up the input of RMFE
*/
void RMFE_GR::set_input(const vector<long>& Input){

    this->input.clear();
    
    long size = Input.size();
    if ( Input.size() != n1*n2*s )
        throw invalid_argument("RMFE_GR::set_message: Input size error!");
    

    for(int i = 0;i<size;i++){
        this->input.push_back(Input[i]);

    }

}


/*
    the composite RMFE map phi
*/
void RMFE_GR::RMFE_GR_PHI(){

    vector<long> Input2; //D*n2*n1*s
    for(int i = 0; i < n2; i++){
        vector<long> Input1(input.begin()+i*n1*s,input.begin()+i*n1*s+n1*s);
        vector<long> result1 = RMFE_GR_PHI1(Input1);
        vector<long> result1_pad;
        result1_pad = PadVectorToLength(result1,r1);
        Input2.insert(Input2.end(),result1_pad.begin(),result1_pad.end());
    }
    RMFE_GR_PHI2(Input2);

}


/*
    the first RMFE map phi1
*/
vector<long> RMFE_GR::RMFE_GR_PHI1(vector<long>& input){

    if(input.size() != s * n1)
        throw invalid_argument("RMFE_GR::RMFE_GR_PHI1: input size error!");
    ZZ q3 = power(p,k);
    ZZ_p::init(q3);
    ZZ_pX F = long2ZZpX(PrimitiveIrred1);
    ZZ_pE::init(F);

    ZZ_pEX h;
    if(power(p,s) + 1 == n1){
        // interpolation points
        vec_ZZ_pE v1 ;
        v1.SetLength(n1 - 1);

        for(int i = 0; i < n1 - 1; i++){
            vector<long> coeff(Interpolation1.begin()+i*s, Interpolation1.begin()+i*s+s);
            v1[i] = long2ZZpE(coeff);
        }
        // input points

        vec_ZZ_pE v2, v2_;
        for (int i = 0; i < n1; i++){
            ZZ_pX f;
            ZZ_pE g;
            for (int j = 0; j < s; j++){
                SetCoeff(f,j,input[i*s+j]);
            }
            conv(g,f);
            v2_.append(g);
        }

        for(int i = 0; i < n1 - 1; i++){
            v2.append(v2_[i]);
        }

        ZZ_pEX tmp;

        SetCoeff(tmp,n1-1,v2_[n1-1]);

        for(int i = 0; i < n1 - 1; i++){
            v2[i] = v2[i] - eval(tmp,v1[i]);
        }
        interpolate_for_GR(h,v1,v2,p,k,s);
        SetCoeff(h,n1-1,v2_[n1-1]);

    }
    else{
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

    interpolate_for_GR(h,v1,v2,p,k,s);
    }

    vector<long> V;
    ZZpEX2long(h,V,s);
    return V;

}

/*
    the second RMFE map phi2
*/
void RMFE_GR::RMFE_GR_PHI2(vector<long>& input){

    if(input.size() != r1 * n2)
        throw invalid_argument("RMFE_GR::RMFE_GR_PHI2: input size error!");
    

    ZZ q3 = power(p,k);
    ZZ_p::init(q3);
    ZZ_pX F = long2ZZpX(PrimitiveIrred2);
    ZZ_pE::init(F);

    ZZ_pEX h;
    if(power(p,r1) + 1 == n2){

        //interpolation points
        vec_ZZ_pE v1 ;
        v1.SetLength(n2 -1);
        for(int i = 0; i < n2 - 1; i++){
            vector<long> coeff(Interpolation2.begin()+i*r1, Interpolation2.begin()+i*r1+r1);
            v1[i] = long2ZZpE(coeff);
        }

        // input points

        vec_ZZ_pE v2_,v2;
        for (int i = 0; i < n2; i++){
            ZZ_pX f;
            ZZ_pE g;
            for (int j = 0; j < r1; j++){
                SetCoeff(f,j,input[i*r1+j]);
            }
            conv(g,f);
            v2_.append(g);
        }
        for(int i = 0; i < n2 - 1; i++){
            v2.append(v2_[i]);
        }

        ZZ_pEX tmp;

        SetCoeff(tmp,n2-1,v2_[n2-1]);

        for(int i = 0; i < n2 - 1; i++){
            v2[i] = v2[i] - eval(tmp,v1[i]);
        }
        interpolate_for_GR(h,v1,v2,p,k,r1);
        SetCoeff(h,n2-1,v2_[n2-1]);
    }
    else{
        //interpolation points
        vec_ZZ_pE v1 ;
        v1.SetLength(n2);
        for(int i = 0; i < n2; i++){
            vector<long> coeff(Interpolation2.begin()+i*r1, Interpolation2.begin()+i*r1+r1);
            v1[i] = long2ZZpE(coeff);
        }

        // input points

        vec_ZZ_pE v2;
        for (int i = 0; i < n2; i++){
            ZZ_pX f;
            ZZ_pE g;
            for (int j = 0; j < r1; j++){
                SetCoeff(f,j,input[i*r1+j]);
            }
            conv(g,f);
            v2.append(g);
        }

        interpolate_for_GR(h,v1,v2,p,k,r1);
    }

    vector<long> h_long;
    if(deg(h) >= 1)
        ZZpEX2long(h,h_long,r1);
    else{
        for(int i=0;i<r1 * n2;i++)
            h_long.push_back(0);
    }


    int r = r2/r1;

    vector<long> res(r2,0);

    for(int i=0;i<n2;i++){
        for(int j = 0; j<r1;j++)
        {
            res[i + j * r] = h_long[r1 * i + j];
        }
    }
    this->output = res;

}

/*
    the composite RMFE map psi
*/
vector<long> RMFE_GR::RMFE_GR_PSI(vector<long> Input){

    ZZ q3 = power(p,k);
    ZZ_p::init(q3);
    ZZ_pX F = long2ZZpX(PrimitiveIrred2);
    ZZ_pE::init(F);

    ZZ_pX H;
    SetCoeff(H,1,1);
    ZZ_pE beta;
    conv(beta,H);

    long r = r2/r1;


    vector<vector<long>> tmp = splitVector(Input,r);
    long numGroups = (Input.size() + r - 1) / r;

    //interpolation points
    vec_ZZ_pE v1 ;
    v1.SetLength(n2);
    for(int i = 0; i < n2; i++){
        vector<long> coeff(Interpolation2.begin()+i*r1, Interpolation2.begin()+i*r1+r1);
        v1[i] = long2ZZpE(coeff);
    }

    vec_ZZ_pE v2;
    // for(int i = 0; i < n2; i++){
    //     v2.append(eval(V,v1[i]));
    // }

    for(int i=0;i<n2;i++){
        ZZ_pE as;
        clear(as);
        for(int j=0;j<r;j++){
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
    return m;

}

void RMFE_GR::get_phi_kernel(vec_ZZ_pE& vec){
    ZZ_pX poly;
    SetCoeff(poly, 0, 0);  // 常数项 0
    SetCoeff(poly, 1, -1); // x^1 项系数 -1
    SetCoeff(poly, 2, 1);  // x^2 项系数 1

    // 将多项式转换为ZZ_pE类型
    ZZ_pE f = conv<ZZ_pE>(poly);

    for(int i =0; i<10; i++){

    ZZ_pX randomPoly;
    random(randomPoly, r2-2);  // randomPoly的次数上限为d-1

    // 将其转换为ZZ_pE类型
    ZZ_pE g = conv<ZZ_pE>(randomPoly);

    vec.append(f*g);

    }


}


