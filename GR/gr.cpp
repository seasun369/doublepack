#include "gr.h"

gr::gr(ZZ p_, long k_, long r_) : p(p_), k(k_), r(r_){
    gr_init_basic_var();
    get_set_T(set_T);
}

void gr::gr_init_basic_var(){
    //gr_prime = conv<ZZ>(2); // p 
    //positive_integer = 64;//k
    ZZ exponent_pk = power(p, k); // p^k
    //gr_degree = 4;//4;//d
    ZZ pd = power(p, r) - 1;// p^d - 1

    ZZ_p::init(p); 
    

    Fps_poly.SetLength(r+1);
    
    int direct_set_poly=1;
    if(direct_set_poly)
    {
        if(r==4)
        {
            SetCoeff(Fps_poly, 0, conv<ZZ_p>(1));//x
            SetCoeff(Fps_poly, 3, conv<ZZ_p>(1));
            SetCoeff(Fps_poly, 4, conv<ZZ_p>(1));  
        }
        else if(r==5)//[1 0 1 0 0 1]
        {
            SetCoeff(Fps_poly, 0, conv<ZZ_p>(1));//x
            SetCoeff(Fps_poly, 2, conv<ZZ_p>(1));         
            SetCoeff(Fps_poly, 5, conv<ZZ_p>(1)); 
        }
        else if(r==6)//[1 1 0 0 1 1 1]
        {
            SetCoeff(Fps_poly, 0, conv<ZZ_p>(1));//x
            SetCoeff(Fps_poly, 1, conv<ZZ_p>(1));
            SetCoeff(Fps_poly, 4, conv<ZZ_p>(1));            
            SetCoeff(Fps_poly, 5, conv<ZZ_p>(1)); 
            SetCoeff(Fps_poly, 6, conv<ZZ_p>(1)); 
        }
        else if(r==10)//[1 0 0 1 1 0 0 0 1 0 1]
        {
            SetCoeff(Fps_poly, 0, conv<ZZ_p>(1));//x
            SetCoeff(Fps_poly, 3, conv<ZZ_p>(1));
            SetCoeff(Fps_poly, 4, conv<ZZ_p>(1));            
            SetCoeff(Fps_poly, 8, conv<ZZ_p>(1)); 
            SetCoeff(Fps_poly, 10, conv<ZZ_p>(1)); 
        }
        else if(r==20)
        {
            SetCoeff(Fps_poly, 0, conv<ZZ_p>(1));
            SetCoeff(Fps_poly, 8, conv<ZZ_p>(1));
            SetCoeff(Fps_poly, 20, conv<ZZ_p>(1));
        }
    }
    /*else
    {
        
    
        //cout << "all_prime_factor of p^d-1 : ";
        vec_long pr_factor = find_all_prime_factor(pr);
        cout << pr_factor << endl; 
        cout << endl;
        while(1) 
        {
            ZZ_pX pr_poly;
            //cout<<"e"<<endl;
            pr_poly.SetLength(conv<long>(pr + 1));//p^r
            cout<<"pr_poly.rep.length() "<<pr_poly.rep.length()<<endl;
            
            SetCoeff(pr_poly, pr_poly.rep.length() - 1, 1);
            SetCoeff(pr_poly, 0, -1);
            cout << "pr_poly "<<pr_poly << endl;   

            vec_pair_ZZ_pX_long poly_factor;
            CanZass(poly_factor, pr_poly);
            //cout<<"e"<<endl;
            cout << "all factors of poly which deg is p^r-1:" << endl;
            cout << "poly_factor: " <<poly_factor << endl;
            
            for(auto p_r : poly_factor) 
            {
                cout << "p_r: " <<p_r << endl;//int i;cin>>i;
                ZZ_pX tmp_poly = p_r.a;
                cout << "tmp_poly: " <<tmp_poly << endl;
                cout << "LeadCoeff(tmp_poly): " <<LeadCoeff(tmp_poly) << endl;
                cout << "tmp_poly.rep.length(): " <<tmp_poly.rep.length() << endl;
                cout << "DetIrredTest(tmp_poly): " <<DetIrredTest(tmp_poly) << endl;
                //int i;cin>>i;
                if(LeadCoeff(tmp_poly) == 1 && tmp_poly.rep.length() == gr_degree+1 && DetIrredTest(tmp_poly) == 1) 
                {
                    if(test_primitive(tmp_poly, pr_factor, pr)) 
                    {
                        Fps_poly = tmp_poly;
                        break;
                    }
                }
            }
            if(Fps_poly.rep.length() != 0) break;
            //initiate_pr_relate(degree_r, pr, pr_factor,gr_prime,gr_degree);          
        }    
    }
    */




    ZZ_pE::init(Fps_poly);
    ZZ_p::init(exponent_pk);
    //cout <<"gr_prime: " << p << endl;  
    //cout <<"positive_integer: " << k << endl;  
    //cout <<"gr_degree: " << d << endl;
    //cout <<"p^k: " << exponent_pk << endl;  
    //cout <<"Fps_poly: " << Fps_poly << endl;
}

void gr::get_set_T(vec_ZZ_pE& set_T){

    //ZZ_pE primtive_element = FindPrimitiveElement(p, k, d);
    //cout <<" get the set T of GR " << endl;
    ZZ pd = power(p, r) - 1;// p^d - 1

    ZZ_pX tmp_poly;
    tmp_poly.SetLength(1);
    set_T.append(conv<ZZ_pE>(tmp_poly));
    //cout <<"tmp_poly: " << tmp_poly << endl;

    tmp_poly.SetLength(2);
    SetCoeff(tmp_poly, 1, 1);
    //cout <<"tmp_poly: " << tmp_poly << endl;

    ZZ_pE tmp_ring(1);
    for(int i = 0; i < pd; i++) 
    {
        set_T.append(tmp_ring);
        tmp_ring *= conv<ZZ_pE>(tmp_poly);
    }
    //cout << "set T: " << endl;
    //cout << set_T << endl;
}