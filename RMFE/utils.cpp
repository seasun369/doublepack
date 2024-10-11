#include "utils.hpp"

/*  
    Find all factors of long n
    eg: n = 8,return [1,2,4,8]
*/
vector<long> FindFactor(long n) {
    vector<long> factors;

    for (long i = 1; i <= n; ++i) {
        if (n % i == 0) {
            factors.push_back(i);
        }
    }

    return factors;
}


/*
    vector<long> -> ZZ_pX
*/
ZZ_pX long2ZZpX(const vector<long>& coefficients) {
    ZZ_pX poly;
    for (size_t i = 0; i < coefficients.size(); i++) {
        SetCoeff(poly, i, to_ZZ_p(coefficients[i]));
    }

    return poly;
}


/*
    vector<long> -> ZZ_pE
*/
ZZ_pE long2ZZpE(const vector<long>& coefficients) {
    ZZ_pX poly;
    for (size_t i = 0; i < coefficients.size(); i++) {
        SetCoeff(poly, i, to_ZZ_p(coefficients[i]));
    }

    ZZ_pE result;
    conv(result, poly);

    return result;
}


/*
    judge if n is power of 2
*/
bool isPowerOfTwo(long n) {
    return n > 0 && (n & (n - 1)) == 0;
}


/*
    vec_ZZ_pE -> vector<long> ; s is the degree of GR
*/
void VeczzpE2Veclong(const vec_ZZ_pE& v1, vector<long>& m, long s) {
    m.clear(); 

    for (long i = 0; i < v1.length(); i++) {
        const ZZ_pE element = v1[i];
        const ZZ_pX polyRep = rep(element);


        for (long j = 0; j < s; j++) {
            long c;
            conv(c,coeff(polyRep, j));
            m.push_back(c);
        }
    }
}


/*
    vector<long> -> string
*/
string Veclong2String(const std::vector<long>& vec) {
    stringstream ss;
    for (size_t i = 0; i < vec.size(); ++i) {

        ss << vec[i]; 
    }
    return ss.str(); 
}


/*
    ZZ_pE -> vector<long> ; s is the degree of GR
*/
void ZzpE2Veclong(const ZZ_pE& F, vector<long>& m, long s) {
    m.clear();

    ZZ_pX polyRep = rep(F);
    for (long i = 0; i < s; i++) {
        long c;
        conv(c,(coeff(polyRep, i)));
        m.push_back(c);
    }
}


/*
    vec_ZZ_pE -> vector<string> ; s is the degree of GR
*/
void VeczzpE2Vecstring(const vec_ZZ_pE& v1, vector<string>& m, long s){
    m.clear();

    for (long i = 0; i < v1.length(); i++){
        vector<long> b;
        ZzpE2Veclong(v1[i],b,s);
        string a = Veclong2String(b);
        m.push_back(a);
    }
}


/*
    ZZ_pX -> vector<long> 
*/
void ZZpX2long(const ZZ_pX& F, vector<long>& Irred) {
    Irred.clear();

    for (long i = 0; i <= deg(F); i++) {
        long c;
        conv(c,(coeff(F, i)));
        Irred.push_back(c);
    }
}


/*
    ZZ_pEX -> vector<long> ; s is the degree of GR
*/
void ZZpEX2long(const ZZ_pEX& v1, vector<long>& m, long s) {
    m.clear(); // 清空原有的数据

    for (long i = 0; i <= deg(v1); i++) {
        const ZZ_pE element = coeff(v1,i);
        const ZZ_pX polyRep = rep(element);


        for (long j = 0; j < s; j++) {
            long c;
            conv(c,coeff(polyRep, j));
            m.push_back(c);
        }
    }
}


/*
    judge if vec is all non-zero
*/
bool allNonZero(const vec_ZZ_pE& vec) {
    for (const auto& elem : vec) {
        if (elem == 0) {
            return false;
        }
    }
    return true;
}


/* 
    output the closest power of 2 that greater than n
*/
long nextPowerOf2(long n) {
    long count = 0;

    if (n && !(n & (n - 1)))
        return n;

    while (n != 0) {
        n >>= 1;
        count++;
    }

    return 1 << count;
}


/*
    print vector<long>
*/
void print(vector<long>& vec){
    cout<<"[";
    for (const auto& element : vec) {
        cout << element << " ";
    }
    cout<<"]"<< endl;
}



// Function to pad a vector to a specified length with zeros
vector<long> PadVectorToLength(const vector<long>& input, size_t targetLength) {
    // Create a vector with the target length, initialized with zeros
    vector<long> paddedVector(targetLength, 0);

    // Determine the number of elements to copy from the input vector
    size_t copyLength = min(input.size(), targetLength);
    // Copy elements from the input vector to the padded vector
    copy(input.begin(), input.begin() + copyLength, paddedVector.begin());

    return paddedVector;
}


// Function to trim trailing zeros from a vector
vector<long> TrimVector(const std::vector<long>& input) {
    size_t lastNonZeroIndex = input.size();

    // Find the index of the last non-zero element
    for (size_t i = input.size(); i > 0; --i) {
        if (input[i - 1] != 0) {
            lastNonZeroIndex = i - 1;
            break;
        }
    }

    // Create a trimmed vector from the input vector up to the last non-zero element
    vector<long> trimmedVector(input.begin(), input.begin() + lastNonZeroIndex + 1);

    return trimmedVector;
}



// Function to split an input vector into segments, pad each segment to a specified length, and concatenate the results
vector<long> SplitAndPadVector(const vector<long>& input, long segmentLength, long numSegments) {
    // Reserve space in the result vector for efficiency
    vector<long> result;
    result.reserve(segmentLength * numSegments);

    // Calculate the number of elements in each segment
    long n = input.size() / numSegments;

    // Iterate over the number of segments
    for (size_t i = 0; i < numSegments; ++i) {
        // Determine the start and end of the current segment
        auto segmentStart = input.begin() + i * n;
        auto segmentEnd = (i + 1) * n <= input.size() ? segmentStart + n : input.end();

        // Create a segment from the input vector
        vector<long> segment(segmentStart, segmentEnd);
        // Pad the segment to the specified length
        vector<long> paddedSegment = PadVectorToLength(segment, segmentLength);
        // Append the padded segment to the result vector
        result.insert(result.end(), paddedSegment.begin(), paddedSegment.end());
    }

    return result;
}


/*
    vector<long> -> ZZ_pEX ; s is the degree of GR
*/
void Long2ZZpEX(const vector<long>& result, ZZ_pEX& V, long s) {

    long index = 0; 

    long n = result.size() / s;
    for (long j = 0; j < n; j++) {
        ZZ_pX polyRep;
        for (long k = 0; k < s; k++) {

            SetCoeff(polyRep, k, to_ZZ_p(result[index++]));
        }
        ZZ_pE a;
        conv(a,polyRep);

        SetCoeff(V, j, a);
    }

}


void Long2ZZpEX2(const vector<long>& result, ZZ_pEX& V, long s,long n) {

    long index = 0; 

    for (long j = 0; j < n; j++) {
        ZZ_pX polyRep;
        for (long k = 0; k < s; k++) {

            SetCoeff(polyRep, k, to_ZZ_p(result[index++]));
        }
        ZZ_pE a;
        conv(a,polyRep);

        SetCoeff(V, j, a);
    }

}



int nearestPerfectSquare(int num) {
    int root = sqrt(num);
    
    while (root * root <= num) {
        root++;
    }

    return root * root;
}


/*
    irreducible polynomial ZZ_pX F -> vector<long>
*/
void fillIrred(const ZZ_pX& F, vector<long>& Irred) {
    Irred.clear();

    for (long i = 0; i <= deg(F); i++) {
        long c;
        conv(c,(coeff(F, i)));
        Irred.push_back(c);
    }
}


/*
    interpolation points vec_ZZ_pE -> vector<long>
*/
void fillInterpolation(const vec_ZZ_pE& v1, vector<long>& Interpolation, long s) {
    Interpolation.clear(); 

    for (long i = 0; i < v1.length(); i++) {
        const ZZ_pE element = v1[i];
        const ZZ_pX polyRep = rep(element);


        for (long j = 0; j < s; j++) {
            long c;
            conv(c,coeff(polyRep, j));
            Interpolation.push_back(c);
        }
    }
}


/*
    Find the primitive polynomial of degree n in Zp[x]
*/
void FindPrimitivePoly(ZZ_pX& g, ZZ p, long n){
    ZZ q2 = power(p,n);
    long q1;
    conv(q1,q2-ZZ(1));

    vector<long> factors = FindFactor(q1);

    bool flag = 0;
    while(!flag){
        ZZ_pX F;
        BuildIrred(F,n);

        ZZ_pX f;
        SetCoeff(f,0,-1);
        flag = 1;
        for(int i=0;i<factors.size()-1;i++){
            if(factors[i]<=n)
                continue;
            
            SetCoeff(f,factors[i],1);
            ZZ_pX q,r;
            DivRem(q,r,f,F);
            if(r == 0){
                flag = 0;
                break;
            }
        }
        if(flag == 1){
            g = F;
        }

    }
    return;
}


/*
    interpolate over galois ring GR(p^l,s) satisfying f(a) = b for all a
*/
void interpolate_for_GR(ZZ_pEX& f, const vec_ZZ_pE& a, const vec_ZZ_pE& b, ZZ p, long l, long s)
{
   long m = a.length();
   if (b.length() != m) LogicError("interpolate: vector length mismatch");

   if (m == 0) {
      clear(f);
      return;
   }
   vec_ZZ_pE prod;
   prod = a;
   ZZ_pE t1, t2;

   long k, i;

   vec_ZZ_pE res;
   res.SetLength(m);

   for (k = 0; k < m; k++) {
      const ZZ_pE& aa = a[k];

      set(t1);
      for (i = k-1; i >= 0; i--) {
         mul(t1, t1, aa);
         add(t1, t1, prod[i]);
      }

      clear(t2);
      for (i = k-1; i >= 0; i--) {
         mul(t2, t2, aa);
         add(t2, t2, res[i]);
      }

      std::cout << "k = " << k << ", aa = " << aa << std::endl;
      std::cout << "Before inversion: t1 = " << t1 << std::endl;
      
      try {
         t1 = Inv(t1, s);
         std::cout << "After inversion: t1 = " << t1 << std::endl;
      } catch (const std::exception& e) {
         std::cerr << "Error in Inv: " << e.what() << std::endl;
         std::cerr << "k = " << k << ", aa = " << aa << ", t1 = " << t1 << std::endl;
         throw;
      }

      sub(t2, b[k], t2);
      mul(t1, t1, t2);

      for (i = 0; i < k; i++) {
         mul(t2, prod[i], t1);
         add(res[i], res[i], t2);
      }

      res[k] = t1;
      if (k < m-1) {
         if (k == 0)
            prod[0]=-prod[0];
         else {
            t1=-a[k];
            add(prod[k], t1, prod[k-1]);
            for (i = k-1; i >= 1; i--) {
               mul(t2, prod[i], t1);
               add(prod[i], t2, prod[i-1]);
            }
            mul(prod[0], prod[0], t1);
         }
      }
   }

   while (m > 0 && IsZero(res[m-1])) m--;
   res.SetLength(m);
   f.rep = res;
}


// 函数：将向量拆分成固定长度的子向量，最后一组不足时补0
vector<vector<long>> splitVector(const vector<long>& input, int groupSize) {
    vector<vector<long>> result;
    int totalSize = input.size();
    int numGroups = (totalSize + groupSize - 1) / groupSize; // 计算总组数

    for (int i = 0; i < numGroups; ++i) {
        vector<long> group;
        for (int j = 0; j < groupSize; ++j) {
            int index = i * groupSize + j;
            if (index < totalSize) {
                group.push_back(input[index]);
            } else {
                group.push_back(0); // 补0
            }
        }
        result.push_back(group);
    }

    return result;
}

/*
    生成 2d-1 x d 的矩阵，第一列为 (b_0, ..., b_{d-1}, 0, ..., 0)，
    后续每列是将前一列向下循环移位。
*/
void generateMatrix(const ZZ_pE& F, long d, vector<vector<long>>& matrix) {
    long rows = 2 * d - 1;
    long cols = d;

    vector<long> vec;

    ZzpE2Veclong(F, vec, d);
    
    // 初始化矩阵
    matrix.resize(rows, vector<long>(cols, 0));

    // 第1列 (vec 填入前 d 个元素，后面填 0)
    for (long i = 0; i < d; i++) {
        matrix[i][0] = vec[i];
    }

    // 每一列都是上一列的循环移位
    for (long j = 1; j < cols; j++) {
        for (long i = 0; i < rows - 1; i++) {
            matrix[i + 1][j] = matrix[i][j - 1];
        }
    }
}

vector<long> multiplyMatrixByVector(const vector<vector<long>>& matrix, const vector<long>& vec) {
    long rows = matrix.size();          // 矩阵的行数
    long cols = matrix[0].size();       // 矩阵的列数（假设是非空矩阵）
    vector<long> result(rows, 0);       // 结果列向量，初始化为0

    // 检查矩阵的列数是否等于向量的大小
    if (cols != vec.size()) {
        cerr << "矩阵列数与向量大小不匹配！" << endl;
        return {};
    }

    // 矩阵乘以列向量
    for (int i = 0; i < rows; ++i) {
        for (int j = 0; j < cols; ++j) {
            result[i] += matrix[i][j] * vec[j];  // 矩阵的第i行乘以向量
        }
    }

    return result;  // 返回结果列向量
}

vector<vector<long>> matrixMultiply(const vector<vector<long>>& A, const vector<vector<long>>& B) {
    int m = A.size();              // A的行数
    int n = A[0].size();           // A的列数，也是B的行数
    int p = B[0].size();           // B的列数

    // 结果矩阵初始化为 m x p 大小，元素全为0
    vector<vector<long>> C(m, vector<long>(p, 0));

    // 进行矩阵乘法
    for (int i = 0; i < m; ++i) {
        for (int j = 0; j < p; ++j) {
            for (int k = 0; k < n; ++k) {
                C[i][j] += A[i][k] * B[k][j];
            }
        }
    }

    return C;
}

/*
    取模2^k
*/
long fresh(long x, int k) {
    while(x < 0){
        x += (1 << k);
    }

    return x & ((1 << k) - 1);
}


/*
    CAFE
*/
vector<vector<long>> compression(const ZZ_pX& F){

    vector<long> f;
    ZZpX2long(F, f);

    long dd = f.size();
    long d = dd - 1;

    //if(c.size() != 2*d-1) LogicError("compression: vector length mismatch");

    if(f.back() != 1) LogicError("compression: Irre error");

    // 存储结果矩阵，每一行表示 x^i 的系数
    vector<vector<long>> result;

    long rows = 2 * d - 1;
    long cols = d;

    result.resize(rows, vector<long>(cols, 0));

    for(int i=0; i<d; i++){
        result[i][i] = 1;
    }

    for(int i=0; i<d; i++){
        result[d][i] = -f[i];
    }

    for(int i=d+1; i<2*d-1; i++){

        for(int j=1; j<d; j++){
            result[i][j] = result[i-1][j-1];
        }

        for(int j=0; j<d; j++){
            result[i][j] = result[i][j] - f[j]*result[i-1][d-1];
        }
    }

    //转置
    vector<vector<long>> transposed(cols, vector<long>(rows));
    
    // 进行转置
    for (int i = 0; i < rows; ++i) {
        for (int j = 0; j < cols; ++j) {
            transposed[j][i] = result[i][j];
        }
    }

    return transposed;
}

/*
vector<vector<long>> expressHighPowers(const vector<long>& coeffs, long d) {
    // 存储结果矩阵，每一行表示 x^i 的系数
    vector<vector<long>> result;

    // 初始化矩阵，前 d-1 行直接为 x^0, x^1, ..., x^{d-1} 的系数表示
    for (long i = 0; i < d; i++) {
        vector<long> row(d, 0);
        row[i] = 1;  // x^i 用自身表示
        result.push_back(row);
    }

    // 开始处理 x^d, x^{d+1}, ..., 这些高次幂需要用 x^0, ..., x^{d-1} 表示
    for (long i = d; i < 2 * d - 1; i++) {
        vector<long> row(d, 0);

        // 使用多项式关系 x^d = -(f_0 + f_1*x + ... + f_{d-2}*x^{d-2})
        for (long j = 0; j < d - 1; j++) {
            row[j] = -coeffs[j];  // 负的多项式系数
        }

        result.push_back(row);
    }

    return result;
}
*/

ZZ_pE ZZpmulZZpE(const ZZ_p& a, const ZZ_pE& b){
    ZZ_pE _a = conv<ZZ_pE>(a);
    return _a*b;
}

vector<vector<ZZ_p>> ZZpEmatrix2ZZpmatrix(vector<vector<ZZ_pE>>& V, long degree, long k){
    long rows = V.size();
    long cols = V[0].size();

    ZZ p = conv<ZZ>(2);  //p 
    ZZ pk = power(p, k);
    ZZ_p::init(pk);

    vector<vector<long>> result(rows*degree, vector<long>(cols*degree, 0));

    for(long i=0; i<rows; i++){
        for(long j=0; j<cols; j++){
            vector<vector<long>> matrix_1; //2d-1 x d
            generateMatrix(V[i][j], degree, matrix_1);
            vector<vector<long>> matrix_2 = compression(ZZ_pE::modulus());
            vector<vector<long>> matrix =  matrixMultiply(matrix_2, matrix_1); //d x d
            for(long m=0; m<degree; m++){
                for(long n=0; n<degree; n++){
                    result[i*degree+m][j*degree+n] = fresh(matrix[m][n],k);
                }
            }
        }
    }

    vector<vector<ZZ_p>> res;
    for(long i=0; i<rows*degree; i++){
        for(long j=0; j<cols*degree; j++){
            res[i][j] = conv<ZZ_p>(result[i][j]);
        }
    }
    return res;
}