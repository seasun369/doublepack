## RMFE over GR

### Environment

RMFE库基于NTL，NTL下载地址 [Download NTL (libntl.org)](https://libntl.org/download.html)

下面介绍在Ubuntu上配置RMFE的运行环境。

下载后在终端运行

```shell
   % gunzip ntl-xxx.tar.gz
   % tar xf ntl-xxx.tar
   % cd ntl-xxx/src
   % ./configure  NTL_GMP_LIP=off NTL_EXCEPTIONS=on  NTL_GF2X_LIB=on

   % sudo make
   % sudo make check
   % sudo make install
```

这将下载NTL到 `/usr/local`

如果不想使用GMP加速，在`configure`传递参数`NTL_GMP_LIP=off` 。GMP库将显著加速NTL库中的函数。

NTL库支持GMP的版本有 3.1, 3.1.1, 4.1.4, 5.1, 6.0, and 6.1。

如果需要高效率地处理GF(2)上的运算，考虑使用`gf2x`库。首先需要在下载`gf2x`库，并在`configure`传递参数 `NTL_GF2X_LIB=on`。

如何编译使用NTL库的程序：假设已经下载NTL至 `/usr/local` 并且需要编译 `foo.cpp` 

```shell
g++ -g -O2 -std=c++11 -pthread -march=native foo.cpp -o foo -lntl -lgmp -lm
```

注意：

- NTL在 `C++11` 模式下编译，必须确保编译器支持 `C++11` 。在编译命令中需要传递 `-std=c++11` 。较新的 GCC 版本（6.1 或更高版本）默认在 `C++14` 模式下编译，因此这可能不是必要的。同样的选项也可以传递给其他类似 GCC 的编译器，例如 CLANG 和 Intel 的 ICC 编译器。
- `-march=native` 选项为针对特定 x86 架构的代码获得最佳性能。同样的选项也可以传递给其他编译器，例如 CLANG 和 Intel 的 ICC 编译器。
- 如果使用 `gf2x`编译,需要添加编译选项 `-lgf2x`  

### Detail

##### HenselLift

我们在`HenselLift.cpp`中实现了HenselLift算法,函数原型为

```c++
void HenselLift(ZZ_pX& g_, const ZZ_pX& f, const ZZ_pX& g, const ZZ p, long n)
```

输入$f,g \in Z_p[X]$,且满足$g|f$ ,函数输出$f,g\_ \in Z_{p^{n+1}}[X]$,且满足$g\_ |f$

##### Inverse

我们在`Inverse.cpp`中实现了伽罗瓦环上的求逆算法,函数原型为

```c++
ZZ_pE Inv(ZZ_pE a, long s)
```

输入元素$a \in GR(p^k,s)$,输出元素$a$的逆元.该算法的复杂度为$O(s^3)$.

##### Primitive Element

我们在`PrimitiveElement.cpp`中实现了求伽罗瓦环$GR(p^k,s)$中阶为$p^d-1$的元素.求得该元素$\zeta$后,我们可以得到该伽罗瓦环的exceptional set $T=\{0,1,\zeta,\ldots \zeta^{2^d-2}\}$.该集合满足两两元素相减可逆,我们可以从该集合中选取插值点.函数原型为

```c++
ZZ_pE FindPrimitiveElement(ZZ p, long k, long s)
```

输出伽罗瓦环$GR(p^k,s)$中阶为$p^d-1$的元素.

##### Utils

我们在`utils.cpp`中实现了一些基本的函数,主要包括

- 伽罗瓦环上的多项式插值算法
- 寻找一定次数本原多项式的算法
- 伽罗瓦环上的元素和多项式与向量之间的转换函数

### Usage

我们在`RMFE.cpp`中实现了伽罗瓦环上的两次级联RMFE算法.

我们定义了`RMFE_GR`类

```c++
class RMFE_GR{
protected:

    ZZ p;
    long s;
    long k;
    long D;
    long n1;
    long n2;
    long r1;
    long r2;
    vector<long> PrimitiveIrred1;
    vector<long> PrimitiveIrred2;

    vector<long> input;
    RMFE_element output;

public:

    vector<long> Interpolation1;
    vector<long> Interpolation2;

    
    RMFE_GR() = default;
    RMFE_GR(ZZ p_, long s_, long k_, long D_, long n1_, long n2_,
            long r1_, long r2_);
    
    void set_input(const vector<long>& Input);

    void RMFE_GR_INIT1();
    void RMFE_GR_INIT2();
    void RMFE_GR_INIT2_cache();

    void RMFE_GR_PHI();
    vector<long> RMFE_GR_PHI1(vector<long>& input);
    void RMFE_GR_PHI2(vector<long>& input);

    vector<long> RMFE_GR_PSI(RMFE_element Input);

    RMFE_element get_result(){ return this->output; }

};
```

该类定义了映射`RMFE_GR_PHI1`:$GR(p^k,s)^{n_1} \to GR(p^k,r_1)$,映射`RMFE_GR_PHI2`:$GR(p^k,r_1)^{n_2} \to GR(p^k,r_2)$,级联后的映射`RMFE_GR_PHI`:$GR(p^k,s)^{n_1n_2} \to GR(p^k,r_2)$,映射`RMFE_GR_PSI`:$GR(p^k,r_2) \to GR(p^k,s)^{n_1n_2}$.

我们用一个例子讲解使用方法.

定义参数,该`RMFE`为$GR(2^{16},1)^{4} \to GR(2^{16},20)$

```c++
    ZZ p (ZZ(2));
    long k = 16;
    long s = 1;
    long D = 2;
    long n1 = 2;
    long n2 = 2;
    long r1 = 5;
    long r2 = 20;
```

给定$GR(2^{16},1)$上的输入向量`Input`,并定义`RMFE_GR `类`a`

```c++
    vector<long> Input = {9,2,4,1};
    RMFE_GR a(p,s,k,D,n1,n2,r1,r2);
```

为`a`设置输入,并运行`RMFE_GR_PHI()`函数,使用`get_result()`得到计算结果

```c++
    a.set_input(Input);
    a.RMFE_GR_PHI();
    RMFE_element result = a.get_result();
```

运行`RMFE_GR_PSI()`函数

```c++
 vector<long> res = a.RMFE_GR_PSI(result);
```

打印`result`和`res`为

```c++
[2 7]
[9 2 4 1 ]
```

第一项代表了映射$\phi$输出的$GR(2^{16},20)$上的元素,第二项为映射$\psi$输出的$GR(2^{16},1)$上的向量.

我们定义另一个`RMFE`映射

```c++
    vector<long> Input2 = {5,9,8,7};
    a.set_input(Input2);
    a.RMFE_GR_PHI();
    RMFE_element result2 = a.get_result();
```

将`result1`和`result2`相乘,并作用$\psi$映射

```c++
    RMFE_element c = result * result2;

    vector<long> res3 = a.RMFE_GR_PSI(c);
    print(res3);
```

打印结果为

```c++
[45 18 32 7 ]
```

正好等于两个输入向量对应元素相乘的结果.



