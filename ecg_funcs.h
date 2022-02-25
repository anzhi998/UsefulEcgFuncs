//
// Created by 须金柯 on 2021/12/7.
//

#ifndef MATLAB2CPP_ECG_FUNCS_H
#define MATLAB2CPP_ECG_FUNCS_H
#include <math.h>
#include <vector>
#include <set>
#include <algorithm>
#include <ctime>
#include <numeric>
#include <assert.h>
#include <thread>

template<class T>
std::vector<int> find(const std::vector<T> &input,T val);
template <class T>
std::vector<double> smooth(const std::vector<T>& input,int span);
template <class T>
void Smooth_rate(std::vector<T> &sig_in, int k, std::vector<T> &sig_out,double rate);
template <class T>

std::vector<int> findPeaks(const std::vector<T> &src, int distance=0,bool parrel=false);
template <class T>

void filter_wq_abs(const std::vector<T> &x, std::vector<T> &y, double* b, double* a, int nfilt,int abs_);
template <class T>

void denoise1(std::vector<T> &sig_in,std::vector<T> &out);
template <class T>
std::vector<T> conv1( const std::vector<T> &src,  const std::vector<T> &kernel);
template <typename T>
std::vector<T> convolution( const std::vector<T> &u, const std::vector<T> &v );


/*
 *  author:xjk
 *  function: use as matlab find(data=x)
 *  input：
 *  input: 	input array
 *  x:	    val to find
 *  return:
 *          index arrary
*/
template<class T>
std::vector<int> find(const std::vector<T> &input,T val)
{
    std::vector<int> ans;
    auto start_iter=input.begin();
    auto ans_iter=std::find(start_iter,input.end(),val);
    while (ans_iter!=input.end())
    {
        ans.push_back(std::distance(input.begin(),ans_iter));
        start_iter=ans_iter+1;
        ans_iter=std::find(start_iter,input.end(),val);
    }
    return ans;
}
/*
 *  author:xjk
 *  function: smooth signal with method 'moving'
 *  input：
 *  input: 	input array
 *  span:	window length
 *  return:
 *          smoothed arrary
*/
template <class T>
std::vector<double> smooth(const std::vector<T>& input,int span)
{
    //using default moving average
    int len=input.size();
    assert(span>=1);
    if(span%2==0) --span;
    std::vector<double> ans(len,0);
    ans[0]=input[0];
    ans[len-1]=input[len-1];
    for (int i = 1; i < len-1; ++i) {
        int win_len=0;
        double temp=0;
        if(i<span/2) win_len=i;
        else if(i>(len-1-span/2)) win_len=len-1-i;
        else win_len=span/2;
        for (int j = i-win_len; j <=i+win_len ; ++j) {
            temp+=input[j];
        }
        ans[i]=temp/static_cast<double>(2*win_len+1);
    }
    return ans;
}
/*
 *  author:xjk
 *  function: smooth signal with windowlengh k and Radio rate
 *  input：
 *  src: 	input array
 *  k:	    windowlength
 *  sigout: output array ( assert(output.size()==input.size()  )
 *  rate:   value_radio
*/
template <class T>
void Smooth_rate(std::vector<T> &sig_in, int k, std::vector<T> &sig_out,double rate)
{
    int end=sig_in.size()-1;
    //用位移运算符替代除法，如果是除以2的次方的话
    int b = (k - 1) >> 1;
    int i=0;//
    //从大for里面拆出来，不用每次都判断属于哪个区间
    T temp=0;
    for (; i <=b; ++i) {
        sig_out[i] = 0;
        //不用每次都从头算，可以用上一个的值
        //本质上就是从头开始累加2i个，比起上次的结果多 2i-1,2i
        if(i==0){
            sig_out[i] += sig_in[0];
        } else{
            sig_out[i] =temp+sig_in[2*i-1]+sig_in[2*i];
        }
        temp=sig_out[i];
        sig_out[i] = sig_out[i] / (i * 2 + 1);
        sig_out[i]*=rate;
    }
    //到这里i=b+1
    for (; i <(end-b); ++i) {
        sig_out[i] = 0;
        //i=b+1,j从1到2b+1
        //i=b+2,j从2到2b+2
        //i=b+3,j从3到2b+3
        //每次减去前一个数（i-b-1),加上一个数(i+b)
        if(i==(b+1)){
            for (int j = i - b; j <= i + b; ++j)
            {
                sig_out[i] += sig_in[j];
            }
        }else{
            sig_out[i]=temp-sig_in[i-b-1]+sig_in[i + b];
        }
        temp=sig_out[i];
        sig_out[i] = sig_out[i] / (2 * b + 1);
        sig_out[i]*=rate;
    }
    //到这里i=end-b
    for (;i<=end;++i){
        sig_out[i] = 0;
        //I=end-b
        //I增加1，2*i-end 增加2
        //也就是后一个比前一个少加了2个数 2*(i-1)-end,2*(i-1)-end+1
        //比如说  i=end-b,j从end-b加到end,i增加1，j从end-b+2加到end
        //就是少加了 2*(i-1)-end=end-b 和 end-b+1 这两个数
        if(i==(end-b)){
            for (int j = 2 * i - end; j <= end; ++j)
            {
                sig_out[i] += sig_in[j];
            }
        }
        else{
            sig_out[i]=temp-sig_in[2*(i-1)-end]-sig_in[2*(i-1)-end+1];
        }
        temp=sig_out[i];
        sig_out[i] = sig_out[i] / (2 * (end - i) + 1);
        sig_out[i]*=rate;
    }
}
/*
 *  author:xjk
 *  function: find peaks in signal fileted by minpeakdistance
 *  input：
 *  src: 	input array
 *  distance:	minpeakdistance
 *  parrel: multi thread computing
 *  output:
 *  peak_index array
*/
template <class T>
std::vector<int> findPeaks(const std::vector<T> &src, int distance,bool parrel)
{
    int length=src.size();
    if(length<=1) return std::vector<int>();
    //we dont need peaks at start and end points
    std::vector<int> sign(length-2,-1);
    std::vector<T> difference(length,0);
    std::vector<int> temp_out;

    const unsigned int nums_kernels=std::thread::hardware_concurrency();
    if(parrel&&nums_kernels>=2){
        //每个线程需要计算的点数
        int len_one=static_cast<int>(src.size()/nums_kernels);
        //线程容器
        std::vector<std::thread> threads(nums_kernels-1);
        //除了主线程外

        auto block_begin=src.begin();
        for (int i = 1; i <nums_kernels; ++i) {
            //第i个线程处理[(i-1)*len,i*len-1]
            //
            auto block_end=block_begin;
            std::advance(block_end,len_one);

        }

        std::for_each(threads.begin(),threads.end(),std::mem_fn(&std::thread::join));
        //主线程计算最后剩下来的点
        int main_start=length-(nums_kernels-1)*len_one;

    }
    else{
        //first-order difference (sign)
        std::adjacent_difference(src.begin(),src.end(),difference.begin());
        difference.erase(difference.begin());
        difference.pop_back();
        for (int i = 0; i < difference.size(); ++i) {
            if(difference[i]>=0) sign[i]=1;
        }
        //second-order difference
        T  diff2 = 0;
        for (int j = 1; j < length-1; ++j)
        {
            int  diff = sign[j] - sign[j - 1];
            diff2 += difference[j-1];
            if ((diff < 0) && diff2 != 0) {
                temp_out.push_back(j);
            }
        }
    }


    if(temp_out.size()==0 || distance==0 ) return temp_out;

    //sort peaks from large to small by src value at peaks
    std::sort(temp_out.begin(),temp_out.end(),[&src](int a,int b){
        return src[a]>=src[b];
    });

    std::vector<int> ans;
    std::set<int> except;

    //Initialize the answer and the collection to be excluded

    ans.push_back(temp_out[0]);
    int left=std::max(0,temp_out[0]-distance);
    int right=std::min(length-1,temp_out[0]+distance);
    for (int i = left;i<=right; ++i) {
        except.insert(i);
    }
    for (int i = 1; i < temp_out.size(); ++i) {
        auto it=except.find(temp_out[i]);
        // if not in except
        // add it to the answer and update except
        if(it==except.end())
        {
            ans.push_back(temp_out[i]);
            int left=std::max(0,temp_out[i]-distance);
            int right=std::min(length-1,temp_out[i]+distance);
            //update except set
            for (int j = left;j<=right; ++j) {
                except.insert(j);
            }
        }
    }
    //sort the ans from small to large by index value
    std::sort(ans.begin(),ans.end());
    return ans;
}
/*
 *  author:xjk
 *  function: filt signal by IIR filter(b,a)
 *  input：
 *  x: 	input array
 *  y:	output array
 *  b:  Filter molecular coefficient
 *  a:  Filter denominator coefficient
 *  nfilt： order of the filter （length of a)
 *  abs_:  1 if need abs result , 0 else
*/
template <class T>
void filter_wq_abs(const std::vector<T> &x, std::vector<T> &y, double* b, double* a, int nfilt,int abs_)
{
    int xlen=x.size();
    double tmp;
    int i, j;
    double EPS = 0.00000001;
    if(y.size()<xlen){
        y.resize(xlen);
    }
    //normalization
    if ((*a - 1.0 > EPS) || (*a - 1.0 < -EPS))  //设nfilt=3,则b0x(n)+b1x(n-1)+b2x(n-2)=a0y(n)+a1y(n-1)+a2y(n-2)；[Y(Z)/X(Z)=b0+b1Z^-1+b2Z^-2/a0+a1Z^-1+a2Z^-2]
    {
        tmp = *a;  //除以a的第一个地址的数即a[0]
        for (i = 0; i < nfilt; i++)
        {
            b[i] /= tmp;
            a[i] /= tmp;
        }
    }

    a[0] = 0.0;  //因为后面循环中减去a0y(n)这一项，所以a0=0
    for (i = 0; i < xlen; i++)
    {
        for (j = 0; i >= j && j < nfilt; j++) //j < nfilt是的a与b的大小在nfit内
        {
            y[i] += (b[j] * x[i - j] - a[j] * y[i - j]);  //例如：y(3)=b(0)x(3)-a(0)y(3)+b(1)x(2)-a(1)y(2)+b(2)x(1)-a(2)y(1)

        }
        if (abs_ == 1)
            y[i] = fabs(y[i]);
        // if(zi&&i<nfilt-1) y[i] += zi[i];   //zi[i],查分方程可以不计，滤波的表达式 1.[y, zf] = filter(b ,a, X) 2.[y, zf] = filter(b, a, X, zi)
    }
    a[0] = 1.0;
}

/*
 *  author:xjk
 *  function: Remove baseline offset
 *  input：
 *  sig_in: 	input array
 *  out:	    output array
*/
template <class T>
void denoise1(std::vector<T> &sig_in,std::vector<T> &out)
{


    auto *sig_out1=new std::vector<T>(sig_in.size(),0);
    auto *baseline=new std::vector<T>(sig_in.size(),0);
    Smooth_rate(sig_in, 5, *sig_out1,1);
    Smooth_rate(*sig_out1,  5, out,1);
    Smooth_rate(out, 150, *baseline,1);
    for (int i = 0; i < sig_in.size(); ++i)
        out.at(i)  -= baseline->at(i);
    delete sig_out1;
    delete baseline;
}
/*
 * convolution and ploynonal multiplication.
 */
template <typename T>
std::vector<T> conv1( const std::vector<T> &signal, const std::vector<T> &filter )
{
    if( signal.size() < filter.size() )
        return convolution( filter, signal );
    else
        return convolution( signal, filter );
}

template <typename T>
std::vector<T> convolution( const std::vector<T> &u, const std::vector<T> &v )
{
    int m = u.size();
    int n = v.size();
    //动态断言
    assert( m >= n );

    int length = m + n - 1;
    std::vector<T> x(length,0);
    //参考matlab conv()的help
    for (int k = 0; k < length; ++k) {
        for (int j = std::max(0,k-n); j <std::min(k+1,m) ; ++j) {
            x[k]+=u[j]*v[k-j];
        }
    }
    return x;
}



#endif //MATLAB2CPP_ECG_FUNCS_H
