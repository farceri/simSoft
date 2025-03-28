//
// Author: Francesco Arceri
// Date: March 28 2025
//

#ifndef DEFS_H_
#define DEFS_H_

#include <vector>
#include <execution>  // For parallel execution policy

// global constants
const double PI = 3.141592653589793238462643383279502884197;
const long MAXDIM = 3;
const long precision = 14;

// FIRE constants
const double alpha0 = 0.2;
const double finc = 1.1;
const double fdec = 0.5;
const double falpha = 0.99;

// WCA cutoff
const double WCAcut = 1.122462048;
// Gravity
const double gravity = 9.81e-02;

// thrust functors
struct square {
    double operator()(const double xi) const { return xi * xi; }
};

struct cube {
    double operator()(const double xi) const { return xi * xi * xi; }
};

struct cosine {
    double operator()(const double xi) const { return std::cos(xi); }
};

struct sine {
    double operator()(const double xi) const { return std::sin(xi); }
};

struct randNum
{
    double a, b;
    mutable std::mt19937 rng;  // Use std::mt19937 for CPU RNG

    randNum(double _a = 0.0, double _b = 1.0) 
        : a(_a), b(_b), rng(std::random_device{}()) {}

    double operator()(const unsigned int n) const
    {
        std::uniform_real_distribution<double> dist(a, b);
        rng.discard(n);
        return dist(rng);
    }
};

struct randInt
{
    int a, b;
    mutable std::mt19937 rng;

    randInt(int _a = 0, int _b = 1) 
        : a(_a), b(_b), rng(std::random_device{}()) {}

    int operator()(const unsigned int n) const
    {
        std::uniform_int_distribution<int> dist(a, b);
        rng.discard(n);
        return dist(rng);
    }
};

struct gaussNum
{
    double a, b;
    mutable std::mt19937 rng;

    gaussNum(double _a = 0.0, double _b = 1.0) 
        : a(_a), b(_b), rng(std::random_device{}()) {}

    double operator()(const unsigned int n) const
    {
        std::normal_distribution<double> dist(a, b);
        rng.discard(n);
        return dist(rng);
    }
};

struct wrappedGaussNum
{
    double a, b;
    mutable std::mt19937 rng;

    wrappedGaussNum(double _a = 0.0, double _b = 1.0) 
        : a(_a), b(_b), rng(std::random_device{}()) {}

    double operator()(const unsigned int n) const
    {
        std::normal_distribution<double> dist(a, b);
        rng.discard(n);
        double wrapNum = dist(rng);
        wrapNum += M_PI;
        wrapNum -= 2.0 * M_PI * std::floor(wrapNum / (2.0 * M_PI));
        wrapNum -= M_PI;
        return wrapNum;
    }
};

// copied from github https://github.com/NVIDIA/thrust/blob/master/examples/strided_range.cu
template <typename Iterator>
class strided_range
{
    public:

    typedef typename thrust::iterator_difference<Iterator>::type difference_type;

    struct stride_functor : public thrust::unary_function<difference_type,difference_type>
    {
        difference_type stride;

        stride_functor(difference_type stride)
            : stride(stride) {}

        __host__ __device__
        difference_type operator()(const difference_type& i) const
        { 
            return stride * i;
        }
    };

    typedef typename thrust::counting_iterator<difference_type>                   CountingIterator;
    typedef typename thrust::transform_iterator<stride_functor, CountingIterator> TransformIterator;
    typedef typename thrust::permutation_iterator<Iterator,TransformIterator>     PermutationIterator;

    // type of the strided_range iterator
    typedef PermutationIterator iterator;

    // construct strided_range for the range [first,last)
    strided_range(Iterator first, Iterator last, difference_type stride)
        : first(first), last(last), stride(stride) {}
   
    iterator begin(void) const
    {
        return PermutationIterator(first, TransformIterator(CountingIterator(0), stride_functor(stride)));
    }

    iterator end(void) const
    {
        return begin() + ((last - first) + (stride - 1)) / stride;
    }
    
    protected:
    Iterator first;
    Iterator last;
    difference_type stride;
};

#endif /* DEFS_H_ */
