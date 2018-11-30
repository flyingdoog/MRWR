/**
 * @file pprclus_mex.cc
 * Implement a PPR clustering scheme.
 * 
 * mex pprclus_mex.cc CXXFLAGS="\$CXXFLAGS -std=c++0x" -largeArrayDims
 */


#include <vector>
#include <queue>
#include <utility> // for pair sorting
#include <assert.h>
#include <limits>
#include <algorithm>
#include<iostream>

#include <unordered_set>
#include <unordered_map>
#define tr1ns std
#ifndef __APPLE__
#define __STDC_UTF_16__ 1
#endif

#include <mex.h>
using namespace std;
/** A replacement for std::queue<int> using a circular buffer array */
class array_queue {
    public:
    std::vector<int> array;
    size_t max_size;
    size_t head, tail;
    size_t cursize;
    array_queue(size_t _max_size)
    : max_size(_max_size), array(_max_size), head(0), tail(0), cursize(0)
    {}
    
    void empty() {
        head = 0;
        tail = 0;
        cursize = 0;
    }
    
    size_t size() {
        return cursize;
    }
    
    void push(int i) {
        assert(size() < max_size);
        array[tail] = i;
        tail ++;
        if (tail == max_size) {
            tail = 0;
        }
        cursize ++;
    }
    
    int front() {
        assert(size() > 0);
        return array[head];
    }
    
    void pop() {
        assert(size() > 0);
        head ++;
        if (head == max_size) {
            head = 0;
        }
        cursize --;
    }
};

struct sparsevec {
  //mxIndex ->int
  typedef tr1ns::unordered_map<mwIndex,double> map_type;
  map_type map;
  /** Get an element and provide a default value when it doesn't exist
   * This command does not insert the element into the vector
   */
  double get(mwIndex index, double default_value=0.0) {
    map_type::iterator it = map.find(index);
    if (it == map.end()) {
      return default_value;
    } else {
      return it->second;
    }
  }

  /** Compute the sum of all the elements
   * Implements compensated summation
   */
  double sum() {
    double s=0.;
    for (map_type::iterator it=map.begin(),itend=map.end();it!=itend;++it) {
      s += it->second;
    }
    return s;
  }
  
  /** Compute the max of the element values 
   * This operation returns the first element if the vector is empty.
   */
  mwIndex max_index() {
    mwIndex index=0;
    double maxval=std::numeric_limits<double>::min();
    for (map_type::iterator it=map.begin(),itend=map.end();it!=itend;++it) {
      if (it->second>maxval) { maxval = it->second; index = it->first; }
    }
    return index;
  }
};

struct sparserow {
    mwSize n, m; 
    mwIndex *ai;
    mwIndex *aj;
    double *a;
};

double sr_degree(sparserow *s, mwIndex u) {
    
  return (s->ai[u+1] - s->ai[u]);
}

double sr_volumn(sparserow *s, mwIndex u) {
    
  double sum = 0;
  for (int i=s->ai[u];i<s->ai[u+1];i++)
    sum+=s->a[i];
  //std:cout<<sum<<endl;
  return sum;
}

template <class Queue>
int compute_local_pagerank(sparserow *s, sparsevec& r, sparsevec& p, 
    double alpha, double epsilon, int max_push_count, Queue& q) 
{
  for (sparsevec::map_type::iterator it=r.map.begin(),itend=r.map.end();
        it!=itend;++it){
    if (it->second > epsilon*sr_volumn(s,it->first)) {
      q.push(it->first);
    }
  }

  // cout<<epsilon<<endl;
  int push_count = 0;
  while (q.size()>0 && push_count < max_push_count) {
    push_count += 1;
    mwIndex u = q.front();
    q.pop();
    double du = sr_volumn(s, u);
    double moving_probability = r.map[u] - 0.5*epsilon*(double)du;
    r.map[u] = 0.5*epsilon*(double)du;
    p.map[u] += (1.-alpha)*moving_probability;
    //cout<<p.map[u]<<endl;
    double neighbor_update = alpha*moving_probability/(double)du;

    for (mwIndex nzi=s->ai[u]; nzi<s->ai[u+1]; nzi++) {
      mwIndex x = s->aj[nzi];
      double dx = sr_volumn(s, x);
      double rxold = r.get(x);
      
      /*
      *Coded by Dongsheng
      */
      double w1 = s->a[nzi];
      double rxnew = rxold + w1*neighbor_update;
      r.map[x] = rxnew;
      if (rxnew > epsilon*dx && rxold <= epsilon*dx) {
        q.push(x);
      }
    }
  }
  
  return (push_count);
}



struct greater2nd {
  template <typename P> bool operator() (const P& p1, const P& p2) {
    return p1.second > p2.second;
  }
};

void cluster_from_sweep(sparserow* G, sparsevec& p, 
      std::vector<mwIndex>& cluster, double *outcond, double* outvolume,
      double *outcut)
{
  // now we have to do the sweep over p in sorted order by value
  typedef std::vector< std::pair<int, double> > vertex_prob_type;
  vertex_prob_type prpairs(p.map.begin(), p.map.end());
  std::sort(prpairs.begin(), prpairs.end(), greater2nd());

  // compute cutsize, volume, and conductance
  std::vector<double> conductance(prpairs.size());
  std::vector<mwIndex> volume(prpairs.size());
  std::vector<mwIndex> cutsize(prpairs.size());

  size_t i=0;
  tr1ns::unordered_map<int,size_t> rank;
  for (vertex_prob_type::iterator it=prpairs.begin(),itend=prpairs.end();
    it!=itend; ++it, ++i) {
    rank[it->first] = i;
  }
  //printf("support=%i\n",prpairs.size());
  // mwIndex total_degree = G->ai[G->m];
  double total_degree=0;
  for(int i=0;i<G->m;i++){
    total_degree+=G->a[i];
  }



  double curcutsize = 0;
  double curvolume = 0;
  i=0;
  for (vertex_prob_type::iterator it=prpairs.begin(),itend=prpairs.end();
    it!=itend; ++it, ++i) {
    mwIndex v = it->first;
    // mwIndex deg = G->ai[v+1]-G->ai[v];

    double  deg = sr_volumn(G,v);
    double change = deg;
    for (mwIndex nzi=G->ai[v]; nzi<G->ai[v+1]; ++nzi) {
      mwIndex nbr = G->aj[nzi];
      if (rank.count(nbr) > 0) {
        if (rank[nbr] < rank[v]) {
          // change -= 2;

          /*
          * Coded by Dongsheng Luo
          */
          int begin=G->ai[v];
          int end = G->ai[v+1];
          for(;begin!=end;begin++){
            if(G->aj[begin]==nbr){
              change -=2*G->a[begin];
            }
          }
          ///////


        }
      }
    }
    curcutsize += change;
    //if (curvolume + deg > target_vol) {
      //break;
    //}
    curvolume += deg;
    volume[i] = curvolume;
    cutsize[i] = curcutsize;
    if (curvolume == 0 || total_degree-curvolume==0) {
      conductance[i] = 1;
    } else {
      conductance[i] = (double)curcutsize/
                        (double)std::min(curvolume,total_degree-curvolume);
    }
    //printf("%5i : cut=%6i vol=%6i prval=%8g cond=%f\n", i, curcutsize, curvolume, it->second, conductance[i]);
  }
  // we stopped the iteration when it finished, or when it hit target_vol
  size_t lastind = i;
  double mincond = std::numeric_limits<double>::max();
  size_t mincondind = 0; // set to zero so that we only add one vertex 
  for (i=0; i<lastind; i++) {
    if (conductance[i] < mincond) {
      mincond = conductance[i];
      mincondind = i;
    }
  }
  //printf("mincond=%f mincondind=%i\n", mincond, mincondind);
  if (lastind == 0) {
    // add a case 
    mincond = 0.0;
  }
  i = 0;
  for (vertex_prob_type::iterator it=prpairs.begin(),itend=prpairs.end();
    it!=itend && i<mincondind+1; ++it, ++i) {
    cluster.push_back(it->first);
  }
  if (outcond) { *outcond = mincond; }
  if (outvolume) { *outvolume = volume[mincondind]; }
  if (outcut) { *outcut = cutsize[mincondind]; }
}

struct local_pagerank_stats {
    double conductance;
    double volume;
    double support;
    double steps;
    double eps;
    double cut;
};

/** Cluster will contain a list of all the vertices in the cluster
 * @param set the set of starting vertices to use
 * @param alpha the value of alpha in the PageRank computation
 * @param target_vol the approximate number of edges in the cluster
 * @param p the pagerank vector
 * @param r the residual vector
 * @param a vector which supports .push_back to add vertices for the cluster
 * @param stats a structure for statistics of the computation
 */
template <class Queue>
int hypercluster_pagerank_multiple(sparserow* G, 
    const std::vector<mwIndex>& set, double alpha, double target_vol, 
    sparsevec& p, sparsevec &r, Queue& q,
    std::vector<mwIndex>& cluster, local_pagerank_stats *stats)
{
  // reset data
  p.map.clear();
  r.map.clear();
  q.empty();
  
  assert(target_vol > 0);
  assert(alpha < 1.0); assert(alpha > 0.0);
  
  //r.map[start] = 1.0;
  double maxdeg = 0;
  for (size_t i=0; i<set.size(); ++i) {
    assert(set[i] >= 0); assert(set[i] < G->n);
    r.map[set[i]] = 1./(double)(set.size());
    //r.map[set[i]] = 1.;
    maxdeg = std::max(maxdeg, sr_volumn(G,set[i]));
  }
  //double pr_eps = 1.0/std::max((double)sr_volumn(G,start)*(double)target_vol, 100.0);
  //double pr_eps = std::min(1.0/std::max(10.*target_vol, 100.0), 
    //1./(double)(set.size()*maxdeg + 1));
  double pr_eps = 1.0/std::max(10.0*target_vol, 100.0);
  if (stats) { stats->eps = pr_eps; }
  
  //printf("find_cluster: target_vol=%7lli alpha=%5.3ld pr_eps=%ld\n", target_vol, alpha, pr_eps);
  
  // calculate an integer number of maxsteps
  double maxsteps = 1./(pr_eps*(1.-alpha));
  maxsteps = std::min(maxsteps, 0.5*(double)std::numeric_limits<int>::max());
      
  int nsteps = compute_local_pagerank(G, r, p, alpha, pr_eps, (int)maxsteps, q);
  if (nsteps == 0) {
    p = r; // just copy over the residual
  }
  int support = r.map.size(); 
  if (stats) { stats->steps = nsteps; }
  if (stats) { stats->support = support; }
  
  //mexPrintf("setsize=%zu, nsteps=%i, support=%i\n", set.size(), nsteps, support);

  // scale the probablities by their degree
  for (sparsevec::map_type::iterator it=p.map.begin(),itend=p.map.end();
    it!=itend;++it) {
    double maxi = sr_volumn(G,it->first);
    if(maxi<1.0)
        maxi = 1.0;
    it->second *= 1.0/maxi;  
  }
  double *outcond = NULL;
  double *outvolume = NULL;
  double *outcut = NULL;
  if (stats) { outcond = &stats->conductance; }
  if (stats) { outvolume = &stats->volume; }
  if (stats) { outcut = &stats->cut; }
  cluster_from_sweep(G, p, cluster, outcond, outvolume, outcut);
  return (0);
}

void pprgrow(sparserow* G, std::vector<mwIndex>& set, double alpha,
    double targetvol, double* fcond, double* fcut,
    double* fvol,sparsevec &p )
{
    sparsevec r;
    std::queue<mwIndex> q;
    local_pagerank_stats stats;
    std::vector<mwIndex> bestclus;
    hypercluster_pagerank_multiple(G, set, alpha, targetvol, 
        p, r, q, bestclus, &stats);
    set = bestclus;
    *fcond = stats.conductance;
    *fcut = stats.cut;
    *fvol = stats.volume;
}

void copy_array_to_index_vector(const mxArray* v, std::vector<mwIndex>& vec)
{
    mxAssert(mxIsDouble(v), "array type is not double");
    size_t n = mxGetNumberOfElements(v);
    double *p = mxGetPr(v);

    vec.resize(n);
    
    for (size_t i=0; i<n; ++i) {
        double elem = p[i];
        mxAssert(elem >= 1, "Only positive integer elements allowed");
        vec[i] = (mwIndex)elem - 1;
    }
}


// USAGE
// [bestset,pindex,pscores,cond,cut,vol] = pprgrow_mex(A,set,targetvol,alpha)
void mexFunction(int nlhs, mxArray* plhs[], int nrhs, const mxArray* prhs[])
{
    mxAssert(nrhs > 2 && nrhs < 5, "2-4 inputs required.");

    const mxArray* mat = prhs[0];
    const mxArray* set = prhs[1];

    mxAssert(mxIsSparse(mat), "Input matrix is not sparse");
    mxAssert(mxGetM(mat) == mxGetN(mat), "Input matrix not square");

    mxArray* cond = mxCreateDoubleMatrix(1,1,mxREAL);
    mxArray* cut = mxCreateDoubleMatrix(1,1,mxREAL);
    mxArray* vol = mxCreateDoubleMatrix(1,1,mxREAL);
    
    if (nlhs > 3) { plhs[3] = cond; }
    if (nlhs > 4) { plhs[4] = cut; }
    if (nlhs > 5) { plhs[5] = vol; }

    //mxAssert(nlhs <= 4, "Too many output arguments");

    double alpha = 0.99;
    if (nrhs >= 4) {
        alpha = mxGetScalar(prhs[3]);
    }
    mxAssert(alpha >= 0. && alpha < 1, "alpha must be 0 <= alpha < 1");

    // use a strange sentinal
    double targetvol = 1000.;
    if (nrhs >= 3) {
        targetvol = mxGetScalar(prhs[2]);
    }

    sparserow r;
    r.m = mxGetM(mat);
    r.n = mxGetN(mat);
    r.ai = mxGetJc(mat);//number of nonzero elements in each colomns
    r.aj = mxGetIr(mat);//nonzero elements index 
    r.a = mxGetPr(mat);//real number
    std::vector< mwIndex > cluster;
    copy_array_to_index_vector( set, cluster );
    sparsevec p;
    pprgrow(&r, cluster, alpha, targetvol,
        mxGetPr(cond), mxGetPr(cut), mxGetPr(vol),p);

    if (nlhs > 0) { 
        mxArray* cassign = mxCreateDoubleMatrix(cluster.size(),1,mxREAL);
        plhs[0] = cassign;

        double *ci = mxGetPr(cassign);
        for (size_t i=0; i<cluster.size(); ++i) {
            ci[i] = (double)(cluster[i] + 1);
        }
    //}

   double thres = 0;//1e-6;
   tr1ns::unordered_map<mwIndex,double> pos;
    typedef tr1ns::unordered_map<mwIndex,double> map_type;
    map_type::iterator it=p.map.begin(), itend=p.map.end();
    for(;it!=itend;it++){
        if(it->second>thres)
            pos[it->first]=it->second;
    }
    //printf("pos size%d",pos.size());
    //Coded by Dongsheng
    //if (nlhs > 2) { 
       mxArray* pindex = mxCreateDoubleMatrix(pos.size(),1,mxREAL);
        mxArray* pscores = mxCreateDoubleMatrix(pos.size(),1,mxREAL);
        plhs[1] = pindex;
        plhs[2] = pscores;
     
        double *pindexi = mxGetPr(pindex);
        double *pscoresi = mxGetPr(pscores);
        
        typedef tr1ns::unordered_map<mwIndex,double> map_type;
        it=pos.begin(), itend=pos.end();

        for (size_t i=0; it!=itend; ++i,++it) {
            pindexi[i] = (double)(it->first + 1);
            pscoresi[i] = (double)(it->second);
        }
    }

}