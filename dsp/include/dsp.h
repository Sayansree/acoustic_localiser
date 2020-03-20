#ifndef DSP_H
#define DSP_H

#include <liquid/liquid.h>
#include <iostream>
#include <vector>
#include <queue>
#include <thread>
#include <chrono>
#include <cmath>
#include <complex>
#include<fstream>
#include<algorithm>
#define PI 3.1415
#define toF 1.0*Fs/time_buffer_size
#define toT 1.0/Fs
#define corr_lowpass 20
//#define resolution_check 2000
#define epsilon 3e3
#define seperation 0.30
#define speed 1500.0
//#define
#define DEBUG
using namespace std;
using namespace std::chrono;

class dsp {

public:
  explicit dsp(int,int,float,int,float,int,int);
  ~dsp();
  bool addToTimeBuffer(float);
  void fft_polar(vector<float>*,vector<float>*);
  int addToSpectralBuffer(float);
  //void corr_freq(dsp*, vector<complex<float>>*);
  //void corr_time( vector<complex<float>>,vector<float>*);
  float fft_power();
  bool check_fft(dsp*);
  void correlation(vector< complex<float> > , vector< complex<float> > ,vector<float>* );
  float get_angle();
  float get_delay();
  float get_phase();
  queue<float> tbuff_right,tbuff_center;



private:

  vector < windowcf >      freq_domain1;//dfft component matrix
  vector< complex<float> > dfft;//fixed time dfft
  int t=0,time_buffer_size,fft_buffer_size,Fs,threshold,fft_threshold,resolution_check;
  float delay,angle,freq;
  const float max_delay =seperation/speed;
  vector<vector<complex<float> >> fft_mat;
  //vector<vector<complex<float>>> ifft_mat;
  //vector< complex<float> > corr_freq;
  //vector<float>  corr_time;

};

#endif
