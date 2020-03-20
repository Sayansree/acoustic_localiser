///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//  Author: Sayansree Paria                                                                                             //
//  email:  sayansreeparia@gmail.com                                                                                    //
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//OBJECTIVE:                                                                                                            //
//  implement  high density discrete time discreste frequency fourier transform dtdft in O(n) using window/ring buffer  //
//  to calculate time delay and phase difference between hydrophone signal to estimate its position                     //
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//NOTE:                                                                                                                 //
//  the library is still under devlopment and requires major bug fix and optimisation                                   //
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
#include <dsp.h>

//pass time domain threshold ,frequency domain power threshold,ping time,small fft buffer size, central freq of small fft ,resolution of small fft(slows the checking time),sampling frequency
//note for phase based system it is recomended to chose small fft window with length integral multiple of the incomming ping and a high speed SAR ADC

dsp::dsp(int threshold,int fft_threshold,float ping_time, int fft_buffer_size,float freq,int resolution_check,int Fs):threshold(threshold),fft_threshold(fft_threshold),
freq_domain1(resolution_check),Fs(Fs),freq(freq),resolution_check(resolution_check),
fft_mat(resolution_check),dfft(resolution_check),fft_buffer_size(fft_buffer_size),time_buffer_size(ping_time*Fs*2),t(0){
  for(int k=0;k<resolution_check;k++)
  {   freq_domain1[k]=(windowcf_create(fft_buffer_size));
    for(int i=0;i<fft_buffer_size;i++)
    {
      fft_mat[k].push_back(complex<float>(cos(2*PI*(freq+(2*k-resolution_check)*epsilon/100.0)*i*toT),-sin(2*PI*(freq+(2*k-resolution_check)*epsilon/100.0)*i*toT)));
    }
  }
}

dsp::~dsp()//cestructor
{

}

//this function is contributed by by Abhijeet Triparty
void dsp::correlation(vector< complex<float> > x1, vector< complex<float> > x2,vector<float>* out)//calculates the circular cross corellation between two signals in frequency domain
{

  std::complex<float>* input1 = x1.data();
  std::complex<float>* output1 = new std::complex<float>[x1.size()];
  fftplan fft1 = fft_create_plan(x1.size(),reinterpret_cast<liquid_float_complex*>(input1),reinterpret_cast<liquid_float_complex*>(output1), LIQUID_FFT_FORWARD, 0);
  fft_execute(fft1);
  fft_destroy_plan(fft1);

  std::complex<float>* input2 = x2.data();
  std::complex<float>* output2 = new std::complex<float>[x2.size()];
  fftplan fft2 = fft_create_plan(x2.size(),reinterpret_cast<liquid_float_complex*>(input2),reinterpret_cast<liquid_float_complex*>(output2), LIQUID_FFT_FORWARD, 0);
  fft_execute(fft2);
  fft_destroy_plan(fft2);

  vector< complex<float> > c;
  for(int i=0;i<x1.size()/2;i++)
  {
    c.push_back((*output1++)*conj((*output2++))*complex<float>((t<corr_lowpass/toF)?1:0,0));//implements an ideal sharp low pass filter
  }
  for(int i=x1.size()/2;i<x1.size();i++)
  {
    c.push_back((*output1++)*conj((*output2++))*complex<float>((t>x1.size()-corr_lowpass/toF)?1:0,0));//implements a ideal sharp low pass filter
  }

  std::complex<float>* input3 = c.data();
  std::complex<float>* output3 = new std::complex<float>[c.size()];
  fftplan ifft = fft_create_plan(c.size(),reinterpret_cast<liquid_float_complex*>(input3),reinterpret_cast<liquid_float_complex*>(output3), LIQUID_FFT_BACKWARD, 0);
  fft_execute(ifft);

  for(int i=0; i<c.size(); i++)
    out->push_back(abs(*output3++));
  rotate(out->begin(),out->begin() + (out->size()/2),out->end());
  }


bool dsp::addToTimeBuffer(const float ft)//store the current signal in the buffer and check if it has crossed a threshold
{
  tbuff_right.push(ft);
  //addToSpectralBuffer(ft);
  if(tbuff_right.size()>time_buffer_size){
    tbuff_center.push(tbuff_right.front());
    tbuff_right.pop();
      }
  if(tbuff_center.size()>time_buffer_size)
        tbuff_center.pop();
return(tbuff_center.size()==time_buffer_size&&abs(tbuff_center.front())>threshold);
}

//note work left on phase based system the fft window would slide  along the ping pulse calculating the signal strength and variance over movement as we reach convergence point the rest signal is discard and phase is published
bool dsp::check_fft(dsp* hyd2)//performs a fast frequency specific fft check on the central buffer if pinger of that freq is detected we proceed to corellation on the central buffer returns true if ping is detected
  {
    #ifdef DEBUG
    std::fstream tfile1,tfile2,tfile3;
    tfile1.open("/home/sayansree/catkin_ws/acoustic_localiser/script/tdata1.txt",std::ios::out);
    tfile2.open("/home/sayansree/catkin_ws/acoustic_localiser/script/tdata2.txt",std::ios::out);
    tfile3.open("/home/sayansree/catkin_ws/acoustic_localiser/script/tcorr2.txt",std::ios::out);
    high_resolution_clock::time_point t2,t1=high_resolution_clock::now();
    #endif
      int k=0;
      vector< complex<float> > x1,x2;


      while(!tbuff_center.empty())
      {

        if(k<fft_buffer_size){
          addToSpectralBuffer(tbuff_center.front());
          hyd2->addToSpectralBuffer(hyd2->tbuff_center.front());
        }else if(k==fft_buffer_size)
        {   if((fft_power()<fft_threshold)&&(hyd2->fft_power()<fft_threshold))
            {
            tbuff_center=queue<float>();
            hyd2->tbuff_center=queue<float>();
          return false;}
          #ifdef DEBUG
          t2=high_resolution_clock::now();
          duration<double> tim = duration_cast<duration<double>>(t2-t1);
          std::cout <<"phase compute time "<<  tim.count()<< '\n';
          #endif
        }
          x1.push_back(complex<float>(tbuff_center.front(),0));
          x2.push_back(complex<float>(hyd2->tbuff_center.front(),0));
          tfile1<<tbuff_center.front()<<" "<<k*toT<<" "<<0<<"\n";
          tfile2<<hyd2->tbuff_center.front()<<" "<<k*toT<<" "<<0<<"\n";
          hyd2->tbuff_center.pop();
          tbuff_center.pop();
          k++;
      }


      vector <float> out;
      correlation(x1,x2,&out);
      int max_index=0,l= out.size();
      for(int k=0;k<l;k++)
      if(out[max_index]<out[k])max_index=k;
      delay=(max_index-l/2)*toT;

      #ifdef DEBUG
    //  cout<<"signal delay"<<1e6*delay<<'\n';
      tfile1.close();
      tfile2.close();
      for(int k=0;k<l;k++)
        tfile3<<out[k]<<" "<<1e6*(k-l/2)*toT<<"\n";
      tfile3.close();
      #endif

      return true;

      //tbuff_center=q2;
  }

float dsp::get_angle() //returns the angular heading assuming object is at distance compared to the hydrophone/microphone spacing
{
  delay=(delay>max_delay)?max_delay:(delay<-max_delay)?-max_delay:delay;
  return asin(delay/max_delay)*180/PI;
}
float dsp::get_delay()//returns the calculated time delay of arrival of ping signal <TDOA>
{
  return delay;
}
float dsp::get_phase()
{
  float mag=0, phase;
  for(complex<float> f: dfft)
    if(mag<abs(f))
    {
      mag=abs(f);
    phase=arg(f);
  }
  return phase*180/PI;
}
void dsp::fft_polar(vector<float> *u,vector<float> *v) //returns the frequency specific fft in polar format u <power distribution of dignal>  and v <corrosponding phase of the signal>
{
  for(complex<float> f: dfft)
  {
    u->push_back(abs(f));
    v->push_back(arg(f));
  }
}
float dsp::fft_power() //calculate relative power distribution across the central frequency over a predefined bandwidth
{
  float s=0;
  for(complex<float> f: dfft)
    s+=abs(f);
  return s;
}
int dsp::addToSpectralBuffer(float fx)  //adds a new value to frequency specific  dft buffer and dft is calculated for the change of datapoint instantly in realtime it would give a O(n)performance
{
  complex<float>  old_val,new_val;
  for(int k=0;k<resolution_check;k++){
    new_val=(complex<float> (fx,0))*fft_mat[k][t];
    windowcf_index(freq_domain1[k],0,reinterpret_cast<liquid_float_complex*>(&old_val));
    windowcf_push(freq_domain1[k],*reinterpret_cast<liquid_float_complex*>(& new_val));
    dfft[k]+=new_val-old_val;
  }
  t=(t+1)%fft_buffer_size;
}

int main()//main function is a sample simulator devloped to test the library
{

  const float pinger_time=10e-3,central_freq=30e3,phase_offset=10.7;

    const int n=512,th=4,fft_threshold=4000,Fs=1000000,time_buffer_size=Fs*pinger_time*2,fft_buffer_size=1000,resolution_check=1000;
    dsp hyd1(th,fft_threshold,pinger_time,fft_buffer_size,central_freq,resolution_check,Fs);
    dsp hyd2(th,fft_threshold,pinger_time,fft_buffer_size,central_freq,resolution_check,Fs);

   float sim_delay=100e-6,phase=30;
   float snr=10;
   cout<<"enter time delay to be simulated";
   cin>>sim_delay;
   cout<<"enter phase delay to be simulated";
   cin>>phase;
   phase*=PI/180;
  std::fstream tfile1,ffile1,tfile2,ffile2,fcorr,tcorr;

  ffile1.open("/home/sayansree/catkin_ws/acoustic_localiser/script/fdata1.txt",std::ios::out);

  ffile2.open("/home/sayansree/catkin_ws/acoustic_localiser/script/fdata2.txt",std::ios::out);
  fcorr.open("/home/sayansree/catkin_ws/acoustic_localiser/script/fcorr.txt",std::ios::out);
  tcorr.open("/home/sayansree/catkin_ws/acoustic_localiser/script/tcorr.txt",std::ios::out);
  high_resolution_clock::time_point t2,t1=high_resolution_clock::now();
    for(int t=0;t<3*time_buffer_size;t++)
    {
      float T=1.0*t/Fs;
      float ft1=(rand()%20/10.0-1);
      float ft2=(rand()%20/10.0-1);
      //hyd1.addToTimeBuffer(10*cos(2*PI*45*T)+6*sin(2*PI*30*T));
      if(T>5e-3&&T<5e-3+pinger_time)
      ft1+=snr*cos(2*PI*30000*T-phase);
      if(T-sim_delay>5e-3&&T<5e-3+sim_delay+pinger_time)
      ft2+=snr*cos(2*PI*30000*T);
      //hyd2.addToBuffer(ft2);
      //);
      if((hyd2.addToTimeBuffer(ft1))||(hyd1.addToTimeBuffer(ft2)))
        { t2=high_resolution_clock::now();
          duration<double> tim = duration_cast<duration<double>>(t2-t1);
            std::cout <<"data push time per datapoint "<<  tim.count()/t<< '\n';
            t1=high_resolution_clock::now();
           if(hyd1.check_fft(&hyd2)){
           cout<<"calculated signal time delay:"<<1e6*hyd1.get_delay()<< " us"<<'\n';
           cout<<"heading angle:"<<hyd1.get_angle()<< " degrees"<<'\n';
           cout<<"phase difference:"<<hyd1.get_phase()-hyd2.get_phase()+phase_offset<< " degrees"<<'\n';
           t2=high_resolution_clock::now();
             duration<double> tim = duration_cast<duration<double>>(t2-t1);
               std::cout <<"angle compute time "<<  tim.count()<< "\n\n";
             }
          //hyd2.check_fft();

          break;
        }

    }

    vector<float> f1,f2,a1,a2;
    vector<complex<float>>fcor;
    vector<float>tcor;
    hyd1.fft_polar(&f1,&a1);
    hyd2.fft_polar(&f2,&a2);


  for(int k=0;k<f1.size();k++)
  {
     ffile1<<f1[k]<<" "<<k<<"\n";
     ffile2<<f2[k]<<" "<<k<<"\n";
  fcorr<<(a1[k]-a2[k])*180/PI+10<<" "<<k<<"\n";
  //tcorr<<a2[k]<<" "<<k<<"\n";
  }
    ffile1.close();

    ffile2.close();
    //tfile2.close();
   fcorr.close();
    tcorr.close();

    //addToBuffer1(i);
  return 0;
}
