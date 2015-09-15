/* to compile this, use
gcc testfft.c /home/tzuchen/fftw-3.1.2/.libs/libfftw3.a -lm
*/



#include </home/tzuchen/fftw/include/fftw3.h>
     
int main(){
         fftw_complex *in, *out;
         fftw_plan p;
         int       i,N = 100;
 

         in  = fftw_malloc(N*sizeof(fftw_complex));
         out = fftw_malloc(N*sizeof(fftw_complex));
   
         for(i=0;i<N;i++){  in[i][0] = i*0.1; in[i][1]= 0; }
       

         in = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * N);
         out = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * N);
         p = fftw_plan_dft_1d(N, in, out, FFTW_FORWARD, FFTW_ESTIMATE);
         
         fftw_execute(p); /* repeat as needed */
         
         fftw_destroy_plan(p);
         fftw_free(in); fftw_free(out);

return 0;

  }
