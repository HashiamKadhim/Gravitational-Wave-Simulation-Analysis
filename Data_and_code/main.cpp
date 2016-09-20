#include <fstream>
#include <rarray>
#include <rarrayio>
#include <complex>
#include <fftw3.h>
#include <cblas.h>
#include <string>
#include <sstream>


int main() {

	// open the file
	std::ifstream f("GWprediction.rat");
	
	rarray<double,1> times_pred;
	rarray<std::complex<double>,1> signal_pred;
	
	// read in the signal
	f >> times_pred;
	f >> signal_pred;
	
	//fourier transform
	int n = signal_pred.extent(0);
	rarray<std::complex<double>,1> four_transf_pred(n);
	fftw_plan plan = fftw_plan_dft_1d(signal_pred.size(),
	(fftw_complex*)signal_pred.data(), (fftw_complex*)four_transf_pred.data(),
	FFTW_FORWARD, FFTW_ESTIMATE);
	fftw_execute(plan);
	fftw_destroy_plan(plan);
	
	//power spectrum
	rarray<double,1>  pow_spec_pred(n);
	for (int ii=0; ii<n;ii++){
		pow_spec_pred[ii]=std::norm(four_transf_pred[ii]); //norm actually gives the L^2 norm squared, which is what we want
	}
	
	
	// In this loop, we do the same for the mock signals, and then find the correlation.
	for (int jj=1; jj<33;jj++){
		
		if(jj < 10){
		
		int mynumb1;
		int mynumb2;
		
		mynumb1=0;
		mynumb2=jj;
		
		std::string mynumbstring;
	
		std::ostringstream convert;   // stream used for the conversion
	
		convert << mynumb1<<mynumb2;
		mynumbstring = convert.str();
		
		std::string generalstr="detection";
		std::string generalstr_end=".rat";
		
		std::string fullstr=generalstr+mynumbstring+generalstr_end;
		
		
		//second file to open
		std::ifstream ff; 
		ff.open(fullstr.c_str());
		// create empty arrays
		rarray<double,1> times;
		rarray<std::complex<double>,1> signal;
		// read in the signal
		ff >> times;
		ff >> signal;
		
		//fourier transform
		rarray<std::complex<double>,1> four_transf(n);
		fftw_plan plan = fftw_plan_dft_1d(signal.size(),
		(fftw_complex*)signal.data(), (fftw_complex*)four_transf.data(),
		FFTW_FORWARD, FFTW_ESTIMATE);
		fftw_execute(plan);
		fftw_destroy_plan(plan);
		
		//power spectrum
		rarray<double,1>  pow_spec(n);
		for (int ii=0; ii<n;ii++){
			pow_spec[ii]=std::norm(four_transf[ii]);
		}
		
		double correlation=(cblas_ddot(n,pow_spec_pred.data(),1,pow_spec.data(),1))/(sqrt(cblas_ddot(n,pow_spec.data(),1,pow_spec.data(),1)*cblas_ddot(n,pow_spec_pred.data(),1,pow_spec_pred.data(),1)));
		
		std::cout<<"correlation of the power spectra of the prediction and  signal  " <<jj<<" is	"<<correlation << "\n";
		
		
		}
		
		else{
				
		int mynumb1;
		
		mynumb1=jj;
		
		std::string mynumbstring;
	
		std::ostringstream convert;   // stream used for the conversion
	
		convert << mynumb1;
		mynumbstring = convert.str();
		
		std::string generalstr="detection";
		std::string generalstr_end=".rat";
		
		std::string fullstr=generalstr+mynumbstring+generalstr_end;
		
		
		//second file to open
		std::ifstream ff; 
		ff.open(fullstr.c_str());
		// create empty arrays
		rarray<double,1> times;
		rarray<std::complex<double>,1> signal;
		// read in the signal
		ff >> times;
		ff >> signal;
		
		//fourier transform
		rarray<std::complex<double>,1> four_transf(n);
		fftw_plan plan = fftw_plan_dft_1d(signal.size(),
		(fftw_complex*)signal.data(), (fftw_complex*)four_transf.data(),
		FFTW_FORWARD, FFTW_ESTIMATE);
		fftw_execute(plan);
		fftw_destroy_plan(plan);
		
		//power spectrum
		rarray<double,1>  pow_spec(n);
		for (int ii=0; ii<n;ii++){
			pow_spec[ii]=std::norm(four_transf[ii]);
		}
		
		double correlation=(cblas_ddot(n,pow_spec_pred.data(),1,pow_spec.data(),1))/(sqrt(cblas_ddot(n,pow_spec.data(),1,pow_spec.data(),1)*cblas_ddot(n,pow_spec_pred.data(),1,pow_spec_pred.data(),1)));
		
		std::cout<<"correlation of the power spectra of the prediction and  signal  " <<jj<<" is	"<<correlation << "\n";
		}
			
	
		
	}
	
	

}
