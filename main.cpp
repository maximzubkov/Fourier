
#include <iostream>
#include <vector>
#include <cstdlib>
#include <string>
#include <complex>
#include <algorithm>
#include <iomanip>
#include <complex>
#include <initializer_list>
#include "fourier.h"
using namespace SubstringMatching;
using namespace FFT;
using namespace std;


int main(){
	std::string v("dokleverwhatevereverkleveruuuwnat");
	std::string c("klever");
	std::vector<size_t> r = SubstringMatching::FindSubstrings(v,c);
	for (const auto& r_elem : r){
		std::cout << r_elem << " ";
	}
	std::cout << std::endl;
	Polynomial<std::complex<double>> p1({{1, 2},{2, 3},{1, 2}});
	Polynomial<std::complex<double>> p2({{1, 2},{2, 3},{1, 2}, {4, 5}});
	Polynomial<std::complex<double>> p3(p1+p2);
	std::cout << "\nsum: " << p3 << "\n";
	Polynomial<std::complex<double>> p4(p1-p2);
	std::cout << "\ndifference: " << p4 << "\n";
	Polynomial<std::complex<double>> p5({{1, 2},{2, 3},{1, 2}});
	std::cout << "\np1 == p5 ? " << (p5 == p1) << "; p1==p2 ? " << (p1 == p2) << "\n";
	Polynomial<std::complex<double>> p6({});
	p6 += p1;
	std::cout << "+=" << p6;
	p6 -= p2;
	std::cout << "-=" << p6;
	std::vector<double> x({1, 4, 5, 9});
	//x = FFT::AddPadding(x, 2);
	x = FFT::AddPadding(x, 7);
	std::cout << "AddPadding: ";
	for (const auto& elem: x){
		std::cout << elem << " ";
	}
	std::cout << std::endl;
	std::complex<double> q = FFT::GetRoot<std::complex<double> >(4);
	std::cout << "4 degree from 1: " << q << "\n";
	std::vector<std::complex<double>> data;
	data.push_back({5,0});
	data.push_back({2,0});
	data.push_back({4,0});
	data.push_back({-1,0});
	std::cout << FFT::NearestPowerOf2((size_t)4);
	std::vector<complex<double> > aaa(FFT::FourierTransform(data));
	// FFT::VectorPrint(FFT::InverseFourierTransform(aaa));
	std::vector<complex<double> > bbb(FFT::FastFourierTransform(data));
	// FFT::VectorPrint(bbb);
	std::vector<complex<double> > ccc(FFT::FastInverseFourierTransform(bbb));
	FFT::VectorPrint(ccc);

	Polynomial<std::complex<double>> alpha({{1, 0}, {2, 0},{1, 0},{0, 0}, {0, 0}});
	Polynomial<std::complex<double>> beta({{1, 0},{2, 0},{1, 0}, {0, 0}});
	Polynomial<std::complex<double>> mul(alpha * beta);
	std::cout << "\nmul: " << mul << "\n";
	alpha = Polynomial<std::complex<double>>({{0, 0}, {10, 0},{11, 0},{12, 0}, {1, 0}});
	beta = Polynomial<std::complex<double>>({{0, 0},{1, 0},{0, 0}, {4, 0}});
	mul = alpha * beta;
	std::cout << "\nmul: " << mul << "\n";
	mul ^= 3;
	std::cout << "degree" << mul;
	// std::vector<std::complex<double>> u (res.GetCoeffs());
	// std::cout << "\n" << res << "\n";
	// // std::vector<std::complex<double>> u = res.get_poly();
	// for (const auto& u_elem : u){
	// 	std::cout << u_elem << " ";
	// }
	return 0;
}