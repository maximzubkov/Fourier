
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

void test_on_FindSubstrings(){
	std::cout << "----------------------------test on FindSubstrings-----------------------------------------------\n";
	std::string text("ive was bor ina new was cos dos plu sss waswaswaswas gae");
	std::string pattern("was");
	std::vector<size_t> r = SubstringMatching::FindSubstrings(text, pattern);
	std::cout << "text: " << text << "\n" << "pattern: " << pattern << "\n";
	for (const auto& r_elem : r){
		std::cout << r_elem << " ";
	}
	std::cout << "\n\n\n\n";
}


void test_on_FindMatches(){
	std::cout << "----------------------------test on FindMatches----------------------------------------------------\n";
	std::string text("ive w*s bor ina *** was cos **s dos plu sss waswaswa*w*s gae ************** hghghghpgshgslghshlgshlgsgslhghlgslhshlgshlsghlsghsglbhjrgivniuglkrh fulhvnsghlavnfklsvndnlvdsfkdnvfdhl");
	std::string pattern("***");
	std::cout << "text: " << text << "\n" << "pattern: " << pattern << "\n";
	std::cout << "expected "<< text.size() << " matchings" << "\n";
	std::vector<size_t> r = SubstringMatching::FindMatches(text, pattern);
	for (const auto& r_elem : r){
		std::cout << r_elem << " ";
	}
	std::cout << "\n\n";

	pattern = "was";
	std::cout << "text: " << text << "\n" << "pattern: " << pattern;
	r = SubstringMatching::FindMatches(text, pattern);
	for (const auto& r_elem : r){
		std::cout << r_elem << " ";
	}
	std::cout << "\n\n\n\n";
}

void test_on_Polynomial(){
	std::cout << "----------------------------test on Polynomial-----------------------------------------------------\n";
	std::cout << "first poly creaded: \n";
	Polynomial<std::complex<double>> p1({{1, 2},{2, 3},{1, 2}});
	std::cout << "\n";
	std::cout << "second poly creaded: \n";
	Polynomial<std::complex<double>> p2({{1, 2},{2, 3},{1, 2}, {4, 5}});
	Polynomial<std::complex<double>> p3(p1+p2);
	std::cout << "\nsum: " << p3 << "\n";
	Polynomial<std::complex<double>> p4(p1-p2);
	std::cout << "\ndifference: " << p4 << "\n";
	std::cout << "another poly creaded: \n";
	Polynomial<std::complex<double>> p5({{1, 2},{2, 3},{1, 2}});
	std::cout << "\np1 == p5 ? " << (p5 == p1) << "; p1==p2 ? " << (p1 == p2) << "\n";
	Polynomial<std::complex<double>> p6({});
	p6 += p1;
	std::cout << "+= operation: " << p6 << "\n";
	p6 -= p2;
	std::cout << "-= operation: " << p6;

	std::vector<std::complex<double>> data;
	data.push_back({5,0});
	data.push_back({2,0});
	data.push_back({4,0});
	data.push_back({-1,0});

	std::cout << "vector data:";
	for (const auto & d: data){
		std::cout << d;
	}
	std::cout << "\n FourierTransform result:";
	std::vector<complex<double> > aaa(FFT::FourierTransform(data));
	FFT::VectorPrint(aaa);
	std::cout << "\n InverseFourierTransform result:";
	std::vector<complex<double> > bbb(FFT::InverseFourierTransform(aaa));
	FFT::VectorPrint(bbb);
	std::cout << "\n FastFourierTransform result:";
	std::vector<complex<double> > ccc(FFT::FastFourierTransform(data));
	FFT::VectorPrint(ccc);
	std::cout << "\n FastInverseFourierTransform result:";
	std::vector<complex<double> > ddd(FFT::FastInverseFourierTransform(bbb));
	FFT::VectorPrint(ddd);

	Polynomial<std::complex<double>> alpha({{1, 0}, {2, 0},{1, 0},{0, 0}, {0, 0}});
	Polynomial<std::complex<double>> beta({{1, 0},{2, 0},{1, 0}, {0, 0}});
	Polynomial<std::complex<double>> mul(alpha * beta);
	
	std::cout << "\nmul: " << mul << "\n";
	alpha = Polynomial<std::complex<double>>({{0, 0}, {10, 0},{11, 0},{12, 0}, {1, 0}});
	beta = Polynomial<std::complex<double>>({{0, 0},{1, 0},{0, 0}, {4, 0}});
	mul *= beta;
	std::cout << "\n *=: " << mul << "\n";
	mul ^= 3;
	std::cout << "3 degree: " << mul;
	std::cout << "3 degree: \n";
	Polynomial<std::complex<double>> deg(alpha^3);
	std::cout << deg;
}
int main(){
	test_on_Polynomial();
	test_on_FindSubstrings();
	test_on_FindMatches();
	// // std::vector<std::complex<double>> u (res.GetCoeffs());
	// // std::cout << "\n" << res << "\n";
	// // std::vector<std::complex<double>> u = res.get_poly();
	// for (const auto& u_elem : u){
	// 	std::cout << u_elem << " ";
	// }
	return 0;
}