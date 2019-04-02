#pragma once 

#include <iostream>
#include <vector>
#include <cstdlib>
#include <string>
#include <complex>
#include <algorithm>
#include <iomanip>
#include <complex>
#include <cmath>
#include <initializer_list>
using namespace std;


// Реализуте пропущенные методы
// в качестве Т будет использоваться std::complex<float> // <double> // <long double>
namespace FFT {
    // Возвращает корень степени degree из 1

    template <typename T>
    void VectorPrint(std::vector<T> vect){
        std:;cout << std::endl;
        for (const auto& elem: vect){
            std::cout << elem << " "; 
        }
        std::cout << std::endl;
    }
    template <typename T>
    T GetRoot(double degree){
        return T(cos(2 * M_PI / degree) , sin(2 * M_PI /degree));
    }

    // Выполняет преобразование фурье квадратичным алгоритмом для вектора произвольной длины
    template <typename T>
    std::vector<T> FourierTransform(const std::vector<T>& data){
        size_t N = data.size();
        T tmp(0);
        // std::vector<std::vector<T> > matrix;
        // std::vector<T> tmp_vect;
        std::vector<T> res;
        for (int i = 0; i < N; i++){
            for (int j = 0; j < N; j++){
                if (i == 0 || j == 0){
                    tmp += (T)data[j];
                    // tmp_vect.push_back((T)1);
                } else {
                    tmp += data[j] * GetRoot<T>((double)N / (j * i));
                    // tmp_vect.push_back(GetRoot<T>((double)N / (j * i)));
                }
            }
            res.push_back(tmp);
            tmp = 0;
            // VectorPrint(tmp_vect);
            // matrix.push_back(tmp_vect);
            // tmp_vect.clear();
            // std::cout << "\n";
        }
        // std::cout << "\n";
        // VectorPrint(res);
        return res;
    }

    // Выполняет обратное преобразование квадратичным алгоритмом для вектора фиксированной длины
    template <typename T>
    std::vector<T> InverseFourierTransform(const std::vector<T>& data){
        size_t N = data.size();
        // std::vector<std::vector<T> > matrix;
        // std::vector<T> tmp_vect;
        T tmp(0);
        std::vector<T> res;
        for (int i = 0; i < N; i++){
            for (int j = 0; j < N; j++){
                if (i == 0 || j == 0){
                    tmp += T(data[j].real()/N, data[j].imag()/N);
                    // tmp_vect.push_back(T(1.0/N, 0));
                } else {
                    T exp = GetRoot<T>(((-1.0) * N / (j * i)));
                    tmp += T(data[j].real()/ N , data[j].imag()/ N) * exp ;
                    // tmp_vect.push_back(T(exp.real()/ N , exp.imag()/ N ));
                }
            }
            res.push_back(tmp);
            tmp = {0, 0};
            // VectorPrint(tmp_vect);
            // matrix.push_back(tmp_vect);
            // tmp_vect.clear();
            // std::cout << "\n";
        }
        // std::cout << "\n";
        // VectorPrint(tmp_vect);
        return res;
    }

    // Добивает вектор в конце нулями до длины expected_length,
    // выбрасывает std::runtime_error если expected_length < data.size()
    template <typename T>
    std::vector<T> AddPadding(const std::vector<T>& data, size_t expected_length){
        try{
            if (expected_length < data.size()){
                 throw std::runtime_error("expected_length < size of data");
            }
        }
        catch(std::string ex){
            std::cout << ex;
            return data;
        }
        std::vector<T> res(data);
        for (int i = data.size(); i < expected_length; i++){
            res.push_back((T)0);
        }
        return res;
    }
    template <typename T>
    size_t NearestPowerOf2(const T size){
        int k = 0; // близжайшя степень двойки 
        uint64_t res = 1ULL;
        while (res < size){
            k++;
            res = res << 1ULL;
        }
        return res;
    }

    // Быстрое преобразование Фурье для вектора длины 2^k
    template <typename T>
    std::vector<T> FastFourierTransform(const std::vector<T>& data){
        if (data.size() == 1){
            return data;
        }
        size_t size = NearestPowerOf2(data.size());

        std::vector<T> extended_data = AddPadding(data, size);
        std::vector<T> a_0;
        std::vector<T> a_1;
        std::vector<T> res(size, 0);
        for (int i = 0; i < size; i += 2){
            a_0.push_back(extended_data[i]);
            a_1.push_back(extended_data[i + 1]);
        }

        std::vector<T> y_0(FastFourierTransform(a_0));
        std::vector<T> y_1(FastFourierTransform(a_1));

        for (int i = 0; i < size / 2; i++){
            res[i] = y_0[i] + GetRoot<T>((1.0) * size / i) * y_1[i];
            res[i + size / 2] = y_0[i] - GetRoot<T>((1.0) * size / i) * y_1[i];
        }

        return res;
    }

    // Обратное быстрое преобразование Фурье для вектора длины 2^k
    template <typename T>
    std::vector<T> FastInverseFourierTransform(const std::vector<T>& data){
        if (data.size() == 1){
            return data;
        }
        size_t size = NearestPowerOf2(data.size());

        std::vector<T> extended_data = AddPadding(data, size);
        std::vector<T> a_0;
        std::vector<T> a_1;
        std::vector<T> res(size, 0);
        for (int i = 0; i < size; i += 2){
            a_0.push_back(extended_data[i]);
            a_1.push_back(extended_data[i + 1]);
        }
    
        std::vector<T> y_0(FastInverseFourierTransform(a_0));
        std::vector<T> y_1(FastInverseFourierTransform(a_1));


        for (int i = 0; i < size / 2; i++){
            res[i] = y_0[i] + GetRoot<T>((-1.0) * size / i) * y_1[i];
            res[i + size / 2] = y_0[i] - GetRoot<T>((-1.0) * size / i) * y_1[i];
            res[i] /= 2;
            res[i + size / 2] /= 2;
        }

        return res;
    }
} // namespace FFT

// Операции над многочленами с помощью ффт
template <typename T>
class Polynomial {
public:
    explicit Polynomial(const std::initializer_list<T>& coefficients):
        coefficients_(coefficients.begin(), coefficients.end()), degree_(coefficients.size()) 
    {
        for(const auto& elem: coefficients_){
            std::cout<<elem;
        }
        std::cout<<std::endl;
    }


    // Чтобы можно было написать Polynomial<std::complex<double>> p({1, 2, 1})
    // И получить представление многочлена 1 + 2x + x^2
    // Polynomial(const std::initializer_list<T>& coefficients)
    //     : coefficients_(coefficients.begin(), coefficients.end()), degree_(coefficients.size()) {
    // }
    Polynomial(const std::vector<T> coefficients):
        coefficients_(coefficients.begin(), coefficients.end()), degree_(coefficients.size())  {}
    Polynomial(const Polynomial&) = default;
    Polynomial(Polynomial&&) noexcept = default;

    // Polynomial& operator=(const Polynomial) = default;
    Polynomial& operator=(Polynomial&&) noexcept = default;

    size_t GetDegree() const {
        return degree_;
    }

    std::vector<T> GetCoeffs() const{
        return coefficients_;
    }

    // Реализуйте следующие методы
    bool operator==(const Polynomial& other){
        std::vector<T> v1 = other.GetCoeffs();
        std::vector<T> v2 = this->GetCoeffs();
        if (v1.size() != v2.size()){
            return false;
        }
        for (int i = 0; i < v1.size(); i++){
            if (v1[i] != v2[i]){
                return false;
            }
        }
        return true;
    }

    // Используя предыдущее
    Polynomial operator+(const Polynomial<T>& other) const{
        std::vector<T> v1 = other.GetCoeffs();
        std::vector<T> v2 = this->GetCoeffs();
        std::vector<T> res = v1.size() < v2.size() ? v2 : v1;
        std::vector<T> v_min = v1.size() < v2.size() ? v1 : v2;
        for (int i = 0; i < v_min.size(); i++){
            res[i] += v_min[i];
        }
        Polynomial<T> sum (res);
        return sum;
    }

    Polynomial operator-(const Polynomial<T>& other) const{
        std::vector<T> v = other.GetCoeffs();
        std::vector<T> res = this->GetCoeffs();
        int max = v.size() < res.size() ? res.size() : v.size();
        for (int i = 0; i < max; i++){
            if (i < res.size()){
                res[i] -= v[i];
            } else {
                res.push_back(-v[i]);
            }
        }
        Polynomial<T> dif (res);
        return dif;
    }

    Polynomial& operator+=(const Polynomial& other){
        *this = *this + other;
        return *this;
    }

    Polynomial& operator-=(const Polynomial& other){
        *this = *this - other;
        return *this;
    }

    Polynomial operator*(const Polynomial& other){
        std::vector<T> a(other.GetCoeffs());
        std::vector<T> b(this->GetCoeffs());
        size_t max = a.size() > b.size() ? a.size() : b.size();
        max = FFT::NearestPowerOf2(max * 2);
        std::vector<T> dft_a = FFT::FastFourierTransform(FFT::AddPadding(a, max));
        std::vector<T> dft_b = FFT::FastFourierTransform(FFT::AddPadding(b, max));
        std::vector<T> inv_res;
        for (int i = 0; i < max; i++){
            inv_res.push_back(dft_a[i] * dft_b[i]);
        }
        return Polynomial(FFT::FastInverseFourierTransform(inv_res));
    }

    Polynomial operator^(size_t pow);
    

    // Возведение в степень pow с помощью комбинации FFT и индийского возведения в степень
    Polynomial& operator^=(size_t pow){
        if (pow == 0){
            std::vector<T> v;
            v.push_back(1);
            Polynomial *e = new Polynomial(v);
            return *e;
        }
        Polynomial tmp = *this;
        Polynomial e = *this;
        while (pow > 0){
            if (pow == 1){
                *this *= e;
            }else {
                if (pow % 2 == 1)
                    *this *= tmp;
            }
            tmp *= tmp;
            pow /= 2;
            std::cout << "pow:  "<< pow << ", poly: "<< *this << "\n";
        }
        return *this;
    }

    // С Использованием быстрого преобразования фурье
    Polynomial& operator*=(const Polynomial& other){
        *this = *this * other;
        return *this;
    }


    friend std::ostream& operator<< (std::ostream &out, const Polynomial& poly){
        std::vector<T> coef(poly.GetCoeffs());
        int i;
        out << coef[0] << " + ";
        for (i = 1; i != static_cast<int>(coef.size()) - 1; i++){
            if (coef[i].real() * coef[i].real() + coef[i].imag() * coef[i].imag() > 1e-4)
            out << coef[i] << "x^" << i << " + ";
        }
        out << coef[i] << "x^" << i << endl;
        return out;
    }
    // И еще один, унарный минус
    friend Polynomial operator-(const Polynomial& other){
        std::vector<T> v = other.GetCoeffs();
        for(int i = v.begin(); i < v.end(); i++){
            v[i] = -v[i];
        }
        return Polynomial(v);
    }

private:
    std::vector<T> coefficients_;
    size_t degree_;
};

// Напишите какие-то тесты, демонстрирующие корректность написанных вами методов



// Задачи, решаемые с помощью умножения многочленов
// Если вы напишете решение, работающее для произольных строк над ascii кодировкой - укажете это и
// возможно, получите небольшой бонусный балл
namespace SubstringMatching {

    // Метод принимает две строки str и pattern, возвращает индексов,
    // указывающих начала подстрок str, совпадающих с pattern
    // Можете считать, что str и pattern состоят из символов 0 и 1
    std::vector<size_t> FindSubstrings(const std::string& str, const std::string& pattern);

    // Аналогично предыдущему, но теперь pattern может содержать символ '?', на месте которого
    // можно поставить любой символ алфавита
    std::vector<size_t> FindMatches(const std::string& str, const std::string& pattern);

} // namespace SubstringMatcher



// Напишите какие-то тесты, демонстрирующие корректность написанных вами методов
// Обратите внимание, что при проведении прямого и обратного преобразования для многочлена с
// целыми коэффициентами, новые коэффициенты могут выйти дробными и трубующими округления.