#include <iostream>
#include <b118/linspace.hpp>
#include <vector>
#include <array>
#include <list>
#include <deque>
#include <forward_list>

template<typename Real>
void example_std_vector(Real a, Real b, std::size_t n);

template<std::size_t N, typename Real>
void example_std_array(Real a, Real b);

template<typename Real>
void example_std_list(Real a, Real b, std::size_t n);

template<typename Real>
void example_std_deque(Real a, Real b, std::size_t n);

template<typename Real>
void example_std_forward_list(Real a, Real b, std::size_t n);

int main() {
    using namespace std;
    using namespace b118;

    cout << "Room b118 Math Library\n";
    cout << "----------------------" << endl;
    cout << "B118_VERSION:       " << B118_VERSION << endl;
    cout << "B118_VERSION_MAJOR: " << B118_VERSION_MAJOR << endl;
    cout << "B118_VERSION_MINOR: " << B118_VERSION_MINOR << endl;
    cout << "B118_VERSION_PATCH: " << B118_VERSION_PATCH << endl;
    cout << "B118_VERSION_TWEAK: " << B118_VERSION_TWEAK << endl;
    cout << "\nlinspace example" << endl;

    size_t constexpr n = 5;
    double const a = -2.1;
    double const b =  9.9;

    example_std_vector(a, b, n);
    example_std_array<n>(a, b);
    example_std_list(a, b, n);
    example_std_deque(a, b, n);
    example_std_forward_list(a, b, n);

    return 0;
}

template<typename Real>
void example_std_vector(Real a, Real b, std::size_t n) {
    using namespace std;
    using namespace b118;

    cout << "\nEXEMPLO STD_VECTOR\n";
    cout << "------------------\n";

    
    // 1. LINSPACE PASSANDO ITERADORES
    cout << "1. LINSPACE PASSANDO ITERADORES\n";
    {
        cout << "\t1.a. sem excluir os pontos:\n";
        vector<Real> x(n);
        cout << "\tx = ";
        linspace(x.begin(), x.end(), a, b);
        for (auto& elem : x)
            cout << elem << ' ';
        cout << endl;
    }
    {
        cout << "\t1.b. excluindo \'a\':\n";
        vector<Real> x(n);
        cout << "\tx = ";
        linspace(x.begin(), x.end(), a, b, exclude::left);
        for (auto& elem : x)
            cout << elem << ' ';
        cout << endl;
    }
    {
        cout << "\t1.c. excluindo \'b\':\n";
        vector<Real> x(n);
        cout << "\tx = ";
        linspace(x.begin(), x.end(), a, b, exclude::right);
        for (auto& elem : x)
            cout << elem << ' ';
        cout << endl;
    }
    {
        cout << "\t1.d. excluindo ambos \'a\' e \'b\':\n";
        vector<Real> x(n);
        cout << "\tx = ";
        linspace(x.begin(), x.end(), a, b, exclude::both);
        for (auto& elem : x)
            cout << elem << ' ';
        cout << endl;
    }

    // 2. LINSPACE PASSANDO CONTAINERS
    cout << "2. LINSPACE PASSANDO CONTAINERS\n";
    {
        cout << "\t2.a. sem excluir os pontos:\n";
        vector<Real> x(n);
        cout << "\tx = ";
        linspace(x, a, b);
        for (auto& elem : x)
            cout << elem << ' ';
        cout << endl;
    }
    {
        cout << "\t2.b. excluindo \'a\':\n";
        vector<Real> x(n);
        cout << "\tx = ";
        linspace(x, a, b, exclude::left);
        for (auto& elem : x)
            cout << elem << ' ';
        cout << endl;
    }
    {
        cout << "\t2.c. excluindo \'b\':\n";
        vector<Real> x(n);
        cout << "\tx = ";
        linspace(x, a, b, exclude::right);
        for (auto& elem : x)
            cout << elem << ' ';
        cout << endl;
    }
    {
        cout << "\t2.d. excluindo ambos \'a\' e \'b\':\n";
        vector<Real> x(n);
        cout << "\tx = ";
        linspace(x, a, b, exclude::both);
        for (auto& elem : x)
            cout << elem << ' ';
        cout << endl;
    }

    // 3. LINSPACE RETORNANDO CONTAINERS
    cout << "3. LINSPACE RETORNANDO CONTAINERS\n";
    {
        cout << "\t3.a. sem excluir os pontos:\n";
        cout << "\tx = ";
        auto x = linspace<vector, Real>(a, b, n);
        for (auto& elem : x)
            cout << elem << ' ';
        cout << endl;
    }
    {
        cout << "\t3.b. excluindo \'a\':\n";
        cout << "\tx = ";
        auto x = linspace<vector, Real>(a, b, n, exclude::left);
        for (auto& elem : x)
            cout << elem << ' ';
        cout << endl;
    }
    {
        cout << "\t3.c. excluindo \'b\':\n";
        cout << "\tx = ";
        auto x = linspace<vector, Real>(a, b, n, exclude::right);
        for (auto& elem : x)
            cout << elem << ' ';
        cout << endl;
    }
    {
        cout << "\t3.d. excluindo ambos \'a\' e \'b\':\n";
        cout << "\tx = ";
        auto x = linspace<vector, Real>(a, b, n, exclude::both);
        for (auto& elem : x)
            cout << elem << ' ';
        cout << endl;
    }
}

template<
    std::size_t N,
    typename Real
>
void example_std_array(Real a, Real b) {
    using namespace std;
    using namespace b118;
    std::size_t constexpr n = N;

    cout << "\nEXEMPLO STD_ARRAY\n";
    cout << "-----------------\n";
    
    // 1. LINSPACE PASSANDO ITERADORES
    cout << "1. LINSPACE PASSANDO ITERADORES\n";
    {
        cout << "\t1.a. sem excluir os pontos:\n";
        array<Real, N> x;
        cout << "\tx = ";
        linspace(x.begin(), x.end(), a, b);
        for (auto& elem : x)
            cout << elem << ' ';
        cout << endl;
    }
    {
        cout << "\t1.b. excluindo \'a\':\n";
        array<Real, N> x;
        cout << "\tx = ";
        linspace(x.begin(), x.end(), a, b, exclude::left);
        for (auto& elem : x)
            cout << elem << ' ';
        cout << endl;
    }
    {
        cout << "\t1.c. excluindo \'b\':\n";
        array<Real, N> x;
        cout << "\tx = ";
        linspace(x.begin(), x.end(), a, b, exclude::right);
        for (auto& elem : x)
            cout << elem << ' ';
        cout << endl;
    }
    {
        cout << "\t1.d. excluindo ambos \'a\' e \'b\':\n";
        array<Real, N> x;
        cout << "\tx = ";
        linspace(x.begin(), x.end(), a, b, exclude::both);
        for (auto& elem : x)
            cout << elem << ' ';
        cout << endl;
    }

    // 2. LINSPACE PASSANDO CONTAINERS
    cout << "2. LINSPACE PASSANDO CONTAINERS\n";
    {
        cout << "\t2.a. sem excluir os pontos:\n";
        array<Real, N> x;
        cout << "\tx = ";
        linspace(x, a, b);
        for (auto& elem : x)
            cout << elem << ' ';
        cout << endl;
    }
    {
        cout << "\t2.b. excluindo \'a\':\n";
        array<Real, N> x;
        cout << "\tx = ";
        linspace(x, a, b, exclude::left);
        for (auto& elem : x)
            cout << elem << ' ';
        cout << endl;
    }
    {
        cout << "\t2.c. excluindo \'b\':\n";
        array<Real, N> x;
        cout << "\tx = ";
        linspace(x, a, b, exclude::right);
        for (auto& elem : x)
            cout << elem << ' ';
        cout << endl;
    }
    {
        cout << "\t2.d. excluindo ambos \'a\' e \'b\':\n";
        array<Real, N> x;
        cout << "\tx = ";
        linspace(x, a, b, exclude::both);
        for (auto& elem : x)
            cout << elem << ' ';
        cout << endl;
    }

    // 3. LINSPACE RETORNANDO CONTAINERS
    cout << "3. LINSPACE RETORNANDO CONTAINERS\n";
    {
        cout << "\t3.a. sem excluir os pontos:\n";
        cout << "\tx = ";
        auto x = linspace<array, Real, N>(a, b);
        for (auto& elem : x)
            cout << elem << ' ';
        cout << endl;
    }
    {
        cout << "\t3.b. excluindo \'a\':\n";
        cout << "\tx = ";
        auto x = linspace<array, Real, N>(a, b, exclude::left);
        for (auto& elem : x)
            cout << elem << ' ';
        cout << endl;
    }
    {
        cout << "\t3.c. excluindo \'b\':\n";
        cout << "\tx = ";
        auto x = linspace<array, Real, N>(a, b, exclude::right);
        for (auto& elem : x)
            cout << elem << ' ';
        cout << endl;
    }
    {
        cout << "\t3.d. excluindo ambos \'a\' e \'b\':\n";
        cout << "\tx = ";
        auto x = linspace<array, Real, N>(a, b, exclude::both);
        for (auto& elem : x)
            cout << elem << ' ';
        cout << endl;
    }
}


template<typename Real>
void example_std_list(Real a, Real b, std::size_t n) {
    using namespace std;
    using namespace b118;

    cout << "\nEXEMPLO STD_LIST\n";
    cout << "----------------\n";

    
    // 1. LINSPACE PASSANDO ITERADORES
    cout << "1. LINSPACE PASSANDO ITERADORES\n";
    {
        cout << "\t1.a. sem excluir os pontos:\n";
        list<Real> x(n);
        cout << "\tx = ";
        linspace(x.begin(), x.end(), a, b);
        for (auto& elem : x)
            cout << elem << ' ';
        cout << endl;
    }
    {
        cout << "\t1.b. excluindo \'a\':\n";
        list<Real> x(n);
        cout << "\tx = ";
        linspace(x.begin(), x.end(), a, b, exclude::left);
        for (auto& elem : x)
            cout << elem << ' ';
        cout << endl;
    }
    {
        cout << "\t1.c. excluindo \'b\':\n";
        list<Real> x(n);
        cout << "\tx = ";
        linspace(x.begin(), x.end(), a, b, exclude::right);
        for (auto& elem : x)
            cout << elem << ' ';
        cout << endl;
    }
    {
        cout << "\t1.d. excluindo ambos \'a\' e \'b\':\n";
        list<Real> x(n);
        cout << "\tx = ";
        linspace(x.begin(), x.end(), a, b, exclude::both);
        for (auto& elem : x)
            cout << elem << ' ';
        cout << endl;
    }

    // 2. LINSPACE PASSANDO CONTAINERS
    cout << "2. LINSPACE PASSANDO CONTAINERS\n";
    {
        cout << "\t2.a. sem excluir os pontos:\n";
        list<Real> x(n);
        cout << "\tx = ";
        linspace(x, a, b);
        for (auto& elem : x)
            cout << elem << ' ';
        cout << endl;
    }
    {
        cout << "\t2.b. excluindo \'a\':\n";
        list<Real> x(n);
        cout << "\tx = ";
        linspace(x, a, b, exclude::left);
        for (auto& elem : x)
            cout << elem << ' ';
        cout << endl;
    }
    {
        cout << "\t2.c. excluindo \'b\':\n";
        list<Real> x(n);
        cout << "\tx = ";
        linspace(x, a, b, exclude::right);
        for (auto& elem : x)
            cout << elem << ' ';
        cout << endl;
    }
    {
        cout << "\t2.d. excluindo ambos \'a\' e \'b\':\n";
        list<Real> x(n);
        cout << "\tx = ";
        linspace(x, a, b, exclude::both);
        for (auto& elem : x)
            cout << elem << ' ';
        cout << endl;
    }

    // 3. LINSPACE RETORNANDO CONTAINERS
    cout << "3. LINSPACE RETORNANDO CONTAINERS\n";
    {
        cout << "\t3.a. sem excluir os pontos:\n";
        cout << "\tx = ";
        auto x = linspace<list, Real>(a, b, n);
        for (auto& elem : x)
            cout << elem << ' ';
        cout << endl;
    }
    {
        cout << "\t3.b. excluindo \'a\':\n";
        cout << "\tx = ";
        auto x = linspace<list, Real>(a, b, n, exclude::left);
        for (auto& elem : x)
            cout << elem << ' ';
        cout << endl;
    }
    {
        cout << "\t3.c. excluindo \'b\':\n";
        cout << "\tx = ";
        auto x = linspace<list, Real>(a, b, n, exclude::right);
        for (auto& elem : x)
            cout << elem << ' ';
        cout << endl;
    }
    {
        cout << "\t3.d. excluindo ambos \'a\' e \'b\':\n";
        cout << "\tx = ";
        auto x = linspace<list, Real>(a, b, n, exclude::both);
        for (auto& elem : x)
            cout << elem << ' ';
        cout << endl;
    }
}

template<typename Real>
void example_std_deque(Real a, Real b, std::size_t n) {
    using namespace std;
    using namespace b118;

    cout << "\nEXEMPLO STD_DEQUE\n";
    cout << "-----------------\n";

    
    // 1. LINSPACE PASSANDO ITERADORES
    cout << "1. LINSPACE PASSANDO ITERADORES\n";
    {
        cout << "\t1.a. sem excluir os pontos:\n";
        deque<Real> x(n);
        cout << "\tx = ";
        linspace(x.begin(), x.end(), a, b);
        for (auto& elem : x)
            cout << elem << ' ';
        cout << endl;
    }
    {
        cout << "\t1.b. excluindo \'a\':\n";
        deque<Real> x(n);
        cout << "\tx = ";
        linspace(x.begin(), x.end(), a, b, exclude::left);
        for (auto& elem : x)
            cout << elem << ' ';
        cout << endl;
    }
    {
        cout << "\t1.c. excluindo \'b\':\n";
        deque<Real> x(n);
        cout << "\tx = ";
        linspace(x.begin(), x.end(), a, b, exclude::right);
        for (auto& elem : x)
            cout << elem << ' ';
        cout << endl;
    }
    {
        cout << "\t1.d. excluindo ambos \'a\' e \'b\':\n";
        deque<Real> x(n);
        cout << "\tx = ";
        linspace(x.begin(), x.end(), a, b, exclude::both);
        for (auto& elem : x)
            cout << elem << ' ';
        cout << endl;
    }

    // 2. LINSPACE PASSANDO CONTAINERS
    cout << "2. LINSPACE PASSANDO CONTAINERS\n";
    {
        cout << "\t2.a. sem excluir os pontos:\n";
        deque<Real> x(n);
        cout << "\tx = ";
        linspace(x, a, b);
        for (auto& elem : x)
            cout << elem << ' ';
        cout << endl;
    }
    {
        cout << "\t2.b. excluindo \'a\':\n";
        deque<Real> x(n);
        cout << "\tx = ";
        linspace(x, a, b, exclude::left);
        for (auto& elem : x)
            cout << elem << ' ';
        cout << endl;
    }
    {
        cout << "\t2.c. excluindo \'b\':\n";
        deque<Real> x(n);
        cout << "\tx = ";
        linspace(x, a, b, exclude::right);
        for (auto& elem : x)
            cout << elem << ' ';
        cout << endl;
    }
    {
        cout << "\t2.d. excluindo ambos \'a\' e \'b\':\n";
        deque<Real> x(n);
        cout << "\tx = ";
        linspace(x, a, b, exclude::both);
        for (auto& elem : x)
            cout << elem << ' ';
        cout << endl;
    }

    // 3. LINSPACE RETORNANDO CONTAINERS
    cout << "3. LINSPACE RETORNANDO CONTAINERS\n";
    {
        cout << "\t3.a. sem excluir os pontos:\n";
        cout << "\tx = ";
        auto x = linspace<deque, Real>(a, b, n);
        for (auto& elem : x)
            cout << elem << ' ';
        cout << endl;
    }
    {
        cout << "\t3.b. excluindo \'a\':\n";
        cout << "\tx = ";
        auto x = linspace<deque, Real>(a, b, n, exclude::left);
        for (auto& elem : x)
            cout << elem << ' ';
        cout << endl;
    }
    {
        cout << "\t3.c. excluindo \'b\':\n";
        cout << "\tx = ";
        auto x = linspace<deque, Real>(a, b, n, exclude::right);
        for (auto& elem : x)
            cout << elem << ' ';
        cout << endl;
    }
    {
        cout << "\t3.d. excluindo ambos \'a\' e \'b\':\n";
        cout << "\tx = ";
        auto x = linspace<deque, Real>(a, b, n, exclude::both);
        for (auto& elem : x)
            cout << elem << ' ';
        cout << endl;
    }
}

template<typename Real>
void example_std_forward_list(Real a, Real b, std::size_t n) {
    using namespace std;
    using namespace b118;

    cout << "\nEXEMPLO STD_FORWARD_LIST\n";
    cout << "------------------------\n";

    
    // 1. LINSPACE PASSANDO ITERADORES
    cout << "1. LINSPACE PASSANDO ITERADORES\n";
    {
        cout << "\t1.a. sem excluir os pontos:\n";
        forward_list<Real> x(n);
        cout << "\tx = ";
        linspace(x.begin(), x.end(), a, b);
        for (auto& elem : x)
            cout << elem << ' ';
        cout << endl;
    }
    {
        cout << "\t1.b. excluindo \'a\':\n";
        forward_list<Real> x(n);
        cout << "\tx = ";
        linspace(x.begin(), x.end(), a, b, exclude::left);
        for (auto& elem : x)
            cout << elem << ' ';
        cout << endl;
    }
    {
        cout << "\t1.c. excluindo \'b\':\n";
        forward_list<Real> x(n);
        cout << "\tx = ";
        linspace(x.begin(), x.end(), a, b, exclude::right);
        for (auto& elem : x)
            cout << elem << ' ';
        cout << endl;
    }
    {
        cout << "\t1.d. excluindo ambos \'a\' e \'b\':\n";
        forward_list<Real> x(n);
        cout << "\tx = ";
        linspace(x.begin(), x.end(), a, b, exclude::both);
        for (auto& elem : x)
            cout << elem << ' ';
        cout << endl;
    }

    // 2. LINSPACE PASSANDO CONTAINERS
    cout << "2. LINSPACE PASSANDO CONTAINERS\n";
    {
        cout << "\t2.a. sem excluir os pontos:\n";
        forward_list<Real> x(n);
        cout << "\tx = ";
        linspace(x, a, b);
        for (auto& elem : x)
            cout << elem << ' ';
        cout << endl;
    }
    {
        cout << "\t2.b. excluindo \'a\':\n";
        forward_list<Real> x(n);
        cout << "\tx = ";
        linspace(x, a, b, exclude::left);
        for (auto& elem : x)
            cout << elem << ' ';
        cout << endl;
    }
    {
        cout << "\t2.c. excluindo \'b\':\n";
        forward_list<Real> x(n);
        cout << "\tx = ";
        linspace(x, a, b, exclude::right);
        for (auto& elem : x)
            cout << elem << ' ';
        cout << endl;
    }
    {
        cout << "\t2.d. excluindo ambos \'a\' e \'b\':\n";
        forward_list<Real> x(n);
        cout << "\tx = ";
        linspace(x, a, b, exclude::both);
        for (auto& elem : x)
            cout << elem << ' ';
        cout << endl;
    }

    // 3. LINSPACE RETORNANDO CONTAINERS
    cout << "3. LINSPACE RETORNANDO CONTAINERS\n";
    {
        cout << "\t3.a. sem excluir os pontos:\n";
        cout << "\tx = ";
        auto x = linspace<forward_list, Real>(a, b, n);
        for (auto& elem : x)
            cout << elem << ' ';
        cout << endl;
    }
    {
        cout << "\t3.b. excluindo \'a\':\n";
        cout << "\tx = ";
        auto x = linspace<forward_list, Real>(a, b, n, exclude::left);
        for (auto& elem : x)
            cout << elem << ' ';
        cout << endl;
    }
    {
        cout << "\t3.c. excluindo \'b\':\n";
        cout << "\tx = ";
        auto x = linspace<forward_list, Real>(a, b, n, exclude::right);
        for (auto& elem : x)
            cout << elem << ' ';
        cout << endl;
    }
    {
        cout << "\t3.d. excluindo ambos \'a\' e \'b\':\n";
        cout << "\tx = ";
        auto x = linspace<forward_list, Real>(a, b, n, exclude::both);
        for (auto& elem : x)
            cout << elem << ' ';
        cout << endl;
    }
}