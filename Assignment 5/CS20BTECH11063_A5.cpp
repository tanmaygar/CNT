#include <iostream>
#include <gmpxx.h>
#include <bits/stdc++.h>

using namespace std;
void print_poly(vector<mpz_class> &a);

// Number theory functions

vector<mpz_class> extendedEuclidean(mpz_class a, mpz_class b)
{
    mpz_class x, y, A, B, u, v;
    // A = max(a, b);
    // B = min(a, b);
    A = a;
    B = b;

    x = 1;
    y = 0;
    u = 0;
    v = 1;

    // Apply the extended euclidean algorithm
    while (B != 0)
    {
        mpz_class q = A / B;
        mpz_class r = A % B;
        mpz_class tmp1 = x - q * u;
        mpz_class tmp2 = y - q * v;
        A = B;
        B = r;
        x = u;
        y = v;
        u = tmp1;
        v = tmp2;
    }

    vector<mpz_class> result = {x, y, A};
    return result;
}

// Function to compute (a^b) mod n and check if a is invertible in mod n
void arithmetic_power_zn(mpz_class n, mpz_class a, mpz_class b)
{
    mpz_class result_a_b = a;
    while (b > 1)
    {
        result_a_b = (result_a_b * a) % n;
        b--;
    }
    cout << result_a_b;

    // Check if a is invertible mod n using the extended Euclidean algorithm
    vector<mpz_class> result = extendedEuclidean(n, a);
    if (result[2] == 1)
    {
        cout << ", true, ";
        // Calculate a^(-1) mod n
        // ax = 1 mod n and then scale x to get the positive value
        cout << (result[1] % n + n) % n;
        // cout << (extendedEuclidean(a, n)[0] % n);
    }
    else
    {
        cout << ", false";
    }
    cout << "\n";
}

// mpz_class arithmetic_power(mpz_class &a, mpz_class &b, mpz_class &n)
// {
//     mpz_class result_a_b = a;
//     while (b > 1)
//     {
//         result_a_b = (result_a_b * a) % n;
//         b--;
//     }
//     return result_a_b;
// }

// Function to compute (a^b) mod n using fast exponentiation
mpz_class arithmetic_power(mpz_class a, mpz_class b, mpz_class n)
{
    mpz_class result = 1;
    while (b > 0)
    {
        if (b % 2 == 1)
            result = (result * a) % n;
        a = (a * a) % n;
        b /= 2;
    }
    return result;
}

// Function to compute (a^b) using fast exponentiation
mpz_class normal_arithmetic_power(mpz_class a, mpz_class b)
{
    mpz_class result = 1;
    while (b > 0)
    {
        if (b % 2 == 1)
            result = (result * a);
        a = (a * a);
        b /= 2;
    }
    return result;
}

// Function to solve Chinese Remainder Theorem
mpz_class crt(mpz_class a, mpz_class m, mpz_class b, mpz_class n)
{
    // x = a (mod m) and x = b (mod n)
    // If m and n are not coprime, then return -1
    vector<mpz_class> gcd_vector = extendedEuclidean(m, n);
    if (gcd_vector[2] != 1)
    {
        return -1;
    }
    // x1 = m^-1 (mod n) and x2 = n^-1 (mod m)
    mpz_class x1 = (gcd_vector[0] % n + n) % n;
    // mpz_class x1 = (extendedEuclidean(m, n)[0] % n);
    mpz_class x2 = (gcd_vector[1] % m + m) % m;
    // mpz_class x2 = (extendedEuclidean(n, m)[0] % m);

    // Compute the solution
    mpz_class result = (m * x1 * b + n * x2 * a) % (m * n);
    return result;
}

// Function returns the smallest k such that b^(2^k) = 1 mod p
mpz_class findk(mpz_class b, mpz_class p)
{
    mpz_class k = 0;
    mpz_class temp = b;
    while (temp != 1)
    {
        temp = (temp * temp) % p;
        k++;
    }
    return k;
}

// Function to find r such that r^((p-1)/2) = -1 mod p
mpz_class findr(mpz_class p)
{
    mpz_class r = 2;
    // mpz_class temp = (p - 1) / 2;
    while (arithmetic_power(r, (p - 1) / 2, p) != p - 1)
    {
        r++;
    }
    return r;
}

// mpz_class findr(mpz_class p)
// {
//     mpz_class r = 2;
//     while (true) {
//         if (arithmetic_power(r, (p - 1) / 2, p) == p - 1) {
//             return r;
//         }
//         r++;
//     }
// }

// Function to find the solution using Tonelli-Shanks algorithm
mpz_class tonelli_shanks(mpz_class a, mpz_class p)
{
    mpz_class m = p - 1;
    mpz_class t = 0;
    while (m % 2 == 0)
    {
        m = m / 2;
        t++;
    }

    mpz_class b = arithmetic_power(a, m, p) % p;
    mpz_class k = findk(b, p);
    if (k == t)
    {
        // cout << "No solution exists\n";
        return 0;
    }

    // mpz_class k_0_case = (m + 1) / 2;
    // cout << m << " " << k << " " << k_0_case << "\n";
    mpz_class x = arithmetic_power(a, (m + 1) / 2, p) % p;
    if (k == 0)
    {
        // cout << "Solution exists\n";
        if (x > p - x)
        {
            return p - x;
        }
        return x;
    }

    mpz_class r = findr(p);
    mpz_class s = arithmetic_power(r, m, p) % p;

    // mpz_class s_power = t - k;
    // mpz_class s_two = 2;
    // mpz_class s_two_power = arithmetic_power(s_two, s_power, p);
    mpz_class S = arithmetic_power(s, normal_arithmetic_power(mpz_class(2), t - k), p) % p;

    while (k > 0)
    {
        b = (b * S) % p;
        // x = (x * S) % p;
        x = (x * arithmetic_power(s, normal_arithmetic_power(mpz_class(2), t - k - 1), p)) % p;
        k = findk(b, p);
        // mpz_class tmp_s_power = t - k - 1;
        // mpz_class tmp_s_two = 2;
        // mpz_class tmp_s_two_power = arithmetic_power(tmp_s_two, tmp_s_power, p);
        S = arithmetic_power(s, normal_arithmetic_power(mpz_class(2), t - k), p) % p;
    }
    // cout << "Normal case\n";

    // return smallest from x and p-x
    if (x > p - x)
    {
        return p - x;
    }
    return x;
}

// Polynomial Functions

// Function to add two polynomials
vector<mpz_class> add_poly(vector<mpz_class> a, vector<mpz_class> b, mpz_class &p)
{
    vector<mpz_class> result(max(a.size(), b.size()), 0);
    for (int i = 0; i < a.size(); i++)
    {
        result[i] = (result[i] + a[i] + p) % p;
    }
    for (int i = 0; i < b.size(); i++)
    {
        result[i] = (result[i] + b[i] + p) % p;
    }
    // Remove highest degree terms if they are 0
    while (result[result.size() - 1] == 0 && result.size() > 0)
    {
        result.pop_back();
    }
    // if(result[result.size() - 1] == 0)
    // {
    //     result.pop_back();
    // }
    if (result.size() == 0)
    {
        result.push_back(0);
    }
    return result;
}

// Function to subtract two polynomials
vector<mpz_class> subtract_poly(vector<mpz_class> a, vector<mpz_class> b, mpz_class &p)
{
    vector<mpz_class> result(max(a.size(), b.size()), 0);
    for (int i = 0; i < a.size(); i++)
    {
        result[i] = (result[i] + a[i] + p) % p;
    }
    for (int i = 0; i < b.size(); i++)
    {
        result[i] = (result[i] - b[i] + p) % p;
    }
    while (result[result.size() - 1] == 0 && result.size() > 0)
    {
        result.pop_back();
    }
    // if(result[result.size() - 1] == 0)
    // {
    //     result.pop_back();
    // }
    if (result.size() == 0)
    {
        result.push_back(0);
    }
    return result;
}

// Function to multiply two polynomials
vector<mpz_class> multiply_poly(vector<mpz_class> a, vector<mpz_class> b, mpz_class &p)
{
    vector<mpz_class> result(a.size() + b.size() - 1, 0);
    for (int i = 0; i < a.size(); i++)
    {
        for (int j = 0; j < b.size(); j++)
        {
            result[i + j] = (result[i + j] + (a[i] * b[j]) + p) % p;
        }
    }
    while (result[result.size() - 1] == 0 && result.size() > 0)
    {
        result.pop_back();
    }
    // if(result[result.size() - 1] == 0)
    // {
    //     result.pop_back();
    // }
    if (result.size() == 0)
    {
        result.push_back(0);
    }
    return result;
}

// Function to divide two polynomials
vector<mpz_class> divide_poly(vector<mpz_class> a, vector<mpz_class> b, mpz_class &p)
{
    vector<mpz_class> remainder = a;
    int degA = a.size() - 1;
    int degB = b.size() - 1;
    // print_poly(remainder);
    // If degreeA is less than degreeB then return 0 quotient
    if(degA < degB)
    {
        vector<mpz_class> quotient = {0};
        return quotient;
    }
    vector<mpz_class> quotient(degA - degB + 1, 0);
    
    // If b is a constant then modify A and return
    if(degB == 0)
    {
        // cout << "degB = 0\n";
        mpz_class b_inv = (extendedEuclidean(b[0], p)[0] % p + p) % p;
        // cout << b_inv << endl;
        for(int i = 0; i < a.size(); i++)
        {
            quotient[i] = (a[i] * b_inv) % p;
        }
        return quotient;
    }

    // long division of polynomial
    while (degA >= degB)
    {
        // add the coefficient of leading term to term
        mpz_class b_inv = (extendedEuclidean(b[degB], p)[0] % p + p) % p;
        mpz_class lead_coeff = (remainder[degA] * b_inv) % p;
        vector<mpz_class> term(degA - degB + 1, 0);
        term[degA - degB] = lead_coeff;
        // update quotient to include the new term
        quotient = add_poly(quotient, term, p);
        // subtract the term from the remainder and update it
        remainder = subtract_poly(remainder, multiply_poly(b, term, p), p);
        degA = remainder.size() - 1;
    }

    return quotient;
}

// Function to compute extended Euclidean algorithm for polynomials
vector<vector<mpz_class>> extendedEuclideanPoly(vector<mpz_class> a, vector<mpz_class> b, mpz_class &p)
{
    vector<vector<mpz_class>> result;
    vector<mpz_class> x, y, A, B, u, v;
    A = a;
    B = b;
    x = {1};
    y = {0};
    u = {0};
    v = {1};
    // cout << "x(x) = "; print_poly(x);
    // cout << "y(x) = "; print_poly(y);
    // cout << "u(x) = "; print_poly(u);
    // cout << "v(x) = "; print_poly(v);
    // vector<mpz_class> x(a.size(), 0), y(a.size(), 0), u(b.size(), 0), v(b.size(), 0);
    // vector<mpz_class> A = a, B = b;

    while (B.size() > 0)
    {
        // cout << "A(x) = "; print_poly(A); cout << "B(x) = "; print_poly(B);
        // print_poly(A);
        vector<mpz_class> q = divide_poly(A, B, p);
        // cout << "A(x) = "; print_poly(A); cout << "B(x) = "; print_poly(B);
        vector<mpz_class> r = subtract_poly(A, multiply_poly(q, B, p), p);
        // cout << "A(x) = "; print_poly(A); cout << "B(x) = "; print_poly(B);
        // cout << "q(x) = "; print_poly(q); cout << "r(x) = "; print_poly(r);
        vector<mpz_class> tmp1 = subtract_poly(x, multiply_poly(q, u, p), p);
        vector<mpz_class> tmp2 = subtract_poly(y, multiply_poly(q, v, p), p);

        A = B;
        B = r;
        x = u;
        y = v;
        u = tmp1;
        v = tmp2;
        // cout << "A(x) = "; print_poly(A);
        // cout << "B(x) = "; print_poly(B);
        // cout << "q(x) = "; print_poly(q);
        // cout << "r(x) = "; print_poly(r);
        // cout << "x(x) = "; print_poly(x);
        // cout << "y(x) = "; print_poly(y);
        // cout << "u(x) = "; print_poly(u);
        // cout << "v(x) = "; print_poly(v);
        // cout << "tmp1(x) = "; print_poly(tmp1);
        // cout << "tmp2(x) = "; print_poly(tmp2);
        if (B.size() == 1 && B[0] == 0)
        {
            break;
        }
        // if(A.size() == 1)
        // {
        //     break;
        // }
    }

    result.push_back(x);
    result.push_back(y);
    result.push_back(A);
    return result;
}

// Function to print polynomial
void print_poly(vector<mpz_class> &a)
{
    int degree = a.size() - 1;
    if (degree == 0)
    {
        cout << a[0] << endl;
        return;
    }
    if(degree == 1)
    {
        // 1 must not be printed before x
        if(a[1] == 1)
        {
            cout << "x";
        }
        else
        {
            cout << a[1] << "*x";
        }
        if(a[0] != 0)
        {
            cout << " + " << a[0] << endl;
        }
        else
        {
            cout << "\n";
        }
        // cout << a[1] << "x + " << a[0] << endl;
        return;
    }
    if(a[degree] == 1)
    {
        cout << "x^" << degree;
    }
    else
    {
        cout << a[degree] << "*x^" << degree;
    }
    for (int i = degree - 1; i > 1; i--)
    {
        if (a[i] != 0)
        {
            if(a[i] == 1)
            {
                cout << " + x^" << i;
            }
            else
            {
                cout << " + " << a[i] << "*x^" << i;
            }
            // cout << " + " << a[i] << "x^" << i;
        }
    }
    if (a[1] != 0)
    {
        if(a[1] == 1)
        {
            cout << " + x";
        }
        else
        {
            cout << " + " << a[1] << "*x";
        }
    }
    if (a[0] != 0)
    {
        cout << " + " << a[0];
    }
    cout << "\n";

    return;
}

void print_poly_2(vector<mpz_class> &a)
{
    int degree = a.size() - 1;
    if (degree == 0)
    {
        cout << a[0] << endl;
        return;
    }
    if(degree == 1)
    {
        // 1 must not be printed before x
        if(a[1] == 1)
        {
            cout << "x";
        }
        else
        {
            cout << a[1] << "*x";
        }
        if(a[0] != 0)
        {
            cout << " + " << a[0];
        }
        else
        {
            // cout << "\n";
        }
        // cout << a[1] << "x + " << a[0] << endl;
        return;
    }
    if(a[degree] == 1)
    {
        cout << "x^" << degree;
    }
    else
    {
        cout << a[degree] << "*x^" << degree;
    }
    for (int i = degree - 1; i > 1; i--)
    {
        if (a[i] != 0)
        {
            if(a[i] == 1)
            {
                cout << " + x^" << i;
            }
            else
            {
                cout << " + " << a[i] << "*x^" << i;
            }
            // cout << " + " << a[i] << "x^" << i;
        }
    }
    if (a[1] != 0)
    {
        if(a[1] == 1)
        {
            cout << " + x";
        }
        else
        {
            cout << " + " << a[1] << "*x";
        }
    }
    if (a[0] != 0)
    {
        cout << " + " << a[0];
    }
    // cout << "\n";

    return;
}

// Function to make the polynomial monic
vector<mpz_class> make_monic(vector<mpz_class> f, mpz_class &p)
{
    vector<mpz_class> result(f.size(), 0);
    mpz_class lead_coeff = f[f.size() - 1];
    mpz_class lead_coeff_inv = (extendedEuclidean(lead_coeff, p)[0] % p + p) % p;
    for (int i = 0; i < f.size(); i++)
    {
        result[i] = (f[i] * lead_coeff_inv) % p;
    }
    return result;
}

// Function to find derivative of a polynomial
vector<mpz_class> derivative_poly(vector<mpz_class> f, mpz_class &p)
{
    vector<mpz_class> result(f.size() - 1, 0);
    for (int i = 1; i < f.size(); i++)
    {
        result[i - 1] = ((f[i] * i) % p + p) % p;
    }

    return make_monic(result, p);
}

// // Function to find square free part of a polynomial
// vector<mpz_class> square_free_part(vector<mpz_class> f, mpz_class &p)
// {
//     vector<mpz_class> g = f;
//     vector<mpz_class> F = f;
//     vector<mpz_class> h;
//     while (g.size() > 1)
//     {
//         h = g;
//         g = derivative_poly(g, p);
//     }
//     // if h belongs to Zp then return f / gcd(f, f')
//     if (h.size() == 1)
//     {
//         vector<mpz_class> gcd = extendedEuclideanPoly(f, derivative_poly(f, p), p)[2];
//         return make_monic(divide_poly(f, gcd, p), p);
//     }
//     F = divide_poly(F, h, p);
//     // if h(x) = H(x^p) then recurse
//     // if(h.size() == 2)
//     // {
//     //     vector<mpz_class> H = {h[0], 0, h[1]};
//     //     return square_free_part(H, p);
//     // }
//     vector<mpz_class> H(h.size(), 0);
//     for (int i = 0; i < h.size(); i++)
//     {
//         mpz_class coeff = h[i];
//         for (int j = 0; j < i; j++)
//         {
//             coeff *= p;
//         }
//         H[i] = coeff;
//     }
//     if (H == h)
//     {
//         h = square_free_part(H, p);
//     }
//     F = divide_poly(multiply_poly(h, F, p), extendedEuclideanPoly(F, derivative_poly(F, p), p)[2], p);
//     F = make_monic(F, p);
//     return F;
// }

// Function to find square free part of a polynomial
vector<mpz_class> square_free_poly(vector<mpz_class> f, mpz_class &p)
{
    vector<mpz_class> g = f;
    vector<mpz_class> h;
    while(true)
    {
        h = extendedEuclideanPoly(g, derivative_poly(g, p), p)[2];
        if(h.size() == 1)
        {
            break;
        }
        g = divide_poly(g, h, p);
    }
    return make_monic(g, p);
}


// // Function to calculate f(x) mod g(x)
vector<mpz_class> modulo_poly(vector<mpz_class> f, vector<mpz_class> g, mpz_class p)
{
    vector<mpz_class> q = divide_poly(f, g, p);
    vector<mpz_class> r = subtract_poly(f, multiply_poly(q, g, p), p);
    return r;
}

// Function to calculate power of polynomial f(x) with modulo p and modulo g(x)
vector<mpz_class> power_poly(vector<mpz_class> f, mpz_class power, vector<mpz_class> g, mpz_class p)
{
    vector<mpz_class> result = {1};
    vector<mpz_class> base = f;
    while (power > 0)
    {
        if (power % 2 == 1)
        {
            vector<mpz_class> tmp = multiply_poly(result, base, p);
            result = modulo_poly(tmp, g, p);
            // print_poly(result);

            // result = extendedEuclideanPoly(multiply_poly(result, base, p), g, p)[2];
        }
        vector<mpz_class> tmp = multiply_poly(base, base, p);
        base = modulo_poly(tmp, g, p);
        // print_poly(base);
        // base = modulo_poly(multiply_poly(base, base, p), g, p);
        // base = extendedEuclideanPoly(multiply_poly(base, base, p), g, p)[2];
        power /= 2;
    }
    return result;
}

vector<vector<mpz_class>> distinct_degree_factor(vector<mpz_class> f, mpz_class p, bool is_square_free = true)
{
    if(!is_square_free)
    {
        vector<mpz_class> square_free_f = square_free_poly(f, p);
    // cout << "f(x) = ";
    // print_poly(f);
    cout << "square free f(x) = ";
    print_poly(square_free_f);
    f = square_free_f;
    }
    vector<vector<mpz_class>> result;
    mpz_class d = 1;
    // h(x) = x
    vector<mpz_class> h = {0, 1};
    while (d < f.size())
    {
        // h := h^p mod f;
        h = power_poly(h, p, f, p);
        // g := gcd(h - x, f);
        vector<mpz_class> g = extendedEuclideanPoly(subtract_poly(h, {0, 1}, p), f, p)[2];
        g = make_monic(g, p);
        result.push_back(g);
        // f = f / g;
        f = divide_poly(f, g, p);
        // if f = 1 then return (g1, . . . , gd);
        if (f.size() == 1)
        {
            return result;
        }
        d++;
    }
    return result;
}

// Function to randomly generate a polynomial of degree less than n
vector<mpz_class> generate_poly(int n, mpz_class p)
{
    vector<mpz_class> result(n + 1, 0);
    for (int i = 0; i < n + 1; i++)
    {
        result[i] = rand() % p;
    }
    // remove leading zeros
    while (result[result.size() - 1] == 0 && result.size() > 0)
    {
        result.pop_back();
    }
    return result;
}

// Function to find irreducible factor of a polynomial
set<vector<mpz_class>> get_irreducible_factor_poly(vector<mpz_class> f, int i, mpz_class p)
{
    // set<vector<mpz_class>> result;
    if (i == f.size() - 1)
    {
        set<vector<mpz_class>> result;
        result.insert(f);
        return result;
    }
    vector<mpz_class> g = generate_poly(f.size() - 1, p);
    g = make_monic(g, p);
    // cout << "g(x) = ";
    // print_poly(g);
    while(g.size() == 1)
    {
        g = generate_poly(i, p);
        g = make_monic(g, p);
        // cout << "g(x) = ";
        // print_poly(g);
    }
    vector<mpz_class> h1 = {1};
    vector<mpz_class> h2 = {1};
    while (h1.size() == 1 || h2.size() == 1)
    {
        // vector<mpz_class> g = generate_poly(i, p);
        // vector<mpz_class> g = generate_poly(i, p);
        // g = make_monic(g, p);
        // cout << "g(x) = ";
        // print_poly(g);
        // vector<mpz_class> inner_g = power_poly(g, (normal_arithmetic_power(p, i) - 1) / 2, f, p) - vector<mpz_class>{1};
        vector<mpz_class> inner_g = subtract_poly(power_poly(g, (normal_arithmetic_power(p, i) - 1) / 2, f, p), {1}, p);
        
        h1 = extendedEuclideanPoly(f, inner_g, p)[2];
        h1 = make_monic(h1, p);
        h2 = divide_poly(f, h1, p);
        h2 = make_monic(h2, p);

        // cout << "h1(x) = ";
        // print_poly(h1);
        // cout << "h2(x) = ";
        // print_poly(h2);
    }
    set<vector<mpz_class>> list_1 = get_irreducible_factor_poly(h1, i, p);
    set<vector<mpz_class>> list_2 = get_irreducible_factor_poly(h2, i, p);
    list_1.insert(list_2.begin(), list_2.end());
    return list_1;
}

int main(int argc, char const *argv[])
{
    // Opening the input file
    // ifstream fin_Q1;
    // fin_Q1.open("input-polygcd1.csv");
    argv[1] = "input-CZ.csv";
    // if(argc != 2)
    // {
    //     cout << "Error: " << argv[0] << " <input-file>\n";
    //     return 1;
    // }
    ifstream fin_Q1;
    fin_Q1.open(argv[1]);

    if (!fin_Q1.is_open())
    {
        cout << "Error in opening file\n";
        return 2;
    }

    // ofstream fout;
    // fout.open("output-gcd.csv");

    // Read each line till end of file
    string line;
    // istringstream input_line(line);
    getline(fin_Q1, line);
    // cout << line << "\n";
    mpz_class p(stoi(line));
    cout << "p = " << p << "\n";

    // run a while loop to read the lines and print the polynomial
    while (getline(fin_Q1, line))
    {
        istringstream input_line(line);
        int degree;
        input_line >> degree;
        // cout << "degree = " << degree << "\n";
        vector<mpz_class> f(degree + 1, 0);
        for (int i = 0; i <= degree; i++)
        {
            getline(input_line, line, ',');
            input_line >> f[i];
            // cout << f[i] << " ";
        }
        // cout << "\n";
        reverse(f.begin(), f.end());
        cout << "f(x) = ";
        print_poly(f);
        // vector<mpz_class> square_free_f_x = square_free_part(f, p);
        // cout << "f(x) = ";
        // print_poly(f);

        vector<vector<mpz_class>> result = distinct_degree_factor(f, p);
        // for (int i = 0; i < result.size(); i++)
        // {
        //     cout << "f_" << i + 1 << "(x) = ";
        //     print_poly(result[i]);
        // }

        // get the irreducible factors
        // vector<vector<mpz_class>> irreducible_factors;
        // for (int i = 0; i < result.size(); i++)
        // {
        //     // if ith result is just a constant polynomial then continue
        //     if (result[i].size() == 1)
        //     {
        //         continue;
        //     }
        //     set<vector<mpz_class>> irreducible_set = get_irreducible_factor_poly(result[i], i + 1, p);
        //     for (auto it = irreducible_set.begin(); it != irreducible_set.end(); it++)
        //     {
        //         irreducible_factors.push_back(*it);
        //     }
        // }
        // // print the irreducible factors
        // for (int i = 0; i < irreducible_factors.size(); i++)
        // {
        //     cout << "f_" << i + 1 << "(x) = ";
        //     print_poly(irreducible_factors[i]);
        // }

        map<int, vector<vector<mpz_class>>> irreducible_factors;
        for (int i = 0; i < result.size(); i++)
        {
            // if ith result is just a constant polynomial then continue
            if (result[i].size() == 1)
            {
                continue;
            }
            set<vector<mpz_class>> irreducible_set = get_irreducible_factor_poly(result[i], i + 1, p);
            for (auto it = irreducible_set.begin(); it != irreducible_set.end(); it++)
            {
                irreducible_factors[i + 1].push_back(*it);
            }
        }
        // print the irreducible factors
        cout << "f(x) = ";
        for (auto it = irreducible_factors.begin(); it != irreducible_factors.end(); it++)
        {
            for (int i = 0; i < it->second.size(); i++)
            {
                if (it != irreducible_factors.begin() || i != 0)
                {
                    cout << " * ";
                }
                cout << "(";
                print_poly_2(it->second[i]);
                cout << ")";
            }
        }
        cout << "\n\n";
        // break;
    }

    // getline(fin_Q1, line);
    // istringstream input_line(line);
    // int degree;
    // input_line >> degree;
    // // cout << "degree = " << degree << "\n";
    // vector<mpz_class> f(degree + 1, 0);
    // for (int i = 0; i <= degree; i++)
    // {
    //     getline(input_line, line, ',');
    //     input_line >> f[i];
    //     // cout << f[i] << " ";
    // }
    // cout << "\n";

    // getline(fin_Q1, line);
    // istringstream input_line2(line);
    // int degree2;
    // input_line2 >> degree2;
    // // cout << "degree2 = " << degree2 << "\n";
    // vector<mpz_class> g(degree2 + 1, 0);
    // for (int i = 0; i <= degree2; i++)
    // {
    //     getline(input_line2, line, ',');
    //     input_line2 >> g[i];
    //     // cout << g[i] << " ";
    // }
    // cout << "\n";

    // reverse the polynomials
    // reverse(f.begin(), f.end());
    // // reverse(g.begin(), g.end());

    // cout << "f(x) = ";
    // print_poly(f);
    // cout << "g(x) = ";
    // print_poly(g);

    // vector<vector<mpz_class>> result = extendedEuclideanPoly(f, g, p);
    // vector<mpz_class> u = result[0];
    // vector<mpz_class> v = result[1];
    // vector<mpz_class> gcd = result[2];

    // // Make the gcd monic
    // mpz_class gcd_lead = gcd[gcd.size() - 1];
    // mpz_class gcd_lead_inv = (extendedEuclidean(gcd_lead, p)[0] % p + p) % p;
    // for (int i = 0; i < gcd.size(); i++)
    // {
    //     gcd[i] = (gcd[i] * gcd_lead_inv) % p;
    // }
    // u = multiply_poly(u, {gcd_lead_inv}, p);
    // v = multiply_poly(v, {gcd_lead_inv}, p);

    // cout << "\ngcd(f(x), g(x)) = ";
    // print_poly(gcd);
    // cout << "u(x) = ";
    // print_poly(u);
    // cout << "v(x) = ";
    // print_poly(v);

    // fout.close();
    fin_Q1.close();

    return 0;
}
