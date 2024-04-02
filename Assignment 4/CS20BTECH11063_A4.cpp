#include <iostream>
#include <gmpxx.h>
#include <bits/stdc++.h>

using namespace std;

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

void print_poly(vector<mpz_class> &a);

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
    while(result[result.size() - 1] == 0 && result.size() > 0)
    {
        result.pop_back();
    }
    // if(result[result.size() - 1] == 0)
    // {
    //     result.pop_back();
    // }
    if(result.size() == 0)
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
    while(result[result.size() - 1] == 0 && result.size() > 0)
    {
        result.pop_back();
    }
    // if(result[result.size() - 1] == 0)
    // {
    //     result.pop_back();
    // }
    if(result.size() == 0)
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
    while(result[result.size() - 1] == 0 && result.size() > 0)
    {
        result.pop_back();
    }
    // if(result[result.size() - 1] == 0)
    // {
    //     result.pop_back();
    // }
    if(result.size() == 0)
    {
        result.push_back(0);
    }
    return result;
}

vector<mpz_class> divide_poly(vector<mpz_class> a, vector<mpz_class> b, mpz_class &p)
{
    vector<mpz_class> remainder = a;
    int degA = a.size() - 1;
    int degB = b.size() - 1;
    // print_poly(remainder);
    if(degA < degB)
    {
        vector<mpz_class> quotient = {0};
        return quotient;
    }
    vector<mpz_class> quotient(degA - degB + 1, 0);
    
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

    while (degA >= degB)
    {
        mpz_class b_inv = (extendedEuclidean(b[degB], p)[0] % p + p) % p;
        mpz_class lead_coeff = (remainder[degA] * b_inv) % p;
        vector<mpz_class> term(degA - degB + 1, 0);
        term[degA - degB] = lead_coeff;
        quotient = add_poly(quotient, term, p);
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
        if(B.size() == 1 && B[0] == 0)
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

int main(int argc, char const *argv[])
{
    // Opening the input file
    // ifstream fin_Q1;
    // fin_Q1.open("input-polygcd1.csv");
    if(argc != 2)
    {
        cout << "Error: " << argv[0] << " <input-file>\n";
        return 1;
    }
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

    getline(fin_Q1, line);
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


    getline(fin_Q1, line);
    istringstream input_line2(line);
    int degree2;
    input_line2 >> degree2;
    // cout << "degree2 = " << degree2 << "\n";
    vector<mpz_class> g(degree2 + 1, 0);
    for (int i = 0; i <= degree2; i++)
    {
        getline(input_line2, line, ',');
        input_line2 >> g[i];
        // cout << g[i] << " ";
    }
    // cout << "\n";

    // reverse the polynomials
    reverse(f.begin(), f.end());
    reverse(g.begin(), g.end());

    cout << "f(x) = ";
    print_poly(f);
    cout << "g(x) = ";
    print_poly(g);
    
    vector<vector<mpz_class>> result = extendedEuclideanPoly(f, g, p);
    vector<mpz_class> u = result[0];
    vector<mpz_class> v = result[1];
    vector<mpz_class> gcd = result[2];

    // Make the gcd monic
    mpz_class gcd_lead = gcd[gcd.size() - 1];
    mpz_class gcd_lead_inv = (extendedEuclidean(gcd_lead, p)[0] % p + p) % p;
    for (int i = 0; i < gcd.size(); i++)
    {
        gcd[i] = (gcd[i] * gcd_lead_inv) % p;
    }
    u = multiply_poly(u, {gcd_lead_inv}, p);
    v = multiply_poly(v, {gcd_lead_inv}, p);

    cout << "\ngcd(f(x), g(x)) = ";
    print_poly(gcd);
    cout << "u(x) = ";
    print_poly(u);
    cout << "v(x) = ";
    print_poly(v);

    // fout.close();
    fin_Q1.close();

    return 0;
}
