#include <iostream>
#include <gmpxx.h>
#include <bits/stdc++.h>

using namespace std;

vector<mpz_class> extendedEuclidean(mpz_class &a, mpz_class &b)
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
void arithmetic_power_zn(mpz_class &n, mpz_class &a, mpz_class &b)
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

// Function to compute (a^b) mod n
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

// Function to compute (a^b)
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
mpz_class crt(mpz_class &a, mpz_class &m, mpz_class &b, mpz_class &n)
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
    mpz_class temp = (p - 1) / 2;
    while (arithmetic_power(r, temp, p) != p - 1)
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
mpz_class tonelli_shanks(mpz_class &a, mpz_class &p)
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

int main()
{
    // Opening the input file

    ifstream fin_Q1;
    fin_Q1.open("inputSquareRoots.csv");

    if (!fin_Q1.is_open())
    {
        cout << "Error in opening file\n";
        return 1;
    }

    // ofstream fout;
    // fout.open("output-gcd.csv");

    // Read each line till end of file
    string line;
    while (getline(fin_Q1, line))
    {
        // cout<<line<<endl;
        // Parse the line from the file into a, b
        istringstream input_line(line);
        string a_str, p_str;
        getline(input_line, a_str, ',');
        getline(input_line, p_str, ',');

        // Convert a_str, b_str to mpz_class
        mpz_class a(a_str);
        mpz_class p(p_str);

        mpz_class ans = tonelli_shanks(a, p);
        cout << ans << "\n";
    }

    // fout.close();
    fin_Q1.close();

    return 0;
}
