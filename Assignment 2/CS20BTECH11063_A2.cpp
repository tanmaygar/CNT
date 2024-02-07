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
    if (extendedEuclidean(n, a)[2] == 1)
    {
        cout << ", true, ";
        // Calculate a^(-1) mod n
        cout << (extendedEuclidean(a, n)[0] % n + n) % n;
        // cout << (extendedEuclidean(a, n)[0] % n);
    }
    else
    {
        cout << ", false";
    }
    cout << "\n";
}

// Function to solve Chinese Remainder Theorem
mpz_class crt(mpz_class &a, mpz_class &m, mpz_class &b, mpz_class &n)
{
    // x = a (mod m) and x = b (mod n)
    // If m and n are not coprime, then return -1
    if (extendedEuclidean(m, n)[2] != 1)
    {
        return -1;
    }
    // x1 = m^-1 (mod n) and x2 = n^-1 (mod m)
    mpz_class x1 = (extendedEuclidean(m, n)[0] % n + n) % n;
    // mpz_class x1 = (extendedEuclidean(m, n)[0] % n);
    mpz_class x2 = (extendedEuclidean(n, m)[0] % m + m) % m;
    // mpz_class x2 = (extendedEuclidean(n, m)[0] % m);

    // Compute the solution
    mpz_class result = (m * x1 * b + n * x2 * a) % (m * n);
    return result;
}

int main()
{
    // Opening the input file
    cout << "---------------------- Question 1 -----------------------\n";
    ifstream fin_Q1;
    fin_Q1.open("testinput-Zn.txt");

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
        string a_str, b_str, n_str;
        getline(input_line, n_str, ',');
        getline(input_line, a_str, ',');
        getline(input_line, b_str, ',');

        // Convert a_str, b_str to mpz_class
        mpz_class n(n_str);
        mpz_class a(a_str);
        mpz_class b(b_str);

        arithmetic_power_zn(n, a, b);
    }

    // fout.close();
    fin_Q1.close();
    cout << "---------------------- Question 1 Ended-----------------------\n";
    cout << "---------------------- Question 2 -----------------------\n";

    ifstream fin_Q2;
    fin_Q2.open("testinput-crt.txt");

    if (!fin_Q2.is_open())
    {
        cout << "Error in opening file\n";
        return 1;
    }

    while (getline(fin_Q2, line))
    {
        istringstream input_line(line);
        string a_str, m_str, b_str, n_str;
        getline(input_line, a_str, ',');
        getline(input_line, m_str, ',');
        getline(input_line, b_str, ',');
        getline(input_line, n_str, ',');

        mpz_class a(a_str);
        mpz_class m(m_str);
        mpz_class b(b_str);
        mpz_class n(n_str);

        cout << crt(a, m, b, n) << endl;
    }

    fin_Q2.close();

    cout << "---------------------- Question 2 Ended-----------------------\n";
    return 0;
}
