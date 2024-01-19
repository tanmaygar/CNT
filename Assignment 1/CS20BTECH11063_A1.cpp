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

int main()
{
    // Opening the input file
    ifstream fin;
    fin.open("input-gcd.csv");

    if (!fin.is_open())
    {
        cout << "Error in opening file\n";
        return 1;
    }

    // ofstream fout;
    // fout.open("output-gcd.csv");

    // Read each line till end of file
    string line;
    while (getline(fin, line))
    {
        // cout<<line<<endl;
        // Parse the line from the file into a, b
        istringstream input_line(line);
        string a_str, b_str;
        getline(input_line, a_str, ',');
        getline(input_line, b_str, ',');

        // Convert a_str, b_str to mpz_class
        mpz_class a(a_str);
        mpz_class b(b_str);

        // Apply the extended euclidean algorithm
        vector<mpz_class> result = extendedEuclidean(a, b);
        mpz_class x = result[0];
        mpz_class y = result[1];
        mpz_class A = result[2];

        // Print the result
        cout << "x = " << x << ", y = " << y << ", c = " << A << "\n";
        // fout << x << "," << y << "," << A << "\n";
        // cout << "A = " << A << ", B = " << B << endl;
        // cout << "a = " << a << ", b = " << b << endl;
    }

    // fout.close();
    fin.close();

    return 0;
}
