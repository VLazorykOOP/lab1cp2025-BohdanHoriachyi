#include <iostream>
#include <fstream>
#include <sstream>
#include <cmath>
#include <map>
#include <string>
#include <vector>
#include <stdexcept>

using namespace std;

map<double, double> Ux = {
    {-5.0, 0.2801}, {-4.5, 0.2093}, {-4.0, 0.6190}, {-3.5, 0.8811}, {-3.0, 1.0422},
    {-2.5, 1.1463}, {-2.0, 1.2176}, {-1.5, 1.2560}, {-1.0, 1.1998}, {-0.5, 1.1209},
    {0.0, 1.0039}, {0.5, 0.8196}, {1.0, 0.5187}, {1.5, 0.0707}, {2.0, 0.4054},
    {2.5, 0.7487}, {3.0, 0.9603}, {3.5, 1.0926}, {4.0, 1.1803}, {4.5, 1.2418},
    {5.0, 1.2338}
};

map<double, double> Tx = {
    {-10.0, 0.7832}, {-9.0, 1.1063}, {-8.0, 1.2486}, {-7.0, 1.1587}, {-6.0, 0.9105},
    {-5.0, 0.2801}, {-4.0, 0.6190}, {-3.0, 1.0422}, {-2.0, 1.2176}, {0.0, 1.0039},
    {0.5, 0.5187}, {1.0, 0.4054}, {3.0, 0.9603}, {4.0, 1.1803}, {5.0, 1.2338},
    {6.0, 1.0761}, {7.0, 0.7068}, {8.0, 0.1450}, {9.0, 0.8533}, {10.0, 1.1347}
};

map<string, double> textX = {
    {"aet", 1.175}, {"bet", 1.278}, {"cet", 1.381}, {"set", 1.484}, {"get", 1.587},
    {"ret", 1.690}, {"het", 1.793}, {"met", 1.896}, {"net", 1.999}, {"qet", 2.102},
    {"tet", 2.205}, {"wet", 2.308}, {"yet", 2.411}, {"iet", 2.514}, {"oet", 2.617},
    {"pet", 2.720}, {"det", 2.823}, {"fet", 2.926}, {"let", 3.029}, {"zet", 3.132},
    {"vet", 3.235}
};

double interpolate(map<double, double>& table, double x) {
    auto it = table.lower_bound(x);
    if (it == table.end()) return prev(it)->second;
    if (it == table.begin()) return it->second;
    if (it->first == x) return it->second;
    auto it2 = it;
    --it;
    double x0 = it->first, y0 = it->second;
    double x1 = it2->first, y1 = it2->second;
    return y0 + (y1 - y0) * (x - x0) / (x1 - x0);
}

double Qnk(double x, double y) {
    return pow(x, 2) + pow(y, 2) + sin(x * y);
}

double Qnr(double x, double y) {
    if (x == 0) return 1;
    if (y == 0 && x != 0) return log(10) * 3 * x;
    if (y == x) return 2 * 10 * x * x;
    if (y == -x) return 3.75 * x * x * 2;
    return 0.0;
}

double Rnk(double x, double y, double z) {
    return Qnk(x, y) + Qnk(y, z) + Qnk(z, x);
}

double Wnk(double x, double y) {
    return Qnk(x, y) + Qnk(y, x) + Qnk(x, x);
}

double Wnr(double x, double y) {
    return interpolate(Tx, x) * interpolate(Tx, y) - interpolate(Ux, x) * interpolate(Ux, y);
}

double Gnk(double x, double y, double z) {
    double sum = x * x + y * y + z * z;
    if (sum < 0.001) return Wnk(z, y);
    return Wnk(x, x) + Wnk(y, x) + Wnk(z, y);
}

double gold(double x, double y, double z) {
    return x + y + z + Gnk(x, y, z);
}

double Tfun(double u, double v, const string& text) {
    double r = 0.0;
    auto it = textX.find(text);
    if (it != textX.end()) {
        r = it->second;
    }
    return r + u * u + v * v - r;
}

double func(double u, double v, const string& text) {
    if (u < 0.5) return Tfun(v, 0, text);
    if (u >= 0.5 && text.empty()) return Tfun(v, u, text);
    if (u >= 0.5 && !text.empty()) return Tfun(u, 0, text);
    return 0.0;
}

int main() {
    double x, y;
    string text;
    cout << "Enter x, y, text: ";
    cin >> x >> y >> text;

    try {
        double result = func(x, y, text);
        cout << "\nResult: " << result << endl;
    }
    catch (const exception& e) {
        cerr << "Exception: " << e.what() << endl;
    }

    return 0;
}
