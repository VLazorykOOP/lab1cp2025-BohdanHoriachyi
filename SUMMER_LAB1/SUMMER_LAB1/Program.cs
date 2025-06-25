using System;
using System.Collections.Generic;
using System.Globalization;
using System.IO;

class Program
{
    static Dictionary<double, double> Ux = new()
    {
        {-5.0, 0.2801}, {-4.5, 0.2093}, {-4.0, 0.6190}, {-3.5, 0.8811},
        {-3.0, 1.0422}, {-2.5, 1.1463}, {-2.0, 1.2176}, {-1.5, 1.2560},
        {-1.0, 1.2908}, {-0.5, 1.2909}, {0.0, 1.2039}, {0.5, 0.8196},
        {1.0, 0.4083}, {1.5, 0.4054}, {2.0, 0.7487}, {2.5, 0.9607},
        {3.0, 1.0926}, {3.5, 1.0962}, {4.0, 1.1803}, {4.5, 1.2418}, {5.0, 1.2338}
    };

    static Dictionary<double, double> Tx = new()
    {
        {-10.0, 0.7832}, {-9.0, 1.1063}, {-8.0, 1.2486}, {-7.0, 1.1587},
        {-6.0, 1.0033}, {-5.0, 0.2801}, {-4.0, 0.6190}, {-3.0, 1.0422},
        {-2.0, 1.2176}, {-1.0, 1.1998}, {0.0, 1.2039}, {1.0, 0.5187},
        {2.0, 0.4095}, {3.0, 0.4054}, {4.0, 1.1803}, {5.0, 1.2338},
        {6.0, 0.7681}, {7.0, 0.7068}, {8.0, 0.1450}, {9.0, 0.8533}, {10.0, 1.1347}
    };

    static Dictionary<string, double> TextX = new()
    {
        {"act", 1.175}, {"bet", 1.278}, {"cet", 1.381}, {"set", 1.484},
        {"get", 1.587}, {"ret", 1.69}, {"het", 1.793}, {"met", 1.896},
        {"net", 1.999}, {"qet", 2.102}, {"tet", 2.205}, {"wet", 2.308},
        {"yet", 2.411}, {"iet", 2.514}, {"pet", 2.617}, {"det", 2.826},
        {"fet", 2.923}, {"let", 3.029}, {"zet", 3.132}, {"vet", 3.235}
    };

    static double FindU(double x) => Ux.ContainsKey(x) ? Ux[x] : 1.0;
    static double FindT(double x) => Tx.ContainsKey(x) ? Tx[x] : 1.0;
    static double FindTextX(string text) => TextX.ContainsKey(text) ? TextX[text] : 1.0;

    static double ApplyAllAlgorithms(double x, double y, double z, string text)
    {
        double Rok = Math.Sin(x + y) + Math.Cos(x * x - y);
        double Qek = x > y && x * y != 0 ? Math.Pow(x + y, 2) + 0.5 * x :
                      x < y ? Math.Pow(x - y, 3) + y : 1;
        double Wek = x * x - z * z + y;

        if (x > 0 && y > 0)
        {
            Rok = Math.Pow(Rok, 2) + Math.Pow(FindTextX(text), 2);
        }

        if (Math.Abs(z) < 2.15)
        {
            Rok = Math.Pow(Rok + 2, 2) + Math.Sin(z);
        }
        else if (Math.Abs(z) < 3.45)
        {
            Rok = Math.Pow(Rok - 1, 3) + Math.Cos(z);
        }
        else
        {
            Rok = Math.Sin(Rok) + Math.Pow(z, 3);
        }

        if (z < 0) Wek = -Wek;

        double U = FindU(x);
        double T = FindT(z);
        double V = FindTextX(text);

        return Rok + Qek + Wek + U + T + V;
    }

    static void ProcessFile(string filename)
    {
        Console.WriteLine($"\nОбробка файлу: {filename}");
        foreach (var line in File.ReadAllLines(filename))
        {
            if (string.IsNullOrWhiteSpace(line)) continue;

            var parts = line.Split(new[] { ' ', '\t' }, StringSplitOptions.RemoveEmptyEntries);
            if (parts.Length != 4) continue;

            double x = double.Parse(parts[0], CultureInfo.InvariantCulture);
            double y = double.Parse(parts[1], CultureInfo.InvariantCulture);
            double z = double.Parse(parts[2], CultureInfo.InvariantCulture);
            string text = parts[3];

            double result = ApplyAllAlgorithms(x, y, z, text);
            Console.WriteLine($"f({x}, {y}, {z}, {text}) = {result:F4}");
        }
    }

    static void Main()
    {
        string[] files = { "dat1.dat", "dat2.dat", "dat3.dat" };
        foreach (var file in files)
        {
            if (File.Exists(file))
                ProcessFile(file);
            else
                Console.WriteLine($"Файл {file} не знайдено!");
        }
    }
}
