// main.cpp -- TC2038 Lab 4.3 · Convex Hulls in Practice
// Equipo: <Integrante 1> (<matrícula>) · <Integrante 2> (<matrícula>)
//
// Banco de pruebas del equipo para las Partes 1 y 2 (no necesita CGAL).
//
// Uso (desde build/):
//   ./main                        corre todos los tests de data/ con graham y jarvis
//   ./main test                   igual que arriba
//   ./main graham < archivo.txt   imprime el hull de un archivo (formato del PDF)
//   ./main jarvis < archivo.txt   idem, y reporta el número de orient() de Jarvis
//
// Los resultados esperados vienen de la tabla de la Parte 1 del PDF.

#include "hull.hpp"
#include <fstream>
#include <iostream>
#include <string>
#include <vector>

#ifndef DATA_DIR
#define DATA_DIR "../data"   // CMake lo sobrescribe con la ruta absoluta
#endif

// ---------------------------------------------------------------- utilidades

static std::vector<Pt> read_points(std::istream& in) {
    int n = 0;
    in >> n;
    std::vector<Pt> P(n);
    for (auto& p : P) in >> p.x >> p.y;
    return P;
}

static void print_hull(const std::vector<Pt>& H) {
    std::cout << H.size() << "\n";
    for (auto& p : H) std::cout << p.x << " " << p.y << "\n";
}

static std::string to_str(const std::vector<Pt>& H) {
    std::string s = std::to_string(H.size()) + ":";
    for (auto& p : H) s += " (" + std::to_string(p.x) + "," + std::to_string(p.y) + ")";
    return s;
}

// ---------------------------------------------------------------- casos de prueba

struct TestCase {
    std::string file;           // archivo dentro de data/
    std::vector<Pt> expected;   // hull esperado según la convención del PDF
    std::string what;           // qué caso degenerado prueba
};

static const std::vector<TestCase> CASES = {
    {"example13.txt",         {{312, 56}, {684, 86}, {720, 170}, {420, 308}, {192, 146}}, "ejemplo de clase"},
    {"square_degenerate.txt", {{0, 0}, {4, 0}, {4, 4}, {0, 4}},                            "duplicados, puntos en aristas, interiores"},
    {"collinear.txt",         {{0, 0}, {3, 3}},                                            "todos colineales + duplicado"},
    {"single.txt",            {{7, -2}},                                                   "n = 1"},
    // TODO (opcional): agreguen sus propios casos aquí, p. ej. n = 0, n = 2, todos iguales...
};

// Corre un algoritmo sobre todos los casos. Regresa cuántos fallaron.
template <class Algo>
static int run_suite(const std::string& name, Algo algo) {
    int fails = 0;
    std::cout << "== " << name << " ==\n";
    for (const auto& tc : CASES) {
        std::ifstream f(std::string(DATA_DIR) + "/" + tc.file);
        if (!f) { std::cout << "  [ERROR] no se pudo abrir " << tc.file << "\n"; fails++; continue; }
        auto P = read_points(f);
        auto H = algo(P);
        bool ok = (H == tc.expected);
        if (!ok) fails++;
        std::cout << "  [" << (ok ? " OK " : "FAIL") << "] " << tc.file << "  (" << tc.what << ")\n";
        if (!ok) {
            std::cout << "         esperado: " << to_str(tc.expected) << "\n"
                      << "         obtenido: " << to_str(H) << "\n";
        }
    }
    return fails;
}

// ---------------------------------------------------------------- main

int main(int argc, char** argv) {
    std::string mode = argc > 1 ? argv[1] : "test";

    if (mode == "graham") {
        print_hull(graham(read_points(std::cin)));
        return 0;
    }
    if (mode == "jarvis") {
        long long tests = 0;
        auto H = jarvis(read_points(std::cin), &tests);
        print_hull(H);
        std::cerr << "orientation tests = " << tests << "\n";   // dato para la pregunta Q2
        return 0;
    }
    if (mode != "test") {
        std::cerr << "uso: " << argv[0] << " [test | graham | jarvis] < input\n";
        return 1;
    }

    // Sanity check rápido de orient() (Tarea 1.1)
    {
        Pt a{0, 0}, b{1, 0}, c{0, 1};
        bool ok = orient(a, b, c) > 0 && orient(a, c, b) < 0 && orient(a, b, {2, 0}) == 0;
        std::cout << "[" << (ok ? " OK " : "FAIL") << "] orient() izquierda/derecha/colineal\n";
        // TODO (Q1.1): agreguen un caso con coordenadas en ±1e9 para el peor caso de overflow
    }

    int fails = 0;
    fails += run_suite("Graham scan (Parte 1)", [](const std::vector<Pt>& P) { return graham(P); });
    fails += run_suite("Jarvis march (Parte 2)", [](const std::vector<Pt>& P) { return jarvis(P); });

    // Graham y Jarvis deben dar exactamente lo mismo (lo pide el PDF)
    std::cout << "\n" << (fails ? "FALLAN " + std::to_string(fails) + " casos" : "TODO OK") << "\n";
    return fails ? 1 : 0;
}
