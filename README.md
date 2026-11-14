# Lab 4.3 · Convex Hulls in Practice (TC2038)

Equipo: **Diego Rodríguez Sánchez (A01742085)** · **Omar Avila Meza (A01743547)**  — Entrega: vie 13 nov 2026, 23:59

## Compilar

```bash
mkdir build && cd build
cmake ..          # Release por defecto
make -j4
./main            # corre los tests de Graham y Jarvis contra data/
```

Sin CGAL solo se compilan `main`, `hull_main` y `app`; con CGAL (`sudo apt install libcgal-dev` / `brew install cgal`) también `cgal_hello`, `stress`, `bench`, `robust` y `hull3d`.

## Reparto sugerido

| Parte | Archivo | Responsable | Estado |
|---|---|---|---|
| 1 · orient + Graham | `hull.hpp`, `hull.cpp` | | ⬜ |
| 2 · Jarvis | `hull.cpp` | | ⬜ |
| 3 · Stress vs CGAL | `stress.cpp` | | ⬜ |
| 4 · Benchmark + plots | `bench.cpp`, `results.csv` | | ⬜ |
| 5 · Perímetro, área, calipers | `app.cpp` | | ⬜ |
| 6 · Robustez | `robust.cpp` (correr y explicar) | | ⬜ |
| Reporte (PDF 3–5 págs) | `report/` | ambos | ⬜ |

## Flujo con git

```bash
git pull                       # antes de empezar
git checkout -b parte-1        # una rama por parte
git add . && git commit -m "Parte 1: orient y graham"
git push -u origin parte-1     # luego PR / merge a main
```

Instrucciones completas en `docs/lab-convex-hull.pdf`.
