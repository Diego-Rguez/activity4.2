# Lab 4.3 · Convex Hulls in Practice — Reporte

**TC2038 Análisis y Diseño de Algoritmos Avanzados**
Equipo: Diego Rodríguez Sánchez (A01742085) · Omar Avila Meza (A01743547)

> Este archivo contiene las partes 0, 1, 2 y 3. Las partes 4, 5 y 6 se agregan en sus secciones.

---

## Parte 0 · Setup

CGAL se instaló con MSYS2 (`pacman -S mingw-w64-ucrt-x86_64-cgal`). Salida de `cgal_hello` con `example13.txt`:

```
5
192 146
312 56
684 86
720 170
420 308
```

Son los 5 vértices esperados en orden antihorario, pero CGAL empieza en `(192,146)` y no en el punto más bajo `(312,56)`.

### Pregunta 0
Para comparar dos implementaciones de forma automática, su salida tiene que ser **idéntica elemento por elemento**. Un mismo polígono convexo se puede escribir de muchas maneras: puede empezar en cualquiera de sus h vértices, recorrerse en sentido horario o antihorario, e incluir o no los puntos colineales que están en medio de una arista. Si no se fija una convención, un `==` marcaría como distintos dos resultados que son correctos. Por eso fijamos: sólo puntos extremos, sin duplicados, orden antihorario y empezando en el punto más bajo (el de menor `y`; en empate, el de menor `x`). Así basta normalizar la salida de CGAL (con un `std::rotate`) para poder comparar miles de casos sin revisarlos a mano.

---

## Parte 1 · Prueba de orientación y Graham scan

### Implementación
- **`orient(a, b, c)`** calcula `(bx − ax)(cy − ay) − (by − ay)(cx − ax)` sólo con `long long`. Regresa el valor completo y no sólo el signo, porque su magnitud (el doble del área del triángulo abc) se usa en los rotating calipers de la Parte 5.
- **`graham(P)`**:
  1. `sort` + `unique` + `erase` con `operator<` para eliminar duplicados. Si quedan menos de 3 puntos, se regresan tal cual (como ya están ordenados, el más bajo queda primero).
  2. Gracias a ese orden, `P[0]` ya es el pivote p0 (el más bajo y, en empate, el de más a la izquierda).
  3. Se ordena el resto por ángulo alrededor de p0 con un comparador basado en `orient`: `a` va antes que `b` si `orient(p0, a, b) > 0`; si son colineales, va primero el más cercano (`dist2`).
  4. Recorrido con una pila `H`: antes de meter `p`, se saca el tope mientras `orient(H[-2], H[-1], p) ≤ 0`. El `≤` elimina también los puntos colineales que quedan en medio de una arista.

Los 4 archivos de prueba (`example13`, `square_degenerate`, `collinear` y `single`) dan exactamente la salida esperada.

### Pregunta 1.1 · ¿Puede desbordarse `orient`?
Con |x|, |y| ≤ 10⁹:

| Término | Valor máximo |
|---|---|
| Cada diferencia, p. ej. `bx − ax` | 10⁹ − (−10⁹) = **2·10⁹** |
| Cada producto, p. ej. `(bx − ax)(cy − ay)` | (2·10⁹)² = **4·10¹⁸** |
| El resultado (resta de dos productos) | 4·10¹⁸ + 4·10¹⁸ = **8·10¹⁸** |

Un `long long` con signo de 64 bits llega hasta 2⁶³ − 1 ≈ **9.22·10¹⁸**. Como 8·10¹⁸ < 9.22·10¹⁸, ni los productos intermedios ni el resultado se desbordan (aunque el margen es pequeño).

Si las coordenadas llegaran a **10¹⁰**, las diferencias serían de hasta 2·10¹⁰ y los productos de hasta 4·10²⁰, lo cual **sí** desborda 64 bits. La solución sería calcular los productos con `__int128` (hasta ~1.7·10³⁸) y comparar sólo el signo. **No** se debe usar `double`, porque perdería exactitud: sólo representa enteros exactos hasta 2⁵³ ≈ 9·10¹⁵.

### Pregunta 1.2 · ¿Qué puntos saca la pila en el ejemplo de clase?
Orden angular alrededor de p0 = (312,56): p0, p1, p2, …, p12. En este ejemplo coincide con la numeración de la clase.

| Al llegar | Se saca (pop) | Por qué |
|---|---|---|
| p3 (720,170) | p2 (636,104) | p1 → p2 → p3 gira a la derecha |
| p5 (576,194) | p4 (612,152) | |
| p6 (492,218) | p5 (576,194) | |
| p9 (402,224) | p8 (378,146), luego p7 (420,200) | dos pops seguidos |
| p10 (420,308) | p9 (402,224), luego p6 (492,218) | dos pops seguidos |
| p12 (192,146) | p11 (288,170) | |

En total hay **8 pops**, en el orden **p2, p4, p5, p8, p7, p9, p6, p11**. Quedan en la pila **p0, p1, p3, p10, p12** = (312,56) (684,86) (720,170) (420,308) (192,146), que es el hull esperado.

**¿Por qué el recorrido es O(n) aunque el `while` saque varios puntos de golpe?** Por análisis amortizado: cada punto entra a la pila **exactamente una vez** y sale **a lo más una vez**. Entonces el número total de pops en todo el recorrido es ≤ n, aunque un solo paso pueda sacar varios. Hay n pushes y ≤ n pops, así que el recorrido cuesta O(n). El costo total de Graham, O(n log n), lo domina el ordenamiento.

---

## Parte 2 · Jarvis march

### Implementación
Mismo preprocesamiento que Graham (quitar duplicados y tratar el caso n < 3). Desde `start = P[0]` (el punto más bajo), en cada vuelta se mete `current` al hull y se busca el siguiente vértice: se toma un candidato `next` y, para cada punto `r ≠ current`, se calcula `orient(current, next, r)`:
- si es **< 0** (r queda a la derecha de current → next), `r` es el nuevo candidato;
- si es **= 0** (colineal) y `r` está **más lejos** de `current`, también se toma `r`, para no reportar puntos en medio de una arista.

El ciclo (`do … while`) termina al regresar a `start`. Cada llamada a `orient` se cuenta y el total se guarda en `*tests`. Las salidas de Jarvis y Graham son idénticas en los 4 archivos de prueba.

### Pregunta 2
Con `example13.txt` (n = 13, h = 5) nuestro Jarvis hace **60 pruebas de orientación**:

```
Get-Content data\example13.txt | .\lab.exe jarvis   →   orientation tests = 60
```

Cada una de las **h** vueltas (una por vértice del hull) compara el candidato contra los otros **n − 1** puntos, así que:

**tests = h · (n − 1) = 5 · 12 = 60**

**Peor caso:** cuando **todos los puntos son extremos** (h = n), por ejemplo puntos sobre una parábola `y = x²` o sobre un círculo. Ahí el algoritmo hace n · (n − 1) pruebas, o sea **Θ(n²)**. Para n = 60 000 (la parábola del benchmark) serían ≈ 3.6·10⁹ pruebas de orientación.

---

## Parte 3 · Stress test contra CGAL

### Implementación
- **`cgal_hull`**: después de `CGAL::convex_hull_2`, se rota el resultado para que empiece en el punto más bajo:
  `std::rotate(H.begin(), std::min_element(H.begin(), H.end()), H.end());`
  (`min_element` usa nuestro `operator<`, que define "más bajo, luego más a la izquierda").
- **Casos difíciles:** además de los puntos aleatorios, según la semilla se agregan
  - `seed % 3 == 1` → **duplicados**: se vuelven a insertar hasta 5 puntos que ya existen;
  - `seed % 3 == 2` → **colineales en una arista**: los puntos (−1000 + 500k, −10⁹) para k = 0..4, que quedan sobre el borde inferior. Sólo los dos extremos deben aparecer en el hull.

Resultado con 3000 entradas (1000 semillas × 3 distribuciones):

```
OK: 0 mismatches
```

### Pregunta 3
**Bugs encontrados:** ninguno. Graham y Jarvis ya manejaban los dos casos degenerados desde su diseño: los duplicados se eliminan con `sort` + `unique` antes de empezar, y los colineales se excluyen con la condición `≤ 0` en la pila de Graham y con la regla "quedarse con el más lejano" en Jarvis.

**¿Por qué es seguro convertir las coordenadas a `double` antes de llamar a CGAL?** Un `double` tiene 53 bits de mantisa, así que representa **exactamente** cualquier entero con |k| ≤ 2⁵³ ≈ 9.007·10¹⁵. Nuestras coordenadas cumplen |x|, |y| ≤ 10⁹ < 2⁵³, por lo que la conversión no pierde información y CGAL recibe exactamente los mismos puntos.

**¿Qué kernel usamos y qué significa?** `CGAL::Exact_predicates_inexact_constructions_kernel` (EPICK):
- **Predicados exactos:** las preguntas que regresan una decisión discreta (como la orientación izquierda/derecha/colineal) siempre dan la respuesta correcta. CGAL las evalúa primero con aritmética de intervalos en `double` y, sólo si el resultado es ambiguo, las recalcula con aritmética exacta (GMP/MPFR).
- **Construcciones inexactas:** las operaciones que crean objetos geométricos nuevos (por ejemplo, el punto de intersección de dos segmentos) se calculan en `double` y pueden tener error de redondeo.

El convex hull sólo usa **predicados** (nunca construye puntos nuevos, sólo selecciona puntos de la entrada), así que con EPICK el resultado es exacto y además rápido.
