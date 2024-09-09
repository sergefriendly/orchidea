# Тестовое задание

  Дан идеальный математический маятник – груз массой $m = 1,0$ кг, прикреплённый к шарниру невесомым нерастяжимым стержнем длиной $L = 5$ м. На маятник действуют сила сопротивления стержня $T$ и сила тяжести $mg(t)$, где $g(t) = 9,81+0,05\sin{2πt}$ – переменное ускорение свободного падения, $t$ – время. Пренебрегая трением в шарнире и сопротивлением воздуха, получаем систему уравнений, описывающую движение груза в декартовой системе координат:  
```math
m\ddot{x}=-\frac{x}{L}T,
```
```math
m\ddot{y}=-\frac{y}{L}T-mg,
```
```math
g=9,81+0,05\,\sin{2πt},
```
```math
x^2+y^2=L^2,
```
где $x$ и $y$ – проекции радиус-вектора груза на оси $OX$ и $OY$ соответственно.  

  Пусть в начальный момент времени груз толкают из точки $x=3$ м, $y=−4$ м вниз и влево со скоростью $1$ м/c. Используя **любой язык программирования, численно решите** эту задачу Коши для дифференциально-алгебраической системы уравнений на интервале $\left[0; 2\right].$ Решите задачу на большем интервале, например, $\left[0; 100\right],$ постройте график $\sqrt{x^2(t)+y^2(t)}$ и проанализируйте его.
  
  **Запрещается** выполнять переход к полярной системе координат или исключать переменную $T$ из системы уравнений: выполнение таких преобразований в автоматическом режиме является нетривиальной задачей, а подобные системы уравнений часто возникают на практике, особенно в случае использования компонентного (физического, объектно-ориентированного) подхода к моделированию сложных динамических систем. Вы должны решить задачу без использования каких-либо особенностей физики процесса.
  
  **Запрещается** применять **готовые** реализации алгоритмов интегрирования задачи Коши – **решатели**. Если необходимо, **рекомендуется** использовать сторонние **библиотеки** для решения **вспомогательных задач**, например, для умножения матриц, решения систем линейных и нелинейных алгебраических уравнений.

# Решение
## Физический анализ явления
<p align="center"><img src="https://github.com/sergefriendly/orchidea/raw/main/drawing.PNG" alt="Чертёж"></p>
По Второму закону Ньютона имеем:  

```math

m \frac{d^2 \vec{r} }{dt^2} = \vec{T} + m\vec{g}.

```

В проекциях на оси $OX$ и $OY$ соответственно:

```math

\begin{cases}
  -m\ddot{x} = -T\,\sin\,\alpha \\
  -m\ddot{y} = T\,\cos\,\alpha - mg
\end{cases}

```
и далее, учитывая условия задачи будем иметь:

```math
\begin{cases}
  m\ddot{x} = \frac{x}{L}T \\
  m\ddot{y} = -\frac{y}{L}T + mg \\
  g=9,81+0,05\,\sin\,2πt \\
  x^2+y^2=L^2
\end{cases}
```
Подводя итоги, решаемости системы алгебро-дифференциальных уравнений численным интегрированием заметим:
1. По начальным условиям $x(0)=3$, $y(0)=-4$.
2. Также по начальным условиям $v(0)=1$ и основываясь на утверждении, что «груз толкают вниз влево», предположим, что $\dot{x}(0) = -\frac{v}{\sqrt{2}}$ и $\dot{y}(0) = -\frac{v}{\sqrt{2}}$.
3. В системе алгебро-дифференциальных уравнений известными постоянными будут велечины $L$, $m$, известными переменными $t$, $g$, а неизвестными переменными — $x$, $y$ и также величина $T$. Получается — 4 уравнения, из которых 2 дифференциальных, в которых 3 переменных в них. Необходимо дополнительное условие для $T$. Сила сопротивления стержня $\vec{T}$, фигурирующая в системе выше, в минимальном положении маятника по модулю полностью равна силе тяжести $mg$, что замечается из приведённого чертежа. В промежуточных же положениях маятника, она меньше по модулю и определяется проекциями $\frac{x}{L}T$ и $\frac{y}{L}T$ что тоже, что и $\frac{x}{L}mg$ и $\frac{y}{L}mg$.

Финальная система алгебро-дифференциальных уравнений будет следующая:

```math
\begin{cases}
  m\ddot{x} = \frac{x}{L}T \\
  m\ddot{y} = -\frac{y}{L}T + mg \\
  T=mg \\
  g=9,81+0,05\,\sin\,2πt \\
  x^2+y^2=L^2 \\
  x(0)=3,\,y(0)=-4,\\
  \dot{x}(0) = -\frac{v}{\sqrt{2}},\,\dot{y}(0) = -\frac{v}{\sqrt{2}} \\
  t\in\left[0;2\right],\quad\left(t\in\left[0;100\right]\right) \\
  v=1,\,L=5,\,m=1
\end{cases}
```

**Замечение**. Сравнивая систему с приведенной в условии, можно найти отличия в знаках у слогаемах в дифференциальных уравнениях. Почему так вышло ответить пока не готов, поэтому перепишем уравнения ещё разок с теми знаками, которые даны по условию:

```math
\begin{cases}
  m\ddot{x} = -\frac{x}{L}T \\
  m\ddot{y} = -\frac{y}{L}T - mg \\
  T=mg \\
  g=9,81+0,05\,\sin\,2πt \\
  x^2+y^2=L^2 \\
  x(0)=3,\,y(0)=-4,\\
  \dot{x}(0) = -\frac{v}{\sqrt{2}},\,\dot{y}(0) = -\frac{v}{\sqrt{2}} \\
  t\in\left[0;2\right],\quad\left(t\in\left[0;100\right]\right) \\
  v=1,\,L=5,\,m=1
\end{cases}
```
Перейдём непосредственно к решению системы.

## Решение методом Эйлера
Решение системы уравнений реализует программа в файле `code_Euler_method.py`, реализующая метод Эйлера, написанная на языке программирования Python. Запуская программу получаем следующие графики:

Решая на интервале $\left(0;2\right)$ получим:
<p align="center"><img src="https://github.com/sergefriendly/orchidea/raw/main/charts_1.png" alt="Чертёж"></p>

Решая на интервале $\left(0;10\right)$ получим:
<p align="center"><img src="https://github.com/sergefriendly/orchidea/raw/main/charts_2.png" alt="Чертёж"></p>
