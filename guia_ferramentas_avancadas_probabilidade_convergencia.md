# Guia de Estudos — Ferramentas Avançadas de Probabilidade e Convergências
### CE309 – Inferência Estatística | Lista 3

> Este guia não resolve a lista. Ele te entrega o mapa, a bússola e os avisos de trecho perigoso. A caminhada é sua.

---

## 1. Conceitos Fundamentais

### 1.1 Duas famílias de ferramentas

A lista inteira gira em torno de duas perguntas diferentes. Saber qual delas está sendo feita já resolve metade do problema.

**Pergunta A — "Quanto no máximo pode valer essa probabilidade?"**
Você conhece pouco sobre a distribuição (só a média, ou média e variância) e quer um **limite** (*bound*) válido para qualquer distribuição. Ferramentas: **Desigualdade de Markov** e **Desigualdade de Chebyshev**. O resultado é uma cota, não um valor exato, e vale para $n$ fixo (não é assintótico).

**Pergunta B — "Para onde essa estatística vai quando $n \to \infty$?"**
Você quer o comportamento limite de uma sequência de variáveis aleatórias. Ferramentas: **Lei dos Grandes Números**, **Teorema Central do Limite**, **Teorema da Aplicação Contínua** e **Teorema de Slutsky**. O resultado é uma aproximação que melhora com $n$.

### 1.2 Modos de convergência (o vocabulário mínimo)

| Modo | Notação | Ideia intuitiva | Onde aparece na lista |
|---|---|---|---|
| Em probabilidade | $X_n \xrightarrow{P} X$ | A chance de $X_n$ estar longe de $X$ some | Consistência de estimadores |
| Em distribuição | $X_n \xrightarrow{d} X$ | A *forma* da distribuição de $X_n$ se estabiliza | TCL, aproximações normais |
| Quase certa | $X_n \xrightarrow{q.c.} X$ | A trajetória converge (mais forte) | Lei Forte dos Grandes Números |

Hierarquia que você precisa ter no reflexo:

$$\xrightarrow{q.c.} \;\Longrightarrow\; \xrightarrow{P} \;\Longrightarrow\; \xrightarrow{d}$$

As setas **não voltam**, com uma exceção útil: se $X_n \xrightarrow{d} c$ com $c$ constante, então $X_n \xrightarrow{P} c$. Essa exceção é o que faz o Teorema de Slutsky funcionar na prática.

### 1.3 O que é "consistência"

Um estimador $\hat{\theta}_n$ é **consistente** para $\theta$ quando $\hat{\theta}_n \xrightarrow{P} \theta$. Na prática, quase sempre você prova consistência de um destes três jeitos:

1. **Direto pela LGN**: o estimador é uma média de v.a.'s iid com esperança igual ao alvo.
2. **Via Chebyshev**: mostra que o viés vai a zero e a variância vai a zero.
3. **Via Aplicação Contínua**: o estimador é uma função contínua de algo que você já sabe que converge.

### 1.4 Distribuição exata vs. distribuição aproximada

Vários itens da lista pedem para "comparar a distribuição aproximada com a exata". Guarde a distinção:

- **Exata**: vale para qualquer $n$, mas normalmente exige suposição forte (por exemplo, normalidade dos dados) ou conta pesada (soma de binomiais, convolução de Poisson).
- **Aproximada (assintótica)**: não exige a forma da distribuição, só momentos finitos, mas só é boa para $n$ grande.

O objetivo pedagógico dos itens computacionais é te fazer ver **a partir de qual $n$** a aproximação fica aceitável e **onde ela erra mais** (spoiler antecipado sem resolver nada: quase sempre nas caudas e quando há assimetria).

---

## 2. Fórmulas e Propriedades

### 2.1 Desigualdade de Markov

Para $X \geq 0$ e $a > 0$:

$$P(X \geq a) \leq \frac{E(X)}{a}$$

Requisitos: **variável não negativa** e esperança finita. Só usa a média. É a ferramenta mais fraca e mais barata.

Versão útil com função crescente não negativa $g$:

$$P(X \geq a) \leq \frac{E[g(X)]}{g(a)}$$

### 2.2 Desigualdade de Chebyshev

Para $X$ com média $\mu$ e variância $\sigma^2 < \infty$, e $k > 0$:

$$P(|X - \mu| \geq k) \leq \frac{\sigma^2}{k^2}$$

Forma equivalente em "número de desvios padrão", com $k = t\sigma$:

$$P(|X - \mu| \geq t\sigma) \leq \frac{1}{t^2}$$

Forma complementar (a que você usa para dimensionar amostra):

$$P(|X - \mu| < k) \geq 1 - \frac{\sigma^2}{k^2}$$

**Aplicada à média amostral** (o caso que mais aparece na lista): como $E(\bar{Y}_n) = \mu$ e $V(\bar{Y}_n) = \sigma^2/n$,

$$P(|\bar{Y}_n - \mu| \geq \varepsilon) \leq \frac{\sigma^2}{n\varepsilon^2}$$

### 2.3 Lei Fraca dos Grandes Números (LFGN)

Se $Y_1, \dots, Y_n$ são iid com $E(Y_i) = \mu$ finita, então

$$\bar{Y}_n \xrightarrow{P} \mu$$

Extensão que vale ouro na lista: se $h(\cdot)$ é uma função tal que $E[h(Y_i)]$ existe, então

$$\frac{1}{n}\sum_{i=1}^{n} h(Y_i) \xrightarrow{P} E[h(Y_1)]$$

É isso que permite tratar médias de quadrados, de desvios ao quadrado, etc.

### 2.4 Teorema Central do Limite (Lindeberg–Lévy)

Se $Y_1, \dots, Y_n$ são iid com média $\mu$ e variância $0 < \sigma^2 < \infty$:

$$\sqrt{n}\,\frac{\bar{Y}_n - \mu}{\sigma} \xrightarrow{d} N(0,1)$$

Formas equivalentes que você vai usar sem pensar:

$$\bar{Y}_n \stackrel{\cdot}{\sim} N\!\left(\mu, \frac{\sigma^2}{n}\right), \qquad \sum_{i=1}^{n} Y_i \stackrel{\cdot}{\sim} N\!\left(n\mu,\; n\sigma^2\right)$$

O símbolo $\stackrel{\cdot}{\sim}$ significa "aproximadamente distribuído como".

### 2.5 Correção de continuidade

Ao aproximar uma variável **discreta** $X$ por uma normal:

$$P(X \geq k) \approx P\!\left(Z \geq \frac{k - 0.5 - \mu}{\sigma}\right), \qquad P(X \leq k) \approx P\!\left(Z \leq \frac{k + 0.5 - \mu}{\sigma}\right)$$

Regra mnemônica: você "estica" o intervalo meia unidade **para incluir** o valor inteiro que está na fronteira.

### 2.6 Teorema da Aplicação Contínua (Continuous Mapping)

Se $g$ é contínua no ponto/limite relevante:

$$X_n \xrightarrow{P} X \;\Longrightarrow\; g(X_n) \xrightarrow{P} g(X)$$
$$X_n \xrightarrow{d} X \;\Longrightarrow\; g(X_n) \xrightarrow{d} g(X)$$

Exemplos de $g$ que aparecem: $g(x)=\sqrt{x}$, $g(x)=1/x$, $g(x)=x^2$, $g(x)=\log x$. Sempre verifique a continuidade **no ponto limite** (por exemplo, $\sqrt{\cdot}$ é contínua em $\sigma^2 > 0$).

### 2.7 Teorema de Slutsky

Se $X_n \xrightarrow{d} X$ e $A_n \xrightarrow{P} a$ (constante), então:

$$X_n + A_n \xrightarrow{d} X + a, \qquad A_n X_n \xrightarrow{d} aX, \qquad \frac{X_n}{A_n} \xrightarrow{d} \frac{X}{a} \;\;(a \neq 0)$$

Este é **o teorema** que resolve a maior parte dos itens teóricos difíceis da lista. Toda vez que você tiver uma razão em que o numerador tem TCL e o denominador é uma quantidade estimada, o roteiro é: separe, mostre que o denominador converge em probabilidade para uma constante, aplique Slutsky.

### 2.8 Método Delta (o "TCL para funções")

Se $\sqrt{n}(T_n - \theta) \xrightarrow{d} N(0, \tau^2)$ e $g$ é diferenciável com $g'(\theta) \neq 0$:

$$\sqrt{n}\left(g(T_n) - g(\theta)\right) \xrightarrow{d} N\!\left(0,\; [g'(\theta)]^2 \tau^2\right)$$

Use quando o exercício pedir a distribuição aproximada de algo como $\sqrt{S^2}$, $1/\bar{Y}$, $\log \bar{Y}$.

### 2.9 Momentos das distribuições que aparecem na lista

| Distribuição | Parametrização | $E(X)$ | $V(X)$ | Propriedade útil |
|---|---|---|---|---|
| Bernoulli($p$) | — | $p$ | $p(1-p)$ | soma de $n$ = Binomial |
| Binomial($n,p$) | — | $np$ | $np(1-p)$ | reprodutiva em $n$ |
| Poisson($\lambda$) | — | $\lambda$ | $\lambda$ | soma de Poisson é Poisson; Poisson($n$) = soma de $n$ Poisson(1) iid |
| Exponencial | taxa $\lambda$ | $1/\lambda$ | $1/\lambda^2$ | soma de $n$ = Gama($n,\lambda$) |
| Exponencial | média $\beta$ | $\beta$ | $\beta^2$ | cuidado com a parametrização |
| $N(\mu,\sigma^2)$ | — | $\mu$ | $\sigma^2$ | $\mu_4 = 3\sigma^4$ |

O fato de que **Poisson($n$) pode ser escrita como soma de $n$ Poisson(1) independentes** é a chave para aplicar o TCL quando o parâmetro cresce em vez do tamanho da amostra.

### 2.10 Momentos centrais e variância de quadrados

Defina $\mu_k = E[(Y - \mu)^k]$. Então:

$$V\!\left[(Y-\mu)^2\right] = E\!\left[(Y-\mu)^4\right] - \left(E\left[(Y-\mu)^2\right]\right)^2 = \mu_4 - \sigma^4$$

Sob normalidade, $\mu_4 = 3\sigma^4$, logo essa variância vira $2\sigma^4$. Guarde essa fórmula: ela é o coração do exercício de distribuição assintótica da variância.

### 2.11 Decomposição da soma de quadrados

Identidade algébrica que aparece em toda demonstração de consistência de $S^2$:

$$\sum_{i=1}^{n}(Y_i - \bar{Y})^2 = \sum_{i=1}^{n}(Y_i - \mu)^2 - n(\bar{Y} - \mu)^2$$

Dividindo por $n$:

$$\frac{1}{n}\sum_{i=1}^{n}(Y_i - \bar{Y})^2 = \underbrace{\frac{1}{n}\sum_{i=1}^{n}(Y_i - \mu)^2}_{\xrightarrow{P}\;\sigma^2 \text{ pela LGN}} - \underbrace{(\bar{Y} - \mu)^2}_{\xrightarrow{P}\;0}$$

E a ponte entre as duas versões da variância amostral:

$$S^2 = \frac{n}{n-1} \cdot \frac{1}{n}\sum_{i=1}^{n}(Y_i-\bar{Y})^2, \qquad \frac{n}{n-1} \to 1$$

### 2.12 Resultados exatos sob normalidade

Se $Y_1,\dots,Y_n \sim N(\mu,\sigma^2)$ iid:

$$\frac{(n-1)S^2}{\sigma^2} \sim \chi^2_{n-1}, \qquad \frac{1}{\sigma^2}\sum_{i=1}^{n}(Y_i-\mu)^2 \sim \chi^2_{n}$$

$$T = \frac{\bar{Y}-\mu}{S/\sqrt{n}} \sim t_{n-1}$$

Lembretes: $E(\chi^2_k) = k$, $V(\chi^2_k) = 2k$, e $t_k \to N(0,1)$ quando $k \to \infty$.

---

## 3. Passo a Passo Lógico (Framework de Resolução)

### Framework A — Problemas de limite superior para probabilidade

**Gatilho:** o enunciado diz "obtenha um limite superior", "no mínimo/no máximo", e te dá **apenas momentos**, sem nomear a distribuição.

1. **Inventário do que você tem.** Só a média? Markov é o caminho. Média **e** variância? Chebyshev.
2. **Verifique a não negatividade** se for usar Markov. Receita, tempo, contagem, quantidade: tudo isso é $\geq 0$. Se a variável puder ser negativa, Markov não se aplica diretamente.
3. **Reescreva o evento no formato canônico da desigualdade.**
   - Markov quer $\{X \geq a\}$.
   - Chebyshev quer $\{|X - \mu| \geq k\}$. Eventos do tipo "pelo menos $A$ ou no máximo $B$" com $A$ e $B$ simétricos em torno de $\mu$ são exatamente isso: identifique o centro $\mu = (A+B)/2$ e a distância $k = (A - B)/2$.
4. **Identifique quem é a variável.** Cuidado: às vezes o evento é sobre a **soma** ou sobre a **contagem**, não sobre a média. Os momentos mudam conforme o objeto.
5. **Aplique a fórmula e interprete.** Se o limite der maior que 1, ele é verdadeiro mas inútil. Se der um número pequeno, comente que é conservador.

**Exemplo de raciocínio (com números fictícios, fora da lista):** se $X \geq 0$ tem média 50 e você quer limitar $P(X \geq 200)$, o formato já é o de Markov, com $a = 200$: o limite é $50/200$. Note como não foi preciso saber nada sobre a forma da distribuição.

### Framework B — Dimensionamento de amostra via Chebyshev

**Gatilho:** "quantas observações são necessárias para que a probabilidade de ... seja de pelo menos $1-\alpha$".

1. Escreva o evento desejado como $\{|\bar{Y}_n - \mu| < \varepsilon\}$ e identifique $\varepsilon$ no enunciado.
2. Use a forma complementar: $P(|\bar{Y}_n - \mu| < \varepsilon) \geq 1 - \dfrac{\sigma^2}{n\varepsilon^2}$.
3. Imponha a exigência sobre o **lado garantido** da desigualdade: $1 - \dfrac{\sigma^2}{n\varepsilon^2} \geq 1-\alpha$.
4. Isole $n$. Vai sobrar algo do tipo $n \geq \dfrac{\sigma^2}{\alpha \varepsilon^2}$.
5. **Arredonde para cima** (é tamanho de amostra, tem que ser inteiro e satisfazer a desigualdade).
6. Opcional e muito valorizado: compare com o $n$ que o TCL exigiria. O de Chebyshev sempre será bem maior, porque não usa informação sobre a forma da distribuição.

### Framework C — Aproximação normal para variáveis discretas

**Gatilho:** binomial, Poisson ou contagem com "qual a probabilidade de pelo menos $k$".

1. **Calcule a probabilidade exata** quando pedido: identifique a distribuição e some as parcelas relevantes, ou use a cauda complementar $P(X \geq k) = 1 - P(X \leq k-1)$.
2. **Monte a aproximação:** obtenha $\mu$ e $\sigma^2$ da variável, padronize.
3. **Aplique a correção de continuidade.** Para $P(X \geq k)$, use $k - 0.5$.
4. **Consulte a tabela da normal** com o $z$ obtido.
5. **Compare os dois números** e comente a qualidade da aproximação. Ela é melhor quando $np$ e $n(1-p)$ são ambos razoavelmente grandes (regra de bolso: $\geq 5$ ou $\geq 10$) e quando $p$ está perto de $0.5$.

### Framework D — "Qual a distribuição aproximada de $\bar{Y}_n$?"

1. Identifique a distribuição de cada $Y_i$ e **anote a parametrização** com cuidado.
2. Calcule $E(Y_i)$ e $V(Y_i)$ pela tabela de momentos.
3. Escreva direto: $\bar{Y}_n \stackrel{\cdot}{\sim} N\!\left(E(Y_i),\, V(Y_i)/n\right)$.
4. Se o enunciado pedir "mostre", apresente a forma padronizada com $\sqrt{n}$ e cite o TCL nominalmente, verificando as hipóteses (iid, variância finita).

### Framework E — Provar convergência em probabilidade

Escolha uma das três rotas conforme a cara do estimador.

**Rota 1 — LGN direta.** O estimador é $\frac{1}{n}\sum h(Y_i)$? Então basta mostrar que $E[h(Y_1)]$ existe e vale o alvo. Fim.

**Rota 2 — Chebyshev.** Escreva $P(|\hat{\theta}_n - \theta| \geq \varepsilon) \leq \dfrac{E[(\hat{\theta}_n-\theta)^2]}{\varepsilon^2} = \dfrac{V(\hat{\theta}_n) + \text{viés}^2}{\varepsilon^2}$ e mostre que o numerador vai a zero. Esta rota dá um argumento autocontido, sem invocar a LGN.

**Rota 3 — Decompor + Aplicação Contínua.** O estimador é combinação de pedaços que convergem? Então:
   - quebre a expressão em pedaços com limites conhecidos;
   - use que soma, produto e razão de convergências em probabilidade convergem para a soma, produto e razão dos limites;
   - se houver uma função contínua envolvendo tudo, aplique o Teorema da Aplicação Contínua no final.

Esta é a rota para qualquer coisa do tipo "para qual valor $\sqrt{\hat{\theta}_n}$ converge?".

### Framework F — Razões com denominador aleatório (o padrão Slutsky)

**Gatilho:** a estatística é uma fração e o denominador contém uma quantidade estimada ou aleatória. Este é o padrão dos exercícios sobre a estatística $t$ e sobre a razão envolvendo Poisson.

O roteiro tem sempre a mesma cara. Dada uma estatística $W_n = \dfrac{N_n}{D_n}$:

1. **Multiplique e divida por uma constante conveniente** para transformar $W_n$ num produto de duas peças:

$$W_n = \underbrace{\frac{N_n}{c_n}}_{\text{peça com TCL}} \times \underbrace{\frac{c_n}{D_n}}_{\text{peça que vai para 1}}$$

onde $c_n$ é o desvio padrão teórico do numerador.

2. **Peça 1:** mostre que $N_n/c_n \xrightarrow{d} N(0,1)$ pelo TCL. Verifique média zero e variância unitária.
3. **Peça 2:** mostre que $D_n/c_n \xrightarrow{P} 1$, tipicamente via LGN mais Aplicação Contínua (a raiz quadrada é contínua e positiva no limite).
4. **Conclua por Slutsky:** produto de convergência em distribuição por convergência em probabilidade a constante.
5. **Diga explicitamente qual é a conclusão para $n$ finito:** a estatística é *aproximadamente* $N(0,1)$, não exatamente.

### Framework G — Ilustração computacional (Monte Carlo)

Todo item que pede "faça uma ilustração computacional" segue o mesmo esqueleto de cinco passos:

1. **Fixe** a distribuição geradora, os parâmetros verdadeiros, o tamanho amostral $n$ e o número de réplicas $M$ (use $M \geq 10.000$ para curvas suaves).
2. **Replique:** para cada uma das $M$ réplicas, gere uma amostra de tamanho $n$ e calcule a estatística de interesse.
3. **Padronize** a estatística exatamente como o teorema manda (subtraia a média teórica, divida pelo desvio padrão teórico), senão a comparação com a $N(0,1)$ não faz sentido.
4. **Compare** de três maneiras complementares:
   - histograma ou densidade empírica com a curva teórica sobreposta;
   - **QQ-plot** contra a distribuição teórica (é aqui que as caudas denunciam o desajuste);
   - função de distribuição empírica contra a teórica, ou um teste de Kolmogorov–Smirnov como métrica resumo.
5. **Varie $n$** (por exemplo 50, 250, 1000) e monte um painel. A narrativa esperada é "a aproximação melhora e o desajuste se concentra nas caudas".

Esqueleto em R:

```r
set.seed(2024)
M <- 10000
ns <- c(50, 250, 1000)

simula <- function(n) {
  replicate(M, {
    y <- rDISTRIBUICAO(n, parametros)   # troque pela geradora do exercício
    estatistica <- ...                  # calcule a estatística
    (estatistica - media_teorica) / desvio_teorico
  })
}

par(mfrow = c(2, length(ns)))
for (n in ns) {
  z <- simula(n)
  hist(z, freq = FALSE, breaks = 40, main = paste("n =", n))
  curve(dnorm(x), add = TRUE, lwd = 2)
}
for (n in ns) {
  z <- simula(n)
  qqnorm(z); qqline(z)
}
```

Esqueleto equivalente em Python:

```python
import numpy as np
import matplotlib.pyplot as plt
from scipy import stats

rng = np.random.default_rng(2024)
M, ns = 10_000, [50, 250, 1000]

def simula(n):
    amostras = rng.DISTRIBUICAO(parametros, size=(M, n))   # troque pela geradora
    est = ...                                              # estatística por linha (axis=1)
    return (est - media_teorica) / desvio_teorico

fig, ax = plt.subplots(2, len(ns), figsize=(14, 7))
for j, n in enumerate(ns):
    z = simula(n)
    ax[0, j].hist(z, bins=40, density=True)
    grid = np.linspace(-4, 4, 200)
    ax[0, j].plot(grid, stats.norm.pdf(grid))
    stats.probplot(z, dist="norm", plot=ax[1, j])
```

Para comparar com a distribuição **exata** (e não com a normal), troque a curva sobreposta pela densidade exata correspondente: qui-quadrado no caso das somas de quadrados, $t_{n-1}$ no caso da estatística $t$, e a massa de probabilidade discreta no caso de contagens.

---

## 4. Dicas e Macetes de Memorização

**"Markov é pobre, Chebyshev é remediado."** Markov só tem a média e por isso só consegue limitar uma cauda de variável positiva. Chebyshev tem média e variância e por isso consegue limitar as duas caudas ao mesmo tempo.

**Chebyshev é Markov disfarçado.** Não decore as duas separadamente. Aplique Markov à variável não negativa $(X-\mu)^2$ com o ponto de corte $k^2$:

$$P(|X-\mu| \geq k) = P\!\left((X-\mu)^2 \geq k^2\right) \leq \frac{E[(X-\mu)^2]}{k^2} = \frac{\sigma^2}{k^2}$$

Se você souber reproduzir essas três linhas, nunca mais erra a fórmula, e ainda ganha pontos em prova quando pedirem a demonstração.

**A regra "1 sobre t ao quadrado".** Chebyshev garante que a chance de estar a mais de 2 desvios da média é no máximo $1/4$, e a mais de 3 desvios é no máximo $1/9$. Compare mentalmente com a normal (cerca de 5% e 0,3%) e você lembra na hora o quanto o limite é frouxo.

**Mnemônico da divisão de trabalho:**
- LGN responde **para onde vai** (o centro).
- TCL responde **como flutua em torno do centro** (a escala e a forma).
- Slutsky responde **o que acontece quando você substitui um parâmetro por um estimador**.
- Aplicação Contínua responde **o que acontece quando você transforma o resultado**.

**O $\sqrt{n}$ nunca some.** Toda vez que escrever um TCL, cheque a "conta de escala": $\bar{Y}$ tem desvio padrão $\sigma/\sqrt{n}$, então multiplicar $(\bar{Y}-\mu)$ por $\sqrt{n}/\sigma$ produz algo com variância 1. Se a sua expressão padronizada não tiver variância 1, tem erro nela.

**Correção de continuidade — a régua mental.** Desenhe a barra do histograma do inteiro $k$: ela ocupa de $k-0.5$ a $k+0.5$. Se o evento **inclui** $k$, a área tem que começar em $k-0.5$. Se o evento **exclui** $k$, começa em $k+0.5$.

**Exponencial: escreva a parametrização antes de qualquer conta.** Duas convenções coexistem, taxa e média. Um erro aqui contamina toda a questão. Deixe explícito na primeira linha: "adotando $\lambda$ como taxa, $E(Y)=1/\lambda$".

**Poisson tem média igual à variância.** Isso simplifica muito as padronizações: o desvio padrão é $\sqrt{\lambda}$.

**"Parâmetro grande = amostra grande".** Quando o TCL parece não se aplicar porque não há amostra, procure a representação de soma. Poisson($n$) é soma de $n$ Poisson(1) iid; Binomial($n,p$) é soma de $n$ Bernoulli($p$) iid; Gama($n,\lambda$) é soma de $n$ Exponenciais($\lambda$) iid. O TCL entra por essa porta.

**A variância de um quadrado precisa do quarto momento.** Sempre que a estatística envolver $(Y_i-\mu)^2$, a variância da parcela depende de $\mu_4$. Anote a fórmula $\mu_4 - \sigma^4$ num cartão. Sob normalidade ela colapsa em $2\sigma^4$, o que é consistente com $V(\chi^2_n) = 2n$.

**Truque de sanidade para resultados assintóticos de variância.** Se você chegou a $\sqrt{n}(\hat\sigma^2 - \sigma^2) \xrightarrow{d} N(0, 2\sigma^4)$ no caso normal, confira contra o resultado exato: $\hat\sigma^2 \sim \frac{\sigma^2}{n}\chi^2_n$ tem variância $\frac{\sigma^4}{n^2} \cdot 2n = \frac{2\sigma^4}{n}$. Bate.

**Fórmula-mestra para o padrão Slutsky.** Decore este encaixe, que serve para a estatística $t$ e para razões com denominador aleatório:

$$\frac{\text{alvo}}{\text{estimado}} = \frac{\text{alvo}}{\text{verdadeiro}} \times \frac{\text{verdadeiro}}{\text{estimado}}$$

A primeira fração tem TCL, a segunda vai para 1.

---

## 5. Zona de Perigo (Armadilhas Comuns)

### 5.1 Nas desigualdades

**Usar Markov em variável que pode ser negativa.** A desigualdade exige $X \geq 0$. Se o enunciado não garante isso, você precisa de outra ferramenta ou de uma transformação não negativa.

**Confundir limite com valor exato.** Markov e Chebyshev entregam uma cota. Escrever "a probabilidade é $0{,}25$" quando o correto é "a probabilidade é no máximo $0{,}25$" custa a questão inteira.

**Esquecer que a variância da média é $\sigma^2/n$.** O erro mais comum da lista inteira. Se o evento é sobre $\bar{Y}_n$, use $\sigma^2/n$. Se é sobre a soma, use $n\sigma^2$. Se é sobre uma única observação, use $\sigma^2$.

**Não centrar corretamente o evento de duas caudas.** "Pelo menos $A$ ou no máximo $B$" só vira $|X-\mu| \geq k$ se $A$ e $B$ forem simétricos em torno da média. Confira: $\mu - B$ deve ser igual a $A - \mu$. Se não forem simétricos, você precisa limitar cada cauda separadamente ou usar a cauda mais distante (o que gera um limite mais frouxo, e isso deve ser dito).

**Reportar um limite maior que 1 sem comentar.** Matematicamente correto, praticamente inútil. Reconheça isso no texto.

**Trocar o sentido da desigualdade ao passar para o complementar.** $P(|X-\mu| < k) \geq 1 - \sigma^2/k^2$. O "$\leq$" vira "$\geq$". Escreva com calma.

### 5.2 No TCL e nas aproximações

**Esquecer a correção de continuidade em variáveis discretas.** Para $n$ moderado o erro é visível. Se o exercício pede para comparar exato e aproximado, a correção normalmente é o que salva a comparação.

**Errar o lado da correção.** $P(X \geq 9)$ vira $P(Z \geq (8.5 - \mu)/\sigma)$, não $9.5$. Pense na barra do histograma.

**Aplicar o TCL sem checar variância finita.** A hipótese $\sigma^2 < \infty$ não é decorativa. Distribuições de cauda pesada quebram o teorema.

**Confundir "$\bar{Y}_n$ é aproximadamente normal" com "$\bar{Y}_n \to$ alguma normal".** Rigorosamente, é a versão padronizada por $\sqrt{n}$ que converge; $\bar{Y}_n$ sozinha converge em probabilidade para a constante $\mu$, e a variância da aproximação encolhe com $n$. Escreva sempre a forma padronizada quando for provar algo.

**Trocar taxa por média na exponencial.** Já mencionado nos macetes, mas repito porque é o erro mais caro em questões desse tipo: um $\lambda$ mal interpretado inverte média e variância de lugar.

**Achar que a aproximação normal é igualmente boa em todo lugar.** Ela é pior nas caudas e pior sob assimetria forte. Nas ilustrações computacionais, é exatamente isso que o QQ-plot vai mostrar, e comentar esse ponto é o que diferencia uma resposta boa de uma mediana.

### 5.3 Em convergências e demonstrações

**Confundir os dois estimadores de variância.** Note bem quem é quem:
- $\hat\sigma^2 = \frac{1}{n}\sum (Y_i - \mu)^2$ usa a média **verdadeira** e é média de iid, então a LGN se aplica de imediato.
- $S^2 = \frac{1}{n-1}\sum (Y_i - \bar{Y})^2$ usa a média **amostral**, e as parcelas $(Y_i-\bar{Y})^2$ **não são independentes**. Aplicar a LGN diretamente aqui é um erro grave. É preciso passar pela identidade da soma de quadrados primeiro.

**Esquecer o fator $n/(n-1)$.** Ele tende a 1 e por isso não muda o limite, mas precisa aparecer no argumento. Omiti-lo sugere que você não percebeu a diferença entre as duas versões.

**Aplicar Aplicação Contínua sem verificar continuidade no ponto limite.** Para $S = \sqrt{S^2}$, a função raiz é contínua em $\sigma^2 > 0$. Se $\sigma^2 = 0$ fosse permitido, o argumento precisaria de cuidado extra. Mencione a condição.

**Achar que $E(S) = \sigma$.** Convergência em probabilidade não implica não viesamento. $S^2$ é não viesado para $\sigma^2$, mas $S$ é viesado para $\sigma$ em amostras finitas, mesmo sendo consistente. Não misture as duas propriedades.

**Usar Slutsky com um denominador que converge para uma variável aleatória.** O teorema exige que a peça multiplicativa convirja em probabilidade para uma **constante**. Se ela converge para algo aleatório, a conclusão não vale.

**Substituir $\sigma$ por $S$ sem justificar.** Escrever "como $S \approx \sigma$, então $t \approx Z$" é hand-waving. A justificativa formal é a decomposição em duas peças mais Slutsky. É pouca linha de escrita e vale muito ponto.

**Confundir o resultado assintótico com o exato na estatística $t$.** Para dados normais, a distribuição é exatamente $t_{n-1}$ para todo $n$. O resultado assintótico $N(0,1)$ vale sem a suposição de normalidade, mas só quando $n$ é grande. São afirmações diferentes com hipóteses diferentes, e o exercício que pede a comparação computacional está justamente testando se você entendeu isso.

**Padronizar errado nas simulações.** Se você simular a estatística e comparar com a $N(0,1)$ sem subtrair a média teórica e dividir pelo desvio teórico corretos, o histograma vai parecer deslocado ou com escala errada, e você vai concluir que o teorema falhou. Antes de duvidar do teorema, confira a padronização.

**Somar variâncias de variáveis dependentes.** A regra $V(A - B) = V(A) + V(B)$ só vale sob independência. Em qualquer questão com duas contagens independentes, diga que está usando a independência.

**Deixar o número de réplicas pequeno.** Com $M = 500$ o histograma fica ruidoso e a comparação vira chute. Suba para dezenas de milhares e fixe a semente para reprodutibilidade.

---

## Checklist final antes de entregar

- [ ] Toda desigualdade está declarada como **limite**, não como valor exato?
- [ ] Verifiquei se a variância que usei é a da variável certa (observação, soma ou média)?
- [ ] Anotei explicitamente a parametrização de cada distribuição?
- [ ] Em todo TCL, a expressão padronizada tem média 0 e variância 1?
- [ ] Usei correção de continuidade sempre que a variável era discreta?
- [ ] Em cada demonstração de convergência, citei nominalmente o teorema usado e verifiquei suas hipóteses?
- [ ] Nas simulações: semente fixada, $M$ grande, padronização correta, e comparação por histograma **e** QQ-plot?
- [ ] Comentei o comportamento conforme $n$ cresce, em vez de só colar os gráficos?

Você já tem todas as ferramentas listadas aqui. A lista é essencialmente uma sequência de aplicações de cinco resultados (Markov, Chebyshev, LGN, TCL e Slutsky) em roupagens diferentes. Identifique qual dos cinco cada questão está pedindo e o resto é álgebra cuidadosa. Bom estudo.
