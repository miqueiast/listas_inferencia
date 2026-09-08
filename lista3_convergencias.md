# CE309 — Inferência Estatística
## Lista 3 — Ferramentas Avançadas de Probabilidade e Convergências
### Resolução comentada e aula passo a passo

*Profa. Amanda M. F. Mendes — Departamento de Estatística*

---

## Sumário

- [Parte 0 — Caixa de ferramentas teórica](#parte-0--caixa-de-ferramentas-teórica)
- [Exercício 1 — Markov](#exercício-1--markov)
- [Exercício 2 — Dimensionamento amostral](#exercício-2--dimensionamento-amostral)
- [Exercício 3 — Chebyshev na binomial](#exercício-3--chebyshev-na-binomial)
- [Exercício 4 — Binomial exata vs. TCL](#exercício-4--binomial-exata-vs-tcl)
- [Exercício 5 — Média amostral de exponenciais](#exercício-5--média-amostral-de-exponenciais)
- [Exercício 6 — Poisson: LGN e TCL](#exercício-6--poisson-lgn-e-tcl)
- [Exercício 7 — Estimador de variância com μ conhecido](#exercício-7--estimador-de-variância-com-μ-conhecido)
- [Exercício 8 — Consistência de S² e de S](#exercício-8--consistência-de-s²-e-de-s)
- [Exercício 9 — A estatística t e o Teorema de Slutsky](#exercício-9--a-estatística-t-e-o-teorema-de-slutsky)
- [Exercício 10 — Diferença de Poissons studentizada](#exercício-10--diferença-de-poissons-studentizada)
- [Apêndice A — Mapa de decisão](#apêndice-a--mapa-de-decisão-qual-ferramenta-usar)
- [Apêndice B — Erros clássicos de prova](#apêndice-b--erros-clássicos-de-prova)

---

# Parte 0 — Caixa de ferramentas teórica

Antes de resolver qualquer exercício, é preciso ter clareza sobre **quais são as ferramentas disponíveis e o que cada uma exige**. Toda esta lista se resolve com os sete resultados abaixo.

## 0.1 Desigualdade de Markov

Seja $X$ uma variável aleatória **não-negativa** com $E(X)<\infty$. Então, para todo $a>0$:

$$P(X\geq a)\;\leq\;\frac{E(X)}{a}$$

**O que ela faz:** limita a probabilidade de uma variável não-negativa assumir valores grandes, usando **apenas a média**.

**Hipóteses:** $X\geq 0$ e média finita. Nada mais. Não precisa de variância, não precisa da distribuição.

**Intuição:** se a média de uma quantidade positiva é 3, ela não pode passar de 12 com frequência alta — porque valores grandes puxariam a média para cima. A desigualdade quantifica esse raciocínio.

**Esboço da prova** (vale escrever em prova se pedido):

$$E(X)=\int_0^\infty x\,f(x)\,dx\;\geq\;\int_a^\infty x\,f(x)\,dx\;\geq\;\int_a^\infty a\,f(x)\,dx = a\,P(X\geq a)$$

O primeiro corte descarta a região $[0,a)$, que contribui com algo não-negativo; o segundo substitui $x$ pelo seu menor valor possível naquela região, que é $a$.

## 0.2 Desigualdade de Chebyshev

Seja $X$ com média $\mu=E(X)$ e variância $\sigma^2<\infty$. Então, para todo $\varepsilon>0$:

$$P\left(|X-\mu|\geq \varepsilon\right)\;\leq\;\frac{\sigma^2}{\varepsilon^2}$$

Ou, na versão em múltiplos do desvio padrão ($\varepsilon = t\sigma$):

$$P\left(|X-\mu|\geq t\sigma\right)\;\leq\;\frac{1}{t^2}$$

**O que ela faz:** limita a probabilidade de $X$ se afastar da sua média, usando média **e** variância.

**De onde vem:** é Markov aplicada a $Y=(X-\mu)^2$, que é não-negativa, com $a=\varepsilon^2$:

$$P\left(|X-\mu|\geq\varepsilon\right)=P\left((X-\mu)^2\geq\varepsilon^2\right)\leq\frac{E\left[(X-\mu)^2\right]}{\varepsilon^2}=\frac{\sigma^2}{\varepsilon^2}$$

Esse é o motivo de Chebyshev vir sempre depois de Markov no curso: ela é um corolário.

**Quando usar:** sempre que o enunciado der só média e variância, e pedir um "limite superior" ou "no mínimo X% de probabilidade". A palavra **limite** é a pista.

**Preço:** a cota é universal, portanto conservadora. Costuma ser várias vezes maior que a probabilidade real.

## 0.3 Os três tipos de convergência

Esta é a distinção conceitual mais importante da lista.

| Tipo | Notação | Definição | O limite é... |
|---|---|---|---|
| Em probabilidade | $X_n\xrightarrow{p}c$ | $\forall\varepsilon>0:\;P(\lvert X_n-c\rvert\geq\varepsilon)\to 0$ | uma **constante** |
| Em distribuição | $X_n\xrightarrow{d}X$ | $F_{X_n}(x)\to F_X(x)$ nos pontos de continuidade | uma **distribuição** |
| Quase certa | $X_n\xrightarrow{q.c.}c$ | $P(\lim X_n=c)=1$ | uma constante (mais forte) |

Relação entre elas:

$$\text{quase certa}\;\Longrightarrow\;\text{em probabilidade}\;\Longrightarrow\;\text{em distribuição}$$

As recíprocas são falsas.

**A leitura prática:**

- Convergência **em probabilidade** responde: *o estimador acerta o alvo?* → consistência.
- Convergência **em distribuição** responde: *qual é o formato da flutuação em torno do alvo?* → construção de intervalos de confiança e testes.

## 0.4 Lei Fraca dos Grandes Números (LGN)

Seja $Y_1,\dots,Y_n$ i.i.d. com $E(Y_i)=\mu$ finita. Então:

$$\bar{Y}_n\;\xrightarrow{p}\;\mu$$

**Hipótese mínima (versão de Khintchine):** apenas média finita. Não é preciso variância finita.

**Versão via Chebyshev** (a que se demonstra em prova): se além disso $\operatorname{Var}(Y_i)=\sigma^2<\infty$, então

$$P\left(|\bar{Y}_n-\mu|\geq\varepsilon\right)\leq\frac{\operatorname{Var}(\bar{Y}_n)}{\varepsilon^2}=\frac{\sigma^2}{n\varepsilon^2}\;\xrightarrow[n\to\infty]{}\;0$$

**O mecanismo é sempre o mesmo:** $\operatorname{Var}(\bar{Y}_n)=\sigma^2/n$ decai como $1/n$; Chebyshev converte essa variância que some em uma probabilidade que some.

**Aplicação essencial (usada o tempo todo nesta lista):** se $g$ é uma função qualquer e $E[g(Y_i)]$ existe, então $g(Y_1),\dots,g(Y_n)$ também são i.i.d., logo

$$\frac{1}{n}\sum_{i=1}^n g(Y_i)\;\xrightarrow{p}\;E[g(Y_1)]$$

Ou seja: **médias amostrais de qualquer transformação convergem para a esperança correspondente.** Metade da lista é aplicação disso com $g(y)=y^2$ ou $g(y)=(y-\mu)^2$.

## 0.5 Teorema Central do Limite (TCL)

Seja $Y_1,\dots,Y_n$ i.i.d. com média $\mu$ e variância $\sigma^2<\infty$ (e $\sigma^2>0$). Então:

$$\frac{\bar{Y}_n-\mu}{\sigma/\sqrt{n}}=\frac{\sqrt{n}\left(\bar{Y}_n-\mu\right)}{\sigma}\;\xrightarrow{d}\;N(0,1)$$

Notação prática: $\bar{Y}_n\stackrel{a}{\sim}N\left(\mu,\sigma^2/n\right)$.

**Hipótese:** variância finita. Isso é mais forte que a LGN exige.

**A sacada do fator $\sqrt{n}$:** pela LGN, $\bar{Y}_n-\mu\to 0$ — o limite é degenerado, não informa nada sobre o formato. Multiplicando por $\sqrt{n}$, a escala é ajustada **exatamente** na medida certa: nem colapsa em zero, nem explode. Sobra uma distribuição não-degenerada, e ela é sempre a normal.

| | LGN | TCL |
|---|---|---|
| Escala | $\bar{Y}_n-\mu$ | $\sqrt{n}(\bar{Y}_n-\mu)$ |
| Limite | $0$ (constante) | $N(0,\sigma^2)$ |
| Diz | o erro some | a que velocidade e com que forma |

**Quão grande deve ser $n$?** Depende da assimetria da população:

| População | $n$ confortável |
|---|---|
| Simétrica e leve (uniforme) | $n\approx 10$ |
| Moderadamente assimétrica | $n\approx 30$ |
| Exponencial (assimetria 2) | $n\approx 50$ |
| Binomial | $np\geq 5$ e $n(1-p)\geq 5$ |
| Poisson | $n\lambda\gtrsim 10$ |

A regra "$n\geq 30$" é referência, não lei.

## 0.6 Teorema da Aplicação Contínua (Continuous Mapping Theorem)

Se $g$ é contínua no ponto $c$, então:

$$X_n\xrightarrow{p}c \quad\Longrightarrow\quad g(X_n)\xrightarrow{p}g(c)$$

E, para convergência em distribuição: se $g$ é contínua, $X_n\xrightarrow{d}X \Rightarrow g(X_n)\xrightarrow{d}g(X)$.

**Por que é tão útil:** transforma um resultado de convergência em vários, sem nenhuma conta. Sabendo que $S^2\xrightarrow{p}\sigma^2$, obtém-se de graça $S\xrightarrow{p}\sigma$, $1/S^2\xrightarrow{p}1/\sigma^2$, $\log S^2\xrightarrow{p}\log\sigma^2$ — bastando checar continuidade no ponto.

**Cuidado clássico:** o teorema fala de **convergência**, não de **esperança**. É falso que $E[g(X_n)]\to g(E(X_n))$. Contra-exemplo: $E(\bar{Y}^2)=\mu^2+\sigma^2/n\neq\mu^2$, embora $\bar{Y}^2\xrightarrow{p}\mu^2$.

## 0.7 Teorema de Slutsky

Se $X_n\xrightarrow{d}X$ e $A_n\xrightarrow{p}a$ (constante), então:

$$X_n+A_n\;\xrightarrow{d}\;X+a \qquad\qquad A_n X_n\;\xrightarrow{d}\;aX \qquad\qquad \frac{X_n}{A_n}\;\xrightarrow{d}\;\frac{X}{a}\;\;(a\neq 0)$$

**O que ele resolve:** o problema central da inferência prática. O TCL entrega

$$\frac{\bar{Y}-\mu}{\sigma/\sqrt{n}}\xrightarrow{d}N(0,1)$$

mas $\sigma$ é **desconhecido**. Na prática usamos $S$ no lugar. Slutsky é o teorema que garante que essa substituição não estraga o limite normal — desde que $S\xrightarrow{p}\sigma$.

**Atenção à hipótese:** o termo que converge em probabilidade tem que convergir para uma **constante**. Se $A_n\xrightarrow{d}A$ aleatório, o teorema não se aplica.

**Nota sobre a versão "$X_n+A_n$":** ela também permite trocar uma quantidade por outra assintoticamente equivalente, escrevendo a diferença como um termo que vai a zero em probabilidade.

## 0.8 Correção de continuidade

Ao aproximar uma variável **discreta** $X$ (binomial, Poisson) por uma **contínua** (normal), trate cada inteiro $k$ como o intervalo $[k-0{,}5;\,k+0{,}5]$:

| Evento discreto | Substituir por |
|---|---|
| $X\geq k$ | $X > k-0{,}5$ |
| $X > k$ | $X > k+0{,}5$ |
| $X\leq k$ | $X < k+0{,}5$ |
| $X < k$ | $X < k-0{,}5$ |
| $X = k$ | $k-0{,}5 < X < k+0{,}5$ |

**Regra mnemônica:** desloque sempre no sentido que **aumenta** a região.

O efeito é grande quando $n$ é pequeno ou o corte está perto da média — exatamente onde a "barra" descartada é mais alta.

---

# Exercício 1 — Markov

> A receita diária de uma pequena empresa tem média de R\$ 3.000,00. Obtenha um limite superior para a probabilidade de a receita de um dia qualquer ser de pelo menos R\$ 12.000,00.

## 1.1 Diagnóstico

O que o enunciado **dá**: a média, e só. Não há variância, não há distribuição.

O que o enunciado **pede**: um limite superior para $P(X\geq 12000)$.

Com apenas a média disponível, existe uma única ferramenta: **Markov**. E ela exige que a variável seja não-negativa — o que é razoável aqui, já que receita não é negativa. Essa hipótese precisa ser **declarada explicitamente** na resposta; é ela que autoriza o uso do teorema.

## 1.2 Resolução

Seja $X$ a receita diária, com $X\geq 0$ e $E(X)=3000$.

Pela desigualdade de Markov, com $a=12000$:

$$P(X\geq 12000)\;\leq\;\frac{E(X)}{12000}=\frac{3000}{12000}=\frac{1}{4}=0{,}25$$

**Resposta:** a probabilidade é de no máximo $25\%$.

## 1.3 Comentários do professor

**Sobre a estrutura da resposta.** Note que $12000 = 4\times 3000$, ou seja, $a=4\mu$. Nesse formato Markov sempre devolve $1/4$:

$$P(X\geq k\mu)\leq\frac{\mu}{k\mu}=\frac{1}{k}$$

Vale memorizar: **a probabilidade de uma variável não-negativa exceder $k$ vezes sua média é no máximo $1/k$.**

**A cota é atingível?** Sim — e isso é o que a torna a melhor cota possível com essa informação. Considere a variável que vale $12000$ com probabilidade $1/4$ e $0$ com probabilidade $3/4$. Sua média é $0{,}25\times 12000=3000$, exatamente como no enunciado, e $P(X\geq 12000)=0{,}25$ exatamente. Ou seja, **existe** uma distribuição compatível com o enunciado que atinge a cota; não é possível melhorá-la sem informação adicional.

**Por que não Chebyshev?** Porque Chebyshev exige a variância, que o enunciado não fornece. Se fornecesse — digamos $\sigma^2 = 1.000.000$ — poderíamos escrever $12000 = 3000 + 9000$ e obter uma cota muito melhor:

$$P(X\geq 12000)\leq P(|X-3000|\geq 9000)\leq\frac{1000000}{9000^2}=0{,}0123$$

De $25\%$ para $1{,}2\%$. É o valor da informação extra: cada momento a mais aperta a cota.

**Interpretação prática.** A cota de Markov é grosseira porque precisa acomodar distribuições muito estranhas. Numa empresa real, receitas diárias oscilam em torno da média com desvio moderado, e a probabilidade verdadeira de quadruplicar seria muito menor que $25\%$. Markov não erra — ela só é cautelosa demais por não saber nada sobre a forma da distribuição.

---

# Exercício 2 — Dimensionamento amostral

> Seja $Y_1,\dots,Y_n$ uma amostra aleatória de uma v.a. com média $\mu$ e variância $\sigma^2=9$. Quantas observações $n$, no mínimo, são necessárias para que a probabilidade da média amostral $\bar{Y}_n$ diferir de $\mu$ por menos de $0{,}5$ seja de pelo menos $95\%$?

## 2.1 Diagnóstico

O enunciado dá média (não especificada) e variância ($\sigma^2=9$), sem distribuição. Pede-se o $n$ que garante:

$$P\left(|\bar{Y}_n-\mu|<0{,}5\right)\geq 0{,}95$$

Como este é o capítulo de desigualdades, a resposta esperada é via **Chebyshev**. Apresento também a versão por TCL, porque a comparação entre as duas é instrutiva e pode ser cobrada.

## 2.2 Passo comum: momentos da média amostral

Por linearidade da esperança:

$$E(\bar{Y}_n)=E\left(\frac{1}{n}\sum_{i=1}^n Y_i\right)=\frac{1}{n}\sum_{i=1}^n E(Y_i)=\frac{n\mu}{n}=\mu$$

Pela independência (variância da soma é a soma das variâncias) e lembrando que constante sai da variância **ao quadrado**:

$$\operatorname{Var}(\bar{Y}_n)=\frac{1}{n^2}\sum_{i=1}^n\operatorname{Var}(Y_i)=\frac{n\sigma^2}{n^2}=\frac{\sigma^2}{n}=\frac{9}{n}$$

**Este é o fato que faz o problema ter solução.** A variância decai como $1/n$; aumentando a amostra, a média amostral se concentra em torno de $\mu$ e a probabilidade pedida cresce. Sem essa dependência em $n$, não haveria o que resolver.

## 2.3 Solução via Chebyshev

Aplicando a desigualdade **à média amostral** (não aos $Y_i$ individuais — atenção), com $\varepsilon=0{,}5$:

$$P\left(|\bar{Y}_n-\mu|\geq 0{,}5\right)\;\leq\;\frac{\operatorname{Var}(\bar{Y}_n)}{(0{,}5)^2}=\frac{9/n}{0{,}25}=\frac{36}{n}$$

Chebyshev limita a **cauda**, mas o exercício pede a **região central**. Passando ao complementar — e a desigualdade **inverte de sentido** ao subtrair de 1:

$$P\left(|\bar{Y}_n-\mu|<0{,}5\right)=1-P\left(|\bar{Y}_n-\mu|\geq 0{,}5\right)\;\geq\;1-\frac{36}{n}$$

Agora a lógica do dimensionamento: temos a cadeia

$$P\left(|\bar{Y}_n-\mu|<0{,}5\right)\;\geq\;1-\frac{36}{n}\;\geq\;0{,}95$$

Se garantirmos a **segunda** desigualdade, a primeira transporta a garantia para a probabilidade verdadeira. Estamos forçando o **pior caso** a ser aceitável. Impondo:

$$1-\frac{36}{n}\geq 0{,}95 \quad\Longleftrightarrow\quad \frac{36}{n}\leq 0{,}05 \quad\Longleftrightarrow\quad n\geq\frac{36}{0{,}05}=720$$

**Resposta (Chebyshev): $n_{\min}=720$.**

## 2.4 Solução via TCL

Pelo Teorema Central do Limite, para $n$ grande:

$$Z=\frac{\bar{Y}_n-\mu}{\sigma/\sqrt{n}}=\frac{\bar{Y}_n-\mu}{3/\sqrt{n}}\;\stackrel{a}{\sim}\;N(0,1)$$

**A sacada algébrica:** o evento $|\bar{Y}_n-\mu|<0{,}5$ não está tabelado. Divida os dois lados pelo erro padrão $3/\sqrt{n}$ — que é uma constante positiva, portanto não altera o evento:

$$|\bar{Y}_n-\mu|<0{,}5 \quad\Longleftrightarrow\quad \frac{|\bar{Y}_n-\mu|}{3/\sqrt{n}}<\frac{0{,}5}{3/\sqrt{n}}=\frac{0{,}5\sqrt{n}}{3}$$

O lado esquerdo é $|Z|$. Chamando $c=\frac{0{,}5\sqrt{n}}{3}$, e usando a simetria da normal padrão:

$$P(|Z|<c)=\Phi(c)-\Phi(-c)=2\Phi(c)-1\;\geq\;0{,}95$$

$$\Phi(c)\geq 0{,}975$$

A probabilidade $0{,}05$ restante se divide igualmente entre as duas caudas ($0{,}025$ de cada lado) — daí o $0{,}975$, e não $0{,}95$. Como $\Phi$ é estritamente crescente:

$$c\geq z_{0{,}975}=1{,}96$$

$$\frac{0{,}5\sqrt{n}}{3}\geq 1{,}96 \quad\Longrightarrow\quad \sqrt{n}\geq\frac{1{,}96\times 3}{0{,}5}=11{,}76 \quad\Longrightarrow\quad n\geq 138{,}30$$

**Resposta (TCL): $n_{\min}=139$.**

Verificação: com $n=139$, $c=1{,}9639>1{,}96$ (satisfaz); com $n=138$, $c=1{,}9568<1{,}96$ (não satisfaz).

## 2.5 Comentários do professor

**Arredondamento.** Em problemas de tamanho amostral, **sempre para cima**, mesmo que a parte decimal seja mínima. $n=138$ daria probabilidade ligeiramente abaixo de $95\%$, violando o "pelo menos".

**A comparação é o aprendizado principal:**

| Método | Hipótese | $n_{\min}$ |
|---|---|---|
| Chebyshev | apenas $\sigma^2<\infty$ | 720 |
| TCL | $n$ grande (aproximação) | 139 |

Chebyshev pede **5 vezes mais dados**. Por quê? Porque ela precisa funcionar para a distribuição mais desfavorável possível com variância 9 — inclusive distribuições de cauda pesadíssima. O TCL usa a informação adicional de que a média amostral de muitas observações é aproximadamente normal, e por isso é muito mais econômico.

**A troca fundamental:** Chebyshev dá uma **garantia** (vale sempre, exatamente); o TCL dá uma **aproximação** (boa, mas aproximação). Escolher entre elas é escolher entre rigor e eficiência.

**Coerência interna do TCL.** Note que a resposta $n=139$ é grande, o que justifica retroativamente o uso do TCL. Se a conta tivesse devolvido $n=4$, o resultado seria autocontraditório e deveríamos desconfiar.

**Como identificar qual usar na prova:** se o enunciado menciona "limite superior/inferior" ou o capítulo é de desigualdades → Chebyshev. Se menciona "aproximadamente", "aproximação" ou TCL → normal. Na dúvida, apresente as duas e comente a diferença.

---

# Exercício 3 — Chebyshev na binomial

> Suponha que uma moeda é lançada 100 vezes. Encontre um limite superior para a probabilidade de o número de caras ser de no mínimo 60 ou no máximo 40.

## 3.1 Modelagem

Cada lançamento é um ensaio de Bernoulli independente, com $p=0{,}5$ (moeda honesta — convenção quando o enunciado não especifica). Seja $X$ o número de caras:

$$X\sim\text{Binomial}(n=100,\;p=0{,}5)$$

$$E(X)=np=50 \qquad\qquad \operatorname{Var}(X)=np(1-p)=100\times 0{,}5\times 0{,}5=25$$

$$\sigma=\sqrt{25}=5$$

## 3.2 A sacada central: reescrever a união como um módulo

Esta é a única passagem realmente engenhosa do exercício. Observe a simetria dos pontos de corte em relação à média:

$$60=50+10 \qquad\qquad 40=50-10$$

Ambos estão **à mesma distância** de $\mu=50$. Portanto:

- $X\geq 60 \iff X-50\geq 10$
- $X\leq 40 \iff X-50\leq -10$

A união dos dois eventos é exatamente "o desvio em relação à média tem módulo pelo menos 10":

$$\{X\geq 60\}\cup\{X\leq 40\}=\{|X-50|\geq 10\}$$

**Por que isso é indispensável:** Chebyshev só fala de eventos da forma $|X-E(X)|\geq \varepsilon$. Sem essa reescrita, a ferramenta não se aplica. O enunciado foi construído de propósito com cortes simétricos — se fossem $\{X\geq 60\}\cup\{X\leq 35\}$, seria preciso usar o corte mais frouxo ou tratar as caudas separadamente.

## 3.3 Aplicação

Com $\varepsilon=10$ (a **distância**, não o valor 60 ou 40):

$$P\left(|X-50|\geq 10\right)\;\leq\;\frac{\operatorname{Var}(X)}{\varepsilon^2}=\frac{25}{100}=0{,}25$$

**Resposta: a probabilidade é de no máximo $25\%$.**

**Forma alternativa, mais rápida.** Como $\sigma=5$, a distância 10 equivale a $2\sigma$. Na versão em múltiplos do desvio padrão:

$$P\left(|X-\mu|\geq 2\sigma\right)\leq\frac{1}{2^2}=\frac{1}{4}$$

Mesma resposta em uma linha. Sempre que os cortes caírem em múltiplos inteiros de $\sigma$, use esta forma.

## 3.4 Quão frouxa é a cota?

Aqui, diferentemente do Exercício 1, **conhecemos a distribuição** — então podemos medir o estrago:

| Método | Probabilidade |
|---|---|
| Cota de Chebyshev | $0{,}2500$ |
| Aproximação normal (com correção) | $\approx 0{,}0357$ |
| Valor exato (binomial) | $\approx 0{,}0352$ |

A cota é cerca de **7 vezes** maior que a verdade.

**Isso não é erro.** Chebyshev usa apenas $\mu$ e $\sigma^2$, e precisa valer para **toda** distribuição com $\sigma^2=25$ — inclusive distribuições muito mais dispersas nas caudas do que a binomial, que é bem-comportada e quase simétrica. Quando a distribuição é conhecida, sempre existirá algo melhor.

**A lição metodológica:** o enunciado pediu "um limite superior", não "a probabilidade". Ele já está sinalizando que aceita uma cota. Se pedisse a probabilidade, a resposta correta seria a soma binomial ou a aproximação normal.

---

# Exercício 4 — Binomial exata vs. TCL

> A probabilidade de um jogador de basquete acertar uma cesta é $p=0{,}5$. Em 20 arremessos, qual é a probabilidade de acertar pelo menos 9? Obtenha a probabilidade exata usando a Binomial e uma aproximação usando o TCL.

## 4.1 Modelagem

Quatro condições caracterizam um experimento binomial, e todas estão presentes: número fixo de ensaios ($n=20$), dois resultados por ensaio (acerta/erra), probabilidade constante ($p=0{,}5$) e independência.

$$X\sim\text{Binomial}(20;\,0{,}5)$$

## 4.2 Parte (a) — Probabilidade exata

### Escolha do complementar

$\{X\geq 9\}$ envolve 12 termos ($k=9,\dots,20$); o complementar $\{X\leq 8\}$ envolve 9 termos. Calcula-se o lado com menos parcelas:

$$P(X\geq 9)=1-P(X\leq 8)=1-\sum_{k=0}^{8}\binom{20}{k}(0{,}5)^k(0{,}5)^{20-k}$$

**Atenção ao ponto de corte:** o complementar de "pelo menos 9" é "no máximo **8**", não 9. Em variáveis discretas essa fronteira é fonte permanente de erro.

### A simplificação que só existe porque $p=0{,}5$

No caso geral, cada termo tem potências diferentes de $p$ e $1-p$. Mas com $p=1-p=0{,}5$:

$$p^k(1-p)^{20-k}=(0{,}5)^k(0{,}5)^{20-k}=(0{,}5)^{20}$$

O fator é **o mesmo para todo $k$** e sai em evidência:

$$P(X\leq 8)=\frac{1}{2^{20}}\sum_{k=0}^{8}\binom{20}{k}=\frac{1}{1\,048\,576}\sum_{k=0}^{8}\binom{20}{k}$$

Somando os coeficientes:

| $k$ | $\binom{20}{k}$ |
|---|---|
| 0 | 1 |
| 1 | 20 |
| 2 | 190 |
| 3 | 1.140 |
| 4 | 4.845 |
| 5 | 15.504 |
| 6 | 38.760 |
| 7 | 77.520 |
| 8 | 125.970 |
| **Soma** | **263.950** |

$$P(X\leq 8)=\frac{263\,950}{1\,048\,576}=0{,}251722$$

$$P(X\geq 9)=1-0{,}251722=0{,}748278\approx 0{,}7483$$

### Verificação pela simetria

Quando $p=0{,}5$ a binomial é simétrica em torno de $\mu=10$, isto é, $P(X=k)=P(X=20-k)$. Logo:

$$P(X\leq 9)=\frac{1-P(X=10)}{2}=\frac{1-0{,}176197}{2}=0{,}411901$$

Conferindo pela soma direta: $\frac{263\,950+167\,960}{1\,048\,576}=0{,}411901$. Confere. **Use este truque em prova para validar a soma dos coeficientes.**

## 4.3 Parte (b) — Aproximação pelo TCL

### Momentos

$$E(X)=np=10 \qquad \operatorname{Var}(X)=np(1-p)=5 \qquad \sigma=\sqrt{5}=2{,}2361$$

Verificação das condições: $np=10\geq 5$ e $n(1-p)=10\geq 5$. A simetria ($p=0{,}5$) torna a aproximação especialmente boa mesmo com $n$ modesto.

### Correção de continuidade

Este é o ponto que mais custa nota na prova. Estamos aproximando uma variável **discreta** por uma **contínua**. Na binomial, $P(X\geq 9)$ inclui integralmente a barra do valor 9; na normal, cortar exatamente em 9 descarta metade dessa barra.

$$\{X\geq 9\}\;\longrightarrow\;\{X>8{,}5\}$$

### Padronização

$$z=\frac{8{,}5-10}{2{,}2361}=\frac{-1{,}5}{2{,}2361}=-0{,}6708$$

$$P(X\geq 9)\approx P(Z\geq -0{,}6708)=\Phi(0{,}6708)\approx 0{,}7488$$

**Uso da tabela.** Se sua tabela é do tipo "área entre 0 e $z$" (a mais comum em livros brasileiros), ela devolve $A(0{,}67)=0{,}2486$. Então:

$$P(Z\geq -0{,}67)=\underbrace{P(-0{,}67\leq Z\leq 0)}_{=A(0{,}67)=0{,}2486}+\underbrace{P(Z\geq 0)}_{0{,}5}=0{,}7486$$

O primeiro pedaço vem da simetria: a área de $-0{,}67$ a $0$ é igual à de $0$ a $0{,}67$. Para $z>0$, a regra geral é $\Phi(z)=0{,}5+A(z)$.

Como saber qual tabela você tem:

| Tipo | Valor em $z=0$ | Valor em $z=3$ | Conversão para $\Phi(z)$, $z>0$ |
|---|---|---|---|
| Área entre 0 e $z$ | $0{,}0000$ | $0{,}4987$ | $0{,}5+A(z)$ |
| Acumulada $\Phi(z)$ | $0{,}5000$ | $0{,}9987$ | já é o valor |
| Cauda superior $Q(z)$ | $0{,}5000$ | $0{,}0013$ | $1-Q(z)$ |

## 4.4 Comparação e a importância da correção

| Método | Resultado | Erro absoluto |
|---|---|---|
| Binomial exata | $0{,}7483$ | — |
| TCL **com** correção | $0{,}7488$ | $0{,}0005$ |
| TCL **sem** correção | $0{,}6726$ | $0{,}0757$ |

Sem correção, o erro é **150 vezes maior**. O efeito é dramático aqui porque $n$ é pequeno e o corte ($k=9$) está muito próximo da média — região onde a barra descartada é justamente a mais alta da distribuição.

Sem a correção teríamos $z=\frac{9-10}{2{,}2361}=-0{,}4472$ e $\Phi(0{,}4472)=0{,}6726$.

**Regra:** ao aproximar discreta por normal, **sempre** aplique a correção.

---

# Exercício 5 — Média amostral de exponenciais

> O tempo de vida de um dispositivo eletrônico segue distribuição exponencial com parâmetro $\lambda=100$. Uma amostra i.i.d. de tamanho $n$ é retirada. Qual é a distribuição aproximada da média amostral?

## 5.1 Alerta de parametrização

Antes de qualquer conta, é preciso resolver uma ambiguidade real do enunciado. A exponencial aparece em duas convenções:

| Convenção | Densidade | $E(Y)$ | $\operatorname{Var}(Y)$ |
|---|---|---|---|
| **Taxa** (padrão em teoria) | $f(y)=\lambda e^{-\lambda y}$ | $1/\lambda$ | $1/\lambda^2$ |
| **Média/escala** | $f(y)=\frac{1}{\lambda}e^{-y/\lambda}$ | $\lambda$ | $\lambda^2$ |

Com a primeira, o dispositivo dura em média $0{,}01$ unidade de tempo — pouco plausível no contexto. Com a segunda, dura 100 unidades, o que faz mais sentido fisicamente.

**Conduta em prova:** declare explicitamente qual convenção adotou e siga coerente. A estrutura da resposta é idêntica; muda só o número. Apresento a taxa como principal, por ser o que "parâmetro $\lambda$" significa formalmente.

## 5.2 Momentos da população

Com $Y_i\sim\text{Exp}(\lambda=100)$ na parametrização por taxa:

$$\mu=E(Y_i)=\frac{1}{\lambda}=\frac{1}{100}=0{,}01$$

$$\sigma^2=\operatorname{Var}(Y_i)=\frac{1}{\lambda^2}=\frac{1}{10\,000}=0{,}0001$$

**Observação interessante:** na exponencial, $\sigma=1/\lambda=\mu$. O desvio padrão é igual à média, ou seja, o coeficiente de variação é exatamente 1. É uma distribuição bastante dispersa em relação ao seu centro.

## 5.3 Momentos da média amostral

$$E(\bar{Y}_n)=\mu=\frac{1}{\lambda}=0{,}01$$

$$\operatorname{Var}(\bar{Y}_n)=\frac{\sigma^2}{n}=\frac{1}{n\lambda^2}=\frac{0{,}0001}{n}$$

Estes dois resultados são **exatos**, valem para todo $n$. O que o TCL acrescenta é apenas o **formato** da distribuição.

## 5.4 Aplicação do TCL

O que o TCL exige é somente que média e variância existam e sejam finitas — a exponencial satisfaz ambas. A forma da densidade não importa.

$$\frac{\bar{Y}_n-1/\lambda}{1/(\lambda\sqrt{n})}=\sqrt{n}\,\lambda\left(\bar{Y}_n-\frac{1}{\lambda}\right)\;\xrightarrow{d}\;N(0,1)$$

**Resposta:**

$$\bar{Y}_n\;\stackrel{a}{\sim}\;N\!\left(\frac{1}{\lambda},\;\frac{1}{n\lambda^2}\right)=N\!\left(0{,}01;\;\frac{0{,}0001}{n}\right)$$

**Se a convenção do curso for $\lambda$ = média** ($E(Y_i)=100$, $\operatorname{Var}(Y_i)=10\,000$):

$$\bar{Y}_n\;\stackrel{a}{\sim}\;N\!\left(100;\;\frac{10\,000}{n}\right)$$

## 5.5 Comentários do professor

**Por que "$\stackrel{a}{\sim}$" e não "$\sim$".** Para $n$ finito, $\bar{Y}_n$ **não** é normal. A normalidade só aparece no limite. Escrever "$\bar{Y}_n\sim N(\cdot)$" sem o "a" é erro conceitual.

**Quanto $n$ é necessário aqui?** A exponencial é fortemente assimétrica à direita — seu coeficiente de assimetria é 2, contra 0 da normal. Pelo teorema de Berry–Esseen, o erro da aproximação normal é da ordem de $1/\sqrt{n}$ com constante proporcional à assimetria. Na prática:

| $n$ | Qualidade |
|---|---|
| 5 | ruim — a cauda direita ainda é visível |
| 30 | aceitável no centro, ruim nas caudas |
| 50+ | boa também nas caudas |

Ou seja, para a exponencial a regra "$n\geq 30$" é o **mínimo**, não o confortável.

**Um detalhe elegante: aqui a distribuição exata é conhecida.** A soma de $n$ exponenciais i.i.d. de mesma taxa é uma Gama (Erlang):

$$\sum_{i=1}^n Y_i\sim\text{Gama}(n,\lambda) \qquad\Longrightarrow\qquad \bar{Y}_n\sim\text{Gama}(n,\;n\lambda)$$

Este é um dos raros casos em que não precisaríamos do TCL. Mas o exercício pede explicitamente a distribuição **aproximada**. E o resultado é consistente: a Gama$(n,n\lambda)$ de fato converge para a normal quando $n\to\infty$, o que confirma o TCL por um caminho independente.

**Ilustração computacional sugerida** (não pedida, mas útil para o estudo):

```r
set.seed(42)
lambda <- 100; M <- 10000
for (n in c(5, 30, 100)) {
  medias <- replicate(M, mean(rexp(n, rate = lambda)))
  z <- (medias - 1/lambda) / (1/(lambda*sqrt(n)))
  hist(z, freq = FALSE, breaks = 50, main = paste("n =", n))
  curve(dnorm(x), add = TRUE, lwd = 2, col = "red")
}
```

Espera-se ver o histograma claramente assimétrico à direita em $n=5$, quase simétrico em $n=30$, e praticamente sobreposto à normal em $n=100$.

---

# Exercício 6 — Poisson: LGN e TCL

> Seja $Y_1,\dots,Y_n$ v.a. i.i.d. da distribuição de Poisson com parâmetro $\lambda$.
> **a)** Mostre que a média amostral converge em probabilidade para $\lambda$ quando $n\to\infty$.
> **b)** Encontre a distribuição aproximada da média amostral.
> **c)** Faça uma ilustração computacional e compare a distribuição empírica com a aproximada.

## 6.1 A peculiaridade da Poisson

$$E(Y_i)=\lambda \qquad\qquad \operatorname{Var}(Y_i)=\lambda$$

**Média e variância coincidem.** Um único parâmetro governa centro e dispersão. Isso simplifica todas as contas: não há um $\sigma^2$ separado para carregar.

*Nota prática:* é justamente por isso que, em análise de dados de contagem, "variância amostral muito maior que a média amostral" é sinal de que a Poisson não serve — o fenômeno da **superdispersão**, que motiva o uso da binomial negativa.

## 6.2 Item (a) — Convergência em probabilidade

### O que precisamos provar

Pela definição:

$$\bar{Y}_n\xrightarrow{p}\lambda \quad\Longleftrightarrow\quad \forall\,\varepsilon>0,\;\;\lim_{n\to\infty}P\left(|\bar{Y}_n-\lambda|\geq\varepsilon\right)=0$$

Em palavras: fixe uma tolerância $\varepsilon$ tão pequena quanto quiser; a chance de a média amostral errar $\lambda$ por mais que isso vai a zero. Não se afirma que $\bar{Y}_n$ "chega" em $\lambda$ — ela continua aleatória para todo $n$ —, mas que sua massa de probabilidade se concentra em torno de $\lambda$.

### Momentos da média amostral

$$E(\bar{Y}_n)=\frac{1}{n}\sum_{i=1}^n E(Y_i)=\frac{n\lambda}{n}=\lambda$$

$$\operatorname{Var}(\bar{Y}_n)=\frac{1}{n^2}\sum_{i=1}^n\operatorname{Var}(Y_i)=\frac{n\lambda}{n^2}=\frac{\lambda}{n}$$

### Demonstração via Chebyshev

Seja $\varepsilon>0$ **arbitrário** (fixe-o no início — a definição exige que o argumento valha para qualquer $\varepsilon$).

Aplicando Chebyshev a $\bar{Y}_n$:

$$P\left(|\bar{Y}_n-\lambda|\geq\varepsilon\right)\;\leq\;\frac{\operatorname{Var}(\bar{Y}_n)}{\varepsilon^2}=\frac{\lambda}{n\varepsilon^2}$$

Agora **esprema** a probabilidade entre zero e algo que vai a zero:

$$0\;\leq\;P\left(|\bar{Y}_n-\lambda|\geq\varepsilon\right)\;\leq\;\frac{\lambda}{n\varepsilon^2}\;\xrightarrow[n\to\infty]{}\;0$$

Na expressão $\frac{\lambda}{n\varepsilon^2}$, tanto $\lambda$ quanto $\varepsilon$ são constantes; só $n$ cresce, e está no denominador. A fração é do tipo $\text{constante}/n\to 0$.

Pelo teorema do confronto, o limite do meio é zero. Como $\varepsilon>0$ era arbitrário:

$$\bar{Y}_n\;\xrightarrow{p}\;\lambda$$

$\blacksquare$

**Alternativa em uma linha:** como $E(Y_i)=\lambda<\infty$, o resultado é consequência imediata da Lei Fraca dos Grandes Números.

### Por que Chebyshev é a ferramenta certa

Chebyshev é a **ponte** entre variância (que sabemos calcular) e probabilidade (que queremos limitar). O ponto crucial é que $\operatorname{Var}(\bar{Y}_n)=\lambda/n$ **depende de $n$** e decai — é isso, e só isso, que faz a cota ir a zero. Se a variância fosse constante em $n$, não haveria convergência alguma.

## 6.3 Item (b) — Distribuição aproximada

Aplicando o TCL com $\mu=\lambda$ e $\sigma^2=\lambda$:

$$\frac{\bar{Y}_n-\lambda}{\sqrt{\lambda/n}}\;\xrightarrow{d}\;N(0,1)$$

$$\bar{Y}_n\;\stackrel{a}{\sim}\;N\!\left(\lambda,\;\frac{\lambda}{n}\right)$$

### A distinção conceitual entre (a) e (b)

Esta é a razão de ser do exercício:

| | Item (a) | Item (b) |
|---|---|---|
| Tipo de convergência | em probabilidade | em distribuição |
| Teorema | LGN | TCL |
| O limite é | uma **constante** ($\lambda$) | uma **distribuição** ($N(0,1)$) |
| Escala | $\bar{Y}_n-\lambda\to 0$ | $\sqrt{n}(\bar{Y}_n-\lambda)\to N(0,\lambda)$ |
| Responde | o erro some? | a que velocidade e com que forma? |

A LGN diz que o erro desaparece. O TCL diz **em que ritmo** e **com que formato**. O fator $\sqrt{n}$ é o ajuste de escala exato: sem ele o desvio colapsa em zero; com ele, sobra uma distribuição não-degenerada.

### Observação: a distribuição exata é conhecida

A soma de Poissons independentes é Poisson com parâmetro somado:

$$\sum_{i=1}^n Y_i\sim\text{Poisson}(n\lambda) \qquad\Longrightarrow\qquad \bar{Y}_n=\frac{1}{n}\,\text{Poisson}(n\lambda)$$

Isso confirma o resultado por outro caminho, já que a Poisson$(n\lambda)$ tende à normal quando $n\lambda\to\infty$.

**Consequência importante:** a qualidade da aproximação depende do **produto $n\lambda$**, não só de $n$. Com $\lambda$ muito pequeno (digamos $0{,}1$), é preciso $n$ bem maior para a assimetria se dissipar. Regra prática: $n\lambda\gtrsim 10$.

## 6.4 Item (c) — Ilustração computacional

### O que o script precisa fazer

A ideia é **construir empiricamente a distribuição amostral** de $\bar{Y}_n$ e confrontá-la com a curva teórica. A lógica em quatro passos:

1. **Fixar** $\lambda$ e um valor de $n$.
2. **Repetir $M$ vezes** (com $M$ grande, tipicamente 10.000): sortear uma amostra de tamanho $n$ da Poisson$(\lambda)$ e guardar a média amostral. Ao final, temos $M$ realizações independentes de $\bar{Y}_n$ — uma amostra da própria distribuição amostral.
3. **Padronizar** cada média: $z_j=\frac{\bar{y}_j-\lambda}{\sqrt{\lambda/n}}$. Padronizar é o que permite comparar diferentes $n$ na mesma escala e sobrepor sempre a mesma $N(0,1)$.
4. **Comparar**: histograma de densidade dos $z_j$ com a curva $\phi(z)$ sobreposta; e um QQ-plot normal, que é mais sensível ao comportamento das caudas.

### Comportamento esperado

| Situação | O que se observa |
|---|---|
| $n$ pequeno, $\lambda$ pequeno | histograma visivelmente assimétrico à direita; degraus discretos aparentes (a média assume valores múltiplos de $1/n$); QQ-plot com curvatura em "S" |
| $n$ moderado | assimetria quase some; centro já bem ajustado, caudas ainda imperfeitas |
| $n$ grande | sobreposição quase perfeita; QQ-plot praticamente sobre a reta |
| $\lambda$ grande | convergência muito mais rápida — a própria Poisson já é quase normal |

O ponto pedagógico central: **é o produto $n\lambda$ que comanda**. Simular $n=10,\lambda=10$ e $n=100,\lambda=1$ produz aproximações de qualidade semelhante.

### Script em R

```r
set.seed(2024)
M <- 10000                 # número de repetições Monte Carlo
lambda <- 2
ns <- c(5, 30, 200)

par(mfrow = c(2, 3))

for (n in ns) {
  # 1) M médias amostrais de tamanho n
  medias <- replicate(M, mean(rpois(n, lambda)))
  
  # 2) padronização pelo TCL
  z <- (medias - lambda) / sqrt(lambda / n)
  
  # 3) histograma empírico vs. densidade N(0,1)
  hist(z, freq = FALSE, breaks = 40, col = "grey85", border = "white",
       main = paste0("n = ", n, "  (n*lambda = ", n*lambda, ")"),
       xlab = "estatística padronizada", xlim = c(-4, 4))
  curve(dnorm(x), add = TRUE, lwd = 2, col = "red")
  
  # diagnóstico numérico
  cat(sprintf("n=%4d | media=%6.3f | var=%6.3f | assimetria=%6.3f\n",
              n, mean(z), var(z), mean((z - mean(z))^3) / sd(z)^3))
}

for (n in ns) {
  medias <- replicate(M, mean(rpois(n, lambda)))
  z <- (medias - lambda) / sqrt(lambda / n)
  qqnorm(z, main = paste("QQ-plot, n =", n), pch = 20, cex = 0.4)
  qqline(z, col = "red", lwd = 2)
}
```

### Script em Python

```python
import numpy as np
import matplotlib.pyplot as plt
from scipy import stats

rng = np.random.default_rng(2024)
M, lam = 10_000, 2
ns = [5, 30, 200]

fig, axes = plt.subplots(2, 3, figsize=(15, 8))

for j, n in enumerate(ns):
    # matriz M x n: cada linha é uma amostra
    amostras = rng.poisson(lam, size=(M, n))
    medias = amostras.mean(axis=1)
    z = (medias - lam) / np.sqrt(lam / n)

    ax = axes[0, j]
    ax.hist(z, bins=40, density=True, color="lightgrey", edgecolor="white")
    grid = np.linspace(-4, 4, 400)
    ax.plot(grid, stats.norm.pdf(grid), "r-", lw=2)
    ax.set_title(f"n = {n}  (n·λ = {n*lam})")

    stats.probplot(z, dist="norm", plot=axes[1, j])
    axes[1, j].set_title(f"QQ-plot, n = {n}")

    print(f"n={n:4d} | média={z.mean():6.3f} | var={z.var():6.3f} "
          f"| assimetria={stats.skew(z):6.3f}")

plt.tight_layout(); plt.show()
```

### Como ler a saída numérica

A média dos $z$ deve ficar próxima de 0 e a variância próxima de 1 **já para $n$ pequeno** — porque esses dois momentos são exatos, não aproximados. O que muda com $n$ é a **assimetria**, que deve cair aproximadamente como $1/\sqrt{n}$:

$$\text{assimetria}(\bar{Y}_n)=\frac{1}{\sqrt{n\lambda}}$$

Com $\lambda=2$: em $n=5$ espera-se $\approx 0{,}32$; em $n=30$, $\approx 0{,}13$; em $n=200$, $\approx 0{,}05$. **Confira esses números na saída** — é a melhor evidência de que a simulação está correta e de que a convergência acontece na velocidade prevista pela teoria.

---

# Exercício 7 — Estimador de variância com μ conhecido

> Seja $Y_1,\dots,Y_n$ i.i.d. com $E(Y_i)=\mu$ e $V(Y_i)=\sigma^2<\infty$, e seja
> $$\hat{\sigma}^2=\frac{1}{n}\sum_{i=1}^n(Y_i-\mu)^2$$
> **a)** Mostre que $\hat{\sigma}^2\xrightarrow{p}\sigma^2$. **b)** Obtenha a distribuição aproximada de $\hat{\sigma}^2$.
> **c)** Ilustre computacionalmente para $n=50,250,1000$. **d)** Compare com a distribuição exata.

## 7.1 A sacada central: reconhecer uma média amostral disfarçada

O exercício parece novo, mas não é. Observe a forma:

$$\hat{\sigma}^2=\frac{1}{n}\sum_{i=1}^n\underbrace{(Y_i-\mu)^2}_{\text{chame de }W_i}=\frac{1}{n}\sum_{i=1}^n W_i=\bar{W}_n$$

**É uma média amostral das $W_i$.** E para médias amostrais de variáveis i.i.d. já temos LGN e TCL prontos.

Toda a dificuldade se reduz a **calcular $E(W_i)$ e $\operatorname{Var}(W_i)$**. Essa reescrita economiza metade da prova.

### Por que as $W_i$ continuam i.i.d.

Cada $W_i=g(Y_i)$ com $g(y)=(y-\mu)^2$ depende **só** de $Y_i$. Funções mensuráveis de variáveis independentes são independentes; aplicar a mesma função a variáveis identicamente distribuídas produz variáveis identicamente distribuídas. Escreva essa frase na prova — é o que autoriza o uso dos teoremas.

**Detalhe decisivo:** isso só funciona porque $\mu$ é uma **constante conhecida**. Se fosse $(Y_i-\bar{Y}_n)^2$, todos os termos compartilhariam $\bar{Y}_n$ e seriam dependentes entre si — é exatamente o que acontece no Exercício 8, que por isso exige outra estratégia.

## 7.2 Momentos de $W_i$

### A média é a própria definição de variância

$$E(W_i)=E\left[(Y_i-\mu)^2\right]=\operatorname{Var}(Y_i)=\sigma^2$$

Isso é literalmente a definição de variância, já que $\mu=E(Y_i)$. Sem nenhuma conta.

**Consequência imediata:** $E(\hat{\sigma}^2)=\sigma^2$, ou seja, **este estimador é não-viesado com divisor $n$**. O familiar "$n-1$" não aparece aqui porque a correção de Bessel existe para compensar o fato de $\bar{Y}$ ser estimado dos dados. Com $\mu$ conhecido, nenhum grau de liberdade é consumido e $n$ é o divisor correto.

### A variância traz o quarto momento

$$\operatorname{Var}(W_i)=E(W_i^2)-\left[E(W_i)\right]^2$$

O termo $E(W_i^2)$ é onde mora a novidade:

$$E(W_i^2)=E\left[\left((Y_i-\mu)^2\right)^2\right]=E\left[(Y_i-\mu)^4\right]=\mu_4$$

Elevar $W$ ao quadrado eleva $(Y-\mu)$ à **quarta potência**. Faz sentido: para saber quão instável é um estimador de dispersão, é preciso saber quão pesadas são as caudas — e caudas são exatamente o que o quarto momento (curtose) mede.

$$\operatorname{Var}(W_i)=\mu_4-\sigma^4$$

**Notação:** $\sigma^4$ significa $(\sigma^2)^2$, o quadrado da variância. Não confunda com "momento de ordem 4".

$$\operatorname{Var}(\hat{\sigma}^2)=\frac{\operatorname{Var}(W_i)}{n}=\frac{\mu_4-\sigma^4}{n}$$

## 7.3 Item (a) — Convergência

### Via LGN (mais econômico)

$W_1,\dots,W_n$ são i.i.d. com $E(W_i)=\sigma^2<\infty$. Pela Lei Fraca dos Grandes Números:

$$\bar{W}_n\xrightarrow{p}E(W_1) \qquad\Longrightarrow\qquad \hat{\sigma}^2\xrightarrow{p}\sigma^2$$

$\blacksquare$

### Via Chebyshev (se o professor pedir "na mão")

Requer a hipótese adicional $\mu_4<\infty$. Seja $\varepsilon>0$ arbitrário:

$$0\;\leq\;P\left(|\hat{\sigma}^2-\sigma^2|\geq\varepsilon\right)\;\leq\;\frac{\operatorname{Var}(\hat{\sigma}^2)}{\varepsilon^2}=\frac{\mu_4-\sigma^4}{n\varepsilon^2}\;\xrightarrow[n\to\infty]{}\;0$$

Com $\mu_4$, $\sigma^4$ e $\varepsilon$ constantes, a fração é $\text{const}/n\to 0$. Teorema do confronto conclui.

## 7.4 Item (b) — Distribuição aproximada

Supondo $\mu_4<\infty$, aplique o TCL **às variáveis $W_i$**:

$$\frac{\bar{W}_n-E(W_i)}{\sqrt{\operatorname{Var}(W_i)/n}}=\frac{\hat{\sigma}^2-\sigma^2}{\sqrt{(\mu_4-\sigma^4)/n}}\;\xrightarrow{d}\;N(0,1)$$

$$\hat{\sigma}^2\;\stackrel{a}{\sim}\;N\!\left(\sigma^2,\;\frac{\mu_4-\sigma^4}{n}\right)$$

Equivalentemente: $\sqrt{n}\left(\hat{\sigma}^2-\sigma^2\right)\xrightarrow{d}N\left(0,\;\mu_4-\sigma^4\right)$.

**Erro típico:** padronizar usando $\sigma^2$ (a variância dos $Y$). Não é isso. O TCL está sendo aplicado à sequência $W$, então a variância que entra é $\operatorname{Var}(W_i)=\mu_4-\sigma^4$. Curiosamente, $\sigma^2$ aparece aqui em dois papéis distintos: como **centro** da distribuição assintótica e, via $\sigma^4$, dentro da **variância** assintótica.

### Caso normal

Para $Y\sim N(\mu,\sigma^2)$, os momentos centrais pares são $E[(Y-\mu)^{2k}]=\sigma^{2k}(2k-1)!!$. Com $k=2$: $(3)!!=3\cdot 1=3$, logo $\mu_4=3\sigma^4$ e:

$$\mu_4-\sigma^4=3\sigma^4-\sigma^4=2\sigma^4 \qquad\Longrightarrow\qquad \hat{\sigma}^2\stackrel{a}{\sim}N\!\left(\sigma^2,\;\frac{2\sigma^4}{n}\right)$$

### Qual hipótese cada item exige — o ponto que vale nota

| Item | Ferramenta | Exige | Comentário |
|---|---|---|---|
| (a) | LGN de Khintchine | $E(W_i)=\sigma^2<\infty$ | sai **de graça** da hipótese do enunciado |
| (a) | Chebyshev | $\mu_4<\infty$ | hipótese extra |
| (b) | TCL | $\mu_4<\infty$ | hipótese extra, **não dada no enunciado** |

O item (a) se resolve com o que o enunciado forneceu. O item (b) exige a existência do quarto momento — algo que o enunciado **não** garante. Mencionar isso explicitamente separa a resposta correta da resposta completa.

*Exemplo de que a ressalva importa:* uma $t$ de Student com 4 graus de liberdade tem variância finita mas $\mu_4=\infty$. Nela, o item (a) vale e o item (b) simplesmente não se aplica.

## 7.5 Item (c) — Ilustração computacional com $n=50,250,1000$

### Lógica do script

1. Escolher uma população com $\mu$ e $\sigma^2$ conhecidos (use a normal para poder comparar com a distribuição exata no item d).
2. Para cada $n\in\{50,250,1000\}$, repetir $M$ vezes: sortear a amostra, calcular $\hat{\sigma}^2=\frac{1}{n}\sum(y_i-\mu)^2$ **usando o $\mu$ verdadeiro** (não $\bar{y}$ — este é o ponto do exercício).
3. Padronizar por $\sqrt{(\mu_4-\sigma^4)/n}$ e comparar com a $N(0,1)$.

### Comportamento esperado

| $n$ | O que se observa |
|---|---|
| 50 | assimetria à direita ainda perceptível — herança da $\chi^2$, que é assimétrica |
| 250 | assimetria bem menor; ajuste bom no centro |
| 1000 | praticamente indistinguível da normal |

**Por que a assimetria persiste mais aqui do que no caso da média?** Porque estamos somando quadrados, que são variáveis já naturalmente assimétricas à direita (limitadas inferiormente por zero, ilimitadas acima). A assimetria de $\hat{\sigma}^2$ no caso normal é $\sqrt{8/n}$ — em $n=50$ isso dá $0{,}40$, ainda visível; em $n=1000$, $0{,}09$, imperceptível.

### Script em R (itens c e d juntos)

```r
set.seed(7)
M <- 20000
mu <- 10; sigma2 <- 4; sigma <- sqrt(sigma2)
mu4 <- 3 * sigma2^2                 # normal: mu_4 = 3*sigma^4
v_assint <- mu4 - sigma2^2          # = 2*sigma^4

par(mfrow = c(1, 3))
for (n in c(50, 250, 1000)) {
  sig2hat <- replicate(M, {
    y <- rnorm(n, mu, sigma)
    mean((y - mu)^2)                # ATENCAO: mu verdadeiro, nao mean(y)
  })
  
  z <- (sig2hat - sigma2) / sqrt(v_assint / n)
  
  hist(z, freq = FALSE, breaks = 60, col = "grey85", border = "white",
       main = paste("n =", n), xlab = "padronizado", xlim = c(-4, 5))
  curve(dnorm(x), add = TRUE, col = "red", lwd = 2)          # aproximada
  
  # item (d): distribuicao EXATA sob normalidade
  # n*sig2hat/sigma2 ~ qui-quadrado com n g.l.
  curve(dchisq(x*sqrt(v_assint/n)/sigma2*n + n, df = n) *
          n*sqrt(v_assint/n)/sigma2,
        add = TRUE, col = "blue", lwd = 2, lty = 2)
  
  legend("topright", c("N(0,1)", "exata (qui-quad.)"),
         col = c("red","blue"), lty = c(1,2), lwd = 2, bty = "n")
  
  cat(sprintf("n=%5d | assimetria empirica=%.3f | teorica=%.3f\n",
              n, mean((z-mean(z))^3)/sd(z)^3, sqrt(8/n)))
}
```

## 7.6 Item (d) — Distribuição exata vs. aproximada

### A distribuição exata sob normalidade

Se $Y_i\sim N(\mu,\sigma^2)$, então $\frac{Y_i-\mu}{\sigma}\sim N(0,1)$ e portanto $\left(\frac{Y_i-\mu}{\sigma}\right)^2\sim\chi^2_1$. Somando $n$ termos independentes:

$$\frac{1}{\sigma^2}\sum_{i=1}^n(Y_i-\mu)^2=\frac{n\hat{\sigma}^2}{\sigma^2}\;\sim\;\chi^2_n$$

**Note: $n$ graus de liberdade, não $n-1$.** O grau de liberdade "perdido" no caso usual é justamente o custo de estimar $\mu$; aqui $\mu$ é conhecido, então nada se perde.

### Verificação de coerência

A $\chi^2_n$ tem média $n$ e variância $2n$. Logo:

$$E(\hat{\sigma}^2)=\frac{\sigma^2}{n}\cdot n=\sigma^2 \qquad\qquad \operatorname{Var}(\hat{\sigma}^2)=\frac{\sigma^4}{n^2}\cdot 2n=\frac{2\sigma^4}{n}$$

**Exatamente o que o TCL previu** para o caso normal. Excelente sinal de que a conta está certa — o resultado assintótico e o exato concordam nos dois primeiros momentos.

### O que a comparação revela

| Aspecto | Exata ($\chi^2_n$ reescalada) | Aproximada ($N$) |
|---|---|---|
| Suporte | $(0,\infty)$ — nunca negativa | $(-\infty,\infty)$ — permite valores negativos! |
| Assimetria | $\sqrt{8/n}>0$ | $0$ |
| Momentos 1 e 2 | idênticos aos da aproximação | idênticos aos da exata |

A aproximação normal comete dois pecados conceituais: atribui probabilidade positiva a $\hat{\sigma}^2<0$ (impossível) e ignora a assimetria. Ambos são irrelevantes para $n$ grande, mas visíveis para $n$ pequeno.

**Leitura prática do resultado.** A variância assintótica $\frac{\mu_4-\sigma^4}{n}$ diz que **estimar variância é mais difícil em populações de cauda pesada**: quanto maior a curtose, maior $\mu_4$, e mais devagar $\hat{\sigma}^2$ se estabiliza. Vale a pena rodar a simulação também com uma população $t_5$ ou uma mistura de normais para ver esse efeito.

---

# Exercício 8 — Consistência de S² e de S

> Seja $Y_1,\dots,Y_n$ amostra aleatória de uma distribuição com média $\mu$ e variância $\sigma^2<\infty$.
> **a)** Mostre que $S^2=\frac{1}{n-1}\sum_{i=1}^n(Y_i-\bar{Y})^2$ é consistente para $\sigma^2$.
> **b)** Determine para qual valor $S=\sqrt{S^2}$ converge em probabilidade.

## 8.1 Por que este exercício é diferente do anterior

No Exercício 7 tínhamos $\frac{1}{n}\sum(Y_i-\mu)^2$, com $\mu$ **constante conhecida**. Cada parcela dependia só de $Y_i$, as parcelas eram i.i.d., e bastou aplicar a LGN.

Aqui aparece $\bar{Y}$ no lugar de $\mu$, e isso muda tudo: $\bar{Y}$ é **função de toda a amostra**. Os termos $(Y_i-\bar{Y})^2$ compartilham $\bar{Y}$ e portanto **não são independentes entre si**. Não se pode aplicar a LGN à soma como ela está escrita.

**É preciso reorganizar a expressão primeiro.**

## 8.2 Item (a) — A identidade da soma de quadrados

O truque é separar o que é média amostral do que não é. Expandindo o quadrado:

$$\sum_{i=1}^n(Y_i-\bar{Y})^2=\sum_{i=1}^n Y_i^2-2\bar{Y}\sum_{i=1}^n Y_i+\sum_{i=1}^n\bar{Y}^2$$

Usando $\sum Y_i=n\bar{Y}$ e notando que $\bar{Y}^2$ é constante na soma (aparece $n$ vezes):

$$=\sum_{i=1}^n Y_i^2-2\bar{Y}(n\bar{Y})+n\bar{Y}^2=\sum_{i=1}^n Y_i^2-2n\bar{Y}^2+n\bar{Y}^2=\sum_{i=1}^n Y_i^2-n\bar{Y}^2$$

Portanto:

$$S^2=\frac{1}{n-1}\left(\sum_{i=1}^n Y_i^2-n\bar{Y}^2\right)=\frac{n}{n-1}\left(\underbrace{\frac{1}{n}\sum_{i=1}^n Y_i^2}_{\overline{Y^2}}-\bar{Y}^2\right)$$

**Por que esta reescrita resolve o problema:** agora $S^2$ está expresso através de **duas médias amostrais legítimas** — $\overline{Y^2}$ (média das $Y_i^2$) e $\bar{Y}$ (média dos $Y_i$) — cada uma sobre variáveis i.i.d. A LGN se aplica a cada uma separadamente.

## 8.3 Convergência de cada parcela

### Primeira parcela

A hipótese $\sigma^2<\infty$ é exatamente o que garante o segundo momento finito:

$$E(Y_i^2)=\operatorname{Var}(Y_i)+\left[E(Y_i)\right]^2=\sigma^2+\mu^2<\infty$$

Como $Y_1^2,\dots,Y_n^2$ são i.i.d. (funções de i.i.d.) com média finita, a LGN dá:

$$\overline{Y^2}=\frac{1}{n}\sum_{i=1}^n Y_i^2\;\xrightarrow{p}\;\sigma^2+\mu^2$$

### Segunda parcela

Pela LGN aplicada aos próprios $Y_i$:

$$\bar{Y}\;\xrightarrow{p}\;\mu$$

Como $g(x)=x^2$ é contínua em $\mu$, o **Teorema da Aplicação Contínua** dá:

$$\bar{Y}^2\;\xrightarrow{p}\;\mu^2$$

**Cuidado importante:** não se pode dizer "$E(\bar{Y}^2)=\mu^2$" — isso é **falso**, pois $E(\bar{Y}^2)=\mu^2+\sigma^2/n$. É a **convergência em probabilidade** que atravessa a função contínua, não a esperança. Confundir os dois é erro conceitual grave.

### Recombinação por Slutsky

$$\overline{Y^2}-\bar{Y}^2\;\xrightarrow{p}\;(\sigma^2+\mu^2)-\mu^2=\sigma^2$$

Note o cancelamento de $\mu^2$: a estimativa da média "desconta" a si mesma, restando exatamente a variância. É a versão amostral da identidade $\operatorname{Var}(Y)=E(Y^2)-[E(Y)]^2$.

### O fator $n/(n-1)$

$$\frac{n}{n-1}=\frac{1}{1-1/n}\;\longrightarrow\;1$$

É uma sequência **determinística**, não aleatória — basta o limite usual do cálculo. Slutsky permite multiplicá-la pelo limite em probabilidade:

$$S^2=\frac{n}{n-1}\left(\overline{Y^2}-\bar{Y}^2\right)\;\xrightarrow{p}\;1\cdot\sigma^2=\sigma^2$$

$\blacksquare$

**O papel duplo da correção de Bessel:**

| Regime | Efeito do $n-1$ |
|---|---|
| Amostra finita | torna $E(S^2)=\sigma^2$ exatamente (não-viesado) |
| Assintótico | irrelevante — $n$ ou $n-1$ dão o mesmo limite |

O estimador com divisor $n$ é viesado e **igualmente consistente**.

## 8.4 Item (b) — Convergência de S

A função $g(x)=\sqrt{x}$ é contínua em $[0,\infty)$, em particular no ponto $x=\sigma^2$. Como $S^2\xrightarrow{p}\sigma^2$, o Teorema da Aplicação Contínua dá imediatamente:

$$S=\sqrt{S^2}\;\xrightarrow{p}\;\sqrt{\sigma^2}=\sigma$$

**Resposta: $S$ converge em probabilidade para o desvio padrão populacional $\sigma$.**

Nenhuma conta nova — só a verificação da continuidade. É por isso que o Teorema da Aplicação Contínua é tão valioso: um resultado gera vários.

*Detalhe técnico:* se $\sigma^2>0$, a raiz é contínua **e diferenciável** no ponto, o que permitiria aplicar o Método Delta e obter também a distribuição assintótica de $S$. Se $\sigma^2=0$ a raiz ainda é contínua (a convergência vale), mas a derivada explode e o Método Delta falharia — caso degenerado sem interesse prático.

## 8.5 Aprofundamento: S é viesado, mas consistente

Este ponto costuma confundir, e vale detalhar.

### As duas perguntas são diferentes

| | Viés | Consistência |
|---|---|---|
| Objeto | a **esperança** do estimador | a **distribuição** do estimador |
| Pergunta | o centro está no lugar certo? | a massa se concentra no alvo? |
| Regime | um $n$ fixo | $n\to\infty$ |

### O que acontece com S

Sabemos que $E(S^2)=\sigma^2$ exatamente. Seria natural esperar $E(S)=\sigma$. **É falso**, pela desigualdade de Jensen: para $g$ estritamente côncava e $X$ não-degenerada, $E[g(X)]<g(E(X))$. Com $g(x)=\sqrt{x}$:

$$E(S)=E\left(\sqrt{S^2}\right)\;<\;\sqrt{E(S^2)}=\sigma$$

$S$ **subestima** $\sigma$ sistematicamente.

### Exemplo numérico

Suponha que $S^2$ assuma quatro valores equiprováveis:

$$S^2\in\{1,\;4,\;9,\;16\} \qquad E(S^2)=7{,}5=\sigma^2 \qquad \sigma=2{,}7386$$

$$S\in\{1,\;2,\;3,\;4\} \qquad E(S)=2{,}5$$

Erro de $0{,}24$ para baixo. A causa é geométrica — veja o que a raiz fez com cada valor:

| $S^2$ | $S$ | Encolhimento |
|---|---|---|
| 1 | 1 | 0 |
| 4 | 2 | $-2$ |
| 9 | 3 | $-6$ |
| 16 | 4 | $-12$ |

A raiz comprime **muito mais** os valores grandes. As amostras que dão $S^2$ alto são achatadas com força; as que dão $S^2$ baixo quase não se movem. O conjunto inteiro é arrastado para baixo.

### O viés é proporcional à dispersão

Repita com um $S^2$ menos disperso, mesma média $7{,}5$:

| Cenário | $\operatorname{Var}(S^2)$ | $E(S)$ | Viés |
|---|---|---|---|
| $\{1,4,9,16\}$ | $32{,}25$ | $2{,}5000$ | $-0{,}2386$ |
| $\{6,7,8,9\}$ | $1{,}25$ | $2{,}7309$ | $-0{,}0077$ |

A variância caiu 26 vezes; o viés caiu 31 vezes. É o que prevê a expansão de Taylor de segunda ordem, com $g''(x)=-\frac{1}{4x^{3/2}}$:

$$E(S)\approx\sigma-\frac{\operatorname{Var}(S^2)}{8\sigma^3}$$

Conferindo o segundo caso: $\frac{1{,}25}{8\times 2{,}7386^3}=0{,}0076$ — praticamente igual ao $0{,}0077$ observado.

### Por que não há contradição

Junte as duas peças:

1. O viés de $S$ é proporcional a $\operatorname{Var}(S^2)$.
2. Consistência de $S^2$ significa precisamente que $\operatorname{Var}(S^2)\to 0$.

**A mesma coisa que garante a consistência é a coisa que faz o viés desaparecer.** Com $n$ pequeno, $S^2$ oscila muito, a raiz distorce muito. Com $n$ grande, $S^2$ fica grudado em $\sigma^2$ e quase não sobra oscilação para distorcer.

No caso normal, $E(S)=c_n\sigma$ com $c_n=\sqrt{\frac{2}{n-1}}\cdot\frac{\Gamma(n/2)}{\Gamma\left(\frac{n-1}{2}\right)}$:

| $n$ | $c_n$ | Subestimação |
|---|---|---|
| 2 | 0,7979 | 20,2% |
| 5 | 0,9400 | 6,0% |
| 10 | 0,9727 | 2,7% |
| 30 | 0,9915 | 0,9% |
| 100 | 0,9975 | 0,25% |

### As quatro combinações possíveis

As duas propriedades são logicamente independentes:

| | Consistente | Inconsistente |
|---|---|---|
| **Não-viesado** | $\bar{Y}$ para $\mu$ | usar só $Y_1$ para estimar $\mu$ |
| **Viesado** | $S$ para $\sigma$ | $\hat\theta=7$ sempre |

O caso curioso é $Y_1$: é não-viesado ($E(Y_1)=\mu$ exatamente), mas nunca melhora — a variância continua $\sigma^2$ por mais dados que se colete. **Não-viesamento sozinho não vale nada.** Por isso, em inferência assintótica, consistência é a propriedade inegociável e o viés é preocupação secundária, relevante apenas em amostras pequenas.

---

# Exercício 9 — A estatística t e o Teorema de Slutsky

> Sejam $Y_1,\dots,Y_n$ i.i.d. com $E(Y_i)=\mu$ e $V(Y_i)=\sigma^2<\infty$.
> **a)** Argumente sobre a afirmação: para $n\to\infty$, $\;t=\dfrac{\bar{Y}-\mu}{S/\sqrt{n}}\to Z\sim N(0,1)$.
> **b)** Ilustre computacionalmente para $n=50,250,1000$. **c)** Compare com a distribuição exata da estatística $t$.

## 9.1 Por que este exercício é o clímax da lista

Tudo o que veio antes converge aqui. O TCL entrega:

$$\frac{\bar{Y}-\mu}{\sigma/\sqrt{n}}\;\xrightarrow{d}\;N(0,1)$$

Mas isso é **inútil na prática**: $\sigma$ é desconhecido. Na vida real substituímos $\sigma$ por $S$, e a pergunta é se o limite normal sobrevive à substituição.

**A resposta é sim, e o teorema que garante isso é Slutsky.** Este exercício é a justificativa teórica de todo intervalo de confiança e teste $t$ que você já usou.

## 9.2 Item (a) — A demonstração

### A sacada: multiplicar e dividir por σ

$$t=\frac{\bar{Y}-\mu}{S/\sqrt{n}}=\frac{\bar{Y}-\mu}{\sigma/\sqrt{n}}\cdot\frac{\sigma/\sqrt{n}}{S/\sqrt{n}}=\underbrace{\frac{\bar{Y}-\mu}{\sigma/\sqrt{n}}}_{X_n}\cdot\underbrace{\frac{\sigma}{S}}_{A_n}$$

O $\sqrt{n}$ cancela no segundo fator. Esta fatoração **separa** o problema em duas peças que já sabemos tratar:

- $X_n$ é a estatística padronizada pelo $\sigma$ verdadeiro — o TCL cuida dela.
- $A_n$ é o "erro" cometido ao usar $S$ no lugar de $\sigma$ — o Exercício 8 cuida dela.

### Peça 1: convergência em distribuição

Como $\sigma^2<\infty$, pelo Teorema Central do Limite:

$$X_n=\frac{\bar{Y}-\mu}{\sigma/\sqrt{n}}\;\xrightarrow{d}\;N(0,1)$$

### Peça 2: convergência em probabilidade

Pelo Exercício 8, $S\xrightarrow{p}\sigma$. Supondo $\sigma>0$, a função $g(x)=\sigma/x$ é contínua no ponto $x=\sigma$, logo pelo **Teorema da Aplicação Contínua**:

$$A_n=\frac{\sigma}{S}\;\xrightarrow{p}\;\frac{\sigma}{\sigma}=1$$

O limite é uma **constante** — condição indispensável para Slutsky.

### Peça 3: Slutsky junta as duas

Com $X_n\xrightarrow{d}Z\sim N(0,1)$ e $A_n\xrightarrow{p}1$:

$$t=A_n X_n\;\xrightarrow{d}\;1\cdot Z=Z\sim N(0,1)$$

**A afirmação é verdadeira.** $\blacksquare$

### Por que Slutsky é indispensável aqui

Sem Slutsky, o argumento não fecha. Note que $X_n$ e $A_n$ são **fortemente dependentes** — ambos são calculados da mesma amostra. Não se pode dizer "o limite do produto é o produto dos limites" apelando a independência, porque não há independência.

O que Slutsky explora é diferente: quando um dos fatores converge para uma **constante**, sua aleatoriedade desaparece assintoticamente. Ele deixa de ser um fator aleatório e vira, no limite, uma multiplicação por número. A dependência com $X_n$ se torna irrelevante.

**Esta é a essência do teorema.** Se $A_n$ convergisse para algo aleatório, o resultado seria falso.

## 9.3 As ressalvas que uma resposta completa deve conter

O enunciado pede para **argumentar**, ou seja, quer a análise crítica — não só o "verdadeiro".

**1. Hipóteses necessárias.** Além de $\sigma^2<\infty$ (dada), é preciso $\sigma^2>0$. Se $\sigma^2=0$ a variável é degenerada, $S=0$ com probabilidade 1, e $t$ nem está definida.

**2. É convergência em distribuição, não em probabilidade.** $t$ não converge para nenhum número — ela continua flutuando para sempre. O que converge é a **forma** da sua distribuição.

**3. Não vale para $n$ finito.** Para qualquer $n$ fixo, $t$ **não** é $N(0,1)$. Ela tem caudas mais pesadas, porque incorpora duas fontes de aleatoriedade: a do numerador ($\bar{Y}$) e a do denominador ($S$).

**4. Sem normalidade, não há distribuição exata conhecida.** Sob $Y_i\sim N(\mu,\sigma^2)$, sabe-se que $t\sim t_{n-1}$ **exatamente**. Sem essa hipótese, a distribuição de $t$ em amostra finita depende da população e em geral não tem forma fechada — só o limite normal está garantido.

**5. Coerência do resultado.** É reconfortante notar que $t_{n-1}\xrightarrow{d}N(0,1)$ quando $n\to\infty$. Os dois caminhos — exato sob normalidade e assintótico geral — apontam para o mesmo lugar.

## 9.4 Item (b) — Ilustração computacional

### Lógica do script

1. Escolher uma população. **Rode com duas**: uma normal (onde a distribuição exata é conhecida) e uma assimétrica, como a exponencial (onde a aproximação sofre). O contraste é o aprendizado principal.
2. Para cada $n\in\{50,250,1000\}$: repetir $M$ vezes o sorteio da amostra e o cálculo de $t=\frac{\bar{y}-\mu}{s/\sqrt{n}}$, usando o $\mu$ verdadeiro (que conhecemos por termos gerado os dados).
3. Sobrepor ao histograma **três** curvas: a $N(0,1)$, a $t_{n-1}$ e (implicitamente) a empírica.

### Comportamento esperado

| População | $n=50$ | $n=250$ | $n=1000$ |
|---|---|---|---|
| Normal | empírica $\equiv t_{49}$ exatamente; $N(0,1)$ já quase idêntica | indistinguíveis | indistinguíveis |
| Exponencial | empírica **assimétrica à esquerda**; ambas as curvas erram nas caudas | assimetria pequena | ajuste bom |

**O detalhe mais interessante da simulação:** com população exponencial, a distribuição empírica de $t$ é assimétrica **à esquerda**, embora a população seja assimétrica à direita. A causa é a correlação entre $\bar{Y}$ e $S$: amostras com média alta tendem a ter $S$ alto também (na exponencial, $\sigma=\mu$), o que reduz o $t$; e vice-versa. Essa correlação inverte o sinal da assimetria. É um efeito que só se enxerga simulando.

### Script em R

```r
set.seed(99)
M <- 20000
ns <- c(50, 250, 1000)

simula_t <- function(n, gerador, mu_verdadeiro) {
  replicate(M, {
    y <- gerador(n)
    (mean(y) - mu_verdadeiro) / (sd(y) / sqrt(n))
  })
}

par(mfrow = c(2, 3))

# ---- populacao NORMAL ----
for (n in ns) {
  tt <- simula_t(n, function(k) rnorm(k, 10, 2), 10)
  hist(tt, freq = FALSE, breaks = 60, xlim = c(-4, 4),
       col = "grey85", border = "white",
       main = paste("Normal, n =", n), xlab = "t")
  curve(dnorm(x), add = TRUE, col = "red", lwd = 2)
  curve(dt(x, df = n - 1), add = TRUE, col = "blue", lty = 2, lwd = 2)
}

# ---- populacao EXPONENCIAL (assimetrica) ----
for (n in ns) {
  tt <- simula_t(n, function(k) rexp(k, rate = 1), 1)
  hist(tt, freq = FALSE, breaks = 60, xlim = c(-4, 4),
       col = "grey85", border = "white",
       main = paste("Exponencial, n =", n), xlab = "t")
  curve(dnorm(x), add = TRUE, col = "red", lwd = 2)
  cat(sprintf("Exp | n=%5d | assimetria de t = %7.3f\n",
              n, mean((tt-mean(tt))^3)/sd(tt)^3))
}
```

### Script em Python

```python
import numpy as np, matplotlib.pyplot as plt
from scipy import stats

rng = np.random.default_rng(99)
M, ns = 20_000, [50, 250, 1000]
grid = np.linspace(-4, 4, 400)

def simula_t(n, amostrador, mu):
    y = amostrador((M, n))
    return (y.mean(1) - mu) / (y.std(1, ddof=1) / np.sqrt(n))

fig, axes = plt.subplots(2, 3, figsize=(15, 8))

for j, n in enumerate(ns):
    t_norm = simula_t(n, lambda s: rng.normal(10, 2, s), 10)
    ax = axes[0, j]
    ax.hist(t_norm, bins=60, range=(-4, 4), density=True,
            color="lightgrey", edgecolor="white")
    ax.plot(grid, stats.norm.pdf(grid), "r-", lw=2, label="N(0,1)")
    ax.plot(grid, stats.t.pdf(grid, df=n-1), "b--", lw=2, label=f"t({n-1})")
    ax.set_title(f"Normal, n = {n}"); ax.legend()

    t_exp = simula_t(n, lambda s: rng.exponential(1, s), 1)
    ax = axes[1, j]
    ax.hist(t_exp, bins=60, range=(-4, 4), density=True,
            color="lightgrey", edgecolor="white")
    ax.plot(grid, stats.norm.pdf(grid), "r-", lw=2)
    ax.set_title(f"Exponencial, n = {n}")
    print(f"Exp | n={n:5d} | assimetria de t = {stats.skew(t_exp):7.3f}")

plt.tight_layout(); plt.show()
```

## 9.5 Item (c) — Exata vs. aproximada

### Sob normalidade: a distribuição exata é $t_{n-1}$

$$Y_i\sim N(\mu,\sigma^2) \quad\Longrightarrow\quad t=\frac{\bar{Y}-\mu}{S/\sqrt{n}}\;\sim\;t_{n-1}\;\;\text{(exatamente, para todo } n)$$

A construção clássica: $\frac{\bar{Y}-\mu}{\sigma/\sqrt{n}}\sim N(0,1)$, $\frac{(n-1)S^2}{\sigma^2}\sim\chi^2_{n-1}$, e — pelo Teorema de Basu / independência entre $\bar{Y}$ e $S^2$ sob normalidade — o quociente

$$\frac{N(0,1)}{\sqrt{\chi^2_{n-1}/(n-1)}}$$

é, por definição, uma $t$ com $n-1$ graus de liberdade.

### Comparação numérica

Quantis de $97{,}5\%$:

| $n$ | $t_{n-1}$ | $N(0,1)$ | Diferença relativa |
|---|---|---|---|
| 5 | 2,776 | 1,960 | 41,6% |
| 10 | 2,262 | 1,960 | 15,4% |
| 30 | 2,045 | 1,960 | 4,3% |
| 50 | 2,010 | 1,960 | 2,5% |
| 250 | 1,970 | 1,960 | 0,5% |
| 1000 | 1,962 | 1,960 | 0,1% |

**Conclusão prática:** com $n=50$ a diferença já é de 2,5%; com $n\geq 250$ é irrelevante. É por isso que na prática se usa $z=1{,}96$ sem constrangimento em amostras grandes — mas se insiste na $t$ em amostras pequenas.

### As duas fontes de erro — o quadro que amarra tudo

Ao usar $N(0,1)$ para $t$ em amostra finita, cometem-se **dois** erros distintos, e é essencial não confundi-los:

| Fonte do erro | Causa | Some quando | Corrigida por |
|---|---|---|---|
| **Estimar $\sigma$ por $S$** | denominador aleatório engorda as caudas | $n$ cresce | usar $t_{n-1}$ |
| **População não-normal** | assimetria e curtose da população | $n$ cresce (via TCL) | nada — só mais dados |

Sob normalidade, apenas o primeiro erro existe, e a $t_{n-1}$ o elimina **exatamente**. Sob população assimétrica, o segundo erro aparece e a $t_{n-1}$ **não ajuda** — ela foi deduzida sob normalidade, e usá-la fora desse contexto não compra nada.

Este é o motivo pelo qual, na simulação com exponencial, tanto a $N(0,1)$ quanto a $t_{n-1}$ erram: o problema ali não é a estimação de $\sigma$, é a assimetria da população.

---

# Exercício 10 — Diferença de Poissons studentizada

> Sejam $X_n$ e $Y_m$ v.a. independentes com distribuição Poisson de parâmetros $n$ e $m$, respectivamente.
> **a)** Mostre que $\;R=\dfrac{(X_n-n)-(Y_m-m)}{\sqrt{X_n+Y_m}}\to Z\sim N(0,1)$.
> **b)** Ilustre computacionalmente e compare com a distribuição empírica de $R$.

## 10.1 Leitura do problema

Repare que $n$ e $m$ aqui **não são tamanhos de amostra** — são os **parâmetros** das Poissons. Temos uma única observação de cada distribuição, e o que cresce são os parâmetros. Assumimos $n,m\to\infty$.

Reconheça a estrutura de $R$:

- **Numerador:** cada variável centrada na sua média ($E(X_n)=n$, $E(Y_m)=m$) e depois subtraídas.
- **Denominador:** $\sqrt{X_n+Y_m}$ — um **estimador** do desvio padrão do numerador, construído com os próprios dados.

Ou seja, **é uma estatística studentizada**, exatamente como a $t$ do Exercício 9. A estratégia de prova será a mesma: TCL no numerador, LGN + Aplicação Contínua no denominador, Slutsky para juntar.

## 10.2 Item (a) — Passo 1: momentos do numerador

Seja $N=(X_n-n)-(Y_m-m)$.

$$E(N)=\underbrace{E(X_n)-n}_{0}-\underbrace{\left[E(Y_m)-m\right]}_{0}=0$$

Pela **independência** entre $X_n$ e $Y_m$ (a variância da diferença é a **soma** das variâncias — o sinal negativo some ao elevar ao quadrado):

$$\operatorname{Var}(N)=\operatorname{Var}(X_n)+\operatorname{Var}(Y_m)=n+m$$

Aqui usamos a propriedade da Poisson: $\operatorname{Var}=\text{média}$.

## 10.3 Passo 2: TCL no numerador

A Poisson tem a propriedade de **infinita divisibilidade**: $X_n\sim\text{Poisson}(n)$ pode ser escrita como soma de $n$ variáveis i.i.d. Poisson$(1)$:

$$X_n\;\stackrel{d}{=}\;\sum_{i=1}^{n}U_i,\qquad U_i\stackrel{iid}{\sim}\text{Poisson}(1)$$

Analogamente $Y_m\stackrel{d}{=}\sum_{j=1}^{m}V_j$ com $V_j\sim\text{Poisson}(1)$. Então:

$$N=\sum_{i=1}^{n}(U_i-1)-\sum_{j=1}^{m}(V_j-1)$$

**$N$ é uma soma de $n+m$ variáveis independentes, cada uma com média 0 e variância 1** (as $-(V_j-1)$ também têm média 0 e variância 1). Aplicando o TCL a essa soma:

$$\frac{N}{\sqrt{n+m}}\;\xrightarrow{d}\;N(0,1) \qquad\text{quando } n+m\to\infty$$

**Por que este passo merece cuidado:** o TCL, na forma que conhecemos, fala de médias de amostras i.i.d. Para aplicá-lo aqui é preciso **exibir** $N$ como soma de parcelas independentes — e é a infinita divisibilidade da Poisson que permite isso. Sem esse argumento, o passo fica sem justificativa.

## 10.4 Passo 3: o denominador converge para 1 (após normalizar)

Precisamos comparar $\sqrt{X_n+Y_m}$ com $\sqrt{n+m}$. Considere:

$$D=\frac{X_n+Y_m}{n+m}$$

Como $X_n+Y_m\sim\text{Poisson}(n+m)$ (soma de Poissons independentes):

$$E(D)=\frac{n+m}{n+m}=1 \qquad\qquad \operatorname{Var}(D)=\frac{n+m}{(n+m)^2}=\frac{1}{n+m}$$

Por Chebyshev, para todo $\varepsilon>0$:

$$P\left(|D-1|\geq\varepsilon\right)\leq\frac{1}{(n+m)\varepsilon^2}\;\xrightarrow[n+m\to\infty]{}\;0$$

$$D\;\xrightarrow{p}\;1$$

Como $g(x)=\sqrt{x}$ é contínua em $x=1$, pelo Teorema da Aplicação Contínua:

$$\sqrt{D}=\sqrt{\frac{X_n+Y_m}{n+m}}\;\xrightarrow{p}\;1$$

## 10.5 Passo 4: Slutsky fecha o argumento

A sacada final é **multiplicar e dividir por $\sqrt{n+m}$**:

$$R=\frac{N}{\sqrt{X_n+Y_m}}=\frac{N}{\sqrt{n+m}}\cdot\frac{\sqrt{n+m}}{\sqrt{X_n+Y_m}}=\underbrace{\frac{N}{\sqrt{n+m}}}_{\xrightarrow{d}\,N(0,1)}\cdot\underbrace{\frac{1}{\sqrt{D}}}_{\xrightarrow{p}\,1}$$

Pelo Teorema de Slutsky:

$$R\;\xrightarrow{d}\;N(0,1)\cdot 1=Z\sim N(0,1)$$

$\blacksquare$

## 10.6 Comentários do professor

**A estrutura da prova é sempre a mesma.** Compare com o Exercício 9:

| | Exercício 9 (estatística $t$) | Exercício 10 |
|---|---|---|
| Numerador padronizado pelo valor teórico | $\dfrac{\bar{Y}-\mu}{\sigma/\sqrt{n}}\xrightarrow{d}N(0,1)$ | $\dfrac{N}{\sqrt{n+m}}\xrightarrow{d}N(0,1)$ |
| Razão entre estimador e valor teórico | $\dfrac{\sigma}{S}\xrightarrow{p}1$ | $\dfrac{\sqrt{n+m}}{\sqrt{X_n+Y_m}}\xrightarrow{p}1$ |
| Junção | Slutsky | Slutsky |

**Guarde este roteiro** — ele resolve qualquer problema de estatística studentizada:

1. Padronize pelo desvio padrão **teórico** e aplique o TCL.
2. Mostre que a razão (estimador ÷ teórico) converge em probabilidade para 1.
3. Aplique Slutsky.

**Por que o denominador usa os dados e não os parâmetros?** Porque essa é a situação realista: numa aplicação, $n$ e $m$ podem ser desconhecidos e $X_n+Y_m$ é o que se observa. O exercício mostra que **substituir a variância teórica pela estimada não estraga o limite normal** — que é a mensagem central de toda esta lista.

**Detalhe de rigor:** o resultado exige $n+m\to\infty$. Se apenas $n\to\infty$ com $m$ fixo, ainda funciona, porque $n+m\to\infty$ do mesmo jeito. Mas se ambos ficam limitados, não há convergência.

**Detalhe técnico adicional:** $R$ não está definida quando $X_n+Y_m=0$. Isso ocorre com probabilidade $e^{-(n+m)}$, que vai a zero exponencialmente rápido — logo não afeta a convergência em distribuição. Vale mencionar em uma linha na prova.

## 10.7 Item (b) — Ilustração computacional

### Lógica do script

1. Fixar valores de $n$ e $m$ (comece pequenos, como $n=m=5$, e aumente).
2. Repetir $M$ vezes: sortear $x\sim\text{Poisson}(n)$ e $y\sim\text{Poisson}(m)$ independentes, calcular $r=\frac{(x-n)-(y-m)}{\sqrt{x+y}}$.
3. Descartar (ou tratar) os casos $x+y=0$.
4. Histograma dos $r$ com a $N(0,1)$ sobreposta, mais QQ-plot.

### Comportamento esperado

| $(n,m)$ | O que se observa |
|---|---|
| $(2,2)$ | histograma com "degraus" bem marcados — $R$ assume poucos valores distintos; ajuste ruim |
| $(10,10)$ | discretização ainda visível, formato já normal no centro |
| $(50,50)$ | bom ajuste; QQ-plot quase reto |
| $(500,500)$ | praticamente indistinguível de $N(0,1)$ |
| $(100,5)$ | assimetria residual — o lado com parâmetro pequeno domina o erro |

**O ponto pedagógico:** como $R$ é razão de variáveis **discretas**, para $n,m$ pequenos ela assume um conjunto discreto de valores e o histograma tem aparência de "pente". A convergência à normal é também a suavização dessa granularidade.

### Script em R

```r
set.seed(11)
M <- 50000
cenarios <- list(c(2,2), c(10,10), c(50,50), c(500,500))

par(mfrow = c(2, 4))

for (cen in cenarios) {
  n <- cen[1]; m <- cen[2]
  x <- rpois(M, n); y <- rpois(M, m)
  ok <- (x + y) > 0                     # evita divisao por zero
  r <- ((x[ok] - n) - (y[ok] - m)) / sqrt(x[ok] + y[ok])
  
  hist(r, freq = FALSE, breaks = 80, xlim = c(-4, 4),
       col = "grey85", border = "white",
       main = paste0("n=", n, ", m=", m), xlab = "R")
  curve(dnorm(x), add = TRUE, col = "red", lwd = 2)
  
  cat(sprintf("n=%4d m=%4d | media=%7.4f | var=%6.4f | assim=%7.4f | descartados=%d\n",
              n, m, mean(r), var(r), mean((r-mean(r))^3)/sd(r)^3, sum(!ok)))
}

for (cen in cenarios) {
  n <- cen[1]; m <- cen[2]
  x <- rpois(M, n); y <- rpois(M, m); ok <- (x+y) > 0
  r <- ((x[ok]-n) - (y[ok]-m)) / sqrt(x[ok]+y[ok])
  qqnorm(r, pch = 20, cex = 0.3, main = paste0("QQ: n=", n, ", m=", m))
  qqline(r, col = "red", lwd = 2)
}
```

### Script em Python

```python
import numpy as np, matplotlib.pyplot as plt
from scipy import stats

rng = np.random.default_rng(11)
M = 50_000
cenarios = [(2,2), (10,10), (50,50), (500,500)]
grid = np.linspace(-4, 4, 400)

fig, axes = plt.subplots(2, 4, figsize=(18, 8))

for j, (n, m) in enumerate(cenarios):
    x = rng.poisson(n, M); y = rng.poisson(m, M)
    ok = (x + y) > 0
    r = ((x[ok] - n) - (y[ok] - m)) / np.sqrt(x[ok] + y[ok])

    ax = axes[0, j]
    ax.hist(r, bins=80, range=(-4, 4), density=True,
            color="lightgrey", edgecolor="white")
    ax.plot(grid, stats.norm.pdf(grid), "r-", lw=2)
    ax.set_title(f"n={n}, m={m}")

    stats.probplot(r, dist="norm", plot=axes[1, j])
    axes[1, j].set_title(f"QQ: n={n}, m={m}")

    print(f"n={n:4d} m={m:4d} | média={r.mean():7.4f} | var={r.var():6.4f} "
          f"| assim={stats.skew(r):7.4f} | descartados={(~ok).sum()}")

plt.tight_layout(); plt.show()
```

### Como validar a simulação

Verifique se a média empírica de $R$ fica próxima de 0 e a variância próxima de 1 conforme $(n,m)$ crescem. Para $(n,m)$ pequenos, espere **variância um pouco maior que 1** — o denominador aleatório engorda as caudas, exatamente como acontece com a estatística $t$. É o mesmo fenômeno em outra roupagem.

---

# Apêndice A — Mapa de decisão: qual ferramenta usar

```
O enunciado pede um LIMITE (superior/inferior) para uma probabilidade?
│
├── SIM ─── Tenho apenas a MÉDIA e a variável é ≥ 0?  ──→ MARKOV
│           Tenho MÉDIA e VARIÂNCIA?                   ──→ CHEBYSHEV
│                                                          (reescreva o evento
│                                                           como |X − μ| ≥ ε)
│
└── NÃO ─── O que se pede?
            │
            ├── "converge para uma CONSTANTE"          ──→ LGN
            │   "é consistente"                            (ou Chebyshev com
            │   "converge em probabilidade"                 Var → 0)
            │
            ├── "distribuição APROXIMADA"              ──→ TCL
            │   "aproximadamente normal"                   (padronize por σ/√n;
            │   probabilidade com n grande                  se discreta, corrija
            │                                               a continuidade)
            │
            ├── já sei que Xₙ →ᵖ c e quero g(Xₙ)       ──→ APLICAÇÃO CONTÍNUA
            │                                               (basta g contínua em c)
            │
            └── tenho um limite em DISTRIBUIÇÃO        ──→ SLUTSKY
                multiplicado/dividido por algo que          (fatore em
                converge em PROBABILIDADE a constante        Xₙ · Aₙ)
```

## Roteiro para estatísticas studentizadas

Sempre que a estatística tiver a forma $\dfrac{\text{algo centrado}}{\text{estimador do desvio padrão}}$:

1. **Fatore:** multiplique e divida pelo desvio padrão **teórico**.
2. **TCL** na primeira parte (padronizada corretamente).
3. **LGN + Aplicação Contínua** na segunda: mostre que (teórico ÷ estimado) $\xrightarrow{p}1$.
4. **Slutsky** para concluir.

Este roteiro resolve os Exercícios 9 e 10 — e praticamente todo problema desse tipo em prova.

---

# Apêndice B — Erros clássicos de prova

| # | Erro | Correção |
|---|---|---|
| 1 | Aplicar Chebyshev aos $Y_i$ quando o evento é sobre $\bar{Y}_n$ | A variância que entra é $\operatorname{Var}(\bar{Y}_n)=\sigma^2/n$, não $\sigma^2$ |
| 2 | Esquecer que Chebyshev limita a **cauda**, e o enunciado pede o **centro** | Passe ao complementar; a desigualdade **inverte** |
| 3 | Usar $\varepsilon$ em vez de $\varepsilon^2$ no denominador | $P(\lvert X-\mu\rvert\geq\varepsilon)\leq\sigma^2/\varepsilon^2$ |
| 4 | Arredondar $n$ para baixo | Tamanho amostral **sempre** arredonda para cima |
| 5 | Confundir o corte do complementar em variável discreta | Complementar de "$\geq 9$" é "$\leq 8$" |
| 6 | Esquecer a correção de continuidade | $\{X\geq k\}\to\{X>k-0{,}5\}$; sem ela o erro pode ser 100× maior |
| 7 | Usar $\Phi(0{,}95)$ em problema bilateral | São $0{,}025$ em cada cauda → $z_{0{,}975}=1{,}96$ |
| 8 | Ler tabela "área entre 0 e $z$" como se fosse $\Phi(z)$ | Para $z>0$: $\Phi(z)=0{,}5+A(z)$ |
| 9 | Escrever $\bar{Y}_n\sim N(\cdot)$ sem o "$a$" | Para $n$ finito não é normal; use $\stackrel{a}{\sim}$ |
| 10 | Aplicar LGN a $\sum(Y_i-\bar{Y})^2$ diretamente | As parcelas **não** são independentes; reescreva como $\overline{Y^2}-\bar{Y}^2$ |
| 11 | Dizer "$E(\bar{Y}^2)=\mu^2$" | Falso: $E(\bar{Y}^2)=\mu^2+\sigma^2/n$. Quem atravessa a função contínua é a **convergência**, não a esperança |
| 12 | Padronizar $\hat\sigma^2$ usando $\sigma^2$ | O TCL foi aplicado às $W_i$; a variância é $\mu_4-\sigma^4$ |
| 13 | Aplicar o TCL sem checar variância finita | O TCL exige $\sigma^2<\infty$; a LGN não |
| 14 | Usar Slutsky com $A_n$ convergindo para algo **aleatório** | O limite de $A_n$ tem que ser **constante** |
| 15 | Confundir viés com inconsistência | São independentes: $S$ é viesado **e** consistente |
| 16 | Não declarar as hipóteses ($X\geq 0$ para Markov, $\mu_4<\infty$ para o TCL em $\hat\sigma^2$, $\sigma^2>0$ para a $t$) | Declarar hipóteses vale nota |

---

# Apêndice C — Fórmulas de bolso

**Desigualdades**

$$P(X\geq a)\leq\frac{E(X)}{a}\;\;(X\geq 0) \qquad\qquad P(|X-\mu|\geq\varepsilon)\leq\frac{\sigma^2}{\varepsilon^2}$$

**Média amostral**

$$E(\bar{Y}_n)=\mu \qquad \operatorname{Var}(\bar{Y}_n)=\frac{\sigma^2}{n} \qquad \bar{Y}_n\stackrel{a}{\sim}N\!\left(\mu,\frac{\sigma^2}{n}\right)$$

**Distribuições da lista**

| Distribuição | $E(Y)$ | $\operatorname{Var}(Y)$ | Assimetria |
|---|---|---|---|
| Bernoulli$(p)$ | $p$ | $p(1-p)$ | $\frac{1-2p}{\sqrt{p(1-p)}}$ |
| Binomial$(n,p)$ | $np$ | $np(1-p)$ | $\frac{1-2p}{\sqrt{np(1-p)}}$ |
| Poisson$(\lambda)$ | $\lambda$ | $\lambda$ | $1/\sqrt{\lambda}$ |
| Exponencial(taxa $\lambda$) | $1/\lambda$ | $1/\lambda^2$ | $2$ |
| Normal$(\mu,\sigma^2)$ | $\mu$ | $\sigma^2$ | $0$ |

**Estimadores de variância**

| Estimador | $E(\cdot)$ | $\operatorname{Var}(\cdot)$ | Exato sob normalidade |
|---|---|---|---|
| $\hat\sigma^2=\frac1n\sum(Y_i-\mu)^2$ | $\sigma^2$ | $\frac{\mu_4-\sigma^4}{n}$ | $\frac{n\hat\sigma^2}{\sigma^2}\sim\chi^2_n$ |
| $S^2=\frac{1}{n-1}\sum(Y_i-\bar Y)^2$ | $\sigma^2$ | $\frac{1}{n}\left(\mu_4-\frac{n-3}{n-1}\sigma^4\right)$ | $\frac{(n-1)S^2}{\sigma^2}\sim\chi^2_{n-1}$ |

**Momentos centrais da normal:** $\mu_3=0$, $\mu_4=3\sigma^4$, $\mu_{2k}=\sigma^{2k}(2k-1)!!$

**Quantis usados com frequência**

| Confiança | $z$ bilateral | $z$ unilateral |
|---|---|---|
| 90% | 1,645 | 1,282 |
| 95% | 1,960 | 1,645 |
| 99% | 2,576 | 2,326 |

