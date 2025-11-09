import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import norm
# Parâmetros
lamb = 5 # lambda da Poisson
n = 100 # tamanho da amostra
num_simulations = 10000 # número de médias amostrais a gerar
# Gerar 10000 médias amostrais de tamanho 100
sample_means = [np.mean(np.random.poisson(lamb, n)) for _ in range(num_simulations)]
# Plotar o histograma das médias
plt.hist(sample_means, bins=30, density=True, alpha=0.6, label='Distribuição Empírica')
# Plotar a densidade da Normal teórica (TLC)
mu_normal = lamb
sigma_normal = np.sqrt(lamb / n)
x = np.linspace(mu_normal - 4*sigma_normal, mu_normal + 4*sigma_normal, 100)
plt.plot(x, norm.pdf(x, mu_normal, sigma_normal), 'r-', lw=2, label='Aproximação Normal')
plt.title(f'Distribuição de Médias Amostrais (n={n}) vs. TLC')
plt.xlabel('Média Amostral')
plt.ylabel('Densidade')
plt.legend()
plt.grid(True)
plt.show()