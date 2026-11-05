# macrorf.opt

Repositório desenvolvido no âmbito das atividades de pesquisa da UFRGS, focado na análise de desempenho e na otimização do pacote `macrorf`, originalmente proposto por Gouilmette Coulombe.

## Descrição do Projeto

O objetivo principal deste trabalho é identificar gargalos computacionais no código original de previsão macroeconômica, reestruturando as rotinas para melhorar a eficiência de processamento. 

## Estrutura do Repositório

O projeto segue a estrutura padrão de desenvolvimento de pacotes em R, organizando-se da seguinte forma:

* **R**
  * `macrorf_core.R`: Implementação original das funções do pacote base, utilizada como referência e validação dos resultados.
  * `optimization.R`: Scripts contendo as reestruturações, vetorizações e melhorias de performance desenvolvidas na pesquisa.
* **man/**: Documentação gerada automaticamente para as funções do pacote.
* **DESCRIPTION**: Metadados, dependências e informações gerais do pacote.
* **NAMESPACE**: Arquivo de gerenciamento de exportações gerado pelo ambiente de desenvolvimento.

## Dados

Os datasets utilizados para os testes e validações empíricas (`FREDMD` e `FREDQD`) originam-se do artigo e do repositório de pesquisa original do autor.
