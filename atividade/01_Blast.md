# Usando o Blast

## O que é o Blast?

O BLAST é um programa que compara uma ou mais sequências fornecidas (sejam elas de proteínas ou nucleotídeos), contra uma base de dados e nos devolve uma lista de quais membros dessa base de dados seriam os mais parecidos com o que nós fornecemos.

## Pra que serve?

É um programa muito útil quando não sabemos o que nós temos em mãos e queremos descobrir, costuma ser eficiente e é uma boa escolha para tentar responder a pergunta que levantamos semana passada: afinal, o que é que nós sequenciamos?

## O que faremos

Aqui nós vamos tentar entender melhor o que é que nós acabamos sequênciando na aula passada. Para isso vamos usar o BLAST para descobrir quais sequências se parecem mais com a nossa amostra.

Nós temos o que se chama um arquivo FASTA na pasta "materiais".

Ele é composto por um cabeçalho que inicia com ">", o que vem depois na mesma linha é tratado como o nome da nossa sequência.
A nossa sequência em si está contida a partir da linha seguinte.

### Buscando no blast

Vamos acessar o site [Blast no NCBI](https://blast.ncbi.nlm.nih.gov/Blast.cgi) e pegar um conjunto de amostras que se pareçam com a nossa amostra misteriosa!

Uma vez no site, nós temos uma sequência de DNA, então nós vamos usar o blastn para esta atividade, portanto clicamos na opção Nucleotide>Nucleotide à esquerda da tela

1. Aqui temos uma janela para colocar a nossa sequência. Da maneira que preferirem, vocês podem fazer ctrl-c ctrl-v da sequência dentro da caixa ou clicar em "upload file" ("carregar arquivo")
2. Nós temos uma sequência de RNA 16S, para os próximos passos da disciplina, será útil sabermos em qual genoma ela está, portanto, vamos escolher a base refseq_reference_genomes
3. Vamos selecionar "blastn" dentre as opções de algorítmos que temos na seção "Program selection"
4. Vamos tentar pegar algumas sequências para termos o suficiente para gerar uma árvore, vá em "+Algorithms Parameters" e troque o valor de "max target sequences" de 100 para 500
5. Aperte BLAST


Uma vez com os resultados podemos observar algumas coisas na página, como os valores de:

Coverage (cobertura): O quanto da minha sequência que é contemplado pelo que eu tenho daquele item na minha base de dados
E-value: qual a probabilidade desse resultado se repetir
Per. ident (Porcentagem de identidade): O quão parecida a minha sequência é em relação a esse correspondente na base de dados

Vamos ver as nossas sequências?

Com todas selecionadas, clicamos em "Download" e selecionamos a opção "Download aligned sequences". 

Isso irá gerar um arquivo de extensão ".txt" para download

Se vocês abrirem, irão reparar que na verdade ele está em formato fasta, o que faz ele poder ser observado no programa **Aliview**


Vamos passar para o próximo passo, mas NÃO FECHE a saída do BLAST
