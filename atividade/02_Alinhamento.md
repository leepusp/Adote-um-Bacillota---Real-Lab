# Preparando o alinhamento 

## Montando o arquivo

Ao abrir o Aliview vamos observar a nossa lista de sequências não alinhadas.

Dando zoom-out (ctrl + rodinha do mouse) vamos observar algumas "bandas" que parecem meio estranhas em relação ao resto.

Aqui percebemos que infelizmente o BLAST não nos deu todas as sequências na mesma orientação. 

Se abrirmos a saída do BLAST do passo anterior e formos na aba de alinhamentos, vamos observar que as sequências estranhas possuem a informação "Plus/Minus"

> O que aconteceu?
> Aqui vemos que as sequências foram baixadas em **orientações diferentes**

Isso *vai* ocasionar problemas na hora de gerar um bom alinhamento.

Uma das opções seria ir uma por uma das sequências que estão invertidas e pedir para o Aliview gerar o "Reverse Complement"

A outra opção é, quando formos alinhar, deixarmos o programa de alinhamento decidir como as sequências deveriam estar orientadas (Vamos escolher essa aqui!)


Nós temos um outro problema! 
O BLAST nos deu **TODOS** os hits do 16S no genoma. Faz sentido, um genoma pode ter mais de uma cópia de 16S, mas aqui queremos apenas 1 hit por genoma


Vamos dar um pulinho no google sheets!

## Filtrando as sequências do blast usando Google Sheets

Importe a tabela "[Código]-Alignment-HitTable.csv" no google sheets

Ela está ordenada de maneira tal a que os melhores resultados aparecem primeiro

Por exemplo:

```
Bacillus_USP_ICB_16s_partial_sequence,NZ_BOQB01000033.1,99.181,733,3,2,1,733,339,1068,0.0,1293
Bacillus_USP_ICB_16s_partial_sequence,NZ_CP040336.1,99.181,733,3,2,1,733,486692,487421,0.0,1293
Bacillus_USP_ICB_16s_partial_sequence,NZ_CP040336.1,99.181,733,3,2,1,733,4390250,4389521,0.0,1293
Bacillus_USP_ICB_16s_partial_sequence,NZ_CP040336.1,99.181,733,3,2,1,733,4593208,4592479,0.0,1293
Bacillus_USP_ICB_16s_partial_sequence,NZ_CP040336.1,99.181,733,3,2,1,733,4819468,4818739,0.0,1293
Bacillus_USP_ICB_16s_partial_sequence,NZ_CP040336.1,99.181,733,3,2,1,733,4852855,4852126,0.0,1293
```
Nós estamos interessados na segunda coluna e queremos selecionar apenas a primeira ocorrência, logo, vamos fazer a seguinte operação:

1. Inserir uma linha acima do primeiro resultado:
2. Na coluna vazia **M** vamos colar a seguinte fórmula na segunda linha:

```
=B2<>B1
```

Ela irá retornar o valor TRUE, nós queremos expandir essa fórmula para todas as linhas de baixo, logo, vamos clicar no canto esquerdo onde está a marcação das colunas, digitar M1:M (seleciona todas as linhas), dar Enter e apertar Ctrl+D

3. Também precisamos de parte do nome que as nossas sequências têm no arquivo fasta. Se vocês repararem, os nomes das nossas sequências sempre começam com: "[IDENTIFICADOR]:[COORDENADAS]", vamos usar isso a nosso favor!

Na coluna N cole a seguinte fórmula:

```
=B1&":"&MIN(I1;J1)&"-"&MAX(I1;J1)
```
E fazemos o mesmo procedimento de expandir essa fórmula para as outras linhas

Ela deve preencher as colunas com este tipo de dado:

```
NZ_BOQB01000033.1:339-1068
NZ_CP040336.1:486692-487421
NZ_CP064875.1:151134-151863
NZ_NWUW01000108.1:192-921
NZ_CM166465.1:28970-29699
```

Agora estamos prontos! 

4. Selecione a tabela e no menú de "Dados" selecione Criar filtro

A primeira linha agora deveria estar com um ícone de seta virada para baixo

Na coluna M clique no ícone e desmarque a opção "FALSE"

5. Selecione a coluna N inteira e aperte "Ctrl+C"

### Vamos voltar para o Aliview

No Aliview, vá no menú "Edit" e encontre a opção "Select sequences from clipboard"

Você verá que apenas algumas sequências do arquivo foram selecionadas, aperte ctrl+c novamente e vá em File>New

Uma nova janela se abrirá, aperte ctrl+V

Temos as nossas sequências do Blast!

Falta apenas a nossa.

Vá em Edit> Add sequences from file e selecione a nossa sequência


## Gerando o alinhamento

Por motivos de simplicidade, vamos utilizar a função do mafft para cuidar das sequências que vieram invertidas

[Servidor do MAFFT](https://mafft.cbrc.jp/alignment/server/index.html)

Ao clicar no link acima abrimos a página do servidor do mafft

Da maneira que preferirem, adicionem o nosso arquivo fasta na janela de "input" do servidor. Ou vocês selecionam todo o conteúdo do arquivo e colam dentro da caixa ou vocês fazem upload do arquivo fasta.

Agora, para resolver o problema das nossas sequências não estarem na mesma orientação, nós vamos marcar **Adjust direction according to the first sequence (accurate enough for most cases)** na seção de **Direction of nucleotide sequences**. 

Certifique-se de fazer isso!

## Inspecionando - Parte 2

Agora, nós vamos fazer alguns passos para garantir que não haja problemas:

1. Baixe o fasta fornecido pelo MAFFT ao clicar em "Fasta Format" na tela
2. Abra o alinhamento no Aliview
3. Vá em "File" e em "Save as Nexus (Ilegal names padded)"
4. Vá em "File" novamente, mas salve desta vez em fasta
5. Dê uma olhada no alinhamento


Uma vez que o alinhamento está pronto, nós vamos baixar o alinhamento do site.
Sempre é recomendado dar uma olhada no alinhamento para ver se não há nenhuma surpresa

Tendo esse alinhamento feito, está na hora de gerarmos a árvore!


