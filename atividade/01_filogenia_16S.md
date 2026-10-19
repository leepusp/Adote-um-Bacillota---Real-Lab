## Objetivo

Olá pessoal, nós vamos montar uma árvore de 16S com o isolado que tivemos  posicionado dentro de um contexto mais amplo de diversidade bacteriana: o próprio gênero Bacillus e parentes próximos (via BLAST), mais dois grupos de referência adicionados manualmente (Clostridium e Escherichia coli/Gram-negativas) para servir de "régua" de comparação.

>Resultado esperado: uma árvore com aproximadamente 700 sequências (500 do BLAST do isolado + 100 de Clostridium + 100 de E. coli + o nosso isolado ambiental).

## Passo 1 — Visualização de cromatogramas no Chromas e montagem do gene no GENEIOUS PRIME

O 1° passo é abir o arquivo (Ab1) no chromas, onde visualizaremos o cromatograma inicial, e iremos interpreta-lo, e avaliando a qualidade do sequenciamento Sanger gerado. No próprio Chromas, temos a opção de copiar a sequência lida (porém ainda não montada), ela devera ser utilizada na identificação inicial via BLASTn.

Após a identificação via BLASTn [Gustavo, você precisa falar como vamos buscar no blastn], utilizamos a 1° sequência de referência para montagem do gene. 

1. Para isso migraremos para o Geneious, após copiar o .FASTA da sequência no programa, e nomear a referência, iremos utilizar a opção "Map to reference" no geneious, ao selecionar os dois cromatogramas, e mapear para a referência selecionada.

2. Após isso, o próprio programa irá fazer o mapeamento da referência, e irá assinar "qual nucleotídeo foi lido", ainda assim, precisamos de uma inspeção manual para averiguar se as bases foram assinaladas corretamente, em alguns casos, podemos corrigir o resultado "na mão", mudando o nucleotído lido, a partir do cromatograma visualizado.

3. Depois de montar as duas reads (Forward/Reverse) no Geneious e gerar o consenso, exportem a sequência consenso em formato FASTA. Essa é a nossa sequência de base que vocês vão usar no BLAST a seguir.

## Passo 2 — BLAST do consenso (top 500)

1. Acessem o BLASTn do NCBI (https://blast.ncbi.nlm.nih.gov/Blast.cgi?PROGRAM=blastn&PAGE_TYPE=BlastSearch&LINK_LOC=blasthome).

2. Colem a sequência consenso do isolado de vocês.

2. Na opção de DATABASE, selecionem "rRNA/ITS databases --> 16S ribosomal RNA sequences (Bacteria and Archaea type strains)" — isso é importante: evita trazer sequências de clones ambientais sem identificação ("uncultured bacterium")
3. Atualizem o "Program selection" para blastn "somewhat similar sequences"
4. Alterem o "+algorithm parameters", e mudem der 100 para 500 o número de respostas
5. Rodem o BLAST.

Baixem os 500 primeiros hits em formato FASTA (botão "Download" → "FASTA (aligned sequences)").

## Passo 3 — Adicionar sequências para comparação

Como o BLAST retorna todos os hits (incluindo duplicatas e variações da mesma espécie), vamos adicionar uma 2° busca para selecionar microrganismos um pouco mais distantes, incluindo as bactérias do gênero Clostridium (também produzem endósporos), mas também um grupo de bactérias gram negativas

1. No NCBI Nucleotide, busquem por uma sequência de referência de 16S de Clostridium botulinum.
2. Rodem um novo BLASTn com essa sequência (mesmo banco: 16S type strains).
3. Baixem os 100 primeiros hits em FASTA.
4. Abram o alinhamento/dataset principal no AliView e usem Sequence → Add Sequences para importar esse FASTA novo pro dataset.
5. Passo 5 — Adicionando grupo de comparação: E. coli (Gram-negativas)
6. Repitam o mesmo processo com uma sequência de referência de 16S de Escherichia coli.
7. Rodem o BLASTn (mesmo banco).

==Atenção: é esperado que apareçam outras Gram-negativas próximas no resultado, como Shigella e Salmonella, isso é normal e faz parte da diversidade do grupo, não se preocupem.==

Baixem os 100 primeiros hits e adicionem ao dataset no AliView, do mesmo jeito do Passo 4.

Ao final deste passo, o dataset deve ter aproximadamente 700 sequências.

# Passo 4 — Alinhamento (MAFFT, servidor Haldane via LEEP)

Vamos usar o modo mais rápido do MAFFT (FFT-NS-1 — sem iteração), suficiente pra esse propósito didático.
Temos que tomar o cuidado das 

## Conectem na Haldane via terminal nas máquinas da LEEP: (CHECAR)

bash
ssh usuario@haldane 

Rodem o alinhamento:

bash
mafft --retree 1 --maxiterate 0 dataset_700.fasta > dataset_700_aligned.fasta ((CHECAR) CAMINHO GENÉRICO, TESTAR)
--retree 1: constrói a árvore-guia uma única vez (sem refazer)
--maxiterate 0: sem rodadas de refinamento iterativo

Isso deve rodar em poucos minutos mesmo com ~700 sequências.

## Passo 5 — Trimagem no AliView
Abram o arquivo dataset_700_aligned.fasta no AliView.
Inspecionem visualmente as pontas do alinhamento (regiões com muitos gaps/baixa qualidade nas extremidades).
Salvem como um novo arquivo FASTA: dataset_700_trimmed.fasta.

# Passo 6 — Árvore (IQ-TREE2, servidor Haldane, modelo JC)

bash
iqtree2 -s dataset_700_trimmed.fasta -m JC -nt AUTO -pre arvore_700 ((CHECAR)CAMINHO GENÉRICO, TESTAR)
-m JC: modelo de Jukes-Cantor (o mais simples — todas as substituições com a mesma probabilidade). É uma simplificação proposital para priorizar velocidade e foco no fluxo de ferramentas; não é o modelo que se usaria numa análise filogenética "para publicação" (ali se usaria seleção de modelo com ModelFinder, geralmente GTR+G).
O resultado final de interesse é o arquivo arvore_700.treefile.


# Passo 7 — Subindo no iTOL

Acessem o iTOL e façam upload do arquivo .treefile.
A partir dessa etapa, começa a parte de interpretação que vocês irão fazer, como vocês vão "marcar" cada árvore e cada ramo
Como vocês terão 4 árvores iguais, a idéia é que nessa etapa vocês formem uma dupla, que irá analizar a mesma árvore, porém um de vocês será responsável por buscar informações a respeito daquele grupo de bactérias sendo analisadas)

1 - Pintar ramo (Selecionar Branches --> Whole Clade --> Selecionar cor) por patogenicidade (bactéria reconhecidamente causa doença em humana --> pintar de vermelho ou outra cor que preferirem; Não causa doença em humana e/ou desconhecido, vocês pintam de azul, por exemplo.)
2 - Colorir ramo por hospedeiro provável (bactéria comensal humana x bactéria ambiental)

Como são muitos ramos, a ideia é que vocês olhem os "grupos" ou também chamados de subclados dentro da árvore, ao invés de analisarem uma a uma, não é o objetivo da atividade.

## Passo 8 Agora que vocês aprenderam como montar um árvore inicial, vamos olhar para um exemplo de uma pesquisa já públicada e vamos ver o padrão de algumas das sequências nossas!

1. Vamos abrir a tabela fornecida pelos autores do artigo [Conservation and Evolution of the Sporulation Gene Set in Diverse Members of the Firmicutes](https://pmc.ncbi.nlm.nih.gov/articles/PMC9210971/) disponível na pasta de materiais, 
2. Essa tabela é grande, para motivos desta disciplina, nós vamos focar nos 112 genes caracterizados como core para as próximas questões, a tabela reduzida também na parte de materiais
   
   
### Encontrando anotações

Infelizmente, não podemos integrar o nosso bacilo ambiental pois não possuímos o genoma dele, aqui cada um de vocês irá escolher uma espécie do resultado do BLAST para tentarmos encontrar os equivalentes desses 112 genes.

#### Precisamos de um genoma

É necessário term


## "" MALU PREENCHER AQUI PASSO A PASSO DAS ATIVIDADES "" ##


# Perguntas guia para a análise

1) Olhando a árvore de vocês em relação à diversidade geral de procariotos (material de apoio com as árvores dos livros): como os grupos próximos do isolado de vocês se posicionam dentro dessa diversidade maior? ##-> isso aqui eles não vão conseguir responder com a árvore de 16S

2) Em relação às sequências mais próximas do isolado: seu Bacillus fica perto de quem? Quais espécies/gêneros aparecem no entorno imediato?
3) Olhando para os grupos de referência que vocês mesmos adicionaram (Clostridium, E. coli e parentes), eles aparecem na árvore de vocês na mesma relação que aparecem no material de apoio (árvore dos Gram-positivos)? Tem alguma coisa que ficou diferente do esperado? Existem mais espécie "entre elas" que não foram representadas" ?
4) Depois de pintar por patogenicidade e por hospedeiro provável: esses dois padrões de cor se agrupam dentro dos mesmos subclados da árvore, ou ficam misturados? Um clado filogeneticamente próximo necessariamente significa bactérias com o mesmo estilo de vida (ambiental vs. comensal, patogênica vs. não)? ## isso aqui, você tem os dados, Gustavo? É um set de informações diferente


## MALU, COLOCAR PERGUNTAR SOBRE A ATIVIDADE DE ENDÓSPOROS
