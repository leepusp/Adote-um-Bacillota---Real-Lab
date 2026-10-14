# Real lab - Bacillota — exercício de filogenia e evolução de caracteres

Este repositório contém material de apoio para a atividade a ser conduzida em sala de aula


## Se você quiser rodar o conteúdo dos notebooks em casa


### Instalação (Linux/macOS; Windows preferencialmente via WSL)

Você irá precisar de um gerenciador de pacotes (você pode instalar tudo de maneira independente, mas eles estão ali para facilitar sua vida), atualmente os mais conhecidos são o `mamba` e o `conda`

Com `mamba` (ou `conda`) instalado, a partir da raiz do repositório:

```bash
mamba env create -f environment.yml 
conda activate bacillota-aula
jupyter lab
```

Se utilizar `conda`, substitua `mamba env create` por `conda env create`.
O arquivo `environment.yml` é uma proposta de dependências; a resolução das versões e a execução dos scripts **ainda precisam ser validadas na máquina local com R**.

## Notebooks (`atividade/`)

1. `01_Coletando_os_16S.ipynb`: junta o 16S do isolado e as referências do grupo às 229 sequências, alinha (MAFFT) e apara (trimAl), com um visualizador de alinhamentos que também abre os pré-computados; termina com as instruções do IQ-TREE.
2. `02_Coletando_proteínas_ribossomais.ipynb`: cobertura dos 47 marcadores ribossomais, efeito do trimAl, visualizador por marcador e concatenação dos 46 usados.
3. `03_Arvores_e_Fritz_D_interativo.ipynb`: árvores anotadas (16S, 16S com o isolado, ribossomais e GTDB), D de Fritz para genes e fenótipo de esporulação, reconstrução do estado ancestral da esporulação e, como extra, do envelope (monodermes/didermes).

O visualizador de alinhamentos fica em `scripts/alinhamentos.py` e é usado pelos notebooks 01 e 02. Ele lista os FASTAs alinhados de `dados/` (pré-computados) e de `resultados/` (gerados em aula).

### Árvore com o isolado

O notebook 01 gera `resultados/16S_com_isolado.fna`, `16S_com_isolado.aln` e `16S_com_isolado.trim.aln`. A árvore do IQ-TREE deve ser salva como `16S_isolado.tree.treefile`, em `dados/` ou `resultados/`. Se algo acontecer, o comando equivalente no terminal é:

```bash
iqtree -s 16S_com_isolado.trim.aln -m MFP -B 1000 -T AUTO --prefix 16S_isolado.tree
```

Para que o alinhamento pré-computado com o isolado apareça no visualizador, coloque-o em `dados/` com extensão `.aln` (por exemplo `dados/16S_isolado.aln`).

Abra o JupyterLab na **raiz do projeto**. As células localizam automaticamente `dados/` a partir do diretório corrente, inclusive se o notebook for aberto dentro de `atividade/`.

### Dados essenciais

- `dados/tree_16S.fna`, `.aln` e `.trim.aln`: 16S original, alinhado e trimado.
- `dados/selribs_df.pkl`, `trimribs_df.pkl` e `ribconcat_df.pkl`: dados ribossomais prontos para o notebook 02.
- `dados/16S_final.tree.treefile`, `ribomarkers_final.tree.treefile` e `bac120_229.tree`: árvores Newick. A última é a GTDB **podada aos mesmos 229 genomas**, para ser leve.
- `dados/tip_meta.tsv`: identificadores e taxonomia
- `dados/gene_presence_absence_final_61.tsv`: 61 genes com 0/1 para 229 assemblies.
- `dados/heatmap_matrix2_spore_categories.tsv`: metadados históricos.
- `dados/additional_gene_calls.tsv`: **fonte exclusiva das presenças adicionais inferidas via orthogroups** (cor amarela no heatmap).

`Spore former` representa a classificação fenotípica de Galperin (**1 esporulante, 2 não esporulante, 0 desconhecido**); `expected` representa a categoria de envelope com base na literatura. Nenhum dos dois é inferido pela matriz híbrida. Ausência de gene aqui significa não detectado pelos critérios usados.

O script R grava PDFs/PNGs e tabelas TSV em `resultados/`, ignorados pelo Git.

### Script R (`scripts/phylogeny.r`)

Os notebooks chamam o script, mas ele também pode ser usado diretamente a partir da raiz do projeto:

```bash
Rscript scripts/phylogeny.r plot      --tree gtdb
Rscript scripts/phylogeny.r fritz     --tree ribo --marker spo0A --permut 999 --scope bacillota
Rscript scripts/phylogeny.r compare   --marker lpxC --scope bacillota
Rscript scripts/phylogeny.r ancestral --tree ribo --trait spore --model auto
Rscript scripts/phylogeny.r plot      --tree isolado
Rscript scripts/phylogeny.r --help
```

- **Enraizamento:** todas as árvores são enraizadas com os genomas fora de Bacillota como outgroup. Se eles não formarem um grupo numa árvore, o script usa midpoint e avisa; a raiz usada aparece no título das figuras e na coluna `rooting` das tabelas.
- **Pontas sem metadados** (o isolado e referências acrescentadas) são mantidas na figura (losango preto, nome em caixa amarela) e ficam de fora do D e da reconstrução ancestral.
- **Caracteres:** os 61 marcadores, `spore` (fenótipo de Galperin: 1 esporulante, 2 não esporulante; 0 = desconhecido, tratado como ausente de informação) e `expected` (envelope, só na reconstrução ancestral).
- **`fritz`:** `--marker` aceita um gene ou `spore`. `--scope bacillota` (padrão) calcula o D só com Bacillota; `--scope all` inclui os outgroups.
- **Figura:** pontos nos nós marcam UFBoot ≥ 95 nas árvores do IQ-TREE e bootstrap ≥ 85 na GTDB (`SUPPORT_MIN` no início do script). Bacillales é sombreado; se não for monofilético, cada clado só de Bacillales recebe sombreado próprio e o script lista as ordens que interrompem o grupo.
- **`compare`:** calcula o D nas três árvores de referência (16S, ribo, GTDB) restritas aos assemblies presentes em todas elas.
- **`ancestral`:** reconstrução marginal com modelo Mk (`ape::ace`), ajustando ER e ARD e escolhendo pelo AIC; a tabela traz uma linha por modelo, com taxas e probabilidades no ancestral de Bacillota e na raiz. `--trait` aceita `spore` (padrão), `expected` ou qualquer um dos 61 marcadores. Como a classificação de Galperin cobre só Firmicutes, a reconstrução de `spore` usa apenas Bacillota.

## Organização

`atividade/`: notebooks dos estudantes. `dados/`: entradas necessárias às três atividades. `scripts/`: análise e visualização em R. `histórico/`: notebooks da seleção original, preservados para referência, **não necessários para executar a atividade**. `preparacao/`: notebook enxuto de curadoria que permite ver por cima como a matriz de presença e ausência foi construída.

Por limitações do git, a tabela das anotações e dos ortólogos necessárias para rodar o notebook não estão disponíveis, mas podem ser solicitadas

## Observações e limites

- A árvore GTDB contém inicialmente 189.801 pontas, podadas para 229 na versão para a aula. Usar a subárvore não equivale a inferir uma árvore nova.
- O destaque de Bacillales depende da taxonomia recuperada e pode ser incompleto para genomas não cobertos; evite interpretar a falta de destaque como evidência taxonômica negativa.
- O notebook 01 **não reproduz a extração Bakta**, pois ela requer os arquivos `.ffn` originais ausentes do repositório. O FASTA das 229 sequências já é fornecido.
- O script R e os notebooks foram testados com dados sintéticos (R 4.3.3, ape 5.7.1, caper 1.0.4, ggtree 4.3.1, ggplot2 4.0.3, ggnewscale 0.5.2, optparse 1.8.2; Python 3.12, Biopython 1.88, pandas 3.0, ipywidgets 8), não com os dados reais nem com o ambiente criado pelo `environment.yml`. Rode os comandos acima com os dados reais e teste `mamba env create` numa máquina limpa antes da aula.
- R usa as mesmas 229 assemblies, mas compara topologias que podem diferir entre os marcadores; D depende da topologia, comprimento de ramos e definição de caracteres.
- A reconstrução ancestral depende do enraizamento, da amostragem de outgroups e do modelo. O modelo Mk trata ganho e perda de membrana externa como eventos comparáveis, o que é uma simplificação forte.
- `fabZ` está presente em praticamente todas as bactérias, monodermes inclusive, e não discrimina envelope; vale tratá-lo como contraexemplo na discussão (pergunta 4 do notebook 03).
