# Real lab - Bacillota — exercício de filogenia e evolução de caracteres

Este repositório contém material de apoio para a atividade a ser conduzida em sala de aula

## Atividade prática — Filogenia de 16S rRNA

Consulte o **[roteiro da atividade de filogenia](atividade/01_filogenia_16S.md)** para preparar as sequências de 16S rRNA com Chromas, Geneious e BLASTn, construir um alinhamento com MAFFT, inspecionar ou filtrar com AliView/trimAl e inferir uma árvore com IQ-TREE2 para visualização no iTOL.

A atividade utiliza o servidor Haldane e contas individuais dos alunos. O conjunto de sequências será montado durante a prática.

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
O arquivo `environment.yml` é uma proposta de dependências

## Notebooks (`atividade/`)

### Árvore com o isolado

## Organização
