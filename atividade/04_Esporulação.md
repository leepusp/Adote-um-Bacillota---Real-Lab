# Análise da distribuição de genes de esporulação em Bacillota

## Explica o Galperin

## Tabela do Gaperin

## Selecionando genomas com BLASTN do gene 16S rRNA

## Usando TBLASTN para identificar genes de esporulação

## Usando Google Sheets para filtrar dados

Vamos selecionar os melhores resultados do tblastn usando o Google Sheets!

Precisamos dar nomes às colunas da planilha que baixamos do BLAST/NCBI.

1. Copie a lista de nomes de colunas abaixo:

```
qaccver,saccver,pident,length,mismatch,gapopen,qstart,qend,sstart,send,evalue,bitscore,positives,qframe,frame
```

2. Insira uma nova primeira linha na planilha e cole a lista copiada na primeira célula.
3. Clique no menu Dados e selecione "Dividir texto em colunas"
4. Certifique-se de que o separador de colunas escolhido foi a vírgula (",")
5. Copie os dados da coluna A, **linhas 1 até a penúltima**, para a coluna P e mude seu nome para _anterior_
6. Na segunda linha da coluna Q, coloque a fórmula ```=A2<>P2```
7. Copie a fórmula em Q2 para todas as outras linhas da coluna Q onde existem valores na linha
8. Clique no botão **filtrar** () e desative o item FALSO
9. Copie as linhas filtradas para outra aba e...
