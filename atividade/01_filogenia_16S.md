# Atividade prática — Filogenia de 16S rRNA

## Objetivo

Vamos construir uma árvore filogenética de 16S rRNA para posicionar o isolado ambiental em um contexto mais amplo de diversidade bacteriana: *Bacillus* e parentes próximos (via BLAST), além de grupos de comparação, como *Clostridium* e *Escherichia coli*/outras Gram-negativas.

> **Resultado esperado:** uma árvore com aproximadamente 650 sequências (500 do BLAST do isolado + 50 de *Clostridium* + 100 de *E. coli* + o nosso isolado ambiental). O total poderá variar devido a resultados repetidos entre buscas.

## Passo 1 — Visualização de cromatogramas no Chromas e montagem do gene no GENEIOUS PRIME

Primeiro, abram os arquivos de sequenciamento Sanger (`.ab1`) no **Chromas/Geneious**. Visualizem os cromatogramas, avaliem a qualidade das bases e copiem uma sequência inicial para uma identificação preliminar por BLASTn.

Após encontrar uma sequência de referência de 16S rRNA, utilizem o **Geneious Prime**:

1. Importem a sequência de referência em **FASTA** e os dois cromatogramas (**Forward/Reverse**).
2. Selecionem os cromatogramas e utilizem **Map to Reference** para mapeá-los à referência.
3. Inspecionem os cromatogramas manualmente e corrijam chamadas de nucleotídeos duvidosas quando houver evidência para isso.
4. Gerem a sequência consenso das leituras e exportem-na em **formato FASTA**.

O consenso é a sequência de base para a busca principal no BLASTn.

## Passo 2 — BLAST do consenso (top 500)

1. Acessem o [BLASTn do NCBI](https://blast.ncbi.nlm.nih.gov/Blast.cgi?PROGRAM=blastn&PAGE_TYPE=BlastSearch&LINK_LOC=blasthome).
2. Colem a sequência consenso do isolado.
3. Em **Database**, selecionem **rRNA/ITS databases → 16S ribosomal RNA sequences (Bacteria and Archaea type strains)**. Esse banco favorece sequências de referência identificadas, reduzindo resultados como *uncultured bacterium*.
4. Em **Program Selection**, selecionem **Somewhat similar sequences (blastn)**.
5. Em **Algorithm parameters**, ajustem o número máximo de sequências retornadas para **500**.
6. Executem o BLASTn e baixem os resultados em **formato FASTA**, preferencialmente a opção **Hit sequences**, mantendo os cabeçalhos com as descrições dos organismos. Se baixarem sequências alinhadas com gaps, eles serão removidos antes do MAFFT.

**Entrada:** consenso em sequência de DNA. **Saída:** arquivo FASTA com os hits de 16S rRNA.

## Passo 3 — Adicionar sequências para comparação

Como o BLAST do isolado retorna principalmente organismos próximos, adicionaremos dois grupos de comparação mais distantes. Algumas sequências podem se repetir entre as buscas. Para isso, vamos voltar na página do BlastN, e clicar no NIH, para acessar ao banco de dados do NCBI.

**Grupo *Clostridium*:**

1. No NCBI Nucleotide, procurem uma sequência de referência de 16S rRNA de *Clostridium botulinum*.
2. Executem um novo BLASTn no mesmo banco de 16S de linhagens-tipo.
3. Baixem os **50 primeiros hits em FASTA**.
4. Abram o conjunto principal no **AliView** e utilizem **Sequence → Add Sequences** para importar esse FASTA.

**Grupo *E. coli* e outras Gram-negativas:**

1. Procurem uma sequência de referência de 16S rRNA de *Escherichia coli*.
2. Repitam a busca no mesmo banco, baixando os **100 primeiros hits em FASTA**.
3. Adicionem as sequências ao conjunto no AliView pelo mesmo procedimento.

**Atenção:** podem aparecer outros gêneros relacionados, como *Shigella* e *Salmonella*. Isso faz parte da diversidade do grupo; não os excluam apenas pelo nome do gênero.

Incluam também o **consenso do isolado ambiental** e salvem o conjunto final como **`database.fasta`**. O conjunto deverá ter aproximadamente 651 sequências; verifiquem eventuais entradas duplicadas.

**Entrada:** arquivos FASTA dos três BLASTs e do isolado. **Saída:** `database.fasta` — FASTA com as sequências reunidas.

## Passo 4 — Preparação do FASTA e alinhamento com MAFFT

### 4.1. Acessar o servidor Haldane

No terminal da LEEP, conectem-se usando a conta fornecida para a atividade:

```bash
ssh -p 25022 user@haldane.icb.usp.br
```

Substituam `user` por `aluno1`, `aluno2`, `aluno3` ou `aluno4`. Depois de conectar:

```bash
conda activate py310
cd ~/adote
```

O arquivo `database.fasta`, montado no AliView, deverá estar nessa pasta do Haldane. Se ele estiver somente no computador da LEEP (e não em uma pasta compartilhada), enviem-no a partir do terminal da LEEP com `scp -P 25022 database.fasta user@haldane.icb.usp.br:~/adote/`, substituindo `user` pelo login individual.

### 4.2. Padronizar os identificadores GenBank

Os nomes das bactérias precisam estar nos cabeçalhos do FASTA **antes do IQ-TREE2**, para aparecerem nas pontas da árvore. O comando abaixo mantém os acessos GenBank, simplifica os nomes e remove gaps importados de alinhamentos anteriores.

```bash
awk '
/^>/ {
    id = $1
    sub(/^>/, "", id)
    gsub(/[^A-Za-z0-9._-]/, "_", id)
    nome = $0
    sub(/^>[^[:space:]]+[[:space:]]*/, "", nome)
    sub(/[[:space:]]+16S[[:space:]]+(ribosomal[[:space:]]+RNA|rRNA).*$/, "", nome)
    gsub(/[^A-Za-z0-9._-]+/, "_", nome)
    gsub(/_+/, "_", nome)
    sub(/^_+/, "", nome)
    sub(/_+$/, "", nome)
    label = (nome == "") ? id : nome "__" id
    copia[label]++
    if (copia[label] > 1) label = label "__copy" copia[label]
    print ">" label
    next
}
{
    gsub(/[-.]/, "")
    gsub(/[[:space:]]/, "")
    if (length($0)) print
}
' database.fasta > bacillus_clostridium_ecoli.fasta
```

**Entrada:** `database.fasta` — FASTA com acesso GenBank seguido da descrição do organismo. **Saída:** `bacillus_clostridium_ecoli.fasta` — FASTA sem gaps e com nomes padronizados.

Exemplo de identificador esperado: `Clostridium_sporogenes_strain_JCM_1416__NR_113245.1`.

**Observação:** se um arquivo do BLAST contiver apenas IDs e não as descrições dos organismos, será necessário obter os nomes nos registros do NCBI antes dessa etapa. Cabeçalhos repetidos recebem o sufixo `__copy2`, mas isso não remove sequências duplicadas.

Confiram as contagens (elas devem ser iguais):

```bash
grep -c '^>' database.fasta
grep -c '^>' bacillus_clostridium_ecoli.fasta
```

### 4.3. MAFFT — alinhamento múltiplo

```bash
mafft --nuc --auto --thread 2 bacillus_clostridium_ecoli.fasta > mafft.fasta
```

**Entrada:** `bacillus_clostridium_ecoli.fasta` — FASTA de DNA não alinhado. **Saída:** `mafft.fasta` — FASTA de DNA alinhado.

O **MAFFT** identifica posições correspondentes entre as sequências, inserindo gaps. `--auto` escolhe uma estratégia e `--thread 2` utiliza duas threads.

## Passo 5 — Inspeção do alinhamento no AliView

Abram **`mafft.fasta`** no AliView disponível nos computadores da LEEP. Observem regiões conservadas, gaps e possíveis problemas no alinhamento. Se fizerem uma filtragem manual, salvem uma cópia como **`manual.fasta`**, preservando `mafft.fasta`. O arquivo final deve continuar sendo um **FASTA alinhado**, com as mesmas sequências e o mesmo comprimento de alinhamento para todas elas.

### trimAl — filtragem automática (opcional)

Como alternativa à edição no AliView, podem executar:

```bash
trimal -in mafft.fasta -out trimal.fasta -automated1
```

**Entrada:** `mafft.fasta` — FASTA alinhado. **Saída:** `trimal.fasta` — FASTA alinhado com colunas filtradas.

O **trimAl** remove automaticamente posições do alinhamento segundo critérios de qualidade. Esta etapa é opcional.

## Passo 6 — Construção da árvore com IQ-TREE2

Neste exemplo, usaremos o modelo evolutivo **GTR** e o alinhamento filtrado pelo trimAl:

```bash
MODELO=GTR
mkdir -p "iqtree/$MODELO"
iqtree2 -s trimal.fasta -st DNA -m "$MODELO" -B 1000 -T AUTO --threads-max 8 -mem 6G -pre "iqtree/$MODELO/arvore"
```

**Entrada:** `trimal.fasta` — FASTA com alinhamento múltiplo de DNA. **Saída principal:** `iqtree/GTR/arvore.treefile` — árvore filogenética em **formato Newick**.

O **IQ-TREE2** infere uma árvore por máxima verossimilhança. `-m GTR` define o modelo evolutivo; `--fast` acelera a busca; `-T AUTO --threads-max 8` ajusta o paralelismo até oito threads; `-pre` define o local dos resultados.

**Se não utilizarem o trimAl**, substituam `-s trimal.fasta` por **`-s manual.fasta`** (edição no AliView) ou **`-s mafft.fasta`** (sem filtragem). Para repetir a análise com outro arquivo de entrada, escolham outra pasta de saída para não sobrescrever os resultados. Mantenham a sessão SSH conectada até o fim da execução.

**Importante:** esta busca rápida não calcula suportes UFBoot nem SH-aLRT. Por isso, a árvore é destinada à exploração didática; a análise filogenética com suporte estatístico requer uma execução adicional.

## Passo 7 — Visualização e interpretação no iTOL

Após terminar o IQ-TREE2, localizem **`iqtree/GTR/arvore.treefile`**. Abram esse arquivo, em **formato Newick**, no [iTOL](https://itol.embl.de/).

- **`.treefile` (Newick):** árvore para visualização, com nomes das bactérias e comprimentos dos ramos.
- **`.iqtree`:** relatório da inferência filogenética.

No iTOL, localizem os principais clados, selecionem **Branches → Whole clade** e escolham cores para representar características biológicas. Como são muitas sequências, concentrem-se nos grupos e subclados, em vez de pintar individualmente todas as pontas.

Como atividade de interpretação, comparem os agrupamentos por bactérias esporulantes e façam a coloração dos ramos passa assinalar bactérias que esporulam e que não esporulam, usem uma cor para cada, grupo.

Por fim, vamos colorir as bactérias de outra forma, assinalem bactérias comensais humanas (Staphylococcus, Escherichia coli, e outras bem caracterizadas) e bactérias ambientais, diferenciem por gênero. 

Organizem-se em duplas para investigar a mesma árvore e discutir diferenças de anotação e interpretação dos clados. Caso vocês não saibam a origem de um grupo específico, organizem-se para pesquisar essas informações na literatura.

---

## Passo 8 - Análise dos genes de esporulação

Agora que vocês aprenderam como montar um árvore inicial, vamos olhar para um exemplo de uma pesquisa já públicada e vamos ver o padrão de algumas das sequências nossas!

1. Vamos abrir a tabela fornecida pelos autores do artigo [Conservation and Evolution of the Sporulation Gene Set in Diverse Members of the Firmicutes](https://pmc.ncbi.nlm.nih.gov/articles/PMC9210971/) disponível na pasta de materiais, 
2. Essa tabela é grande, para motivos desta disciplina, nós vamos focar nos 112 genes caracterizados como core para as próximas questões, a tabela reduzida também na parte de materiais (Galperin_reduced.xlsx)
   
   
### Encontrando anotações

Infelizmente, não podemos integrar o nosso bacilo ambiental pois não possuímos o genoma dele, aqui cada um de vocês irá escolher uma espécie do resultado do BLAST para tentarmos encontrar os equivalentes desses 112 genes.

#### Precisamos de um genoma

É necessário termos um genoma para ver se algum dos genomas que nós podemos recuperar através do blast possui esse cog
Vamos primeiro selecionar o genoma

Da saída do BLAST que conseguimos anteriormente, vamos ver quais sequências podem servir para esse objetivo. Estamos em 4. Vamos dividir um genoma para cada

Os top 4 resultados no BLAST feito para elaborar a aula são:

- Bacillus wiedmannii strain FSL W8-0169 
- Bacillus proteolyticus strain MCCC 1A00365
- Bacillus toyonensis strain BCT-7112
- Bacillus sanguinis strain BML-BC004

Vamos tentar encontrar um genoma para essas espécies e strain!

##### Encontrando um genoma

1. Acesse o site da base de dados do ncbi [Databses](https://www.ncbi.nlm.nih.gov/datasets/genome/)
2. Digite primeiro o nome da espécie desejada no campo de busca
3. Filtre pelas palavras chaves após "strain" (exemplo: FSL W8-0169)
   1. Se não existir o genoma da strain, vamos pegar o genoma da espécie!
4. Clique no link para o assembly, role a página que você encontrou até a região de "Chromossomes"
   1. Aqui, se existir apenas um ou dois, copie os identificadores do nucleotídeo, caso existam mais, recomendamos que você aperte o botão de Download e salve a tabela

##### Encontrando o nosso gene

Entremos novamente no BLAST, mas desta vez, iremos selecionar tblastn

Carregue o fasta galperin_cog_queries.faa

Marque "Align with one or more sequences"

No campo que se abriu, cole os identificadores de núcleotídeo (na tabela está como "RefSeq seq accession") anotados anteriormente

Abra os parâmetros extras +Alignment parameters mude o parâmetro de resultados de 50 para 10 e o threshold de 0.05 para 0.005. Estamos sendo mais estritos.

Aqui baixaremos a HitsTable.csv e faremos um processo muito parecido com o que fizemos da última vez!


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
6. Na segunda linha da coluna Q, coloque a fórmula ```=A1<>P2```
7. Copie a fórmula em Q2 para todas as outras linhas da coluna Q onde existem valores na linha
8. Clique no botão **filtrar** () e desative o item FALSO
9. Copie as linhas filtradas para outra aba e:
   1.  Aperte ctrl+H, selecione "expressões regulares" e cole `|*.+$`. Aperte enter. Apenas sobrarão os identificadores do COG
   2.  Vamos aproveitar e substituir o 'COG5820' e o 'COG5821' por 'COG5820+21'
   3.  Agora vá na tabela reduzida e copie as colunas para conseguir os COGs na ordem certa
   4.  Cole os valores de COG na coluna B1
   5.  Na célula B2, digite `=TRANSPOSE(B1:DI1)`
   6.  Na célula C2, digite `=IF(COUNTIF(A:A;B2)>0;"'+'";"'-'")` e arraste para baixo
   7.  Cole os resultados da Coluna C em uma nova aba e faça a mesma operação de transposição, só que agora teremos uma linha
10. Essa linha pode ser colada diretamente na tabela reduzida do Galperin

### Embelezando a tabela

Vamos marcar a tabela!

Uma vez que tenhamos unido os dados da turma em uma tabela só, vamos fazer o seguinte:

1. Vamos remover as aspas, no ctrl+H substitua ' por vazio
2. Selecione toda a tabela e vá em "Formatar">"Formatação condicional". Na aba que se abre à sua direita, selecione "texto é exatamente" e crie duas regras com as cores de sua escolha para + e para -

# Perguntas guia para a análise

1) Após concluir cada etapa, desde a montagem do gene, até a interpretação final das árvores, qual foi a etapa mais desafiadora para vocês ? 

2) Olhando a árvore de vocês em relação à diversidade geral de procariotos (material de apoio com as árvores dos livros): como os grupos próximos do isolado de vocês se posicionam dentro dessa diversidade maior ? 

3) Em relação às sequências mais próximas do isolado: seu Bacillus fica perto de quem? Quais espécies/gêneros aparecem no entorno imediato ?

4) Olhando para os grupos de referência que vocês mesmos adicionaram (Clostridium, E. coli e parentes), eles aparecem na árvore de vocês na mesma relação que aparecem no material de apoio (árvore dos Gram-positivos)? Por que o gênero clostridium está separado do clado de Bacillus, embora ambas esporulam e são bactérias ambientais ? qual seria a hipótese de vocês ? 

5) Depois de pintar por bactérias esporulantes e por hospedeiro provável: esses dois padrões de cor se agrupam dentro dos mesmos subclados da árvore, ou ficam misturados? Um clado filogeneticamente próximo necessariamente significa bactérias com o mesmo estilo de vida (ambiental vs. comensal, esporulante vs. não)? 

6) Será que os genomas que selecionamos são produtores de espóros? Porque?
   A) Existe algum gene compartilhado pelos produtores?
