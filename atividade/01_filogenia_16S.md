# Atividade prática — Filogenia de 16S rRNA

## Objetivo

Vamos construir uma árvore filogenética de 16S rRNA para posicionar o isolado ambiental em um contexto mais amplo de diversidade bacteriana: *Bacillus* e parentes próximos (via BLAST), além de grupos de comparação, como *Clostridium* e *Escherichia coli*/outras Gram-negativas.

> **Resultado esperado:** uma árvore com aproximadamente 700 sequências (500 do BLAST do isolado + 100 de *Clostridium* + 100 de *E. coli* + o nosso isolado ambiental). O total poderá variar devido a resultados repetidos entre buscas.

## Passo 1 — Visualização de cromatogramas no Chromas e montagem do gene no GENEIOUS PRIME

Primeiro, abram os arquivos de sequenciamento Sanger (`.ab1`) no **Chromas**. Visualizem os cromatogramas, avaliem a qualidade das bases e copiem uma sequência inicial para uma identificação preliminar por BLASTn.

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

Como o BLAST do isolado retorna principalmente organismos próximos, adicionaremos dois grupos de comparação mais distantes. Algumas sequências podem se repetir entre as buscas.

**Grupo *Clostridium*:**

1. No NCBI Nucleotide, procurem uma sequência de referência de 16S rRNA de *Clostridium botulinum*.
2. Executem um novo BLASTn no mesmo banco de 16S de linhagens-tipo.
3. Baixem os **100 primeiros hits em FASTA**.
4. Abram o conjunto principal no **AliView** e utilizem **Sequence → Add Sequences** para importar esse FASTA.

**Grupo *E. coli* e outras Gram-negativas:**

1. Procurem uma sequência de referência de 16S rRNA de *Escherichia coli*.
2. Repitam a busca no mesmo banco, baixando os **100 primeiros hits em FASTA**.
3. Adicionem as sequências ao conjunto no AliView pelo mesmo procedimento.

**Atenção:** podem aparecer outros gêneros relacionados, como *Shigella* e *Salmonella*. Isso faz parte da diversidade do grupo; não os excluam apenas pelo nome do gênero.

Incluam também o **consenso do isolado ambiental** e salvem o conjunto final como **`database.fasta`**. O conjunto deverá ter aproximadamente 700 sequências; verifiquem eventuais entradas duplicadas.

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
cd ~/database
```

O arquivo `database.fasta`, montado no AliView, deverá estar nessa pasta do Haldane. Se ele estiver somente no computador da LEEP (e não em uma pasta compartilhada), enviem-no a partir do terminal da LEEP com `scp -P 25022 database.fasta user@haldane.icb.usp.br:~/database/`, substituindo `user` pelo login individual.

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
iqtree2 -s trimal.fasta -st DNA -m "$MODELO" --fast -T AUTO --threads-max 8 -mem 6G -pre "iqtree/$MODELO/arvore"
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

Como atividade de interpretação, comparem os agrupamentos por patogenicidade e por hospedeiro/ambiente. Organizem-se em duplas para investigar a mesma árvore e discutir diferenças de anotação e interpretação dos clados.

---

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
