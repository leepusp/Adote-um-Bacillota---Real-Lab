# Preparando o alinhamento 

## Montando o arquivo

O arquivo que baixamos do blast **não inclui** a nossa sequência de 16S.
Se queremos ver como ela se posiciona na árvore, precisamos incluí-la.

Também faltaria na nossa árvore algumas sequências de bacilos conhecidos. 
Copie e cole a nossa sequência e os conteúdos do arquivo "16S_referencias.fasta" dentro do arquivo que nós baixamos do BLAST

## Inspecionando o arquivo

Ao abrir o Aliview vamos observar a nossa lista de sequências não alinhadas.

Dando zoom-out (ctrl + rodinha do mouse) vamos observar algumas "bandas" que parecem meio estranhas em relação ao resto.

Aqui percebemos que infelizmente o BLAST não nos deu todas as sequências na mesma orientação. 

Se abrirmos a saída do BLAST do passo anterior e formos na aba de alinhamentos, vamos observar que as sequências estranhas possuem a informação "Plus/Minus"

