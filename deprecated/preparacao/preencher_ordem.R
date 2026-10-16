#!/usr/bin/env Rscript
# Preenche a coluna `order` de dados/tip_meta.tsv a partir da taxonomia GTDB.
#
# Fontes, em ordem de prioridade:
#   1. --taxonomy: arquivo bac120_taxonomy*.tsv(.gz) da GTDB (accession<TAB>d__;p__;c__;o__;...).
#   2. Rótulos dos nós internos da árvore GTDB podada (ex.: '100.0:o__Bacillales; f__Bacillaceae').
#      A ordem de cada ponta é a do ancestral mais próximo cujo rótulo traz "o__".
#
# Só preenche ordens ausentes; ordens já presentes são conferidas e as
# divergências são listadas, nunca sobrescritas. Grava uma cópia de segurança.
#
# Uso (na raiz do projeto):
#   Rscript preparacao/preencher_ordem.R
#   Rscript preparacao/preencher_ordem.R --taxonomy bac120_taxonomy_r226.tsv.gz

suppressPackageStartupMessages({ library(optparse); library(ape); library(phangorn) })

opt <- parse_args(OptionParser(option_list=list(
  make_option('--meta', default='dados/tip_meta.tsv'),
  make_option('--tree', default='dados/bac120_229.tree'),
  make_option('--taxonomy', default=NULL, help='bac120_taxonomy da mesma release da árvore'),
  make_option('--dry-run', action='store_true', default=FALSE, dest='dry_run',
              help='Só mostra o que seria preenchido')
)))

meta <- read.delim(opt$meta, check.names=FALSE, stringsAsFactors=FALSE, na.strings=c('', 'NA'))
if (!'order' %in% names(meta)) meta$order <- NA_character_
meta$order <- sub('^o__', '', meta$order)
source_col <- rep(NA_character_, nrow(meta))

order_from_string <- function(x) {
  o <- regmatches(x, regexpr("o__[^;']*", x))
  ifelse(length(o), trimws(sub('^o__', '', o)), NA_character_)
}

# 1. Arquivo de taxonomia da GTDB
gtdb_order <- setNames(rep(NA_character_, nrow(meta)), meta$assembly)
if (!is.null(opt$taxonomy)) {
  tax <- read.delim(opt$taxonomy, header=FALSE, col.names=c('accession', 'taxonomy'),
                    stringsAsFactors=FALSE)
  hit <- match(meta$assembly, tax$accession)
  gtdb_order[] <- vapply(tax$taxonomy[hit], function(x) if (is.na(x)) NA_character_ else order_from_string(x), '')
  cat('Taxonomia GTDB: ', sum(!is.na(hit)), ' de ', nrow(meta), ' assemblies encontrados\n', sep='')
}

# 2. Rótulos da árvore GTDB
if (file.exists(opt$tree)) {
  tr <- read.tree(opt$tree)
  labs <- tr$node.label
  n_tax <- sum(grepl('o__', labs))
  cat('Árvore GTDB: ', n_tax, ' nós internos com rótulo de ordem\n', sep='')
  if (n_tax) {
    node_order <- vapply(labs, function(x) if (grepl('o__', x)) order_from_string(x) else NA_character_, '')
    anc <- Ancestors(tr, seq_len(Ntip(tr)), 'all')   # do pai até a raiz
    tip_order <- vapply(anc, function(a) {
      o <- node_order[a - Ntip(tr)]; o <- o[!is.na(o)]
      if (length(o)) o[1] else NA_character_
    }, '')
    names(tip_order) <- tr$tip.label
    from_tree <- tip_order[meta$assembly]
    use <- is.na(gtdb_order) & !is.na(from_tree)
    gtdb_order[use] <- from_tree[use]
  }
}

# Conferência e preenchimento
known <- !is.na(meta$order) & !is.na(gtdb_order)
conflict <- known & meta$order != gtdb_order
if (any(conflict)) {
  cat('\nDivergências (mantida a ordem existente):\n')
  print(data.frame(assembly=meta$assembly[conflict], tip_meta=meta$order[conflict],
                   gtdb=gtdb_order[conflict]), row.names=FALSE)
}
fill <- is.na(meta$order) & !is.na(gtdb_order)
cat('\nOrdens ausentes: ', sum(is.na(meta$order)), '; preenchidas agora: ', sum(fill), '\n', sep='')
if (any(fill)) print(data.frame(assembly=meta$assembly[fill], genus=meta$genus[fill],
                                order=gtdb_order[fill]), row.names=FALSE)
still <- is.na(meta$order) & !fill
if (any(still)) cat('Continuam sem ordem: ', paste(meta$assembly[still], collapse=', '), '\n', sep='')

if (!opt$dry_run && any(fill)) {
  backup <- paste0(opt$meta, '.bak')
  file.copy(opt$meta, backup, overwrite=FALSE)
  meta$order[fill] <- gtdb_order[fill]
  write.table(meta, opt$meta, sep='\t', quote=FALSE, row.names=FALSE, na='')
  cat('\nGravado: ', opt$meta, ' (original em ', backup, ')\n', sep='')
}
