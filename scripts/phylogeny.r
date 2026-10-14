#!/usr/bin/env Rscript
# Bacillota: árvores anotadas, D de Fritz e reconstrução de estado ancestral.
#
# Uso:
#   Rscript scripts/phylogeny.r plot      --tree gtdb
#   Rscript scripts/phylogeny.r fritz     --tree ribo --marker spo0A --permut 199 --scope bacillota
#   Rscript scripts/phylogeny.r compare   --marker lpxC --scope bacillota
#   Rscript scripts/phylogeny.r ancestral --tree ribo --trait spore
#   Rscript scripts/phylogeny.r plot      --tree isolado
#
# Árvores: 16S, ribo e gtdb ficam em dados/. "isolado" é a árvore 16S com as
# sequências da turma (dados/ ou resultados/16S_isolado.tree.treefile); pontas
# sem metadados (o isolado, referências extras) aparecem destacadas na figura e
# ficam fora das estatísticas.
#
# Caracteres: os 61 marcadores, "spore" (fenótipo de Galperin: esporulante vs.
# não esporulante; desconhecido fica de fora) e "expected" (envelope, literatura).
#
# Todas as saídas vão para resultados/. Os Newicks em dados/ nunca são alterados.

suppressPackageStartupMessages({
  library(optparse)
  library(ape); library(phangorn); library(ggtree); library(ggplot2)
  library(dplyr); library(forcats); library(RColorBrewer); library(ggnewscale)
})

# ---------------------------------------------------------------------------
# Argumentos
# ---------------------------------------------------------------------------
option_list <- list(
  make_option('--project', default='.', help='Raiz do projeto [%default]'),
  make_option('--tree', default='gtdb', help='16S, ribo ou gtdb [%default]'),
  make_option('--marker', default='spo0A',
              help='Caráter para o D de Fritz: um marcador ou spore [%default]'),
  make_option('--permut', type='integer', default=199L, help='Permutações do D [%default]'),
  make_option('--scope', default='bacillota',
              help='bacillota (só Bacillota) ou all (inclui outgroups) [%default]'),
  make_option('--trait', default='spore',
              help='Caráter para a reconstrução ancestral: spore, expected ou um marcador [%default]'),
  make_option('--model', default='auto', help='auto, ER ou ARD [%default]'),
  make_option('--genes', default=NULL, help='Matriz de genes alternativa')
)
parsed <- parse_args(OptionParser(usage='%prog plot|fritz|compare|ancestral [opções]',
                                  option_list=option_list),
                     positional_arguments=1)
command <- parsed$args[1]
opt <- parsed$options
if (!command %in% c('plot','fritz','compare','ancestral')) stop('Comando desconhecido: ', command)
if (!opt$scope %in% c('bacillota','all')) stop('--scope deve ser bacillota ou all')
if (!opt$model %in% c('auto','ER','ARD')) stop('--model deve ser auto, ER ou ARD')

project <- normalizePath(opt$project, mustWork=TRUE)
f <- function(...) file.path(project, ...)
trees <- c('16S'='16S_final.tree.treefile', 'ribo'='ribomarkers_final.tree.treefile',
           'gtdb'='bac120_229.tree', 'isolado'='16S_isolado.tree.treefile')
main_trees <- c('16S','ribo','gtdb')   # as três comparadas por `compare`
tree_kind <- opt$tree
if (!tree_kind %in% names(trees) && command != 'compare') stop('Árvore desconhecida: ', tree_kind)
INGROUP <- 'Bacillota'
SHADE_ORDER <- 'Bacillales'
# Suporte mínimo para o ponto preto no nó: UFBoot nas árvores do IQ-TREE,
# bootstrap convencional na GTDB.
SUPPORT_MIN <- c('16S'=95, 'ribo'=95, 'isolado'=95, 'gtdb'=85)
# Pontas sem metadados (isolado da turma): sem vermelho, que já indica
# esporulação e presenças por orthogroup.
EXTRA_POINT <- 'black'; EXTRA_FILL <- '#ffe066'

# ---------------------------------------------------------------------------
# Leitura e validação dos dados
# ---------------------------------------------------------------------------
spore_file <- f('dados','heatmap_matrix2_spore_categories.tsv')
meta_file <- f('dados','tip_meta.tsv')
gene_file <- if (is.null(opt$genes)) f('dados','gene_presence_absence_final_61.tsv') else opt$genes
additional_file <- f('dados','additional_gene_calls.tsv')
output_dir <- f('resultados'); dir.create(output_dir, recursive=TRUE, showWarnings=FALSE)

required <- c(meta_file, spore_file, gene_file, additional_file)
missing <- required[!file.exists(required)]
if (length(missing)) stop('Arquivos ausentes:\n', paste(missing, collapse='\n'))
meta <- read.delim(meta_file, check.names=FALSE, stringsAsFactors=FALSE)
old <- read.delim(spore_file, check.names=FALSE, stringsAsFactors=FALSE)
genes <- read.delim(gene_file, check.names=FALSE, stringsAsFactors=FALSE)
additional <- read.delim(additional_file, check.names=FALSE, stringsAsFactors=FALSE)

for (name in c('label','assembly','phylum','class','genus','species'))
  if (!name %in% names(meta)) stop('Falta coluna em tip_meta.tsv: ', name)
for (name in c('assembly','label','Spore former','expected'))
  if (!name %in% names(old)) stop('Falta coluna de metadados: ', name)
if (!'assembly' %in% names(genes)) stop('Matriz gênica sem assembly')
if (anyDuplicated(meta$assembly) || anyDuplicated(meta$label) ||
    anyDuplicated(old$assembly) || anyDuplicated(genes$assembly)) stop('Identificadores duplicados')
if (!setequal(old$assembly, genes$assembly)) stop('Genomas diferentes entre matrizes histórica e gênica')
if (!all(genes$assembly %in% meta$assembly)) stop('Há assemblies na matriz gênica ausentes de tip_meta.tsv')
# A árvore 16S usa `label`; os labels das duas tabelas precisam concordar.
label_meta <- meta$label[match(old$assembly, meta$assembly)]
bad_labels <- old$assembly[old$label != label_meta]
if (length(bad_labels))
  stop('label difere entre tip_meta.tsv e a matriz histórica para: ', paste(head(bad_labels, 10), collapse=', '))

markers <- setdiff(names(genes), 'assembly')
if (!length(markers)) stop('Nenhum marcador na matriz')
for (g in markers) {
  vals <- suppressWarnings(as.numeric(genes[[g]]))
  if (any(is.na(vals)) || !all(vals %in% c(0,1))) stop('Marcador não binário ou ausente: ', g)
  genes[[g]] <- vals
}
meta$class <- sub('^c__', '', meta$class)
if (!'order' %in% names(meta)) meta$order <- NA_character_
meta$order <- sub('^o__', '', meta$order)
meta$order[!is.na(meta$order) & meta$order == ''] <- NA_character_
if (anyNA(meta$order)) message('Taxonomia: ordem indisponível para ', sum(is.na(meta$order)),
                               ' genomas; não inferimos a partir de gênero/classe.')

# Presenças acrescentadas apenas pelos orthogroups (vermelho no heatmap).
required_extra <- c('assembly','marker','bakta_direct','orthogroup_inferred')
if (!all(required_extra %in% names(additional))) stop('Colunas ausentes em additional_gene_calls.tsv')
if (anyDuplicated(additional[, c('assembly','marker')])) stop('Chaves assembly+marker duplicadas em additional_gene_calls.tsv')
if (any(!additional$assembly %in% genes$assembly) || any(!additional$marker %in% markers))
  stop('Presenças adicionais com assembly ou marcador desconhecido')
if (any(additional$bakta_direct != 0 | additional$orthogroup_inferred != 1))
  stop('A tabela de adições contém linhas que não são exclusivas de orthogroup')
extra_vals <- mapply(function(i, g) genes[[g]][i], match(additional$assembly, genes$assembly), additional$marker)
if (any(is.na(extra_vals)) || any(extra_vals != 1))
  stop('Presenças adicionais não correspondem a presenças na matriz final')
extra_keys <- paste(additional$assembly, additional$marker, sep='\t')

om_names <- c('bamA','fabZ','kdsA','kdsB','kdsC','kdsD',
              'lptA','lptB','lpxA','lpxB','lpxC','lpxD','lpxK','ompH','waaA')
om_cols <- intersect(om_names, markers)
spo_cols <- setdiff(markers, om_cols)

# ---------------------------------------------------------------------------
# Árvores
# ---------------------------------------------------------------------------
id_field <- function(kind) if (kind %in% c('16S','isolado')) 'label' else 'assembly'
tree_path <- function(kind) {
  candidates <- c(f('dados', trees[[kind]]), f('resultados', trees[[kind]]))
  hit <- candidates[file.exists(candidates)]
  if (!length(hit)) stop('Árvore ausente: ', trees[[kind]], ' (procurada em dados/ e resultados/)')
  hit[1]
}

# Caracteres binários/discretos por assembly. Devolve um vetor nomeado por
# assembly; NA = sem informação.
trait_values <- function(trait, labels=c('0','1')) {
  if (trait == 'spore') {
    v <- old[['Spore former']]
    out <- ifelse(v == 1, labels[2], ifelse(v == 2, labels[1], NA))
  } else if (trait == 'expected') {
    out <- old$expected
    out[!is.na(out) & trimws(out) == ''] <- NA
  } else if (trait %in% markers) {
    out <- ifelse(genes[[trait]][match(old$assembly, genes$assembly)] == 1, labels[2], labels[1])
  } else stop('Caráter desconhecido: ', trait, ' (use spore, expected ou um dos marcadores)')
  setNames(out, old$assembly)
}
ingroup_ids <- function(id) meta[[id]][meta$phylum == INGROUP]

# Enraizamento único para as três árvores: todos os genomas fora de Bacillota
# são o outgroup. O enraizamento é aceito se esse outgroup formar um grupo em
# alguma raiz possível (is.monophyletic com reroot=TRUE). Caso contrário, usa
# midpoint e isso fica registrado na saída, porque a raiz é decisiva para
# qualquer conclusão sobre estado ancestral.
root_tree <- function(t, id, kind) {
  # Pontas sem metadados (isolado da turma) não entram na decisão.
  out <- intersect(meta[[id]][meta$phylum != INGROUP], t$tip.label)
  if (!length(out)) return(list(tree=t, rooting='original (sem outgroup na amostra)'))
  if (ape::is.monophyletic(t, out, reroot=TRUE)) {
    t <- ape::root(t, outgroup=out, resolve.root=TRUE, edgelabel=TRUE)
    return(list(tree=t, rooting=sprintf('outgroup (%d genomas fora de %s)', length(out), INGROUP)))
  }
  t <- phangorn::midpoint(t)
  cat('ATENÇÃO [', kind, ']: os genomas fora de ', INGROUP, ' não formam um grupo nesta árvore; ',
      'enraizada por midpoint.\n', sep='')
  list(tree=t, rooting='midpoint (outgroup não monofilético)')
}

read_prepared_tree <- function(kind) {
  t <- read.tree(tree_path(kind))
  if (anyDuplicated(t$tip.label)) stop('Pontas duplicadas na árvore ', kind)
  id <- id_field(kind)
  available <- meta[[id]][meta$assembly %in% genes$assembly]
  tips <- intersect(t$tip.label, available)
  extra <- setdiff(t$tip.label, tips)
  if (length(tips) < 4) stop('Poucas pontas compatíveis na árvore ', kind, ': ', length(tips))
  message(kind, ': ', length(tips), ' pontas com metadados; ', length(extra), ' sem metadados')
  if (length(extra))
    cat('Pontas sem metadados (destacadas, fora das estatísticas): ',
        paste(head(extra, 15), collapse=', '), if (length(extra) > 15) ', ...', '\n', sep='')
  # A subárvore GTDB podada preserva nós com um único descendente; eles não
  # representam bifurcações entre as amostras. collapse.singles soma os
  # comprimentos e preserva as distâncias entre pontas.
  if (ape::has.singles(t)) {
    t <- ape::collapse.singles(t)
    message(kind, ': removidos nós com um único descendente (comprimentos somados)')
  }
  r <- root_tree(t, id, kind)
  list(tree=r$tree, id=id, rooting=r$rooting, extra=extra)
}

# Limpeza para métodos que simulam ou otimizam na árvore (caper, ace):
# ramos internos de comprimento zero viram politomias. Isso não inventa
# suporte nem comprimentos.
clean_for_models <- function(tr, kind) {
  tr$node.label <- NULL
  if (is.null(tr$edge.length) || anyNA(tr$edge.length) ||
      any(!is.finite(tr$edge.length)) || any(tr$edge.length < 0))
    stop('A árvore ', kind, ' não possui comprimentos de ramo válidos')
  if (ape::has.singles(tr)) tr <- ape::collapse.singles(tr)
  is_tip_edge <- tr$edge[,2] <= ape::Ntip(tr)
  if (any(tr$edge.length[is_tip_edge] <= 0))
    stop('A árvore ', kind, ' contém ramos terminais de comprimento zero')
  n_zero <- sum(tr$edge.length[!is_tip_edge] == 0)
  if (n_zero > 0L) {
    tr <- ape::di2multi(tr, tol=1e-12)
    if (ape::has.singles(tr)) tr <- ape::collapse.singles(tr)
    message(kind, ': colapsados ', n_zero, ' ramos internos de comprimento zero (politomias)')
  }
  # di2multi pode transformar a raiz bifurcante numa politomia basal; root.edge=0
  # informa ao ape/caper que ela continua sendo a raiz.
  if (!ape::is.rooted(tr)) tr$root.edge <- 0
  if (!ape::is.rooted(tr)) stop('Enraizamento inválido: ', kind)
  tr
}

# ---------------------------------------------------------------------------
# Figura: árvore + heatmap
# ---------------------------------------------------------------------------
build_plot_data <- function() {
  df <- merge(old[, c('assembly','label','Spore former','expected')], genes, by='assembly', sort=FALSE)
  # 0 não detectado; 1 presente; 2 presença acrescentada APENAS via orthogroup.
  for (g in markers) {
    v <- df[[g]]
    extra <- paste(df$assembly, g, sep='\t') %in% extra_keys
    df[[g]] <- ifelse(v == 0, '0', ifelse(extra, '2', '1'))
  }
  df
}

arrow_label <- function(p, node, txt, xoff, yoff, col='black') {
  xmax <- max(p$data$x, na.rm=TRUE)
  nd <- p$data[p$data$node == node, ]
  p + annotate('segment', x=nd$x - xoff*xmax, xend=nd$x, y=nd$y + yoff, yend=nd$y,
               colour=col, linewidth=.4,
               arrow=grid::arrow(length=grid::unit(4,'pt'), type='closed')) +
    annotate('label', x=nd$x - xoff*xmax, y=nd$y + yoff, label=txt, hjust=1, size=2.8,
             fontface='italic', colour=col, fill=scales::alpha('white', .85))
}

# Maiores clados formados só por pontas de `allowed` e com pelo menos uma de
# `required`. Um grupo monofilético dá um único clado; um parafilético, vários.
pure_clades <- function(phy, allowed, required=allowed) {
  desc <- phangorn::Descendants(phy, seq_len(Ntip(phy) + Nnode(phy)), 'tips')
  pure <- vapply(desc, function(d) all(phy$tip.label[d] %in% allowed), logical(1))
  parent <- phy$edge[match(seq_along(pure), phy$edge[, 2]), 1]
  keep <- which(pure & (is.na(parent) | !pure[parent]))
  keep <- keep[vapply(keep, function(k) any(phy$tip.label[desc[[k]]] %in% required), logical(1))]
  lapply(keep, function(k) phy$tip.label[desc[[k]]])
}

mark_taxa <- function(p, phy, ids, kind, extra=character(), order_dy=3, order_dx=.15, dy=6, dx=.08) {
  # A monofilia é avaliada sem as pontas sem metadados: o isolado da turma cai
  # dentro de Bacillales/Bacillus e não deve apagar o destaque.
  known <- if (length(extra)) drop.tip(phy, extra) else phy
  ord <- intersect(ids$id[!is.na(ids$order) & ids$order == SHADE_ORDER], phy$tip.label)
  # Genomas sem ordem atribuída não contam como evidência contra Bacillales:
  # podem ficar dentro do sombreado, e a figura informa quantos são.
  no_order <- intersect(ids$id[is.na(ids$order)], known$tip.label)
  if (!length(ord)) {
    found <- sort(unique(na.omit(ids$order[ids$id %in% phy$tip.label])))
    message(kind, ': nenhum genoma com ordem "', SHADE_ORDER, '". Valores na coluna order: ',
            paste(head(found, 20), collapse=', '))
  } else {
    clades <- pure_clades(known, c(ord, no_order), ord)
    nodes <- vapply(clades, function(tips)
      if (length(tips) == 1) match(tips, phy$tip.label) else getMRCA(phy, tips), numeric(1))
    sizes <- lengths(clades)
    n_unknown <- sum(unlist(clades) %in% no_order)
    p <- p + geom_hilight(data=data.frame(node=nodes), mapping=aes(node=node),
                          fill='#fee8c8', alpha=.6, extend=0)
    n <- length(p$layers)
    p$layers <- c(p$layers[n], p$layers[-n])   # sombreado por baixo da árvore
    label <- if (length(clades) == 1) SHADE_ORDER else sprintf('%s (%d clados)', SHADE_ORDER, length(clades))
    if (n_unknown) label <- sprintf('%s\n(inclui %d sem ordem atribuída)', label, n_unknown)
    p <- arrow_label(p, nodes[which.max(sizes)], label, order_dx, order_dy)
    cat(kind, ': ', SHADE_ORDER, ' = ', length(ord), ' genomas em ', length(clades),
        ' clado(s); ', n_unknown, ' genomas sem ordem dentro do sombreado.\n', sep='')
    if (n_unknown)
      cat('  Gêneros sem ordem dentro do sombreado: ',
          paste(sort(unique(ids$genus[match(intersect(unlist(clades), no_order), ids$id)])), collapse=', '),
          '\n  Para preencher a coluna order: Rscript preparacao/preencher_ordem.R\n', sep='')
    inside <- setdiff(known$tip.label[phangorn::Descendants(known, getMRCA(known, ord), 'tips')[[1]]],
                      c(ord, no_order))
    if (length(inside)) {
      intr <- table(ids$order[match(inside, ids$id)])
      cat('  Genomas de outras ordens dentro do ancestral comum de ', SHADE_ORDER, ': ',
          paste(sprintf('%s (%d)', names(intr), as.integer(intr)), collapse=', '), '\n', sep='')
    }
  }
  bac <- intersect(ids$id[ids$genus == 'Bacillus'], phy$tip.label)
  if (length(bac) > 1 && is.monophyletic(known, bac))
    p <- arrow_label(p, getMRCA(phy, bac), 'Bacillus', dx, dy)
  if (length(bac)) p <- p + geom_tiplab(aes(subset=label %in% !!bac, label=species),
                                        size=1.6, fontface='italic', offset=.005,
                                        colour='goldenrod4')
  p
}

base_tree_plot <- function(tr, id, kind) {
  md <- meta; md$label <- md[[id]]
  # As cores de classe só valem para Bacillota; o agrupamento em "Other"
  # considera apenas as classes de Bacillota.
  md$class_top <- NA_character_
  ing <- md$phylum == INGROUP
  md$class_top[ing] <- as.character(fct_lump_n(md$class[ing], n=5, other_level='Other'))
  lv <- setdiff(sort(unique(md$class_top[ing])), 'Other')
  class_colors <- setNames(c(brewer.pal(8, 'Set2')[seq_along(lv)], 'grey55'), c(lv, 'Other'))
  outphyla <- sort(setdiff(unique(md$phylum), INGROUP))
  shapes <- setNames(rep(c(15,16,17,18,8,0,1,2,5,6), length.out=length(outphyla)), outphyla)
  p <- ggtree(tr, linewidth=.2) %<+% md
  # Rótulos: IQ-TREE "SH-aLRT/UFBoot" ou "UFBoot"; GTDB "'100.0:c__Bacilli'".
  lab <- gsub("'", '', p$data$label)
  p$data$support <- suppressWarnings(as.numeric(sub('.*/', '', sub(':.*', '', lab))))
  min_support <- SUPPORT_MIN[[kind]]
  p + geom_tippoint(aes(subset=!is.na(phylum) & phylum == INGROUP, colour=class_top), size=1.5) +
    geom_tippoint(aes(subset=!is.na(phylum) & phylum != INGROUP, shape=phylum), colour='#7b3294', size=1.8) +
    geom_nodepoint(aes(subset=!isTip & support >= !!min_support), size=.7, colour='grey20') +
    scale_colour_manual(values=class_colors, na.translate=FALSE, name='Classe (Bacillota)') +
    scale_shape_manual(values=shapes, name='Outgroup (filo)')
}

support_note <- function(kind) {
  if (kind == 'gtdb') sprintf('Pontos pretos: bootstrap >= %g (GTDB). Não é UFBoot; o limiar é diferente do das outras árvores.', SUPPORT_MIN[[kind]])
  else sprintf('Pontos pretos: UFBoot >= %g (IQ-TREE).', SUPPORT_MIN[[kind]])
}

plot_one <- function(kind) {
  obj <- read_prepared_tree(kind); tr <- obj$tree; id <- obj$id
  ids <- meta %>% transmute(id=.data[[id]], genus, order)
  p <- base_tree_plot(tr, id, kind)
  p <- mark_taxa(p, tr, ids, kind, obj$extra, order_dy=if (id == 'label') 28 else 3,
                 order_dx=if (id == 'label') .05 else .15)
  extra_off <- .01 * max(p$data$x, na.rm=TRUE)
  if (length(obj$extra))
    p <- p + geom_tippoint(aes(subset=label %in% !!obj$extra), colour=EXTRA_POINT, size=2.6, shape=18) +
      geom_label(data=function(d) d[d$isTip & d$label %in% obj$extra, ],
                 aes(x=x + !!extra_off, y=y, label=label), inherit.aes=FALSE,
                 hjust=0, size=2.4, fontface='bold', colour='black', fill=EXTRA_FILL,
                 label.padding=grid::unit(1.5, 'pt'))
  raw <- build_plot_data()
  raw$tip <- meta[[id]][match(raw$assembly, meta$assembly)]
  raw <- raw[raw$tip %in% tr$tip.label, , drop=FALSE]
  make_block <- function(cols) {
    m <- raw[, cols, drop=FALSE]
    rownames(m) <- raw$tip
    m[] <- lapply(m, as.character)
    m
  }
  xmax <- max(p$data$x, na.rm=TRUE)
  hm <- function(plot, mat, off, width) gheatmap(plot, mat, offset=off, width=width,
    color=NA, colnames_position='top', colnames_angle=90, colnames_offset_y=1,
    font.size=2.5, hjust=0)
  # Quatro blocos: envelope (literatura), fenótipo (Galperin), OM/LPS, esporulação.
  gap <- xmax*.03; status_w <- .035; gene_w <- .025
  off_env <- gap
  off_sp <- off_env + status_w*xmax
  off_om <- off_sp + status_w*xmax + gap
  off_spo <- off_om + length(om_cols)*gene_w*xmax + gap
  # gheatmap já define uma escala de fill (a substituição é intencional) e avisa
  # sobre pontas sem linha na matriz (o isolado, que de fato não tem dados).
  suppressWarnings(suppressMessages({
    p <- p + new_scale_fill()
    p <- hm(p, make_block('expected'), off_env, status_w) +
      scale_fill_manual(values=c(diderm='#72874E', monoderm='#476F84'), na.value='grey92',
                        breaks=c('diderm','monoderm'), name='Envelope (literatura)')
    p <- p + new_scale_fill()
    p <- hm(p, make_block('Spore former'), off_sp, status_w) +
      scale_fill_manual(values=c('0'='grey92','1'='#A92222','2'='#BE7977'), na.value='grey92',
        breaks=c('1','2','0'), labels=c('esporulante','não esporulante','desconhecido'),
        name='Galperin')
    p <- p + new_scale_fill()
    p <- hm(p, make_block(om_cols), off_om, length(om_cols)*gene_w) +
      scale_fill_manual(values=c('0'='white','1'='black','2'='#fca103'), na.value='grey92',
        breaks=c('1','2'), labels=c('presente (Bakta)','acrescentado (orthogroup)'),
        name='Membrana externa/LPS')
    p <- p + new_scale_fill()
    p <- hm(p, make_block(spo_cols), off_spo, length(spo_cols)*gene_w) +
      scale_fill_manual(values=c('0'='white','1'='#2166ac','2'='#fca103'), na.value='grey92',
        breaks=c('1','2'), labels=c('presente (Bakta)','acrescentado (orthogroup)'),
        name='Genes de esporulação')
  }))
  # Espaço extra para os rótulos verticais dos 61 genes.
  p <- p + hexpand(.12, direction=-1) + expand_limits(y=c(0, Ntip(tr) + 22)) +
    labs(title=sprintf('Árvore %s | raiz: %s', kind, obj$rooting), caption=support_note(kind)) +
    theme(plot.margin=margin(t=30, r=42, b=18, l=14, unit='mm'),
          legend.text=element_text(size=7), legend.title=element_text(size=8),
          legend.key.size=grid::unit(4, 'mm'))
  outfile <- f('resultados', paste0('heatmap_', kind))
  plot_width <- max(23, 10 + length(markers)*0.24)
  ggsave(paste0(outfile, '.png'), p, width=plot_width, height=19, dpi=165, bg='white', limitsize=FALSE)
  ggsave(paste0(outfile, '.pdf'), p, width=plot_width, height=19, bg='white', limitsize=FALSE)
  cat('Raiz: ', obj$rooting, '\n', sep='')
  cat('SAVED\t', paste0(outfile, '.png'), '\n', sep='')
}

# ---------------------------------------------------------------------------
# D de Fritz e Purvis
# ---------------------------------------------------------------------------
# scope = 'bacillota': só pontas de Bacillota. Os outgroups não têm genes de
# esporulação e quase sempre têm os de LPS; incluí-los infla o agrupamento com
# o contraste trivial "Bacillota vs. resto".
# restrict: vetor de assemblies (usado por compare para igualar as amostras).
fritz_one <- function(kind, marker, permut, scope, restrict=NULL) {
  obj <- read_prepared_tree(kind); tr <- obj$tree
  tv <- trait_values(marker)
  df <- data.frame(assembly=names(tv), trait=as.numeric(tv),
                   tip=meta[[obj$id]][match(names(tv), meta$assembly)],
                   phylum=meta$phylum[match(names(tv), meta$assembly)],
                   stringsAsFactors=FALSE)
  df <- df[!is.na(df$trait) & !is.na(df$tip) & df$tip %in% tr$tip.label, , drop=FALSE]
  if (scope == 'bacillota') df <- df[df$phylum == INGROUP, , drop=FALSE]
  if (!is.null(restrict)) df <- df[df$assembly %in% restrict, , drop=FALSE]
  n1 <- sum(df$trait == 1); n0 <- sum(df$trait == 0)
  result <- data.frame(tree=kind, marker=marker, scope=scope, n=nrow(df), present=n1, absent=n0,
                       D=NA_real_, Pval1=NA_real_, Pval0=NA_real_, rooting=obj$rooting, status='OK')
  if (nrow(df) < 4 || n1 < 2 || n0 < 2) { result$status <- 'Variação insuficiente'; return(result) }
  tr <- clean_for_models(keep.tip(tr, df$tip), kind)
  v <- data.frame(tip=df$tip, trait=df$trait)
  x <- tryCatch(caper::phylo.d(data=v, phy=tr, names.col=tip, binvar=trait, permut=permut),
                error=function(e) e)
  if (inherits(x, 'error')) { result$status <- conditionMessage(x); return(result) }
  result$D <- x$DEstimate; result$Pval1 <- x$Pval1; result$Pval0 <- x$Pval0
  result
}

# Assemblies presentes nas três árvores (16S traduzido de label para assembly).
common_assemblies <- function() {
  sets <- lapply(main_trees, function(k) {
    tips <- read.tree(tree_path(k))$tip.label
    id <- id_field(k)
    meta$assembly[meta[[id]] %in% tips]
  })
  Reduce(intersect, sets)
}

# ---------------------------------------------------------------------------
# Reconstrução de estado ancestral (ape::ace, modelo Mk)
# ---------------------------------------------------------------------------
# Usa a árvore inteira (outgroups incluídos), porque são eles que informam o
# estado na base de Bacillota. A probabilidade no ancestral comum de Bacillota
# depende do enraizamento, da amostragem de outgroups e do modelo.
ancestral_one <- function(kind, trait, model_choice) {
  obj <- read_prepared_tree(kind); tr <- obj$tree; id <- obj$id
  if (trait == 'spore') {
    vals <- trait_values('spore', c('não esporulante', 'esporulante'))
    states_label <- 'Esporulação (Galperin)'
  } else if (trait == 'expected') {
    vals <- trait_values('expected'); states_label <- 'Envelope (literatura)'
  } else {
    vals <- trait_values(trait, c('ausente', 'presente')); states_label <- trait
  }
  tip_of <- meta[[id]][match(names(vals), meta$assembly)]
  keep <- !is.na(vals) & tip_of %in% tr$tip.label
  if (sum(!keep & tip_of %in% tr$tip.label))
    message(kind, ': ', sum(!keep & tip_of %in% tr$tip.label), ' pontas sem estado conhecido foram removidas')
  x <- setNames(factor(vals[keep]), tip_of[keep])
  if (nlevels(x) < 2) stop('O caráter ', trait, ' não varia nesta árvore')
  tr <- clean_for_models(keep.tip(tr, names(x)), kind)
  # ace exige árvore dicotômica: politomias viram ramos de comprimento zero,
  # o que não muda a verossimilhança de um modelo Mk.
  tr <- multi2di(tr, random=FALSE)
  x <- x[tr$tip.label]

  fit <- function(m) tryCatch(suppressWarnings(ace(x, tr, type='discrete', model=m, marginal=TRUE)), error=function(e) e)
  fits <- list()
  for (m in if (model_choice == 'auto') c('ER','ARD') else model_choice) fits[[m]] <- fit(m)
  ok <- Filter(function(z) !inherits(z, 'error'), fits)
  if (!length(ok)) stop('ace falhou: ', conditionMessage(fits[[1]]))
  aics <- sapply(ok, AIC)
  best <- names(which.min(aics))
  a <- ok[[best]]

  lik <- a$lik.anc
  rownames(lik) <- Ntip(tr) + seq_len(Nnode(tr))
  ing <- intersect(ingroup_ids(id), tr$tip.label)
  mrca <- if (length(ing) > 1 && is.monophyletic(tr, ing)) getMRCA(tr, ing) else NA
  root_node <- Ntip(tr) + 1
  # Uma linha por modelo: comparar ER e ARD faz parte da interpretação. O ARD
  # pode estimar uma das taxas como ~0 (mudança irreversível), o que deixa as
  # probabilidades ancestrais próximas de 0 ou 1.
  res <- do.call(rbind, lapply(names(fits), function(m) {
    z <- fits[[m]]
    row <- data.frame(tree=kind, trait=trait, n=Ntip(tr), rooting=obj$rooting, model=m,
                      chosen=(m == best), check.names=FALSE)
    if (inherits(z, 'error')) {
      row$logLik <- NA; row$AIC <- NA; row$status <- conditionMessage(z); return(row)
    }
    row$logLik <- z$loglik; row$AIC <- AIC(z)
    im <- z$index.matrix; lv <- levels(x); tx <- character()
    for (i in seq_along(lv)) for (j in seq_along(lv)) if (i != j && !is.na(im[i, j]))
      tx <- c(tx, sprintf('%s->%s=%.3g', lv[i], lv[j], z$rates[im[i, j]]))
    row$rates <- paste(unique(tx), collapse='; ')
    l <- z$lik.anc; rownames(l) <- Ntip(tr) + seq_len(Nnode(tr))
    for (s in colnames(l)) {
      row[[paste0('P_', s, '_ancestral_', INGROUP)]] <- if (is.na(mrca)) NA else l[as.character(mrca), s]
      row[[paste0('P_', s, '_raiz')]] <- l[as.character(root_node), s]
    }
    row$status <- 'OK'
    row
  }))
  if (is.na(mrca)) res$nota <- paste(INGROUP, 'não monofilético nesta árvore; sem nó ancestral definido')

  # Figura: pontas pelo estado observado; nós internos pela probabilidade do
  # estado focal (diderm, ou "presente" para marcadores).
  focal <- intersect(c('diderm', 'presente', 'esporulante'), colnames(lik))
  focal <- if (length(focal)) focal[1] else colnames(lik)[2]
  nd <- data.frame(node=as.integer(rownames(lik)), p=lik[, focal])
  tipdf <- data.frame(label=names(x), state=as.character(x))
  md <- meta; md$label <- md[[id]]
  p <- ggtree(tr, linewidth=.25) %<+% md %<+% tipdf
  p$data <- left_join(p$data, nd, by='node')
  st_cols <- setNames(c('#1b7837', 'grey70', '#2166ac', '#b2182b')[seq_len(nlevels(x))], levels(x))
  if (all(c('diderm','monoderm') %in% levels(x))) st_cols[c('diderm','monoderm')] <- c('#1b7837','grey70')
  if (all(c('presente','ausente') %in% levels(x))) st_cols[c('presente','ausente')] <- c('#1b7837','grey70')
  if (all(c('esporulante','não esporulante') %in% levels(x)))
    st_cols[c('esporulante','não esporulante')] <- c('#1b7837','grey70')
  p <- p + geom_tiplab(aes(label=species), size=1.3, offset=.01*max(p$data$x, na.rm=TRUE),
                       fontface='italic', colour='grey30') +
    hexpand(.25) +
    geom_tippoint(aes(colour=state), size=1.4) +
    scale_colour_manual(values=st_cols, name=states_label) +
    new_scale_colour() +
    geom_nodepoint(aes(subset=!isTip, colour=p), size=1.6) +
    scale_colour_gradient(low='grey85', high='#1b7837', limits=c(0,1),
                          name=paste0('P(', focal, ')\nno nó'))
  if (!is.na(mrca)) {
    p <- p + geom_point2(aes(subset=node == mrca), shape=21, size=4, stroke=1, colour='black', fill=NA)
    ndm <- p$data[p$data$node == mrca, ]
    p <- p + annotate('label', x=ndm$x + .015*max(p$data$x, na.rm=TRUE), y=ndm$y,
                      label=sprintf('Ancestral de %s\nP(%s) = %.2f', INGROUP, focal,
                                    lik[as.character(mrca), focal]),
                      hjust=0, size=2.8, fill=scales::alpha('white', .85))
  }
  p <- p + labs(title=sprintf('Estado ancestral de "%s": árvore %s', trait, kind),
                subtitle=sprintf('Modelo %s (Mk, ape::ace, reconstrução marginal); raiz: %s',
                                 best, obj$rooting))
  outfile <- f('resultados', paste0('ancestral_', kind, '_', trait))
  ggsave(paste0(outfile, '.png'), p, width=10, height=14, dpi=150, bg='white')
  ggsave(paste0(outfile, '.pdf'), p, width=10, height=14, bg='white')
  cat('SAVED\t', paste0(outfile, '.png'), '\n', sep='')
  res
}

# ---------------------------------------------------------------------------
# Execução
# ---------------------------------------------------------------------------
set.seed(42)
write_result <- function(ans, name) {
  dest <- f('resultados', name)
  write.table(ans, dest, sep='\t', row.names=FALSE, quote=FALSE, na='NA')
  print(ans, row.names=FALSE)
  cat('SAVED\t', dest, '\n', sep='')
}

if (command == 'plot') {
  plot_one(tree_kind)
} else if (command %in% c('fritz','compare')) {
  suppressPackageStartupMessages(library(caper))
  if (!opt$marker %in% c(markers, 'spore')) stop('Marcador desconhecido: ', opt$marker)
  if (is.na(opt$permut) || opt$permut < 10) stop('permut deve ser pelo menos 10')
  if (command == 'fritz') {
    ans <- fritz_one(tree_kind, opt$marker, opt$permut, opt$scope)
    write_result(ans, paste0('fritz_', tree_kind, '_', opt$marker, '_', opt$scope, '.tsv'))
  } else {
    common <- common_assemblies()
    cat('Comparação restrita aos ', length(common), ' assemblies presentes nas três árvores.\n', sep='')
    ans <- do.call(rbind, lapply(main_trees, function(k)
      fritz_one(k, opt$marker, opt$permut, opt$scope, restrict=common)))
    write_result(ans, paste0('fritz_comparacao_', opt$marker, '_', opt$scope, '.tsv'))
  }
} else if (command == 'ancestral') {
  ans <- ancestral_one(tree_kind, opt$trait, opt$model)
  write_result(ans, paste0('ancestral_', tree_kind, '_', opt$trait, '.tsv'))
}
