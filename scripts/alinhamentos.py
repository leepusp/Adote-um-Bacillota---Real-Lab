"""Visualização de alinhamentos para os notebooks da atividade.

Uso nos notebooks:

    from alinhamentos import visualizador, arquivos_de_alinhamento
    display(visualizador(lambda: arquivos_de_alinhamento(DATA_DIR, RESULT_DIR)))

`visualizador` recebe uma função que devolve {nome exibido: fonte}. A fonte pode
ser um caminho para um FASTA alinhado ou um DataFrame com a coluna `sequence`
(e, opcionalmente, `assembly`). A função é chamada de novo no botão
"Atualizar lista", para que alinhamentos recém-gerados em resultados/ apareçam.
"""
from pathlib import Path

import numpy as np
import matplotlib.pyplot as plt
import ipywidgets as widgets
from Bio import AlignIO
from IPython.display import display, clear_output
from matplotlib.colors import ListedColormap
from matplotlib.patches import Patch

ALIGNMENT_SUFFIXES = ('.aln', '.fasta', '.fas', '.fa', '.afa')

# Cores dos resíduos. Nucleotídeos: uma cor por base. Aminoácidos: grupos
# físico-químicos no estilo do Clustal X.
NUC_GROUPS = [('A', 'A', '#4daf4a'), ('C', 'C', '#377eb8'), ('G', 'G', '#ff7f00'),
              ('T', 'TU', '#e41a1c')]
AA_GROUPS = [('hidrofóbicos', 'AILMFWV', '#80a0f0'), ('positivos', 'KR', '#f01505'),
             ('negativos', 'DE', '#c048c0'), ('polares', 'NQST', '#15c015'),
             ('cisteína', 'C', '#f08080'), ('glicina', 'G', '#f09048'),
             ('prolina', 'P', '#c0c000'), ('aromáticos', 'HY', '#15a4a4')]
GAP_COLOR, OTHER_COLOR = '#ffffff', '#bdbdbd'


def arquivos_de_alinhamento(*folders):
    """Lista FASTAs alinhados (mesmo comprimento) nas pastas indicadas."""
    found = {}
    for folder in folders:
        for path in sorted(Path(folder).glob('*')):
            if path.suffix.lower() in ALIGNMENT_SUFFIXES and path.is_file():
                found[f'{Path(folder).name}/{path.name}'] = path
    return found


def carregar(source):
    """Devolve (ids, matriz de caracteres) de um FASTA alinhado ou DataFrame."""
    if isinstance(source, (str, Path)):
        try:
            aln = AlignIO.read(source, 'fasta')
        except ValueError as e:
            raise ValueError(f'{Path(source).name} não é um alinhamento '
                             f'(as sequências precisam ter o mesmo comprimento): {e}')
        ids = [rec.id for rec in aln]
        matrix = np.array(aln)
    else:
        df = source.dropna(subset=['sequence'])
        seqs = df['sequence'].astype(str).str.upper().tolist()
        if len({len(s) for s in seqs}) != 1:
            raise ValueError('As sequências não têm o mesmo comprimento.')
        ids = df['assembly'].astype(str).tolist() if 'assembly' in df else [str(i) for i in df.index]
        matrix = np.array([list(s) for s in seqs])
    return ids, np.char.upper(matrix.astype(str))


def _alphabet(matrix):
    letters = set(np.unique(matrix)) - {'-', '.', 'N', 'X', '?'}
    return 'nuc' if letters <= set('ACGTU') else 'aa'


def _codes(matrix, groups):
    """Converte a matriz de letras em índices de cor (0 = gap, último = outro)."""
    codes = np.full(matrix.shape, len(groups) + 1, dtype=np.uint8)
    codes[np.isin(matrix, ['-', '.'])] = 0
    for i, (_, letters, _) in enumerate(groups, start=1):
        codes[np.isin(matrix, list(letters))] = i
    return codes


def _order(ids, highlight):
    """Coloca as sequências destacadas no topo, mantendo a ordem do arquivo."""
    if not highlight:
        return list(range(len(ids))), 0
    hit = [i for i, s in enumerate(ids) if highlight.lower() in s.lower()]
    rest = [i for i in range(len(ids)) if i not in set(hit)]
    return hit + rest, len(hit)


def mapa_de_gaps(ax, ids, matrix, title, highlight=''):
    order, n_hit = _order(ids, highlight)
    gaps = np.isin(matrix[order], ['-', '.'])
    ax.imshow(gaps, aspect='auto', interpolation='nearest', vmin=0, vmax=1,
              cmap=ListedColormap(['#2878B5', '#F4A340']))
    if n_hit:
        ax.axhline(n_hit - 0.5, color='black', lw=1.5)
        if n_hit <= 30:
            ax.set_yticks(range(n_hit), [ids[i] for i in order[:n_hit]], fontsize=7, color='#b2182b')
    ax.set(xlabel='Posição no alinhamento', ylabel='Sequências', title=title)
    ax.legend(handles=[Patch(color='#2878B5', label='Resíduo'), Patch(color='#F4A340', label='Gap')],
              loc='upper left', bbox_to_anchor=(1.01, 1))


def residuos(ax, ids, matrix, title, start, width, highlight='', max_rows=80):
    order, n_hit = _order(ids, highlight)
    order = order[:max_rows]
    window = matrix[order, start:start + width]
    groups = NUC_GROUPS if _alphabet(matrix) == 'nuc' else AA_GROUPS
    colors = [GAP_COLOR] + [c for *_, c in groups] + [OTHER_COLOR]
    ax.imshow(_codes(window, groups), aspect='auto', interpolation='nearest',
              cmap=ListedColormap(colors), vmin=0, vmax=len(colors) - 1)
    if window.shape[1] <= 120 and window.shape[0] <= 60:
        for (r, c), letter in np.ndenumerate(window):
            ax.text(c, r, letter, ha='center', va='center', fontsize=6, family='monospace')
    ax.set_yticks(range(len(order)), [ids[i] for i in order], fontsize=6)
    for tick, i in zip(ax.get_yticklabels(), order):
        if i in order[:n_hit]:
            tick.set_color('#b2182b')
    xt = np.linspace(0, window.shape[1] - 1, min(10, window.shape[1])).astype(int)
    ax.set_xticks(xt, xt + start + 1)
    ax.set(xlabel='Posição no alinhamento', title=title)
    ax.legend(handles=[Patch(facecolor=c, edgecolor='grey', label=n) for n, _, c in groups] +
                      [Patch(facecolor=GAP_COLOR, edgecolor='grey', label='gap')],
              loc='upper left', bbox_to_anchor=(1.01, 1), fontsize=8)


def visualizador(sources):
    """Widget: escolha um alinhamento e veja o mapa de gaps ou os resíduos."""
    style = {'description_width': '110px'}
    source_select = widgets.Dropdown(description='Alinhamento:', style=style,
                                     layout=widgets.Layout(width='520px'))
    refresh = widgets.Button(description='Atualizar lista', icon='refresh')
    mode = widgets.ToggleButtons(options=[('Mapa de gaps', 'gaps'), ('Resíduos (janela)', 'res')],
                                 value='gaps', style={'button_width': '150px'})
    start = widgets.BoundedIntText(value=1, min=1, max=10**6, description='Início:',
                                   style=style, layout=widgets.Layout(width='220px'))
    width = widgets.IntSlider(value=80, min=20, max=400, step=10, description='Largura:',
                              style=style, layout=widgets.Layout(width='420px'))
    highlight = widgets.Text(value='', placeholder='ex.: AMB_grupo01', description='Destacar:',
                             style=style, layout=widgets.Layout(width='420px'))
    show = widgets.Button(description='Mostrar alinhamento', icon='eye', button_style='primary')
    out = widgets.Output()
    window_box = widgets.HBox([start, width])
    current = {}

    def fill(_=None):
        current.clear()
        current.update(sources())
        previous = source_select.value
        names = list(current) or ['(nenhum alinhamento encontrado)']
        source_select.options = names
        source_select.value = previous if previous in names else names[0]

    def toggle(_=None):
        window_box.layout.display = '' if mode.value == 'res' else 'none'

    def draw(_):
        with out:
            clear_output(wait=True)
            name = source_select.value
            if name not in current:
                print('Nenhum alinhamento disponível.'); return
            try:
                ids, matrix = carregar(current[name])
            except Exception as e:
                print('Erro:', e); return
            n, length = matrix.shape
            print(f'{name}: {n} sequências × {length} posições; '
                  f'{100 * np.isin(matrix, ["-", "."]).mean():.1f}% de gaps')
            if highlight.value and not any(highlight.value.lower() in i.lower() for i in ids):
                print(f'Nenhuma sequência contém "{highlight.value}".')
            if mode.value == 'gaps':
                fig, ax = plt.subplots(figsize=(15, 6))
                mapa_de_gaps(ax, ids, matrix, name, highlight.value)
            else:
                s = min(start.value, length) - 1
                rows = min(n, 80)
                fig, ax = plt.subplots(figsize=(15, max(4, rows * 0.16)))
                residuos(ax, ids, matrix, f'{name} — posições {s + 1}–{min(s + width.value, length)}',
                         s, width.value, highlight.value)
                if n > rows:
                    print(f'Mostrando as primeiras {rows} sequências (as destacadas vêm primeiro).')
            plt.tight_layout()
            plt.show()

    refresh.on_click(fill)
    mode.observe(toggle, names='value')
    show.on_click(draw)
    fill(); toggle()
    return widgets.VBox([widgets.HBox([source_select, refresh]), mode, window_box,
                         highlight, show, out])
