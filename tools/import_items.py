#!/usr/bin/env python3
"""問題管理表の形式（1問1フォルダの item.json ＋ problem.svg / solution.svg）を content/ に取り込む

    python3 tools/import_items.py <content フォルダ> [--sheet CELESTIA問題管理表.xlsx] [--verify <verify フォルダ>]

  ・<content フォルダ>/<ID>/item.json を読み、content/problems/<分野>/<ID>/problem.json を書く。
    問題ID は item.json の id をそのまま使う（変えない）。同じ ID がすでにあれば上書きする。
  ・単元は item.json の field（単元名）で taxonomy.json から引く（表記ゆれは ALIASES で吸収）。
  ・難易度 1〜5 → level：1=基礎(0)、2=標準(1)、3=国公立2次(3)、4=難関(4)、5=最難関(5)。
  ・問題文の enumerate を小問に、答え "(1) … (2) …" を小問ごとの答えに分ける（入力採点はしない：自己採点）。
  ・参考元（source・src_no）はアプリに出さないので取り込まない。
  ・管理表があれば、作成日・作成方法・検証3項目（二重解答／SymPy等／人間確認）を quality と provenance に写す。
  ・--verify があれば検証スクリプト <ID>.py を verify/ にコピーする。
  ・SVG の id は図ごとに接頭辞を付けて、同じページに並んでも衝突しないようにする。
"""
import json
import re
import shutil
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from content_lib import CONTENT, ROOT  # noqa: E402

ALIASES = {'剛体のつり合い': '剛体の運動', '気体の状態変化と熱力学第1法則': '気体の状態変化と熱力学第一法則'}
LEVEL = {1: 0, 2: 1, 3: 3, 4: 4, 5: 5}

# タグ・テーマ・核に含まれる語 → 概念（上から順に、最大3つ）
KEYWORDS = [
    (r'放物|投射|投げ', 'projectile-motion'), (r'等加速度', 'uniform-acceleration'), (r'相対運動|相対速度|相対加速度', 'kinematics-components'),
    (r'万有引力|ケプラー|人工衛星|面積速度|宇宙速度|楕円軌道', 'universal-gravitation'),
    (r'単振動|ばね振り子|振動の中心', 'simple-harmonic-motion'), (r'振り子', 'simple-pendulum'), (r'円運動|向心|遠心力', 'circular-motion'),
    (r'反発|衝突', 'restitution'), (r'運動量|力積|重心', 'momentum-conservation'),
    (r'力学的エネルギー|エネルギー保存', 'mechanical-energy'), (r'仕事', 'work-energy'),
    (r'モーメント|剛体|回転|転がり|重心のつり合い', 'torque'), (r'摩擦', 'static-friction'), (r'垂直抗力', 'normal-force'),
    (r'運動方程式|慣性力', 'equation-of-motion'), (r'つり合い', 'equilibrium'),
    (r'比熱|熱量|熱容量', 'heat-temperature'), (r'分子運動|気体分子|平均自由行程', 'kinetic-theory'),
    (r'熱効率|熱機関|サイクル|カルノー', 'heat-engine'), (r'熱力学第一|内部エネルギー|断熱|定圧|定積|等温|モル比熱', 'first-law'),
    (r'状態方程式|ボイル|シャルル|理想気体', 'ideal-gas'),
    (r'ドップラー', 'doppler'), (r'うなり', 'beats'), (r'気柱|弦|共鳴|固有振動', 'string-resonance'), (r'定常波|重ね合わせ', 'superposition'),
    (r'屈折|全反射|プリズム', 'refraction-light'), (r'レンズ|凹面鏡|凸面鏡|結像|写像', 'lens-imaging'),
    (r'ヤング|回折格子|薄膜|ニュートンリング|くさび|光の干渉|光路差', 'light-interference'), (r'干渉', 'interference'),
    (r'波長|波の式|位相|進行波|波形', 'wave-basics'),
    (r'コンデンサー|電気容量|誘電体|極板', 'capacitor'), (r'電場|電位|クーロン|静電気力|ガウス|導体', 'coulomb-field'),
    (r'ジュール|電力|消費電力', 'electric-power'), (r'キルヒホッフ|直流回路|電流計|電圧計|ホイートストン|抵抗の', 'series-parallel'),
    (r'ローレンツ|サイクロトロン|磁場中の荷電|らせん|ホール効果|電流が受ける力', 'lorentz-force'),
    (r'交流|インピーダンス|リアクタンス|変圧器|振動回路|共振', 'ac-circuit'),
    (r'電磁誘導|誘導起電力|自己誘導|相互誘導|レンツ|渦電流|インダクタンス', 'electromagnetic-induction'),
    (r'光電|光子|コンプトン|X線', 'photon'), (r'物質波|ド・ブロイ|ブラッグ|電子線', 'matter-wave'),
    (r'ボーア|エネルギー準位|水素原子|スペクトル', 'bohr-model'), (r'崩壊|半減期|放射', 'nuclear-decay'),
    (r'核反応|結合エネルギー|質量欠損|核分裂|核融合|E=mc', 'mass-energy'),
]
# 語が1つも当たらないときの、単元ごとの概念
UNIT_DEFAULT = {
    'MEC-01': 'uniform-acceleration', 'MEC-02': 'equilibrium', 'MEC-03': 'equation-of-motion', 'MEC-04': 'work-energy',
    'MEC-05': 'momentum-conservation', 'MEC-06': 'simple-harmonic-motion', 'MEC-07': 'universal-gravitation', 'MEC-08': 'torque',
    'MSY-01': 'equation-of-motion', 'THM-01': 'heat-temperature', 'THM-02': 'kinetic-theory', 'THM-03': 'first-law',
    'WAV-01': 'wave-basics', 'WAV-02': 'string-resonance', 'WAV-03': 'light-interference', 'ELM-01': 'coulomb-field',
    'ELM-02': 'capacitor', 'ELM-03': 'series-parallel', 'ELM-04': 'lorentz-force', 'ELM-05': 'electromagnetic-induction',
    'ELM-06': 'ac-circuit', 'ESY-01': 'electromagnetic-induction', 'ATM-01': 'photon', 'ATM-02': 'bohr-model',
}


def outside_math(s):
    """数式の外の位置だけを True にした配列（$…$・\\[…\\]・\\(…\\)）"""
    mask = [True] * len(s)
    for m in re.finditer(r'\$\$.*?\$\$|\\\[.*?\\\]|\\\(.*?\\\)|\$[^$]*\$', s, flags=re.S):
        for i in range(m.start(), m.end()):
            mask[i] = False
    return mask


def top_envs(text, env='enumerate'):
    """入れ子を考えて、いちばん外側の \\begin{env}…\\end{env} の (開始, 本文開始, 本文終了, 終了) を返す"""
    out, depth, start, body = [], 0, 0, 0
    for m in re.finditer(r'\\(begin|end)\{%s\}(\[[^\]]*\])?' % env, text):
        if m.group(1) == 'begin':
            if depth == 0:
                start, body = m.start(), m.end()
            depth += 1
        else:
            depth -= 1
            if depth == 0:
                out.append((start, body, m.start(), m.end()))
    return out


def top_items(body):
    """入れ子の enumerate/itemize の中の \\item は分けない"""
    cuts, depth = [], 0
    for m in re.finditer(r'\\(begin|end)\{(?:enumerate|itemize)\}|\\item\b', body):
        if m.group(1) == 'begin':
            depth += 1
        elif m.group(1) == 'end':
            depth -= 1
        elif depth == 0:
            cuts.append(m)
    return [body[m.end():(cuts[k + 1].start() if k + 1 < len(cuts) else len(body))].strip() for k, m in enumerate(cuts)]


def split_statement(text):
    """問題文 → blocks（text / parts）。外側の enumerate の項目を小問にする。\\setcounter{enumi}{k} で番号を続ける"""
    blocks, num, pos = [], 0, 0
    SC = r'\\setcounter\{enumi\}\{(\d+)\}'
    for (st, bs, be, en) in top_envs(text):
        before, body = text[pos:st], text[bs:be]
        head = body.split('\\item', 1)[0]
        sc = re.findall(SC, before) + re.findall(SC, head)
        before = re.sub(SC, '', before).strip()
        if before:
            blocks.append({'type': 'text', 'tex': before})
        if sc:
            num = int(sc[-1])
        parts = []
        for it in top_items(re.sub(SC, '', body)):
            lab = re.match(r'\[([^\]]*)\]\s*', it)   # \item[(3)] は小問 3、\item[I.] は見出し（文章として入れる）
            if lab and not re.fullmatch(r'\(?\d+\)?\.?', lab.group(1).strip()):
                if parts:
                    blocks.append({'type': 'parts', 'items': parts})
                    parts = []
                blocks.append({'type': 'text', 'tex': f'\\textbf{{{lab.group(1).strip()}}}　' + it[lab.end():]})
                continue
            if lab:
                num = int(re.sub(r'\D', '', lab.group(1)))
                it = it[lab.end():]
            else:
                num += 1
            parts.append({'id': str(num), 'tex': it})
        if parts:
            blocks.append({'type': 'parts', 'items': parts})
        pos = en
    rest = re.sub(SC, '', text[pos:]).strip()
    if rest:
        blocks.append({'type': 'text', 'tex': rest})
    return blocks


def split_answer(ans, ids):
    """"(1) … (2) …" → 小問ごとの答え。小問番号は先頭か空白の直後にあるものだけを区切りとみなす。
    番号が小問と合わなければ 1つにまとめる"""
    mask = outside_math(ans)
    cands = [m for m in re.finditer(r'(?:^|(?<=[\s　]))\((\d+)\)', ans) if mask[m.start()]]
    cuts, k = [], 0
    for want in ids:
        while k < len(cands) and cands[k].group(1) != want:
            k += 1
        if k == len(cands):
            break
        cuts.append(cands[k])
        k += 1
    if ids and len(cuts) == len(ids) and cuts[0].start() == len(ans) - len(ans.lstrip()):
        out = []
        for j, m in enumerate(cuts):
            end = cuts[j + 1].start() if j + 1 < len(cuts) else len(ans)
            out.append({'id': m.group(1), 'tex': ans[m.end():end].strip().strip('　').strip()})
        return out, True
    return [{'id': 'main', 'tex': ans.strip()}], False


def prefix_svg_ids(svg, pre):
    ids = set(re.findall(r'\bid="([^"]+)"', svg))
    for i in sorted(ids, key=len, reverse=True):
        svg = re.sub(r'\bid="%s"' % re.escape(i), f'id="{pre}-{i}"', svg)
        svg = svg.replace(f'url(#{i})', f'url(#{pre}-{i})').replace(f'href="#{i}"', f'href="#{pre}-{i}"')
    return svg


def concepts_for(item, unit):
    text = ' '.join([item.get('tags', ''), item.get('theme', ''), item.get('core', '')])
    out = []
    for pat, cid in KEYWORDS:
        if re.search(pat, text) and cid not in out:
            out.append(cid)
    if unit.startswith('WAV-03') and 'interference' in out and 'light-interference' not in out:
        out[out.index('interference')] = 'light-interference'
    d = UNIT_DEFAULT[unit]
    if d in out:   # 単元の代表概念を先頭に
        out.remove(d)
    return ([d] + out)[:3]


def load_sheet(path):
    if not path:
        return {}
    import openpyxl
    wb = openpyxl.load_workbook(path, read_only=True, data_only=True)
    ws = wb['問題一覧']
    rows = list(ws.iter_rows(values_only=True))
    head = rows[0]
    col = {h: i for i, h in enumerate(head) if h}
    out = {}
    for r in rows[1:]:
        pid = r[col['ID']]
        if not pid or str(r[col.get('メモ', 0)] or '').startswith('記入例'):
            continue
        out[str(pid).strip()] = {h: r[i] for h, i in col.items()}
    return out


def main(argv):
    if not argv or argv[0] in ('-h', '--help'):
        print(__doc__)
        return 0
    src = Path(argv[0])
    sheet = verify = None
    if '--sheet' in argv:
        sheet = argv[argv.index('--sheet') + 1]
    if '--verify' in argv:
        verify = Path(argv[argv.index('--verify') + 1])
    S = load_sheet(sheet)
    tax = json.loads((CONTENT / 'taxonomy.json').read_text(encoding='utf-8'))
    units = {u['name']: u for u in tax['units']}
    concepts = {c['id'] for c in json.loads((CONTENT / 'concepts.json').read_text(encoding='utf-8'))['concepts']}
    n, merged, problems = 0, [], []
    for d in sorted(p for p in src.iterdir() if (p / 'item.json').exists()):
        it = json.loads((d / 'item.json').read_text(encoding='utf-8'))
        pid = it['id']
        if not re.fullmatch(r'[A-Z]{3}-[A-Z]{3}-\d{3,}', pid) or pid != d.name:
            raise SystemExit(f'{d.name}: 問題ID {pid!r} が規則（ABC-DEF-001）に合わないか、フォルダ名と違います')
        uname = ALIASES.get(it['field'], it['field'])
        u = units.get(uname) or next((v for k, v in units.items() if ALIASES.get(k) == uname), None)
        if not u:
            raise SystemExit(f'{pid}: 単元「{it["field"]}」が taxonomy.json にありません')
        blocks = split_statement(it['problem'])
        ids = [x['id'] for b in blocks if b['type'] == 'parts' for x in b['items']]
        parts, ok = split_answer(it['answer'], ids)
        if not ok and ids:   # 答えを小問に分けられない → 小問は本文のまま（text）にする
            merged.append(pid)
            blocks = [{'type': 'text', 'tex': it['problem']}]
        figs, sol_blocks = [], []
        out = CONTENT / 'problems' / u['field'] / pid
        if out.exists():
            shutil.rmtree(out)
        out.mkdir(parents=True)
        if it.get('problem_fig'):
            svg = prefix_svg_ids((d / it['problem_fig']).read_text(encoding='utf-8'), 'q')
            (out / 'q-fig.svg').write_text(svg, encoding='utf-8')
            figs.append({'id': 'q-fig', 'role': 'question', 'kind': 'diagram', 'format': 'svg', 'src': 'q-fig.svg',
                         'alt': f'問題図：{it["theme"]}', 'caption': ''})
            k = 1 if blocks and blocks[0]['type'] == 'text' else 0
            blocks.insert(k, {'type': 'figure', 'ref': 'q-fig'})
        if it.get('solution_fig'):
            svg = prefix_svg_ids((d / it['solution_fig']).read_text(encoding='utf-8'), 'a')
            (out / 'a-fig.svg').write_text(svg, encoding='utf-8')
            figs.append({'id': 'a-fig', 'role': 'answer', 'kind': 'diagram', 'format': 'svg', 'src': 'a-fig.svg',
                         'alt': f'解答図：{it["theme"]}', 'caption': '解答図'})
            sol_blocks.append({'type': 'figure', 'ref': 'a-fig'})
        sol_blocks.append({'type': 'text', 'tex': it['solution']})
        row = S.get(pid, {})
        done = lambda k: str(row.get(k) or '').strip() == '済'   # noqa: E731
        created = row.get('作成日')
        created = created.strftime('%Y-%m-%d') if hasattr(created, 'strftime') else '2026-11-07'
        cs = [c for c in concepts_for(it, u['key']) if c in concepts]
        syn = u['field'] in ('MSY', 'ESY')
        p = {
            'id': pid, 'schema': 'celestia.problem/1', 'field': u['field'], 'unit': u['key'], 'topic': it['theme'], 'course': '物理',
            'level': LEVEL[int(it['difficulty'])], 'minutes': int(it.get('time_min') or 20),
            'format': ['derivation', 'synthesis'] if syn else ['derivation'], 'examStyle': 'national-2nd',
            'statement': {'blocks': blocks}, 'figures': figs,
            'answer': {'parts': parts},
            'solutions': [{'id': 'main', 'kind': 'main', 'blocks': sol_blocks}],
            'idea': {'tex': it['core']} if it.get('core') else None,
            'concepts': cs, 'laws': [], 'misconceptions': [],
            'tags': [t.strip() for t in re.split(r'[,，、]', it.get('tags') or '') if t.strip()],
            'provenance': {'created': created, 'method': str(row.get('作成方法') or 'AI生成') + '（問題管理表から取り込み）'},
            'quality': {'status': 'published',
                        'checks': {'schema': True, 'units': False, 'figures': False,
                                   'physics': bool(it.get('verified')) or done('検証:SymPy等'), 'human': done('検証:人間確認')},
                        'verified': all(done(k) for k in ('検証:二重解答', '検証:SymPy等', '検証:人間確認'))},
        }
        if p['idea'] is None:
            del p['idea']
        (out / 'problem.json').write_text(json.dumps(p, ensure_ascii=False, indent=1) + '\n', encoding='utf-8')
        if verify and (verify / f'{pid}.py').exists():
            (ROOT / 'verify').mkdir(exist_ok=True)
            shutil.copy(verify / f'{pid}.py', ROOT / 'verify' / f'{pid}.py')
        problems.append(pid)
        n += 1
    print(f'{n}問を取り込みました')
    if merged:
        print(f'答えを小問に分けられず、1つにまとめた問題 {len(merged)}問：' + ' '.join(merged))
    return 0


if __name__ == '__main__':
    sys.exit(main(sys.argv[1:]))
