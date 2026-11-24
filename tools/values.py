#!/usr/bin/env python3
"""頻出度（content/values.json）の点検と更新

    python3 tools/values.py                     点検：すべての公開問題に今の問題文に対する判定があるか（無い・古いものがあれば終了コード 1）
    python3 tools/values.py --stale             判定が無い・古い（問題文が変わった）問題を、判定に必要な情報つきで一覧
    python3 tools/values.py --set ID=3 ID=1 …   判定を書く（今の問題文のハッシュを付ける）
    python3 tools/values.py --import FILE       「問題ID 頻出度」の行を並べたファイルからまとめて書く
    python3 tools/values.py --stats             頻出度の分布（分野別）

頻出度 v：1=低・2=中・3=高。大学入試でその問題の中心の考え方・手法がどれだけ繰り返し出るか（判定の基準は docs/data-spec.md）。
content/values.json：{ 問題ID: { "h": 問題文のハッシュ（content_lib.statement_hash）, "v": 1〜3 } }（Lemma と同じ形）
problem.json に "value" を書いた問題はそれが優先され、values.json の判定は使われない（無くてもよい）。
問題文を直すとハッシュが変わり、その判定は使われなくなる（build_data.py は v を出さない＝中の扱い）。
--stale で見つけて、判定し直して --set で書く。文の手直しだけで判定が変わらないときも --set で書き直す。
"""
import json
import sys
from collections import Counter
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from content_lib import CONTENT, load_problems, statement_hash  # noqa: E402

PATH = CONTENT / 'values.json'
NAME = {1: '低', 2: '中', 3: '高'}


def load():
    return json.loads(PATH.read_text(encoding='utf-8')) if PATH.exists() else {}


def save(values):
    PATH.write_text(json.dumps(dict(sorted(values.items())), ensure_ascii=False, indent=0) + '\n', encoding='utf-8')


def published():
    return [p for _, p in load_problems() if p['quality']['status'] == 'published']


def status(values, probs):
    """(判定のない問題, 古い判定, 存在しない問題の判定, 値の誤り)"""
    ids = {p['id'] for p in probs}
    missing, stale = [], []
    for p in probs:
        if p.get('value'):
            continue
        j = values.get(p['id'])
        if not j:
            missing.append(p)
        elif j.get('h') != statement_hash(p):
            stale.append(p)
    orphan = sorted(set(values) - ids)
    bad = sorted(i for i, j in values.items() if j.get('v') not in (1, 2, 3))
    return missing, stale, orphan, bad


def effective(values, p):
    """build_data.py と同じ決め方：problem.json の value、なければ問題文が同じときの判定、なければ None（中）"""
    if p.get('value'):
        return p['value']
    j = values.get(p['id'])
    return j['v'] if j and j.get('h') == statement_hash(p) else None


def set_many(pairs):
    probs = {p['id']: p for p in published()}
    values = load()
    for i, v in pairs:
        if i not in probs:
            sys.exit(f'問題 {i} がありません')
        if v not in (1, 2, 3):
            sys.exit(f'{i}: 頻出度は 1〜3 です（{v}）')
        values[i] = {'h': statement_hash(probs[i]), 'v': v}
    save(values)
    print(f'{len(pairs)}問の頻出度を書きました → {PATH.relative_to(CONTENT.parent)}')


def main():
    args = sys.argv[1:]
    if args and args[0] == '--set':
        set_many([(a.split('=')[0], int(a.split('=')[1])) for a in args[1:]])
        return
    if args and args[0] == '--import':
        pairs = []
        for line in Path(args[1]).read_text(encoding='utf-8').splitlines():
            line = line.split('#')[0].strip()
            if line:
                i, v = line.split()[:2]
                pairs.append((i, int(v)))
        set_many(pairs)
        return
    probs = published()
    values = load()
    if args and args[0] == '--stats':
        by = {}
        for p in probs:
            by.setdefault(p['field'], Counter())[effective(values, p) or 2] += 1
        tot = Counter()
        for f, c in by.items():
            tot += c
            print(f'{f}：' + '・'.join(f'{NAME[k]} {c[k]}' for k in (3, 2, 1)))
        print('全体：' + '・'.join(f'{NAME[k]} {tot[k]}' for k in (3, 2, 1)))
        return
    missing, stale, orphan, bad = status(values, probs)
    if args and args[0] == '--stale':
        for p in missing + stale:
            j = values.get(p['id'])
            print('\t'.join([p['id'], '判定なし' if not j else f'古い（前の判定 {j["v"]}）', f'L{p["level"]}', p.get('topic', ''),
                             '・'.join(p.get('tags', [])), '・'.join(p.get('laws', []))]))
        return
    overridden = sum(1 for p in probs if p.get('value'))
    print(f'頻出度：公開 {len(probs)}問のうち 判定済み {len(probs) - len(missing) - len(stale) - overridden}問・problem.json で指定 {overridden}問、'
          f'判定なし {len(missing)}問・問題文が変わって古い {len(stale)}問、使われない判定 {len(orphan)}件、値の誤り {len(bad)}件')
    for p in (missing + stale)[:20]:
        print('  要判定', p['id'], p.get('topic', ''))
    for i in orphan[:10]:
        print('  問題が無い判定', i)
    for i in bad[:10]:
        print('  値の誤り', i, values[i])
    if missing or stale:
        print('  → python3 tools/values.py --stale で一覧し、判定して python3 tools/values.py --set ID=頻出度 で書く')
    sys.exit(1 if (missing or stale or orphan or bad) else 0)


if __name__ == '__main__':
    main()
