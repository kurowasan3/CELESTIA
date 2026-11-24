// ブラウザでの通し操作（Playwright・Chromium、iPhone 幅 390×844）
//   node tests/e2e.mjs
// 主な画面、図（解答図は答え合わせまで読み込まない）、一部のデータだけ古いときの警告、オフライン、同期を確かめる。
import fs from 'fs';
import http from 'http';
import path from 'path';
import { fileURLToPath } from 'url';
import { createMockSupabase } from './mock-supabase.mjs';

import { launchBrowser } from './browser.mjs';
const root = path.resolve(path.dirname(fileURLToPath(import.meta.url)), '..');
const TYPES = { '.html': 'text/html', '.js': 'text/javascript', '.mjs': 'text/javascript', '.css': 'text/css', '.json': 'application/json',
  '.png': 'image/png', '.webmanifest': 'application/manifest+json', '.woff2': 'font/woff2', '.svg': 'image/svg+xml' };
let fakeDataVersion = null;
const mock = createMockSupabase({ key: 'sb_publishable_e2e' });   // 同期のテスト用（このサーバーの /mock/ で応える）   // Service Worker のテスト用：index.json の dataVersion を差し替える
const server = http.createServer((req, res) => {
  let p = decodeURIComponent(new URL(req.url, 'http://x').pathname);
  if (p.endsWith('/')) p += 'index.html';
  if (p.startsWith('/mock/')) {
    let body = ''; req.on('data', c => { body += c; });
    req.on('end', async () => {
      const r = await mock.handle({ method: req.method, url: req.url, headers: req.headers, body: body || null });
      res.writeHead(r.status, { 'Content-Type': 'application/json' }); res.end(r.body);
    });
    return;
  }
  if (p === '/js/config.js') {
    res.writeHead(200, { 'Content-Type': 'text/javascript' });
    const v = fs.readFileSync(path.join(root, 'js/config.js'), 'utf8').match(/FILE_VERSION = '([^']+)'/)[1];
    return res.end(`export const FILE_VERSION = '${v}';
export const SUPABASE_URL = 'http://localhost:${server.address().port}/mock';
export const SUPABASE_KEY = 'sb_publishable_e2e';`);
  }
  if (fakeDataVersion && p === '/data/index.json') {
    const d = JSON.parse(fs.readFileSync(path.join(root, p), 'utf8')); d.dataVersion = fakeDataVersion;
    res.writeHead(200, { 'Content-Type': 'application/json' }); return res.end(JSON.stringify(d));
  }
  const f = path.join(root, p);
  if (!f.startsWith(root) || !fs.existsSync(f)) { res.writeHead(404); return res.end(); }
  res.writeHead(200, { 'Content-Type': TYPES[path.extname(f)] || 'application/octet-stream' });
  fs.createReadStream(f).pipe(res);
});
await new Promise(r => server.listen(0, r));
const U = `http://localhost:${server.address().port}/`;

let fail = 0, pass = 0;
const ok = (cond, msg) => { if (cond) pass++; else { fail++; console.log('FAIL', msg); } };
const browser = await launchBrowser();

async function newPage({ state = null, raw = null, route = null, sw = false } = {}) {
  const ctx = await browser.newContext({ viewport: { width: 390, height: 844 }, acceptDownloads: true, serviceWorkers: sw ? 'allow' : 'block' });
  if (state || raw) await ctx.addInitScript((s) => { if (!sessionStorage.getItem('seeded')) { sessionStorage.setItem('seeded', '1'); localStorage.setItem('celestia.v1', s); } }, raw || JSON.stringify(state));
  const page = await ctx.newPage();
  page.errors = [];
  page.on('pageerror', e => page.errors.push(String(e)));
  page.on('console', m => { if (m.type() === 'error') page.errors.push(m.text()); });
  page.on('dialog', d => d.accept());
  page.requested = [];
  page.on('request', r => page.requested.push(new URL(r.url()).pathname));
  if (route) await route(page);
  return page;
}
const next = (page) => page.locator('button:has-text("次")').first().click();

try {
  const fixture = JSON.parse(fs.readFileSync(path.join(root, 'tests/fixtures/state-v1.json'), 'utf8'));
  const LS = (page) => page.evaluate(() => JSON.parse(localStorage.getItem('celestia.v1')));
  // 1. 起動・タブ・図つきの問題・解答図は答え合わせまで出ない・入力採点
  {
    const page = await newPage();
    await page.goto(U + '#today'); await page.waitForTimeout(1200);
    const tabs = await page.locator('.tabbar a span').allInnerTexts();
    ok(tabs.join() === '今日,シミュレーション,概念,定着度,設定', `タブが物理版の構成 ${tabs}`);
    ok(await page.locator('#lv button').count() === 6, '難易度ボタンが6つ');
    ok(await page.locator('.banner').count() === 0, '起動時に警告が出ない');
    await page.evaluate("__celestia.open(['MEC-INC-001', 'ELM-DCC-001'])"); await page.waitForTimeout(1200);
    ok(await page.locator('.question figure.pfig svg').count() === 1, '問題図が本文の中に SVG で入る');
    ok(await page.locator('.question figure.pfig[aria-label], .question figure.pfig [role="img"]').count() === 1, '問題図に代替テキスト');
    ok(await page.locator('figure.fig-answer').count() === 0, '答え合わせ前は解答図が無い');
    ok(!page.requested.some(u => u.includes('a-forces')), '答え合わせ前は解答図のファイルを読み込まない');
    await page.locator('#conf button').nth(2).click();
    const inputs = page.locator('#inputs input');
    ok(await inputs.count() === 3, '小問3つに入力欄');
    await inputs.nth(0).fill('mgcos(θ)'); await inputs.nth(1).fill('m g sinθ'); await inputs.nth(2).fill('μ>=tanθ');
    await page.click('#judge'); await page.waitForTimeout(800);
    ok((await page.locator('.verdict').innerText()).includes('完答'), 'MEC-INC-001 を正しく答えて完答');
    await page.waitForTimeout(500);
    ok(await page.locator('figure.fig-answer svg').count() === 1, '答え合わせのあとに解答図が出る');
    ok(page.requested.some(u => u.includes('a-forces')), '解答図のファイルは答え合わせのあとに読む');
    ok((await page.locator('figure.fig-answer figcaption').innerText()).includes('分解'), '解答図の説明');
    ok((await page.locator('.concept-box').innerText()).includes('垂直抗力'), '解説の下に確かめる概念・よくある誤解');
    await next(page); await page.waitForTimeout(1000);
    ok(await page.locator('.question figure.pfig svg').count() === 1, '回路図の問題');
    ok((await page.locator('#inputs .unit').allInnerTexts()).join() === 'Ω,A,A', '単位 [Ω] [A] が入力欄の横に出る');
    await page.locator('#conf button').nth(1).click();
    await page.locator('#inputs input').nth(0).fill('6'); await page.locator('#inputs input').nth(1).fill('3'); await page.locator('#inputs input').nth(2).fill('4/3');
    await page.click('#judge'); await page.waitForTimeout(400);
    ok((await page.locator('#feedback').innerText()).includes('2つ正解'), '1つだけ違うと「3つのうち2つ正解」');
    const st = await LS(page);
    ok(st.pstat['MEC-INC-001'] && st.log.some(l => l.p === 'MEC-INC-001'), '結果が問題ID で記録される');
    ok(page.errors.length === 0, `ページのエラーなし ${page.errors}`);
  }
  // 2. 自己採点の問題・シミュレーション・概念・設定・定着度
  {
    const page = await newPage();
    await page.goto(U + '#today'); await page.waitForTimeout(1200);
    await page.evaluate("__celestia.open(['MEC-CON-001'])"); await page.waitForTimeout(1000);
    await page.locator('#conf button').nth(1).click();
    await page.click('#reveal'); await page.waitForTimeout(500);
    const order = await page.evaluate(() => { const a = document.querySelector('#solbox'), b = document.querySelector('#self'); return !!(a && b && (a.compareDocumentPosition(b) & Node.DOCUMENT_POSITION_FOLLOWING)); });
    ok(order, '自己採点の問題で、解答・解説が「どうだった？」より前に出る');
    await page.goto(U + '#sim'); await page.waitForTimeout(800);
    ok(await page.locator('.sim-card').count() === 24, 'シミュレーション一覧（登録簿の24件）');
    await page.goto(U + '#concepts'); await page.waitForTimeout(800);
    ok(await page.locator('details.concept').count() === 0 && await page.locator('.view a.law, #view a.law, a.law').count() === 68, '「概念」タブは法則の一覧だけ（68件）');
    await page.goto(U + '#concept=normal-force'); await page.waitForTimeout(800);
    ok(page.url().endsWith('#concepts') && await page.locator('a.law').count() === 68, '古い #concept= のリンクは法則の一覧へ');
    await page.goto(U + '#drill'); await page.waitForTimeout(500);
    ok(await page.locator('.drill-card').count() === 0, '計算トレーニングは無い');
    await page.goto(U + '#settings'); await page.waitForTimeout(600);
    await page.locator('.segment.theme button:has-text("ダーク")').click();
    ok(await page.evaluate(() => document.documentElement.dataset.theme) === 'dark', 'ダークモード');
    await page.goto(U + '#map'); await page.waitForTimeout(600);
    ok((await page.locator('#view').innerText()).includes('定着度'), '定着度画面');
    ok(page.errors.length === 0, `ページのエラーなし ${page.errors}`);
  }
  // 3. 学習記録のある端末：CSV・続きに記録
  {
    const page = await newPage({ state: fixture });
    await page.goto(U + '#settings'); await page.waitForTimeout(1200);
    const dl = page.waitForEvent('download');
    await page.click('#csvall');
    const csv = fs.readFileSync(await (await dl).path(), 'utf8');
    ok(csv.split('\n').filter(l => l.trim()).length === fixture.log.length + 1, `CSV に記録 ${fixture.log.length}件がすべて出る`);
    await page.evaluate("__celestia.open(['MEC-SHM-001'])"); await page.waitForTimeout(1000);
    await page.locator('#conf button').nth(1).click();
    const inp = page.locator('#inputs input');
    await inp.nth(0).fill('2π√(m/k)'); await inp.nth(1).fill('A sqrt(k/m)'); await inp.nth(2).fill('√2');
    await page.click('#judge'); await page.waitForTimeout(600);
    ok((await page.locator('.verdict').innerText()).includes('完答'), '添字（v_max）や √ を含む答えを採点できる');
    const st = await LS(page);
    ok(Object.keys(fixture.pstat).every(k => st.pstat[k]), '前の復習状態がすべて残る');
    ok(st.log.length === fixture.log.length + 1, '学習記録が続きに足される');
    ok(st.predictions.length === fixture.predictions.length && st.misc['MC-SHM-PERIOD-LINEAR-M'], '予測・誤概念の記録が残る');
    ok(page.errors.length === 0, `ページのエラーなし ${page.errors}`);
  }
  // 4. 一部の data/ だけ古い（f-MEC-02.json の版が違う）
  {
    const page = await newPage({ route: (p) => p.route('**/data/f-MEC-02.json', async (r) => {
      const res = await r.fetch(); const d = await res.json(); d._meta = { dataVersion: 'old000000000' };
      await r.fulfill({ response: res, body: JSON.stringify(d) });
    }) });
    await page.goto(U + '#today'); await page.waitForTimeout(1200);
    await page.evaluate("__celestia.open(['MEC-INC-001'])"); await page.waitForTimeout(800);
    const b = page.locator('.banner');
    ok(await b.count() === 1 && (await b.innerText()).includes('data/f-MEC-02.json'), '古い本文ファイルを読むと警告が出る');
    ok(await page.locator('#inputs input').count() > 0, '警告が出ても問題は解ける');
  }
  // 5. 検索・問題セット・共有リンク
  {
    const page = await newPage({ state: fixture });
    await page.goto(U + '#today'); await page.waitForTimeout(1500);
    ok(await page.locator('#sets .set-card').count() >= 1, '今日の画面に問題セット');
    await page.fill('#sq', '垂直抗力'); await page.press('#sq', 'Enter'); await page.waitForTimeout(1200);
    ok(await page.locator('.res-row').count() >= 2, '概念名で検索できる');
    await page.goto(U + '#set=col-rikigaku-hajime'); await page.waitForTimeout(1500);
    ok(page.url().endsWith('#today') && await page.locator('#conf button').count() === 3, '#set= のリンクで問題セットが始まる');
    await page.goto(U + '#ids=MEC-PRJ-001'); await page.waitForTimeout(1500);
    ok((await page.locator('.chip').allInnerTexts()).some(t => t.includes('力学・運動の表し方')), '#ids= のリンクで指定の問題が出る');
    ok(page.errors.length === 0, `ページのエラーなし ${page.errors}`);
  }
  // 6. Service Worker：自分のキャッシュだけを管理し、同じドメインの他のアプリ（Lemma など）のキャッシュは消さない。オフラインで図も出る
  {
    const page = await newPage({ sw: true });
    await page.goto(U + '#today');
    await page.evaluate(async () => { await (await caches.open('lemma-v0.9.0')).put('/x', new Response('lemma')); await (await caches.open('butsuriba-v0.4.0')).put('/y', new Response('b')); });
    await page.waitForFunction(() => navigator.serviceWorker.controller, null, { timeout: 15000 }).catch(() => {});
    await page.reload(); await page.waitForTimeout(1200);
    await page.evaluate("__celestia.open(['MEC-INC-001'])"); await page.waitForTimeout(1500);
    const v = JSON.parse(fs.readFileSync(path.join(root, 'data/index.json'), 'utf8')).dataVersion;
    const keys = await page.evaluate(async () => caches.keys());
    ok(keys.includes('celestia-data-' + v) && keys.includes('celestia-v0.6.0'), `コード用とデータ用のキャッシュ ${keys}`);
    ok(keys.includes('lemma-v0.9.0') && keys.includes('butsuriba-v0.4.0'), '同じドメインの他のアプリのキャッシュを消さない');
    await page.context().setOffline(true);
    await page.goto(U + '#today'); await page.waitForTimeout(1200);
    await page.evaluate("__celestia.open(['MEC-INC-001'])"); await page.waitForTimeout(1500);
    ok(await page.locator('#inputs input').count() === 3 && await page.locator('.question figure.pfig svg').count() === 1, 'オフラインでも開いた問題が図つきで解ける');
    await page.goto(U + '#sim=mech.pendulum'); await page.waitForTimeout(1500);
    ok(await page.locator('.sp-canvas').count() === 1 && await page.locator('.timeplot canvas').count() >= 2, 'オフラインでもシミュレーションが動く');
    await page.context().setOffline(false);
    fakeDataVersion = 'e2enewversion';
    await page.goto(U + '#today'); await page.waitForTimeout(1500);
    await page.evaluate(async () => fetch('data/index.json', { cache: 'no-cache' }));
    await page.waitForTimeout(500);
    const keys2 = await page.evaluate(async () => caches.keys());
    ok(keys2.includes('celestia-data-e2enewversion') && !keys2.includes('celestia-data-' + v) && keys2.includes('lemma-v0.9.0'), `データの版が変わると自分の古いデータ用キャッシュだけを消す ${keys2}`);
    fakeDataVersion = null;
  }
  // 8. シミュレーション・グラフの図・法則の導出
  {
    const page = await newPage();
    await page.goto(U + '#sim'); await page.waitForTimeout(1200);
    ok(await page.locator('a.sim-card').count() === 24, 'シミュレーション一覧に24本（力学7・波動5・光3・電磁気5・熱2・原子2）');
    for (const id of ['mech.projectile', 'mech.incline', 'mech.spring', 'mech.pendulum', 'mech.collision', 'mech.circular', 'mech.gravitation']) {
      await page.goto(U + '#sim=' + id); await page.waitForTimeout(900);
      const ok1 = await page.locator('.sp-canvas').count() === 1 && await page.locator('.timeplot canvas').count() >= 2;
      await page.click('.sp-play'); await page.waitForTimeout(700); await page.click('.sp-play');
      const t = parseFloat(await page.locator('.sp-time').innerText());
      ok(ok1 && t > 0, `${id}：アニメーションとグラフが出て、再生すると時刻が進む（${t} s）`);
    }
    // 他分野：時間で動くものは再生で時刻が進み、時間によらないもの（静的）は再生バーがなく分布のグラフが出る
    for (const id of ['wave.traveling', 'wave.standing', 'wave.beats', 'wave.doppler', 'wave.interference', 'elm.rc', 'elm.lorentz', 'elm.induction', 'elm.rlc', 'thermo.molecules', 'thermo.pv', 'atom.decay']) {
      await page.goto(U + '#sim=' + id); await page.waitForTimeout(800);
      await page.click('.sp-play'); await page.waitForTimeout(600); await page.click('.sp-play');
      const t = parseFloat(await page.locator('.sp-time').innerText());
      ok(t > 0 && await page.locator('.timeplot canvas').count() >= 1, `${id}：再生すると時刻が進み、グラフが出る`);
    }
    for (const id of ['optics.refraction', 'optics.lens', 'optics.double-slit', 'elm.charges', 'atom.photoelectric']) {
      await page.goto(U + '#sim=' + id); await page.waitForTimeout(800);
      ok(await page.locator('.sp-bar:visible').count() === 0 && await page.locator('.timeplot canvas').count() >= 1 && await page.locator('.sp-canvas').count() === 1, `${id}：時間によらない（再生バーなし）で、図と分布のグラフが出る`);
    }
    await page.goto(U + '#sim=optics.refraction'); await page.waitForTimeout(800);
    await page.locator('input[data-k="theta1"]').fill('60'); await page.locator('input[data-k="theta1"]').dispatchEvent('input'); await page.waitForTimeout(300);
    ok((await page.locator('.sp-warn').innerText()).includes('全反射'), 'n₁=1.5→n₂=1.0 で入射角 60° にすると全反射の表示');
    await page.goto(U + '#sim=mech.spring'); await page.waitForTimeout(900);
    const before = await page.locator('.sp-measures td').first().innerText();
    await page.locator('input[data-k="m"]').fill('2'); await page.locator('input[data-k="m"]').dispatchEvent('input'); await page.waitForTimeout(300);
    const after = await page.locator('.sp-measures td').first().innerText();
    ok(parseFloat(after) / parseFloat(before) > 1.40 && parseFloat(after) / parseFloat(before) < 1.43, `質量を2倍にすると周期が √2 倍（${before} → ${after}）`);
    await page.goto(U + '#today'); await page.waitForTimeout(800);
    await page.evaluate("__celestia.open(['MEC-PRJ-003'])"); await page.waitForTimeout(1500);
    ok(await page.locator('.question figure.pfig .timeplot canvas').count() === 1, 'グラフ読解の問題に、モデルから描いた v-t グラフ（sim-plot）');
    const inp = page.locator('#inputs input');
    await page.locator('#conf button').nth(1).click();
    for (const [i, v] of ['10', '1', '10', '5', '20√3'].entries()) await inp.nth(i).fill(v);
    await page.click('#judge'); await page.waitForTimeout(800);
    ok((await page.locator('.verdict').innerText()).includes('完答'), 'グラフ読解の問題を単位つきで完答');
    await page.locator('.sim-open').click(); await page.waitForTimeout(1200);
    ok(await page.locator('.sim-inline .sp-canvas').count() === 1 && await page.locator('.sim-inline input[data-k="g"]').count() === 1, '解いたあとに「シミュレーションで確かめる」で問題の条件のシミュレーションが開く');
    ok((await page.locator('.concept-box').innerText()).includes('使う法則'), '解説の下に使う法則へのリンク');
    // 波形（時刻ごとに形が変わるグラフ）の図・全分野の問題
    await page.goto(U + '#today'); await page.waitForTimeout(800);
    await page.evaluate("__celestia.open(['WAV-BAS-001', 'WAV-REF-001'])"); await page.waitForTimeout(1500);
    ok(await page.locator('.question figure.pfig .timeplot canvas').count() === 1, '波形のグラフ（y-x）の図が問題に出る');
    await page.locator('#conf button').nth(1).click();
    for (const [i, v] of ['0.2', '2', '1', '2'].entries()) await page.locator('#inputs input').nth(i).fill(v);
    await page.click('#judge'); await page.waitForTimeout(700);
    ok((await page.locator('.verdict').innerText()).includes('完答'), '波形のグラフを読む問題を完答');
    await next(page); await page.waitForTimeout(1000);
    ok(await page.locator('.question figure.pfig svg').count() === 1, '屈折の光線の図');
    await page.evaluate("__celestia.open(['ELM-ACC-001'])"); await page.waitForTimeout(1200);
    await page.locator('#conf button').nth(1).click();
    for (const [i, v] of ['2πfL', '1/(2πfC)', 'sqrt(R^2+(2πfL-1/(2πfC))^2)', '1/(2π√(LC))', 'V_0/R'].entries()) await page.locator('#inputs input').nth(i).fill(v);
    await page.click('#judge'); await page.waitForTimeout(700);
    ok((await page.locator('.verdict').innerText()).includes('完答'), '交流の問題：C（電気容量）を含む式を正しく採点（積分定数として扱わない）');
    await page.goto(U + '#law=centripetal-acceleration'); await page.waitForTimeout(1000);
    ok(await page.locator('ol.derivation li.open').count() === 1, '法則のページ：導出は1段目だけ開いている');
    await page.click('#dnext'); await page.waitForTimeout(200);
    ok(await page.locator('ol.derivation li.open').count() === 2, '「次の段を見る」で1段ずつ開く');
    await page.click('#dall'); await page.waitForTimeout(200);
    ok(await page.locator('ol.derivation li:not(.open)').count() === 0 && await page.locator('.law-page .katex').count() > 5, 'すべて表示で導出が全部出て、数式が描ける');
    await page.goto(U + '#law=newton-2nd'); await page.waitForTimeout(800);
    ok((await page.locator('.law-page').innerText()).includes('この法則の根拠') && await page.locator('ol.derivation').count() === 0, '出発点となる法則は導出ではなく根拠を出す');
    await page.goto(U + '#concepts'); await page.waitForTimeout(800);
    await page.locator('a.law').first().click(); await page.waitForTimeout(800);
    ok(page.url().includes('#law='), '法則の一覧から法則のページへ');
    ok(page.errors.length === 0, `ページのエラーなし ${page.errors}`);
  }
  // 7. アカウントと同期：2台の端末でログインすると記録が合わさる（にせ Supabase）
  {
    const login = async (page) => {
      await page.goto(U + '#settings'); await page.waitForTimeout(1500);
      await page.fill('#acct-email', 'me@example.com'); await page.click('#acct-send'); await page.waitForTimeout(500);
      await page.fill('#acct-code', '123456'); await page.click('#acct-verify'); await page.waitForTimeout(1500);
    };
    const A = await newPage({ state: { ...fixture, log: fixture.log.slice(0, 25) } });
    await login(A);
    ok((await A.locator('#acct').innerText()).includes('Lemma と共通'), 'アカウントが Lemma と共通であることを表示');
    const B = await newPage({ state: { ...fixture, log: fixture.log.slice(20) } });
    await login(B);
    ok(await B.locator('.chooser [data-k="merge"]').count() === 1, 'ログイン時：端末にもアカウントにも記録があるときは、合わせるかを聞く');
    await B.click('.chooser [data-k="merge"]'); await B.waitForTimeout(1500);
    const nB = (await LS(B)).log.length;
    ok(nB === fixture.log.length, `端末B：両方の端末の記録が合わさる（${nB}件）`);
    ok(await B.evaluate(() => !!localStorage.getItem('celestia.auth') && !localStorage.getItem('lemma.auth')), 'ログイン情報は celestia.auth に保存（Lemma と別）');
    ok([...mock.rows.values()][0] && mock.rows.size === 1, 'サーバーの表 celestia_state に1行');
    const real = (es) => es.filter(e => !/status of 403/.test(e));
    ok(real(A.errors).length === 0 && real(B.errors).length === 0, `ページのエラーなし ${A.errors} ${B.errors}`);
  }
  // 9. v0.6（Lemma v0.28 から）：今日の学習・中断と再開・類題・自信度の一言・ブックマーク・検索の絞り込み・頻出度
  {
    const page = await newPage();
    await page.goto(U + '#today'); await page.waitForTimeout(1200);
    ok((await page.locator('h1').innerText()).includes('今日の学習') && await page.locator('#todaygo').count() === 1, '今日の画面に「今日の演習を始める」');
    ok(await page.locator('.first-run #firstgo').count() === 1 && await page.locator('#resume').count() === 0, '初めての利用では「まず8問解いてみる」、続きは無い');
    await page.click('#todaygo'); await page.waitForTimeout(1500);
    ok(await page.locator('#conf button').count() === 3 && (await page.locator('.why').innerText()).includes('新しい問題'), '今日の演習が始まり、出題の理由が出る');
    ok((await page.locator('.chip.val').allInnerTexts()).some(t => /頻出度 [高中低]/.test(t)), '頻出度のチップ');
    const first = await page.evaluate(() => JSON.parse(localStorage.getItem('celestia.session.v1')).steps[0].pid);
    await page.click('#quit'); await page.waitForTimeout(500);
    await page.goto(U + '#settings'); await page.waitForTimeout(800);
    await page.click('#amode [data-m="self"]'); await page.waitForTimeout(300);
    ok(await page.evaluate(() => JSON.parse(localStorage.getItem('celestia.v1')).settings.answer) === 'self', '設定：解答の形式を自己採点にできる');
    await page.goto(U + '#today'); await page.waitForTimeout(1200);
    ok(await page.locator('#resume').count() === 1, 'まだ答えていなくても、始めた演習は「前回の続き」に出る');
    await page.click('#resume'); await page.waitForTimeout(1500);
    const cur = await page.evaluate(() => JSON.parse(localStorage.getItem('celestia.session.v1')));
    ok(cur.i === 0 && cur.steps[0].pid === first, '再開すると1問目から');
    ok(await page.locator('#inputs').count() === 0 && await page.locator('.pen').count() === 1, '自己採点の形式では入力欄が出ない');
    await page.locator('#conf button').nth(0).click();
    await page.click('#reveal'); await page.waitForTimeout(500);
    ok((await page.locator('#self button').allInnerTexts()).join() === '不正解,一部正解,完答', '自己採点のボタンは 不正解→一部正解→完答 の順');
    await page.locator('#self button[data-g="full"]').click(); await page.waitForTimeout(500);
    const g = await page.evaluate(() => JSON.parse(localStorage.getItem('celestia.v1')).log.at(-1));
    ok(g.c === 'low' && g.ok && g.g === 2 && g.md === 'self', `自信なしで正解 → 早めに復習（grade ${g.g}）・形式が記録に残る`);
    ok((await page.locator('.conf-note').innerText()).includes('早めに復習'), '自信なしで正解したときの一言');
    ok(await page.locator('.ansrow.big').count() === 1 && await page.locator('.sol .solsec').count() >= 1, '答えの枠と、小問ごとの解説の枠');
    const nSteps = await page.evaluate(() => JSON.parse(localStorage.getItem('celestia.session.v1')).steps.length);
    if (await page.locator('#similar').count()) {
      await page.click('#similar'); await page.waitForTimeout(1200);
      const s2 = await page.evaluate(() => JSON.parse(localStorage.getItem('celestia.session.v1')));
      ok(s2.steps.length === nSteps + 1 && s2.i === 1 && /^類題/.test(s2.steps[1].why), '類題を1問追加すると次の問題に入る');
    } else ok(false, '類題のボタンが出る');
    await page.click('#quit'); await page.waitForTimeout(800);
    ok((await page.locator('#resume').innerText()).includes('2問目から'), '中断すると「前回の続き」が2問目から');
    await page.click('#resume'); await page.waitForTimeout(1500);
    ok((await page.locator('.why').innerText()).includes('類題'), '前回の続きから類題を再開できる');
    await page.locator('#conf button').nth(2).click();
    await page.click('#reveal'); await page.waitForTimeout(400);
    await page.locator('#self button[data-g="none"]').click(); await page.waitForTimeout(400);
    ok(await page.locator('.conf-note.warn').count() === 1, '確実なのに不正解のときは思い込みの確認を促す');
    // ブックマーク
    const pid = await page.evaluate(() => JSON.parse(localStorage.getItem('celestia.v1')).log.at(-1).p);
    await page.click(`[data-bm="${pid}"]`); await page.waitForTimeout(200);
    ok(await page.locator(`[data-bm="${pid}"].on`).count() === 1, 'ブックマーク：問題画面の ☆ で付けられる');
    await page.click('#quit'); await page.waitForTimeout(300);
    await page.goto(U + '#today'); await page.reload(); await page.waitForTimeout(1500);
    ok((await page.locator('#bmset').innerText()).includes('1問'), 'ブックマーク：ホームに「★ ブックマーク 1問」');
    await page.click('#bmset'); await page.waitForTimeout(1200);
    ok(await page.locator(`[data-bm="${pid}"].on`).count() === 1, 'ブックマーク：そこから出題できる');
    await page.click(`[data-bm="${pid}"]`); await page.waitForTimeout(200);
    const bm = await page.evaluate((p) => JSON.parse(localStorage.getItem('celestia.v1')).bookmarks[p], pid);
    ok(bm && bm.on === false && bm.t > 0, 'ブックマーク：外すと on:false の印が残る（同期でよみがえらない）');
    await page.click(`[data-bm="${pid}"]`); await page.waitForTimeout(200);
    await page.click('#quit'); await page.waitForTimeout(300);
    // 検索の絞り込み・並べ方
    await page.goto(U + '#search=' + encodeURIComponent('力')); await page.waitForTimeout(1500);
    const allN = await page.locator('.res-row').count();
    await page.selectOption('#sfilter', 'wrong'); await page.waitForTimeout(1000);
    const wrongIds = await page.$$eval('.res-row', els => els.map(x => x.dataset.pid));
    const wrongAll = await page.evaluate(() => { const m = new Map(); for (const l of JSON.parse(localStorage.getItem('celestia.v1')).log) m.set(l.p, l.g); return [...m].filter(([, g]) => g < 3).map(([p]) => p); });
    ok(wrongIds.every(id => wrongAll.includes(id)), `検索で「間違えた」に絞れる（${wrongIds.length}件）`);
    await page.selectOption('#sfilter', 'marked'); await page.waitForTimeout(1000);
    const mk = await page.$$eval('.res-row', els => els.map(x => x.dataset.pid));
    ok(mk.every(id => id === pid), `検索で「ブックマーク」に絞れる（${mk.length}件）`);
    await page.selectOption('#sfilter', 'fresh'); await page.waitForTimeout(1000);
    const fr = await page.$$eval('.res-row', els => els.map(x => x.dataset.pid));
    const solvedIds = await page.evaluate(() => Object.keys(JSON.parse(localStorage.getItem('celestia.v1')).problems));
    ok(fr.length > 0 && fr.every(id => !solvedIds.includes(id)), '検索で「未着手」に絞れる');
    await page.selectOption('#sfilter', 'due'); await page.waitForTimeout(800);
    ok(await page.locator('#sres').innerText().then(t => t.includes('復習待ち')), '検索で「復習待ち」に絞れる（今は無い）');
    await page.selectOption('#sfilter', 'all'); await page.waitForTimeout(800);
    await page.selectOption('#ssort', 'hard'); await page.waitForTimeout(1000);
    const lvs = await page.$$eval('.res-row .chip[class*="lv"]', els => els.map(x => +x.className.match(/lv(\d)/)[1]));
    ok(lvs.length > 3 && lvs.every((v, i) => !i || lvs[i - 1] >= v) && await page.locator('.res-row').count() === Math.min(allN, 30), '検索結果を「難しい順」に並べられる');
    await page.selectOption('#ssort', 'easy'); await page.waitForTimeout(800);
    const lvs2 = await page.$$eval('.res-row .chip[class*="lv"]', els => els.map(x => +x.className.match(/lv(\d)/)[1]));
    ok(lvs2.every((v, i) => !i || lvs2[i - 1] <= v), '「やさしい順」');
    await page.selectOption('#ssort', 'rel'); await page.waitForTimeout(500);
    await page.goto(U + '#map'); await page.waitForTimeout(800);
    ok(await page.locator('.next-acts').count() === 1 && (await page.locator('.trend').innerText()).includes('最近の伸び'), '定着度：次にやること・最近の伸び');
    ok(page.errors.length === 0, `ページのエラーなし ${page.errors}`);
  }
  // 10. 記録の保全（v0.6）：v0.5 の記録の移行・書き出しと読み込み・記録の消去・壊れた記録・新しい版・2つのタブ・ログイン時の確認
  {
    mock.users.clear(); mock.rows.clear();
    const v05 = JSON.parse(fs.readFileSync(path.join(root, 'tests/fixtures/state-v0.5.json'), 'utf8'));
    const P = await newPage({ state: v05 });
    await P.goto(U + '#today'); await P.waitForTimeout(1500);
    ok(await P.evaluate(() => JSON.parse(localStorage.getItem('celestia.v1')).schema) === 1, '起動しただけでは保存しない（v0.5 の記録はそのまま。画面の中では移行した記録を使う）');
    await P.locator('#lv button').first().click(); await P.locator('#lv button').first().click(); await P.waitForTimeout(300);   // 何か操作すると保存する
    const mig = await P.evaluate(() => ({ st: JSON.parse(localStorage.getItem('celestia.v1')), snaps: JSON.parse(localStorage.getItem('celestia.snap.v1') || '[]') }));
    ok(await P.locator('.banner').count() === 0, 'v0.5 の記録で起動しても警告は出ない');
    ok(mig.st.schema === 5 && mig.st.log.length === v05.log.length && JSON.stringify(mig.st.pstat) === JSON.stringify(v05.pstat) && mig.st.predictions.length === 2 && mig.st.misc['MC-PHOTO-INTENSITY'].n === 2 && mig.st.simPrefs['mech.pendulum'].l === 1.5, 'v0.5 の記録（schema 1）を移行し、解答の記録・復習状態・予測・誤概念・シミュレーションの条件を失わない');
    ok(mig.st.practice.fields.includes('MSY-01') && !mig.st.practice.fields.includes('SYN-01'), '移行で、選んだ「分野横断」を「力学総合」「電磁気学総合」に読み替える');
    ok(mig.snaps.length === 1 && /移行の前/.test(mig.snaps[0].reason) && JSON.parse(mig.snaps[0].raw).schema === 1, '移行の前の記録を控え（celestia.snap.v1）に残す');
    ok(!(await P.locator('.first-run').count()) && (await P.locator('.goal').count()) === 1, '記録のある端末では今日の目標が出る');
    await P.goto(U + '#sim=mech.pendulum'); await P.waitForTimeout(1200);
    ok((await P.locator('input[data-k="l"]').inputValue()) === '1.5', 'シミュレーションは前回の条件（simPrefs）で開く');
    // 書き出し（形式・版・確認用の値つき）
    await P.goto(U + '#settings'); await P.waitForTimeout(1000);
    ok((await P.locator('.integrity summary').innerText()).includes('記録の点検') && await P.locator('.snap-back').count() === 1, '設定：記録の点検と自動の控え');
    const [exp] = await Promise.all([P.waitForEvent('download'), P.click('#exp')]);
    const bk = JSON.parse(fs.readFileSync(await exp.path(), 'utf8'));
    ok(bk.format === 'celestia-backup' && bk.version === 1 && bk.schema === 5 && bk.app === '0.6.0' && /^fnv1a32:/.test(bk.hash) && bk.state.log.length === v05.log.length && bk.counts.predictions === 2, '書き出し：形式・版・作成時刻・件数・確認用の値つきのバックアップ');
    await P.setInputFiles('#imp', { name: 'x.json', mimeType: 'application/json', buffer: Buffer.from('{"hello":1}') });
    await P.waitForTimeout(400);
    ok((await P.locator('.toast').innerText()).includes('CELESTIA の学習記録のファイルではありません'), '読み込み：別のファイルは断る');
    const tampered = JSON.parse(JSON.stringify(bk)); tampered.state.log.pop();
    await P.setInputFiles('#imp', { name: 't.json', mimeType: 'application/json', buffer: Buffer.from(JSON.stringify(tampered)) });
    await P.waitForTimeout(400);
    ok((await P.locator('.toast').innerText()).includes('壊れているか'), '読み込み：中身が欠けたバックアップは断る');
    await P.setInputFiles('#imp', { name: 'b2.json', mimeType: 'application/json', buffer: Buffer.from(JSON.stringify(bk)) });
    await P.waitForTimeout(500);
    await P.click('.chooser [data-k="merge"]'); await P.waitForTimeout(600);
    ok((await P.locator('.chooser').innerText()).includes('読み込みました'), '読み込み：結果（件数・点検）を表示する');
    await P.click('.chooser [data-k="ok"]'); await P.waitForTimeout(400);
    ok((await P.evaluate(() => JSON.parse(localStorage.getItem('celestia.v1')).log.length)) === v05.log.length, '読み込み：同じバックアップを合わせても記録は重複しない');
    // 置き換え（v0.5 の書き出し＝記録そのもの、の旧形式）
    const small = { ...v05, log: v05.log.slice(0, 10) };
    await P.goto(U + '#settings'); await P.waitForTimeout(800);
    await P.setInputFiles('#imp', { name: 'old.json', mimeType: 'application/json', buffer: Buffer.from(JSON.stringify(small)) });
    await P.waitForTimeout(500);
    await P.click('.chooser [data-k="replace"]'); await P.waitForTimeout(600);
    await P.click('.chooser [data-k="ok"]'); await P.waitForTimeout(400);
    const rs = await P.evaluate(() => JSON.parse(localStorage.getItem('celestia.v1')));
    ok(rs.log.length === 10 && rs.resetAt > 0 && rs.schema === 5, '読み込み：v0.5 の書き出し（旧形式）で置き換えると、その記録になり消した印が付く');
    // 控えに戻す
    await P.goto(U + '#settings'); await P.waitForTimeout(800);
    await P.locator('.snap-back').first().click(); await P.waitForTimeout(400);
    await P.click('.chooser [data-k="replace"]'); await P.waitForTimeout(600);
    await P.click('.chooser [data-k="ok"]'); await P.waitForTimeout(400);
    ok((await P.evaluate(() => JSON.parse(localStorage.getItem('celestia.v1')).log.length)) === v05.log.length, '控えの時点（置き換える前）に戻せる');
    // 記録をすべて消す（ログインしていない）
    await P.goto(U + '#settings'); await P.waitForTimeout(800);
    await P.click('#reset'); await P.waitForTimeout(300);
    ok(await P.locator('.chooser [data-k="all"]').count() === 1 && await P.locator('.chooser [data-k="device"]').count() === 0, '消去：ログインしていなければ「消す」だけ');
    await P.click('.chooser [data-k="all"]'); await P.waitForTimeout(800);
    const wiped = await P.evaluate(() => JSON.parse(localStorage.getItem('celestia.v1')));
    ok(wiped.log.length === 0 && wiped.resetAt > 0 && wiped.settings.theme === 'dark' && wiped.simPrefs['mech.pendulum'], '消去：記録は消え、設定・シミュレーションの条件は残り、消した印が付く');
    // ログインして同期 → 別の端末でもログイン → 「記録をすべて消す」（アカウントも）でほかの端末の記録も消える
    const login = async (page, email) => {
      await page.goto(U + '#settings'); await page.waitForTimeout(1200);
      await page.fill('#acct-email', email); await page.click('#acct-send'); await page.waitForTimeout(400);
      await page.fill('#acct-code', '123456'); await page.click('#acct-verify'); await page.waitForTimeout(1500);
    };
    const A = await newPage({ state: { ...v05, settings: { ...v05.settings, answer: 'self' } } });
    await login(A, 'reset@example.com');
    const srvN = () => { const r = [...mock.rows.values()][0]; return r ? r.state.log.length : -1; };
    ok(srvN() === v05.log.length && [...mock.rows.values()][0].state.schema === 5, `同期：アカウントに移行後の記録が上がる（${srvN()}件）`);
    const Q = await newPage({ state: { ...v05, log: v05.log.slice(0, 10) } });
    await login(Q, 'reset@example.com');
    if (await Q.locator('.chooser').count()) { await Q.click('.chooser [data-k="merge"]'); await Q.waitForTimeout(1500); }
    await A.goto(U + '#settings'); await A.waitForTimeout(1000);
    await A.click('#reset'); await A.waitForTimeout(300);
    ok(await A.locator('.chooser [data-k="all"]').count() === 1 && await A.locator('.chooser [data-k="device"]').count() === 1, '消去：ログイン中は「アカウントも」「この端末だけ」を選べる');
    await A.click('.chooser [data-k="all"]'); await A.waitForTimeout(1500);
    ok(srvN() === 0, `消去：アカウントの記録も消える（${srvN()}件）`);
    await Q.goto(U + '#settings'); await Q.waitForTimeout(800);
    await Q.click('#acct-now'); await Q.waitForTimeout(1500);
    const nQ = await Q.evaluate(() => JSON.parse(localStorage.getItem('celestia.v1')).log.length);
    ok(nQ === 0 && srvN() === 0, `消去：ほかの端末の古い記録も、同期で消える（端末 ${nQ}件・サーバー ${srvN()}件）`);
    // ログイン時に「アカウントの記録だけを使う」：端末の記録は書き出してから消える
    await A.evaluate("__celestia.open(['MEC-CON-001'])"); await A.waitForTimeout(1000);
    await A.locator('#conf button').nth(1).click(); await A.click('#reveal'); await A.waitForTimeout(400);
    await A.locator('#self button[data-g="full"]').click(); await A.waitForTimeout(400);
    await A.goto(U + '#settings'); await A.waitForTimeout(800);
    await A.click('#acct-now'); await A.waitForTimeout(1500);
    const R = await newPage({ state: { ...v05, log: v05.log.slice(0, 5) } });
    await R.goto(U + '#settings'); await R.waitForTimeout(1200);
    await R.fill('#acct-email', 'reset@example.com'); await R.click('#acct-send'); await R.waitForTimeout(400);
    await R.fill('#acct-code', '123456'); await R.click('#acct-verify'); await R.waitForTimeout(1200);
    const [dl] = await Promise.all([R.waitForEvent('download', { timeout: 5000 }).catch(() => null), R.click('.chooser [data-k="remote"]')]);
    await R.waitForTimeout(1500);
    const nR = await R.evaluate(() => JSON.parse(localStorage.getItem('celestia.v1')).log.map(l => l.p));
    ok(dl && JSON.stringify(nR) === '["MEC-CON-001"]', `ログイン時「アカウントの記録だけ」：端末の記録を書き出してから、アカウントの記録になる（${nR}）`);
    await R.goto(U + '#settings'); await R.waitForTimeout(800);
    await R.click('#acct-out'); await R.waitForTimeout(300);
    await R.click('.chooser [data-k="wipe"]'); await R.waitForTimeout(1500);
    const nR2 = await R.evaluate(() => JSON.parse(localStorage.getItem('celestia.v1')).log.length);
    ok(nR2 === 0 && await R.locator('#acct-email').count() === 1 && srvN() === 1, 'ログアウト時に端末の記録も消せる（サーバーの記録は残る）');
    // 壊れた記録：上書きせずに退避する
    const BR = await newPage({ raw: '{"log":[{"p":"MEC-CON-001","t":1790000000000,"g":3}], こわれた' });
    await BR.goto(U + '#today'); await BR.waitForTimeout(1200);
    ok((await BR.locator('.banner').innerText()).includes('読めませんでした'), '壊れた記録：起動時に知らせる');
    const brk = await BR.evaluate(() => Object.keys(localStorage).filter(k => k.startsWith('celestia.v1.broken.')).map(k => JSON.parse(localStorage.getItem(k)).raw));
    ok(brk.length === 1 && brk[0].includes('こわれた'), '壊れた記録：元の文字列を消さずに退避する（celestia.v1.broken.）');
    // 2つのタブで解いても、どちらの記録も残る
    const T1 = await newPage({ state: { settings: { size: 6, answer: 'self' } } });
    await T1.goto(U + '#today'); await T1.waitForTimeout(1000);
    const T2 = await T1.context().newPage(); T2.errors = [];
    await T2.goto(U + '#today'); await T2.waitForTimeout(1000);
    const solve = async (pg, pid) => {
      await pg.evaluate(`__celestia.open(['${pid}'])`); await pg.waitForTimeout(1000);
      await pg.locator('#conf button').nth(1).click(); await pg.click('#reveal'); await pg.waitForTimeout(400);
      await pg.locator('#self button[data-g="full"]').click(); await pg.waitForTimeout(300);
    };
    await solve(T1, 'MEC-INC-001');
    await solve(T2, 'MEC-CON-001');
    const both = await T2.evaluate(() => JSON.parse(localStorage.getItem('celestia.v1')).log.map(l => l.p));
    ok(both.includes('MEC-INC-001') && both.includes('MEC-CON-001'), `2つのタブ：あとから保存したタブが、先のタブの記録を消さない（${both}）`);
    // 新しい版（schema 99）の記録は保存しない
    const newerRaw = JSON.stringify({ ...v05, schema: 99, settings: { ...v05.settings, answer: 'self' }, future: { x: 1 } });
    const NW = await newPage({ raw: newerRaw });
    await NW.goto(U + '#today'); await NW.waitForTimeout(1200);
    ok((await NW.locator('.banner').innerText()).includes('保存も同期もしません'), '新しい schema：起動時に保存しないことを知らせる');
    await solve(NW, 'MEC-CON-001');
    ok(await NW.evaluate(() => localStorage.getItem('celestia.v1')) === newerRaw, '新しい schema：解いても保存された記録は書き換わらない');
    // 同じドメインの Lemma の記録・キャッシュに触れない
    const LM = await newPage();
    await LM.addInitScript(() => { if (!sessionStorage.getItem('lm')) { sessionStorage.setItem('lm', '1'); localStorage.setItem('lemma.v2', '{"lemma":1}'); localStorage.setItem('lemma.auth', 'x'); } });
    await LM.goto(U + '#today'); await LM.waitForTimeout(1000);
    await solve(LM, 'MEC-CON-001');
    ok(await LM.evaluate(() => localStorage.getItem('lemma.v2') === '{"lemma":1}' && localStorage.getItem('lemma.auth') === 'x'), '同じドメインの Lemma の記録（lemma.v2・lemma.auth）に触れない');
    const errs = [...P.errors, ...A.errors, ...Q.errors, ...R.errors, ...BR.errors, ...T1.errors, ...T2.errors, ...NW.errors, ...LM.errors].filter(e => !/status of 40[13]/.test(e));
    ok(errs.length === 0, `ページのエラーなし ${errs}`);
  }
} finally {
  await browser.close();
  server.close();
}
console.log(`e2e: ${pass} 件 OK${fail ? `、${fail} 件失敗` : ''}`);
process.exitCode = fail ? 1 : 0;
