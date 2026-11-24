// 解答の形式（設定）を切り替えて、記述式・選択式・自己採点の画面と記録を確かめる（Playwright・Chromium。Lemma v0.28 の tests/e2e_modes.mjs を CELESTIA 向けにしたもの）
//   node tests/e2e_modes.mjs [問題ID]      SHOTS=フォルダ を付けるとスクリーンショットを保存
// data に選択肢（ch）がある問題があればそれを使う。まだ無ければ、テストの中だけで MEC-INC-001 に選択肢を付けて配信する（data は変えない）。
import http from 'http'; import fs from 'fs'; import path from 'path';
import { fileURLToPath } from 'url';
import { launchBrowser } from './browser.mjs';
const root = path.resolve(path.dirname(fileURLToPath(import.meta.url)), '..');
const out = process.env.SHOTS || null;
const SHOT = !!out;
let failed = 0, passed = 0;
const idx = JSON.parse(fs.readFileSync(path.join(root, 'data/index.json'), 'utf8'));
const real = idx.problems.find(e => e.ch);
const ID = process.argv[2] || (real ? real.id : 'MEC-INC-001');
// テスト用の選択肢（入力欄と同じ並び。null の欄は入力のまま）：data にまだ選択肢が無いときだけ使う
const FAKE = !real && !process.argv[2] ? { 'MEC-INC-001': ['mg\\cos\\theta', 'mg\\sin\\theta', 'mg', 'mg\\tan\\theta'] } : null;
const fakeCh = FAKE && [FAKE['MEC-INC-001'], ['mg\\sin\\theta', 'mg\\cos\\theta', 'mg', '\\mu mg'], null];
const T = { '.html': 'text/html', '.js': 'text/javascript', '.css': 'text/css', '.json': 'application/json', '.woff2': 'font/woff2', '.svg': 'image/svg+xml', '.webmanifest': 'application/manifest+json', '.png': 'image/png' };
const srv = http.createServer((q, r) => {
  let p = decodeURIComponent(new URL(q.url, 'http://x').pathname); if (p.endsWith('/')) p += 'index.html';
  const f = path.join(root, p); if (!f.startsWith(root) || !fs.existsSync(f)) { r.writeHead(404); return r.end(); }
  if (fakeCh && /\/data\/(index|f-[A-Z]{3}-\d\d)\.json$/.test(p)) {
    const d = JSON.parse(fs.readFileSync(f, 'utf8'));
    if (p.endsWith('index.json')) { const e = d.problems.find(x => x.id === ID); if (e) e.ch = 1; }
    else if (d[ID]) d[ID].ch = fakeCh;
    r.writeHead(200, { 'Content-Type': 'application/json' }); return r.end(JSON.stringify(d));
  }
  r.writeHead(200, { 'Content-Type': T[path.extname(f)] || 'application/octet-stream' }); fs.createReadStream(f).pipe(r);
});
await new Promise(r => srv.listen(0, r));
const U = `http://localhost:${srv.address().port}/`;
const b = await launchBrowser();
const page = await b.newPage({ viewport: { width: 390, height: 844 }, deviceScaleFactor: 2 });
const errs = []; page.on('pageerror', e => errs.push(e.message)); page.on('console', m => m.type() === 'error' && errs.push(m.text()));
page.on('dialog', d => d.accept());
const ok = (c, m) => { if (!c) { failed++; console.log('✗', m); } else passed++; };
const open = async () => { await page.goto(U + '#ids=' + ID); await page.waitForTimeout(1800); };
await page.goto(U + '#settings'); await page.waitForTimeout(1200);
ok(await page.locator('#amode button').count() === 3, '設定に3つの形式');
ok((await page.locator('#amode button.on').innerText()) === '記述式', '既定は記述式');
await page.locator('#amode button[data-m="choice"]').click(); await page.waitForTimeout(300);
SHOT && await page.screenshot({ path: out + '/1-settings.png' });
// 選択式
await open();
const nOpts = await page.locator('#inputs .choice').count();
ok(nOpts >= 4, `選択式：選択肢ボタン ${nOpts}個`);
ok((await page.locator('.chip').allInnerTexts()).includes('選択式'), 'チップが選択式');
ok(await page.locator('#inputs .choice .katex').count() === nOpts, '選択肢は数式として描画される');
if (fakeCh) ok(await page.locator('#inputs input').count() === 1, '選択肢の無い欄（null）は入力のまま');
SHOT && await page.screenshot({ path: out + '/2-choice.png', fullPage: true });
const keys = [...new Set(await page.locator('#inputs .choice').evaluateAll(bs => bs.map(x => x.dataset.k)))];
// 最初の欄は正解（data-i=0）、ほかは誤答（data-i=1）
for (const [j, k] of keys.entries()) await page.click(`#inputs .choice[data-k="${k}"][data-i="${j === 0 ? 0 : 1}"]`);
for (const inp of await page.locator('#inputs input').all()) await inp.fill(fakeCh ? 'μ>=tanθ' : '1');
await page.locator('#conf button').nth(2).click();
await page.click('#judge'); await page.waitForTimeout(800);
ok(await page.locator('#inputs .choice.right').count() === keys.length, '判定後：正解の選択肢に印');
ok(await page.locator('#inputs .choice.wrong').count() === keys.length - 1, '判定後：選んだ誤答に印');
ok(await page.locator('#judge').count() === 0, '選択式はやり直しなし');
ok(await page.locator('.field.pick .dispute').count() === 0, '選択欄に「判定がおかしい？」は出ない');
SHOT && await page.screenshot({ path: out + '/3-choice-judged.png', fullPage: true });
if (await page.locator('#self button').count()) await page.locator('#self button[data-g="full"]').click();
await page.waitForTimeout(300);
let last = await page.evaluate(() => JSON.parse(localStorage.getItem('celestia.v1')).log.at(-1));
ok(last.md === 'choice' && (keys.length > 1 ? last.sg !== 'full' : last.sg === 'full'), `選択式の記録（${last.md}・${last.sg}）`);
// 選択式で全部正解・「確実」でも、復習の間隔は 3 まで
await open();
for (const k of keys) await page.click(`#inputs .choice[data-k="${k}"][data-i="0"]`);
for (const inp of await page.locator('#inputs input').all()) await inp.fill(fakeCh ? 'μ>=tanθ' : '1');
await page.locator('#conf button').nth(2).click();
await page.click('#judge'); await page.waitForTimeout(800);
if (await page.locator('#self button').count()) await page.locator('#self button[data-g="full"]').click();
await page.waitForTimeout(300);
last = await page.evaluate(() => JSON.parse(localStorage.getItem('celestia.v1')).log.at(-1));
ok(fakeCh ? last.g === 3 && last.sg === 'full' : last.g <= 3, `選択式の正解は「確実」でも grade 3 まで（${last.g}）`);
// 選択肢の無い問題は入力に戻る
await page.goto(U + '#ids=MEC-PRJ-003'); await page.waitForTimeout(1500);
const hasCh = idx.problems.find(e => e.id === 'MEC-PRJ-003').ch;
if (!hasCh) ok((await page.locator('.pen').innerText()).includes('選択肢がない') && await page.locator('#inputs input').count() > 0, '選択式でも選択肢の無い問題は入力で答える');
// 自己採点
await page.goto(U + '#settings'); await page.waitForTimeout(600);
await page.locator('#amode button[data-m="self"]').click(); await page.waitForTimeout(200);
await open();
ok(await page.locator('#inputs').count() === 0 && await page.locator('.pen').count() === 1, '自己採点：入力欄なし');
ok((await page.locator('.chip').allInnerTexts()).includes('自己採点'), 'チップが自己採点');
await page.locator('#conf button').nth(1).click(); await page.click('#reveal'); await page.waitForTimeout(500);
ok(await page.locator('#self button').count() === 3, '自己採点のボタン');
ok(await page.locator('#solbox').count() === 1, '自己採点：先に解答・解説（解答図を含む）を見せる');
SHOT && await page.screenshot({ path: out + '/4-self.png', fullPage: true });
await page.locator('#self button[data-g="part"]').click(); await page.waitForTimeout(300);
// 記述式
await page.goto(U + '#settings'); await page.waitForTimeout(600);
await page.locator('#amode button[data-m="input"]').click(); await page.waitForTimeout(200);
await open();
ok(await page.locator('#inputs input').count() > 0 && await page.locator('#inputs .choice').count() === 0, '記述式：入力欄');
ok((await page.locator('.chip').allInnerTexts()).includes('入力して採点'), 'チップが入力して採点');
SHOT && await page.screenshot({ path: out + '/5-input.png' });
const log = await page.evaluate(() => JSON.parse(localStorage.getItem('celestia.v1') || '{}').log?.map(l => l.md));
ok(JSON.stringify(log) === '["choice","choice","self"]', '記録に形式が残る ' + JSON.stringify(log));
ok(!errs.length, 'エラーなし ' + errs.join(' / '));
await b.close(); srv.close();
console.log(`解答の形式：${passed} 件 OK${failed ? `、${failed} 件失敗` : ''}${fakeCh ? '（data に選択肢がまだ無いので、テストの中だけで MEC-INC-001 に選択肢を付けた）' : ''}`);
process.exit(failed ? 1 : 0);
