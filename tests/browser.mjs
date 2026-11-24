// Playwright の Chromium を起動する（e2e のテストで共通）。
//   Playwright の同梱ブラウザが無い環境（npx playwright install をしていない）では、
//   CHROME_PATH か、よくある場所の Chrome／Chromium を使う。
import fs from 'fs';
import { createRequire } from 'module';
const require = createRequire(import.meta.url);
export const { chromium } = require('playwright');
export async function launchBrowser(opts = {}) {
  try { return await chromium.launch(opts); }
  catch (e) {
    const p = [process.env.CHROME_PATH, '/usr/bin/google-chrome', '/usr/bin/google-chrome-stable', '/usr/bin/chromium', '/usr/bin/chromium-browser'].find(x => x && fs.existsSync(x));
    if (!p) throw e;
    return chromium.launch({ ...opts, executablePath: p });
  }
}
