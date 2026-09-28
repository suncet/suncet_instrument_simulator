import assert from 'node:assert/strict';
import {readFile, writeFile, mkdir, mkdtemp, cp, rm} from 'node:fs/promises';
import {createServer} from 'node:http';
import {fileURLToPath, pathToFileURL} from 'node:url';
import {tmpdir} from 'node:os';
import {test} from 'node:test';
import path from 'node:path';
const {chromium} = await import(process.env.PLAYWRIGHT_MODULE || 'playwright');
const site = path.resolve(fileURLToPath(new URL('../../../docs/colortables/', import.meta.url)));
const manifest = JSON.parse(await readFile(path.join(site, 'manifest.json')));
const screenshotDir = process.env.GALLERY_SCREENSHOTS;

test('gallery works locally and public votes persist, aggregate, undo, export, and recover', async () => {
  assert.equal(await readFile(path.join(site,'gallery.js'),'utf8'), await readFile(new URL('../gallery.js', import.meta.url),'utf8'));
  const localRoot = await mkdtemp(path.join(tmpdir(),'suncet-gallery-test-'));
  await cp(site,localRoot,{recursive:true});
  await writeFile(path.join(localRoot,'config.js'),'window.SUNCET_VOTING = {};');
  const ledger = new Map();
  let origin, signups = 0, failNext = false, closed = false;
  const server = createServer(async (req, res) => {
    const url = new URL(req.url, origin);
    const json = (body, status = 200) => { res.writeHead(status, {'Content-Type':'application/json'}); res.end(JSON.stringify(body)); };
    try {
      if (url.pathname === '/config.js') {
        res.setHeader('Content-Type','text/javascript');
        return res.end(`window.SUNCET_VOTING=${JSON.stringify({supabaseUrl:origin,publishableKey:'sb_publishable_browser_test',turnstileSiteKey:''})};`);
      }
      if (url.pathname === '/auth/v1/signup') {
        signups++;
        const id = `00000000-0000-0000-0000-${String(signups).padStart(12,'0')}`;
        const claims = {sub:id,exp:Math.floor(Date.now()/1000)+3600,role:'authenticated',aud:'authenticated'};
        const token = `${Buffer.from(JSON.stringify({alg:'HS256',typ:'JWT'})).toString('base64url')}.${Buffer.from(JSON.stringify(claims)).toString('base64url')}.test`;
        return json({access_token:token,token_type:'bearer',expires_in:3600,refresh_token:'refresh-'+id,user:{id,aud:'authenticated',role:'authenticated',is_anonymous:true,app_metadata:{},user_metadata:{},created_at:new Date().toISOString()}});
      }
      const token = req.headers.authorization?.split(' ')[1];
      const user = token?.includes('.') ? JSON.parse(Buffer.from(token.split('.')[1],'base64url')).sub : null;
      if (url.pathname === '/rest/v1/colortable_studies') return json({is_open:!closed});
      if (url.pathname === '/rest/v1/colortable_favorites') return json([...(ledger.get(user)||[])].map(palette_slug=>({palette_slug})));
      if (url.pathname === '/rest/v1/rpc/colortable_counts') return json(manifest.options.map(item=>({palette_slug:item.slug,favorites:[...ledger.values()].filter(set=>set.has(item.slug)).length})));
      if (url.pathname === '/rest/v1/rpc/set_colortable_favorite') {
        let body = ''; for await (const chunk of req) body += chunk;
        if (failNext) { failNext=false; return json({message:'Temporary failure'},500); }
        if (closed || !user) return json({message:'Voting unavailable'},403);
        const input=JSON.parse(body),votes=ledger.get(user)||new Set();
        input.p_favorite ? votes.add(input.p_palette) : votes.delete(input.p_palette);
        ledger.set(user,votes);
        return json(input.p_favorite);
      }
      const file = path.resolve(site, '.' + decodeURIComponent(url.pathname === '/' ? '/index.html' : url.pathname));
      if (!file.startsWith(site + path.sep)) return json({},403);
      const mime={'.js':'text/javascript','.html':'text/html','.png':'image/png','.webp':'image/webp','.json':'application/json','.csv':'text/csv'};
      res.setHeader('Content-Type',mime[path.extname(file)]||'text/plain');
      res.end(await readFile(file));
    } catch { res.writeHead(404); res.end(); }
  });
  await new Promise(resolve=>server.listen(0,'127.0.0.1',resolve));
  origin=`http://127.0.0.1:${server.address().port}`;
  const browser = await chromium.launch({headless:true,executablePath:process.env.CHROME_EXECUTABLE});
  try {
    const context = await browser.newContext({viewport:{width:1440,height:1000}});
    const page = await context.newPage(), errors=[];
    page.on('pageerror',error=>errors.push(error.message));
    await page.goto(pathToFileURL(path.join(localRoot,'index.html')).href);
    await page.waitForFunction(()=>document.querySelector('#voteStatus').textContent.includes('Local preview'));
    assert.equal(await page.locator('.option').count(),40);
    assert.equal(await page.title(),'SunCET | Color Voting');
    const socialImage='https://suncet.github.io/suncet_instrument_simulator/social-preview-v1.png';
    assert.equal(await page.locator('meta[property="og:image"]').getAttribute('content'),socialImage);
    assert.equal(await page.locator('meta[name="twitter:image"]').getAttribute('content'),socialImage);
    assert.equal(await page.locator('meta[name="twitter:card"]').getAttribute('content'),'summary_large_image');
    const socialBytes=await readFile(path.join(site,'social-preview-v1.png'));
    assert.equal(socialBytes.readUInt32BE(16),1200);
    assert.equal(socialBytes.readUInt32BE(20),630);
    assert.equal(await page.locator('h1').textContent(),'SunCET | Color table options');
    assert.equal(await page.locator('#subtitle').count(),0);
    assert.equal(await page.locator('#order option[value="newest"]').count(),0);
    assert.equal(manifest.options.some(item=>[24,29].includes(item.id)),false);
    assert.equal(manifest.options.find(item=>item.slug==='tequila-sunrise').id,38);
    assert.deepEqual(manifest.options.filter(item=>item.slug.startsWith('euvi')).map(item=>item.source),['euvi171','euvi195','euvi284','euvi304']);
    assert.equal(manifest.options.filter(item=>item.family==='SunCET branding').every(item=>item.title.startsWith('SunCET branding /')),true);
    await page.getByRole('button',{name:'Favorite Poster / Blue to rose',exact:true}).click();
    await page.reload();
    assert.equal(await page.locator('#favoritesCount').textContent(),'1');
    assert.equal(await page.getByRole('tab',{name:'Results',exact:true}).isDisabled(),true);
    await page.goto(origin);
    await page.waitForFunction(()=>document.querySelector('#voteStatus').textContent==='Vote (with the star) for as many as you like');
    assert.equal(signups,0);
    assert.equal(await page.locator('#favoritesCount').textContent(),'0');
    await page.selectOption('#order','number');
    assert.match(await page.locator('.caption h2').first().textContent(),/Inferno/);
    await page.selectOption('#family','SunCET NASA Poster');
    assert.equal(await page.locator('.option').count(),5);
    const vote = page.getByRole('button',{name:'Favorite Poster / Blue to rose',exact:true});
    await vote.click();
    await page.waitForFunction(()=>document.querySelector('#voteStatus').textContent==='Favorite counted.');
    assert.equal(signups,1);
    await page.reload();
    await page.waitForFunction(()=>document.querySelector('#favoritesCount').textContent==='1');
    await page.getByRole('tab',{name:'Results',exact:true}).click();
    await page.waitForFunction(()=>document.querySelector('#resultsStatus').textContent.startsWith('Updated'));
    assert.match(await page.locator('#rankings tr').first().textContent(),/Blue to rose.*1/);
    const secondContext = await browser.newContext(), second = await secondContext.newPage();
    await second.goto(origin);
    await second.waitForFunction(()=>document.querySelector('#voteStatus').textContent==='Vote (with the star) for as many as you like');
    await second.getByRole('button',{name:'Favorite Poster / Blue to rose',exact:true}).click();
    await second.waitForFunction(()=>document.querySelector('#voteStatus').textContent==='Favorite counted.');
    assert.equal(signups,2);
    await page.getByRole('button',{name:'Refresh',exact:true}).click();
    await page.waitForFunction(()=>document.querySelector('#rankings .votes').textContent==='2');
    const downloadPromise=page.waitForEvent('download');
    await page.getByRole('button',{name:'Export CSV',exact:true}).click();
    const download=await downloadPromise;
    const csv=await readFile(await download.path(),'utf8');
    assert.match(csv,/"33","Poster \/ Blue to rose","SunCET NASA Poster","2"/);
    assert.equal(csv.trim().split('\r\n').length,41);
    if(screenshotDir){await mkdir(screenshotDir,{recursive:true});await page.screenshot({path:path.join(screenshotDir,'results-desktop.png')});}
    await page.getByRole('tab',{name:'Gallery',exact:true}).click();
    await page.getByRole('button',{name:'Favorite Poster / Blue to rose',exact:true}).click();
    await page.waitForFunction(()=>document.querySelector('#voteStatus').textContent==='Favorite removed.');
    assert.equal([...ledger.values()].filter(set=>set.has('poster-blue-rose')).length,1);
    failNext=true;
    await page.getByRole('button',{name:'Favorite Poster / Blue + lilac',exact:true}).click();
    await page.waitForFunction(()=>document.querySelector('#voteStatus').textContent.includes('not confirmed'));
    assert.equal(await page.locator('#favoritesCount').textContent(),'0');
    await page.getByRole('button',{name:'Retry',exact:true}).click();
    await page.waitForFunction(()=>document.querySelector('#voteStatus').textContent==='Vote (with the star) for as many as you like');
    await page.selectOption('#family','SunCET NASA Poster');
    await page.selectOption('#stretch','asinh');
    await page.getByRole('button',{name:'Compare Poster / Blue to rose',exact:true}).click();
    await page.locator('#selectedImage').evaluate(image=>image.decode());
    assert.match(await page.locator('#selectedImage').getAttribute('src'),/^asinh\//);
    assert.equal(await page.locator('#previewStretch').inputValue(),'asinh');
    await page.selectOption('#previewStretch','current');
    assert.equal(await page.locator('#stretch').inputValue(),'current');
    for (const id of ['selectedImage','referenceImage']) assert.match(await page.locator('#'+id).getAttribute('src'),/^current\//);
    assert.match(await page.locator('#fullImage').getAttribute('href'),/^current\//);
    assert.match(await page.locator('#referenceCaption').textContent(),/fourth root/);
    await page.getByRole('button',{name:'Next palette',exact:true}).click();
    assert.equal(await page.locator('#previewStretch').inputValue(),'current');
    await page.selectOption('#previewStretch','asinh');
    await page.getByRole('button',{name:'Single',exact:true}).click();
    assert.equal(await page.locator('#compareImages').evaluate(node=>node.classList.contains('single')),true);
    await page.selectOption('#previewStretch','current');
    assert.match(await page.locator('#selectedImage').getAttribute('src'),/^current\//);
    await page.getByRole('button',{name:'Compare',exact:true}).click();
    for (const width of [1440,390,320]) {
      await page.setViewportSize({width,height:900});
      await page.selectOption('#previewStretch','asinh');
      await page.locator('#selectedImage').evaluate(image=>image.decode());
      await page.locator('#referenceImage').evaluate(image=>image.decode());
      assert.equal(await page.locator('#preview').evaluate(node=>node.scrollWidth<=node.clientWidth),true);
      assert.equal(await page.locator('#previewStretch').evaluate(node=>{const box=node.getBoundingClientRect();return box.left>=0&&box.right<=innerWidth&&box.top>=0&&box.bottom<=innerHeight;}),true);
      if(screenshotDir)await page.screenshot({path:path.join(screenshotDir,`preview-${width}.png`)});
    }
    await page.keyboard.press('Escape');
    await page.locator('#poster-reference img').evaluate(image=>image.decode());
    assert.equal(await page.locator('#poster-reference img').evaluate(image=>image.naturalWidth>0),true);
    for(const width of [1920,1440,390,320]) {
      await page.setViewportSize({width,height:900});
      await page.evaluate(()=>scrollTo(0,0));
      await page.evaluate(()=>Promise.all([...document.querySelectorAll('.grid img')].map(image=>image.decode())));
      assert.equal(await page.evaluate(()=>document.documentElement.scrollWidth>innerWidth),false);
      if(screenshotDir)await page.screenshot({path:path.join(screenshotDir,`gallery-${width}.png`)});
    }
    await page.selectOption('#family','All');
    for(const width of [1440,390,320]) {
      await page.setViewportSize({width,height:900});
      await page.evaluate(()=>scrollTo(0,1200));
      const scrollBefore=await page.evaluate(()=>scrollY);
      assert.equal(await page.locator('.toolbar').evaluate(node=>Math.round(node.getBoundingClientRect().top)),0);
      await page.selectOption('#stretch','current');
      await page.selectOption('#stretch','asinh');
      assert.ok(Math.abs(await page.evaluate(()=>scrollY)-scrollBefore)<2);
      assert.equal(await page.locator('#stretch').evaluate(node=>{const box=node.getBoundingClientRect();return box.top>=0&&box.bottom<innerHeight;}),true);
      if(screenshotDir)await page.screenshot({path:path.join(screenshotDir,`sticky-${width}.png`)});
    }
    closed=true;
    await page.reload();
    await page.waitForFunction(()=>document.querySelector('#voteStatus').textContent.includes('Voting has closed'));
    assert.equal(await page.getByRole('button',{name:'Favorite Poster / Blue to rose',exact:true}).isDisabled(),true);
    assert.equal(await page.getByRole('tab',{name:'Results',exact:true}).isEnabled(),true);
    assert.deepEqual(errors,[]);
    await secondContext.close();
  } finally {
    await browser.close();
    await new Promise(resolve=>server.close(resolve));
    await rm(localRoot,{recursive:true,force:true});
  }
});
