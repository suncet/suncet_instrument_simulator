(() => {
  'use strict';
  const study = window.SUNCET_STUDY;
  const config = window.SUNCET_VOTING || {};
  const items = study.options;
  const bySlug = new Map(items.map(item => [item.slug, item]));
  const $ = id => document.getElementById(id);
  const configured = Boolean(config.supabaseUrl && config.publishableKey);
  const local = location.protocol === 'file:' && !configured;
  const localKey = 'suncet-colortable-favorites';
  const pending = new Set();
  let favorites = new Set(), selected = items[0], visible = items;
  let client, ready = local, sessionPromise, countRequest = 0, counts = null;
  let currentView = 'gallery', initialising = false, favoriteRevision = 0;

  function readStored(key, fallback) {
    try { return JSON.parse(localStorage.getItem(key)) ?? fallback; } catch { return fallback; }
  }
  function store(key, value) {
    try { localStorage.setItem(key, JSON.stringify(value)); } catch { /* Storage can be disabled. */ }
  }
  const saved = readStored(localKey, []);
  if (local && Array.isArray(saved)) favorites = new Set(saved.filter(slug => bySlug.has(slug)));
  const orderKey = `${study.study_id}-order`;
  let order = readStored(orderKey, []);
  order = Array.isArray(order) ? [...new Set(order)].filter(slug => bySlug.has(slug)) : [];
  const newSlugs = items.map(item => item.slug).filter(slug => !order.includes(slug));
  for (let i = newSlugs.length - 1; i > 0; i--) {
    const j = Math.floor(Math.random() * (i + 1));
    [newSlugs[i], newSlugs[j]] = [newSlugs[j], newSlugs[i]];
  }
  order.push(...newSlugs);
  store(orderKey, order);

  function status(message, error = false) {
    $('voteStatus').textContent = message;
    $('statusLine').dataset.error = String(error);
  }
  function element(tag, className, text) {
    const node = document.createElement(tag);
    if (className) node.className = className;
    if (text !== undefined) node.textContent = text;
    return node;
  }
  function matches(item) {
    return ($('family').value === 'All' || item.family === $('family').value || item.id === 1) &&
      (!$('favoritesOnly').checked || favorites.has(item.slug));
  }
  function favoriteButton(item) {
    const button = element('button', 'favorite' + (favorites.has(item.slug) ? ' active' : ''),
      favorites.has(item.slug) ? '\u2605' : '\u2606');
    button.title = `Favorite ${item.title}`;
    button.setAttribute('aria-label', button.title);
    button.setAttribute('aria-pressed', favorites.has(item.slug));
    button.disabled = !ready || pending.has(item.slug);
    button.onclick = () => toggleFavorite(item.slug);
    return button;
  }
  function render() {
    const stretch = $('stretch').value;
    visible = items.filter(matches).sort((a, b) => {
      if ($('order').value === 'number') return a.id - b.id;
      return order.indexOf(a.slug) - order.indexOf(b.slug);
    });
    $('grid').replaceChildren();
    $('count').textContent = `${visible.length} / ${items.length} color tables`;
    $('favoritesCount').textContent = favorites.size;
    $('empty').hidden = visible.length > 0;
    $('overview').href = study.overviews?.[stretch] || `overview-${stretch}.png`;
    for (const item of visible) {
      const figure = element('figure', 'option');
      const button = element('button', 'image-button');
      button.title = `Compare ${item.title}`;
      button.setAttribute('aria-label', button.title);
      button.onclick = () => openPreview(item);
      const img = element('img');
      img.src = item.thumbnails?.[stretch] || item.images[stretch];
      img.alt = item.title;
      img.width = study.shape[1];
      img.height = study.shape[0];
      img.loading = 'lazy';
      button.append(img);
      const caption = element('figcaption', 'caption');
      caption.append(element('span', 'number', String(item.id).padStart(2, '0')),
        element('h2', '', item.title), favoriteButton(item));
      const ramp = element('img', 'ramp');
      ramp.src = item.ramp;
      ramp.alt = 'Dark-to-bright color ramp';
      figure.append(button, caption, element('p', 'family', item.family), ramp);
      $('grid').append(figure);
    }
    renderResults();
    if ($('preview').open) updatePreview();
  }

  const scripts = new Map();
  function loadScript(src) {
    if (!scripts.has(src)) {
      scripts.set(src, new Promise((resolve, reject) => {
        const script = document.createElement('script');
        script.src = src;
        script.onload = resolve;
        script.onerror = () => { scripts.delete(src); script.remove(); reject(new Error('A required service could not be reached.')); };
        document.head.append(script);
      }));
    }
    return scripts.get(src);
  }
  async function verifyVisitor() {
    if (!config.turnstileSiteKey) return undefined;
    await loadScript('https://challenges.cloudflare.com/turnstile/v0/api.js?render=explicit');
    return new Promise((resolve, reject) => {
      let widget, complete = false;
      $('verificationStatus').textContent = '';
      const finish = (error, token) => {
        if (complete) return;
        complete = true;
        $('verification').removeEventListener('close', cancelled);
        if (widget !== undefined) window.turnstile.remove(widget);
        $('verification').close();
        error ? reject(error) : resolve(token);
      };
      const cancelled = () => finish(new Error('Voting cancelled. No favorite was changed.'));
      $('verification').addEventListener('close', cancelled);
      $('verification').showModal();
      widget = window.turnstile.render('#captcha', {
        sitekey: config.turnstileSiteKey, theme: 'dark',
        callback: token => finish(null, token),
        'error-callback': () => { $('verificationStatus').textContent = 'Verification could not load. Please cancel and try again.'; },
        'expired-callback': () => { $('verificationStatus').textContent = 'Verification expired. Please try again.'; }
      });
    });
  }
  async function ensureSession() {
    if (!sessionPromise) {
      sessionPromise = (async () => {
        const {data, error} = await client.auth.getSession();
        if (error) throw error;
        if (data.session) return data.session;
        const captchaToken = await verifyVisitor();
        const result = await client.auth.signInAnonymously({options: {captchaToken}});
        if (result.error) throw result.error;
        return result.data.session;
      })().finally(() => { sessionPromise = null; });
    }
    return sessionPromise;
  }
  async function readFavorites() {
    const revision = favoriteRevision;
    const {data, error} = await client.auth.getSession();
    if (error) throw error;
    if (!data.session) { favorites = new Set(); return; }
    const result = await client.from('colortable_favorites').select('palette_slug').eq('study_id', study.study_id);
    if (result.error) throw result.error;
    if (revision === favoriteRevision) {
      favorites = new Set(result.data.map(row => row.palette_slug).filter(slug => bySlug.has(slug)));
    }
  }
  async function initialiseVoting() {
    if (initialising) return;
    $('retry').hidden = true;
    if (!configured) {
      status(local ? 'Local preview / Favorites saved on this device.' : 'Voting has not opened yet.');
      $('resultsTab').disabled = true;
      $('resultsTab').title = 'Results are not available yet';
      return;
    }
    initialising = true;
    ready = false;
    status('Connecting to voting...');
    try {
      if (!config.publishableKey.startsWith('sb_publishable_')) throw new Error('Voting configuration needs a public publishable key.');
      if (!client) {
        await loadScript('vendor/supabase-2.117.2.js');
        client = window.supabase.createClient(config.supabaseUrl, config.publishableKey, {
          auth: {persistSession: true, autoRefreshToken: true, detectSessionInUrl: false}
        });
      }
      const result = await client.from('colortable_studies').select('is_open').eq('id', study.study_id).single();
      if (result.error) throw result.error;
      await readFavorites();
      ready = result.data.is_open;
      status(ready ? 'Vote (with the star) for as many as you like' : 'Voting has closed / Results remain available.');
    } catch (error) {
      status('Voting is unavailable. Your votes have not been changed. Please retry.', true);
      $('retry').hidden = false;
      console.warn('Voting connection:', error.message);
    } finally { initialising = false; render(); }
  }
  async function toggleFavorite(slug) {
    if (!ready || pending.has(slug)) return;
    const desired = !favorites.has(slug);
    favoriteRevision++;
    if (local) {
      desired ? favorites.add(slug) : favorites.delete(slug);
      store(localKey, [...favorites]);
      render();
      return;
    }
    pending.add(slug);
    render();
    try {
      await ensureSession();
      const {error} = await client.rpc('set_colortable_favorite', {
        p_study: study.study_id, p_palette: slug, p_favorite: desired
      });
      if (error) throw error;
      desired ? favorites.add(slug) : favorites.delete(slug);
      status(desired ? 'Favorite counted.' : 'Favorite removed.');
      if (currentView === 'results') await loadResults();
    } catch (error) {
      // A response can be lost after a committed write; reconcile before another toggle.
      ready = false;
      status('Vote not confirmed. Retry to check your saved favorites.', true);
      $('retry').hidden = false;
      console.warn('Voting request:', error.message);
    } finally { pending.delete(slug); render(); }
  }

  function ranked() {
    return [...items].sort((a, b) => (counts?.get(b.slug) || 0) - (counts?.get(a.slug) || 0) || a.id - b.id);
  }
  function renderResults() {
    $('rankings').replaceChildren();
    if (!counts) return;
    let rank = 0, previousCount = null;
    ranked().forEach((item, index) => {
      const total = counts.get(item.slug) || 0;
      if (total !== previousCount) rank = index + 1;
      previousCount = total;
      if (!matches(item)) return;
      const row = element('tr');
      const cell = element('td');
      const wrapper = element('div', 'palette-cell');
      const image = element('img');
      image.src = item.thumbnails?.[$('stretch').value] || item.images[$('stretch').value];
      image.alt = '';
      image.loading = 'lazy';
      const button = element('button', '', `${String(item.id).padStart(2, '0')} / ${item.title}`);
      button.append(element('small', 'family-name', item.family));
      button.onclick = () => openPreview(item);
      wrapper.append(image, button, favoriteButton(item));
      cell.append(wrapper);
      row.append(element('td', '', String(rank)), cell, element('td', 'votes', String(total)));
      $('rankings').append(row);
    });
  }
  async function loadResults() {
    if (!client) return;
    const request = ++countRequest;
    $('resultsStatus').textContent = 'Loading totals...';
    $('refreshResults').disabled = true;
    try {
      const result = await client.rpc('colortable_counts', {p_study: study.study_id});
      if (result.error) throw result.error;
      if (request !== countRequest) return;
      counts = new Map(result.data.filter(row => bySlug.has(row.palette_slug)).map(row => [row.palette_slug, Number(row.favorites)]));
      $('resultsStatus').textContent = `Updated ${new Date().toLocaleTimeString()} / ${[...counts.values()].reduce((a,b) => a+b,0)} favorites`;
      $('exportResults').disabled = false;
      renderResults();
    } catch {
      if (request !== countRequest) return;
      $('resultsStatus').textContent = counts ? 'Could not refresh. Showing previous totals.' : 'Totals unavailable. Please retry.';
    } finally { if (request === countRequest) $('refreshResults').disabled = false; }
  }
  function showView(view) {
    currentView = view;
    $('galleryPanel').hidden = view !== 'gallery';
    $('resultsPanel').hidden = view !== 'results';
    $('galleryTab').setAttribute('aria-selected', view === 'gallery');
    $('resultsTab').setAttribute('aria-selected', view === 'results');
    $('orderLabel').hidden = view !== 'gallery';
    if (view === 'results') loadResults();
  }
  function exportResults() {
    if (!counts) return;
    const quote = value => `"${String(value).replaceAll('"', '""')}"`;
    const rows = [['number', 'palette', 'family', 'favorites'],
      ...ranked().map(item => [item.id, item.title, item.family, counts.get(item.slug) || 0])];
    const blob = new Blob([rows.map(row => row.map(quote).join(',')).join('\r\n') + '\r\n'], {type: 'text/csv;charset=utf-8'});
    const url = URL.createObjectURL(blob);
    const link = element('a');
    link.href = url;
    link.download = 'suncet-colortable-results.csv';
    link.click();
    setTimeout(() => URL.revokeObjectURL(url), 1000);
  }
  function updatePreview() {
    const stretch = $('stretch').value;
    $('previewStretch').value = stretch;
    $('referenceCaption').textContent = `01 / Current Inferno / ${stretch === 'current' ? 'fourth root' : 'asinh'}`;
    $('previewTitle').textContent = `${String(selected.id).padStart(2,'0')} / ${selected.title}`;
    $('referenceImage').src = items[0].images[stretch];
    $('selectedImage').src = selected.images[stretch];
    $('selectedImage').alt = selected.title;
    $('selectedCaption').textContent = `${selected.title} / ${stretch === 'current' ? 'fourth root' : 'asinh'}`;
    $('description').textContent = selected.description;
    $('fullImage').href = selected.images[stretch];
    $('lutLink').href = selected.lut;
    $('detailFavorite').textContent = favorites.has(selected.slug) ? 'Remove favorite' : 'Add to favorites';
    $('detailFavorite').disabled = !ready || pending.has(selected.slug);
  }
  function openPreview(item) { selected = item; updatePreview(); $('preview').showModal(); }
  function step(delta) {
    const list = currentView === 'results' ? ranked().filter(matches) : visible;
    if (!list.length) return;
    selected = list[(list.indexOf(selected) + delta + list.length) % list.length];
    updatePreview();
  }
  for (const family of new Set(items.map(item => item.family))) $('family').append(element('option', '', family));
  for (const id of ['family', 'stretch', 'order', 'favoritesOnly']) $(id).onchange = render;
  $('previewStretch').oninput = $('previewStretch').onchange = () => {
    if ($('stretch').value === $('previewStretch').value) return;
    $('stretch').value = $('previewStretch').value;
    // Request the full-resolution comparison before refreshing background thumbnails.
    updatePreview();
    render();
  };
  $('galleryTab').onclick = () => showView('gallery');
  $('resultsTab').onclick = () => showView('results');
  $('retry').onclick = initialiseVoting;
  $('refreshResults').onclick = loadResults;
  $('exportResults').onclick = exportResults;
  $('close').onclick = () => $('preview').close();
  $('next').onclick = () => step(1);
  $('previous').onclick = () => step(-1);
  $('detailFavorite').onclick = () => toggleFavorite(selected.slug);
  $('cancelVerification').onclick = () => $('verification').close();
  for (const [id, single] of [['singleMode', true], ['compareMode', false]]) {
    $(id).onclick = () => {
      $('compareImages').classList.toggle('single', single);
      $('singleMode').setAttribute('aria-pressed', single);
      $('compareMode').setAttribute('aria-pressed', !single);
    };
  }
  $('preview').addEventListener('click', event => {
    if (event.target !== $('preview')) return;
    const rect = $('preview').getBoundingClientRect();
    if (event.clientX < rect.left || event.clientX > rect.right || event.clientY < rect.top || event.clientY > rect.bottom) $('preview').close();
  });
  document.addEventListener('keydown', event => {
    if (!$('preview').open || $('verification').open) return;
    if (event.key === 'ArrowRight') { event.preventDefault(); step(1); }
    if (event.key === 'ArrowLeft') { event.preventDefault(); step(-1); }
  });
  setInterval(() => { if (currentView === 'results' && !document.hidden) loadResults(); }, 30000);
  window.addEventListener('focus', async () => {
    if (!client || !ready || pending.size) return;
    try { await readFavorites(); render(); } catch { /* Explicit retry handles connectivity failures. */ }
  });
  render();
  initialiseVoting();
})();
