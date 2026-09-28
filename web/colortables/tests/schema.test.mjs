import assert from 'node:assert/strict';
import {readFile} from 'node:fs/promises';
import {test} from 'node:test';
const {PGlite} = await import(process.env.PGLITE_MODULE || '@electric-sql/pglite');
const source = new URL('../', import.meta.url);

test('votes enforce ownership, idempotency, private records, closure, and rate limits', async () => {
  const db = new PGlite();
  try {
    await db.exec(`
      create role anon; create role authenticated;
      create schema auth;
      create table auth.users(id uuid primary key);
      create function auth.uid() returns uuid language sql stable as
        $$ select nullif(current_setting('request.jwt.claim.sub', true), '')::uuid $$;
      grant usage on schema public, auth to anon, authenticated;
      grant execute on function auth.uid() to anon, authenticated;
      insert into auth.users values
        ('00000000-0000-0000-0000-000000000001'), ('00000000-0000-0000-0000-000000000002');
    `);
    const schema = await readFile(new URL('schema.sql', source), 'utf8');
    const catalog = await readFile(new URL('catalog.sql', source), 'utf8');
    await db.exec(schema);
    await db.exec("insert into public.colortable_studies values ('suncet-frame300-v1','Study',true); insert into public.colortable_palettes(study_id,slug,number,title) values ('suncet-frame300-v1','cividis',29,'Cividis'); insert into public.colortable_favorites(study_id,palette_slug,user_id) values ('suncet-frame300-v1','cividis','00000000-0000-0000-0000-000000000001');");
    await db.exec(catalog);
    await db.exec(schema);
    await db.exec(catalog);
    const count = async slug => Number((await db.query(
      "select favorites from public.colortable_counts('suncet-frame300-v1') where palette_slug=$1", [slug])).rows[0].favorites);
    const vote = (slug, state) => db.query(
      "select public.set_colortable_favorite('suncet-frame300-v1',$1,$2)", [slug, state]);
    const user = async n => {
      await db.exec('reset role');
      await db.query("select set_config('request.jwt.claim.sub',$1,false)", [`00000000-0000-0000-0000-${String(n).padStart(12,'0')}`]);
      await db.exec('set role authenticated');
    };
    await db.exec('set role anon');
    assert.equal((await db.query("select * from public.colortable_counts('suncet-frame300-v1')")).rows.length,36);
    assert.equal((await db.query("select * from public.colortable_counts('suncet-frame300-v1') where palette_slug='cividis'")).rows.length,0);
    assert.equal(await count('poster-blue-rose'),0);
    await assert.rejects(db.query('select * from public.colortable_favorites'), /permission denied/);
    await assert.rejects(vote('poster-blue-rose',true), /permission denied/);
    await user(1);
    assert.equal((await db.query("select * from public.colortable_favorites where palette_slug='cividis'")).rows.length,1);
    await assert.rejects(vote('cividis',true), /Unknown palette/);
    await vote('tequila-sunrise',true);
    assert.equal(await count('tequila-sunrise'),1);
    await vote('tequila-sunrise',false);
    await vote('poster-blue-rose',true);
    await vote('poster-blue-rose',true);
    assert.equal(await count('poster-blue-rose'),1);
    await user(2);
    assert.equal((await db.query('select * from public.colortable_favorites')).rows.length,0);
    await vote('poster-blue-rose',false);
    assert.equal(await count('poster-blue-rose'),1);
    await vote('poster-blue-rose',true);
    assert.equal(await count('poster-blue-rose'),2);
    await assert.rejects(db.query('delete from public.colortable_favorites'), /permission denied/);
    await assert.rejects(db.query("insert into public.colortable_favorites values ('suncet-frame300-v1','coronal-ice','00000000-0000-0000-0000-000000000001',now())"), /permission denied/);
    await assert.rejects(db.query('select * from public.colortable_vote_limits'), /permission denied/);
    await assert.rejects(vote('not-a-palette',true), /Unknown palette/);
    await assert.rejects(vote('poster-blue-rose',null), /favorite state/);
    await vote('poster-blue-rose',false);
    assert.equal(await count('poster-blue-rose'),1);
    await db.exec("reset role; update public.colortable_studies set is_open=false where id='suncet-frame300-v1'; set role authenticated;");
    await assert.rejects(vote('coronal-ice',true), /Voting is closed/);
    assert.equal(await count('poster-blue-rose'),1);
    await db.exec("reset role; update public.colortable_studies set is_open=true; update public.colortable_vote_limits set requests=90,window_start=now(); set role authenticated;");
    await assert.rejects(vote('coronal-ice',true), /Too many vote changes/);
    assert.equal(await count('coronal-ice'),0);
    await db.exec("reset role; update public.colortable_vote_limits set window_start=now()-interval '2 minutes'; set role authenticated;");
    await vote('coronal-ice',true);
    assert.equal(await count('coronal-ice'),1);
    await db.query("select set_config('request.jwt.claim.sub','',false)");
    await assert.rejects(vote('coronal-ice',true), /Sign in anonymously/);
  } finally { await db.close(); }
});
