# SunCET public color table voting

The static gallery in `docs/colortables` is ready for GitHub Pages. The editable
template is `exploration_notebooks/colortable_gallery.html`; behavior lives in
`web/colortables/gallery.js`. Generated PNGs need no Python or FITS files at runtime.

## Supabase setup

1. In the project's SQL Editor, run `schema.sql`, then `catalog.sql` from this folder.
   Both can be rerun without clearing votes. The catalog has 40 active palette IDs.
2. Under Authentication > Sign In / Providers, enable anonymous sign-ins.
3. Set the project URL and `sb_publishable_...` key in `docs/colortables/config.js`.
   These values are public. Never use a secret, service-role key, or database password.
4. For public sharing, enable Cloudflare Turnstile under Authentication > Attack
   Protection and set its secret there. Put only the Turnstile site key in `config.js`.
   Allow `suncet.github.io` in Turnstile's hostname list. Without this optional
   provider configured, Supabase's signup rate limit and the SQL vote-change limit
   still apply; anonymous visitors can create new identities after clearing storage.

The Supabase Data API must be enabled with the `public` schema exposed. Anonymous
sign-ins create private browser sessions; a session is created on the first vote,
not when someone merely browses the gallery. No email address is requested.

Official references:
- https://supabase.com/docs/guides/auth/auth-anonymous
- https://supabase.com/docs/guides/getting-started/api-keys
- https://supabase.com/docs/guides/database/functions

## GitHub Pages

Enable Settings > Pages > Source: GitHub Actions in the repository. The Pages
workflow publishes `docs/colortables` on relevant pushes to `main`, or on manual
dispatch. The expected URL is https://suncet.github.io/suncet_instrument_simulator/.
The site is not live until these files are pushed and that workflow succeeds.

The committed browser configuration can be overridden with repository Actions
variables `SUPABASE_URL`, `SUPABASE_PUBLISHABLE_KEY`, and `TURNSTILE_SITE_KEY`.
The workflow rejects private keys. Do not publish the local `output/` directory:
its manifest and notes retain local source paths for reproducibility.

## Voting behavior

- Gallery order is shuffled once per browser. Number order is also available.
- Palette IDs remain stable. Options 24 and 29 are retired; their historical votes
  are retained but excluded from totals. Tequila Sunrise is option 38, and the
  SunPy STEREO/EUVI 171, 195, 284, and 304 tables are options 39-42.
- The SunCET NASA Poster reference appears below the gallery.
- Rankings appear only in Results; votes apply to palettes, independent of stretch.
- Each authenticated browser identity has at most one favorite for each palette.
- SQL derives identity from `auth.uid()`, validates the palette and open study,
  and serializes a 90-changes-per-minute limit per identity.
- Direct inserts, edits, and deletes by browser roles are denied. The vote function
  is the only mutation route; retries are idempotent.
- Visitors can read only their own vote rows. Public totals contain no voter IDs.
- Results refresh when opened and every 30 seconds while visible. CSV export contains
  only aggregate counts. A failed request never appears as a confirmed vote.
- `file:` previews with an empty configuration retain local favorites, including
  those from the original gallery. Local favorites are not silently cast as public votes.

To close voting without losing results, run:

```sql
update public.colortable_studies set is_open = false where id = 'suncet-frame300-v1';
```

Do not delete anonymous Auth users to tidy the dashboard: their votes cascade when
their accounts are deleted. These are informal preference counts, not verified people.

## Regenerate and test

```sh
MPLBACKEND=Agg python -m exploration_notebooks.colortable_options \
  --fits '/path/to/frame300.fits' --public-output docs/colortables \
  --poster '/path/to/SunCET Poster v3.png'
cd web/colortables
npm install
npx playwright install chromium
npm test
```

The published gallery now uses the pipeline's provisional exposure-normalized
Level 1 frame 300 (DN/s), with effective exposures of 0.07 s inside and 11.25 s
outside. This includes stack rejection and bit-shift normalization, not the full
Level 1 calibration chain. Both stretches use 43.2777 DN/s to the image maximum,
without a radial filter: log10 and asinh with 200 DN/s softening (local experiment
options 1 and 7). The `current` asset directory now contains log10 images.
The manifest records the source checksum, exposure metadata, and exact formulas.
The voting generator requires a Level 1 DN/s input.
Generation removes local filesystem paths from public metadata and versions image
URLs by content hash. Replacing images does not change the study ID or palette
slugs and requires no database updates. Existing votes include earlier imagery.
Configure `CHROME_EXECUTABLE` to use installed Chrome for browser tests. PGlite runs
the SQL access-policy tests against disposable PostgreSQL with mocked Supabase Auth
roles; no production votes are created by these tests.
