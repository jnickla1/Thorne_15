"""Build and check the dashboard's translation files.

site/dist/i18n/en-GB.json is the source locale. This script:
1. refreshes its keys for page text from site/dist/index.html (elements marked data-i18n /
   data-i18n-attr), so the HTML and the English strings cannot drift apart;
2. adds reference keys for text that comes with the data (method labels, notes, full names,
   warming definitions) under "dataset.*", so translators can translate them too;
3. writes site/dist/i18n/en-GB.js (the English bundle, loaded synchronously) and
   site/dist/i18n/locales.js (the registry of every i18n/<code>.json present);
4. checks that every key the scripts use exists in en-GB, and reports, for each other locale,
   missing keys, unknown keys and placeholder mismatches.

Run from the climate-dashboard directory after editing index.html, the scripts or a locale:
python3 source/build_i18n.py            (exit status 1 if a script key is missing from en-GB)
"""
import argparse
import json
import re
import sys
from html.parser import HTMLParser
from pathlib import Path

VOID = {'area', 'base', 'br', 'col', 'embed', 'hr', 'img', 'input', 'link', 'meta', 'source', 'track', 'wbr'}
SCRIPTS = ['app.js', 'tour.js', 'rate-history.js', 'member-ribbon.js', 'i18n.js']
PLACEHOLDER = re.compile(r'(?:^|[^%])\{(\w+)\}')


class Extractor(HTMLParser):
    """Collect the inner HTML of data-i18n elements and the data-i18n-attr values."""

    def __init__(self, source):
        super().__init__(convert_charrefs=False)
        self.source = source
        self.lines = [0]
        for line in source.splitlines(keepends=True):
            self.lines.append(self.lines[-1] + len(line))
        self.stack, self.found, self.attrs = [], {}, {}

    def char_offset(self):
        line, col = self.getpos()
        return self.lines[line - 1] + col

    def handle_starttag(self, tag, attrs):
        a = dict(attrs)
        end = self.char_offset() + len(self.get_starttag_text())
        if 'data-i18n-attr' in a:
            for pair in a['data-i18n-attr'].split(';'):
                attr, key = pair.split(':')
                self.attrs[key] = a[attr]
        if tag in VOID:
            return
        self.stack.append((tag, a.get('data-i18n'), end))

    def handle_endtag(self, tag):
        while self.stack:
            open_tag, key, start = self.stack.pop()
            if open_tag == tag:
                if key:
                    self.found[key] = self.source[start:self.char_offset()]
                return


def placeholders(text):
    return sorted(set(PLACEHOLDER.findall(text)))


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    root = Path(__file__).resolve().parents[1]
    parser.add_argument('--dist', type=Path, default=root/'site/dist')
    parser.add_argument('--data', type=Path, default=root/'data')
    args = parser.parse_args()
    dist, i18n = args.dist, args.dist/'i18n'
    source_path = i18n/'en-GB.json'
    en = json.loads(source_path.read_text())

    ex = Extractor((dist/'index.html').read_text())
    ex.feed(ex.source)
    page = {**ex.found, **ex.attrs}

    data_keys = {}
    for name in ['rpi_real', 'rpi_anthro']:
        d = json.loads((args.data/f'{name}.json').read_text())
        data_keys[f"dataset.definition.{d['definition']}"] = d['definition']
        for m, v in d['methods'].items():
            data_keys[f'dataset.method.{m}.label'] = v['label']
            if v.get('note'):
                data_keys[f'dataset.method.{m}.note'] = v['note']
    refs = json.loads((args.data/'method_references.json').read_text())
    for m, v in refs['methods'].items():
        if v.get('full_name'):
            data_keys[f'dataset.method.{m}.fullName'] = v['full_name']

    strings = en['strings']
    # Drop stale generated keys, then refresh them; script keys are kept as edited.
    page_prefixes = tuple(sorted({k.split('.')[0] + '.' for k in page}))
    previous_page = set(en.get('_generated_page_keys', []))
    for k in list(strings):
        if k in previous_page and k not in page:
            del strings[k]
        if k.startswith('dataset.') and k not in data_keys:
            del strings[k]
    clashes = [k for k in page if k in strings and k not in previous_page and strings[k] != page[k]]
    if clashes:
        print('Page keys clash with script keys of different text:', clashes)
        return 1
    strings.update(page)
    strings.update(data_keys)
    en['strings'] = dict(sorted(strings.items()))
    en['_generated_page_keys'] = sorted(page)
    source_path.write_text(json.dumps(en, indent=2, ensure_ascii=False) + '\n')

    used, prefixes = set(), set()
    for script in SCRIPTS:
        text = (dist/script).read_text()
        used |= set(re.findall(r"\bt\(\s*'([\w.]+)'", text))
        prefixes |= set(re.findall(r"\bt\(\s*'([\w.]+\.)'\s*\+", text))
    literal_used = {k for k in used if not k.endswith('.')}
    missing = sorted(k for k in literal_used if k not in en['strings'])
    missing += sorted(p + '*' for p in prefixes if not any(k.startswith(p) for k in en['strings']))
    status = 0
    if missing:
        print('Keys used by the scripts but missing from en-GB.json:', missing)
        status = 1

    locales = {'en-GB': {k: en['_meta'][k] for k in ['name', 'dir', 'version'] if k in en['_meta']}}
    for path in sorted(i18n.glob('*.json')):
        code = path.stem
        if code == 'en-GB':
            continue
        loc = json.loads(path.read_text())
        meta = loc.get('_meta', {})
        locales[code] = {'name': meta.get('name', code), 'dir': meta.get('dir', 'ltr'), 'version': meta.get('version', '1')}
        ls = loc.get('strings', {})
        miss = [k for k in en['strings'] if k not in ls and not k.startswith('dataset.')]
        extra = [k for k in ls if k not in en['strings']]
        bad = [k for k in ls if k in en['strings'] and placeholders(ls[k]) != placeholders(en['strings'][k])]
        print(f'{code}: {len(ls)} strings; {len(miss)} still in English; {len(extra)} unknown keys; {len(bad)} placeholder mismatches')
        for label, keys in [('unknown', extra), ('placeholder mismatch', bad)]:
            for k in keys:
                print(f'  {label}: {k}')

    bundle = {'name': en['_meta']['name'], 'dir': en['_meta']['dir'], 'strings': en['strings']}
    (i18n/'en-GB.js').write_text('/* Generated by source/build_i18n.py from en-GB.json; do not edit. */\n'
                                 'window.I18N_BUNDLES=window.I18N_BUNDLES||{};window.I18N_BUNDLES["en-GB"]='
                                 + json.dumps(bundle, ensure_ascii=False, separators=(',', ':')) + ';\n')
    (i18n/'locales.js').write_text('/* Generated by source/build_i18n.py from the i18n/*.json files present; do not edit. */\n'
                                   'window.I18N_LOCALES=' + json.dumps(locales, ensure_ascii=False) + ';\n')
    print(f"en-GB: {len(en['strings'])} strings ({len(page)} from index.html, {len(data_keys)} dataset references); "
          f"scripts use {len(literal_used)} keys; locales: {', '.join(locales)}")
    return status


if __name__ == '__main__':
    sys.exit(main())
